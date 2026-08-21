# -*- coding: utf-8 -*-
"""
baseline_manager.py - offline PubMed master database builder (pipeline Step 0).

Downloads the complete NLM PubMed Baseline (and optional daily update files)
over FTP and compiles them into a single local SQLite "master" database of
PMIDs, publication dates, and MeSH annotations that the rest of the pipeline
queries entirely offline.

Parsing is parallelized across CPU cores with a MapReduce sharding pattern: each
worker parses a chunk of gzipped XML into its own temporary SQLite shard, and
the shards are merged into the master database. To stay safe on local disks,
network-attached storage, and cloud-synced folders, the build runs in a fast
local workspace and copies back with SHA-256 checksum verification and periodic
atomic checkpoints, so an interrupted run can resume.
"""

import os
import sys
import json
import ftplib
import gzip
import sqlite3
import time
import shutil
import tempfile
import hashlib
import subprocess
import multiprocessing
import xml.etree.ElementTree as ET
from pathlib import Path
from uuid import uuid4


def _keep_on_device(path):
    """Best-effort: mark a file 'always keep on this device' so OneDrive will not
    dehydrate it into a placeholder.

    A dehydrated placeholder reads unreliably - a read-only open can return a
    stale/empty view and a copy can come back malformed, especially under disk
    pressure when OneDrive aggressively frees space. Windows/OneDrive only; a
    harmless no-op elsewhere.
    """
    if sys.platform != "win32":
        return
    try:
        subprocess.run(["attrib", "+P", "-U", str(path)],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=15)
    except Exception:
        pass

# ==========================================
# MODULE-LEVEL FUNCTIONS (Required for OS-Agnostic Multiprocessing)
# ==========================================

def _force_os_sync():
    """Forces the OS to flush write buffers to the physical disk (OS-Agnostic)."""
    if hasattr(os, 'sync'):
        os.sync()

def get_file_hash(filepath):
    """Calculates the SHA-256 digital fingerprint of a file in low-memory blocks."""
    hasher = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(8192 * 1024):
            hasher.update(chunk)
    return hasher.hexdigest()

def verified_safe_transfer(src_path: Path, dest_path: Path, max_retries=5):
    """Copies a file and strictly verifies the copy is identical byte-for-byte."""
    print(f"      [Sync] Profiling local database fingerprint...", end=" ", flush=True)
    src_size = os.path.getsize(src_path)
    src_hash = get_file_hash(src_path)
    print("Done.")

    tmp_dest = Path(str(dest_path) + ".tmp")

    for attempt in range(1, max_retries + 1):
        print(f"      [Sync] Attempt {attempt}/{max_retries} - Transferring to target storage...", end=" ", flush=True)
        try:
            if tmp_dest.exists():
                tmp_dest.unlink()

            shutil.copy2(src_path, tmp_dest)
            _force_os_sync()

            dest_size = os.path.getsize(tmp_dest)
            if dest_size != src_size:
                print(f"FAIL (Size mismatch: {dest_size} != {src_size})")
                continue

            dest_hash = get_file_hash(tmp_dest)
            if src_hash == dest_hash:
                os.replace(tmp_dest, dest_path)
                # Keep OneDrive from dehydrating the freshly written DB into a
                # placeholder, which reads back stale/empty/malformed.
                _keep_on_device(dest_path)
                print("SUCCESS (Cryptographic Match Verified).")
                return True
            else:
                print(f"FAIL (Hash mismatch detected).")
                continue

        except Exception as e:
            print(f"ERROR ({e})")
            time.sleep(5)

    raise Exception("\n[CRITICAL ERROR] Failed to transfer a healthy database after maximum retries.")

def _extract_pub_date(elem) -> str:
    """Safely extracts and formats the publication date into YYYY-MM-DD."""
    # 1. Try to get the official PubMed indexing date
    pub_date_node = elem.find('.//PubMedPubDate[@PubStatus="pubmed"]')
    if pub_date_node is not None:
        y = pub_date_node.findtext('Year')
        m = pub_date_node.findtext('Month')
        d = pub_date_node.findtext('Day')
        if y:
            m = m.zfill(2) if m else '01'
            d = d.zfill(2) if d else '01'
            return f"{y}-{m}-{d}"

    # 2. Fallback to Journal Issue PubDate
    pub_date_node = elem.find('.//PubDate')
    if pub_date_node is not None:
        y = pub_date_node.findtext('Year')
        m_str = pub_date_node.findtext('Month')
        d = pub_date_node.findtext('Day')

        if y:
            m = '01'
            if m_str:
                if m_str.isdigit():
                    m = m_str.zfill(2)
                else:
                    m_map = {
                        'jan':'01', 'feb':'02', 'mar':'03', 'apr':'04',
                        'may':'05', 'jun':'06', 'jul':'07', 'aug':'08',
                        'sep':'09', 'oct':'10', 'nov':'11', 'dec':'12'
                    }
                    m = m_map.get(m_str.lower()[:3], '01')
            d = d.zfill(2) if d else '01'
            return f"{y}-{m}-{d}"

    # 3. Default fallback if absolutely no date is found
    return "1900-01-01"

def build_local_shard(job):
    """Executes on parallel CPU cores to parse XML chunks into local SQLite shards.

    `job` is (chunk, workspace). The workspace travels with the job because these
    run in spawned processes on Windows, where a module-level setting in the
    parent is not inherited.
    """
    local_file_chunk, workspace = job
    shard_path = Path(workspace) / f"shard_{uuid4().hex}.db"
    conn = sqlite3.connect(shard_path)

    conn.execute("PRAGMA journal_mode = OFF;")
    conn.execute("PRAGMA synchronous = OFF;")
    # UPDATED SCHEMA: Now includes pub_date
    conn.execute("CREATE TABLE shard_data (pmid INTEGER, pub_date TEXT, mesh_terms TEXT, source_file TEXT)")

    processed_names = []
    total_articles = 0
    cursor = conn.cursor()

    for filepath in local_file_chunk:
        batch = []
        source_filename = Path(filepath).name
        try:
            with gzip.open(filepath, 'rb') as f:
                # Track the root to drop parsed siblings; otherwise each ~30k-article
                # baseline file accumulates in memory until the whole file is done.
                context = ET.iterparse(f, events=('start', 'end'))
                _, root = next(context)
                for event, elem in context:
                    if event != 'end' or elem.tag != 'PubmedArticle':
                        continue

                    pmid_node = elem.find('.//PMID')
                    mesh_list = elem.find('.//MeshHeadingList')

                    if pmid_node is not None and mesh_list is not None:
                        pmid = int(pmid_node.text)
                        pub_date = _extract_pub_date(elem)
                        terms = [f"{'*' if d.get('MajorTopicYN') == 'Y' else ''}{d.text}"
                                 for d in mesh_list.findall('.//DescriptorName')]

                        if terms:
                            batch.append((pmid, pub_date, ";".join(terms), source_filename))

                    elem.clear()
                    root.clear()

            if batch:
                cursor.executemany("INSERT INTO shard_data (pmid, pub_date, mesh_terms, source_file) VALUES (?, ?, ?, ?)", batch)
                total_articles += len(batch)

            processed_names.append(source_filename)
        except Exception as e:
            print(f"\n  [!] Error parsing {source_filename}: {e}")

    conn.commit()
    conn.close()
    return str(shard_path), processed_names, total_articles


# ==========================================
# MAIN MANAGER CLASS
# ==========================================

class PubMedBaselineManager:
    """Downloads the NLM PubMed XML archives and compiles them into the local master database."""

    def __init__(self, raw_data_dir: Path, master_db_path: Path, workspace_dir=None):
        """Set up the raw/download directories, fast local workspace, and FTP endpoints.

        `workspace_dir` is where the build does its scratch work. It defaults to
        the system temp directory, which is the fastest disk on most machines but
        is also the system drive - and the workspace holds a full working copy of
        the master database, so on a small C: a build can run it out of space.
        Pointing it at a roomier drive is the remedy.
        """
        self.raw_data_dir = Path(raw_data_dir)
        self.master_db_path = Path(master_db_path)

        self.baseline_dir = self.raw_data_dir / "pubmed_baseline"
        self.updates_dir = self.raw_data_dir / "pubmed_updates"

        self.baseline_dir.mkdir(parents=True, exist_ok=True)
        self.updates_dir.mkdir(parents=True, exist_ok=True)

        # High-Speed Local Workspace (OS-Agnostic). The named subfolder is always
        # appended, so a configured location is never itself treated as ours -
        # the uninstaller removes this folder whole.
        _ws_base = Path(workspace_dir) if workspace_dir else Path(tempfile.gettempdir())
        self.local_workspace = _ws_base / "mesh_etl_workspace"
        self.local_workspace.mkdir(exist_ok=True)
        self.local_db_path = self.local_workspace / "local_active_master.db"
        self.local_xml_dir = self.local_workspace / "xml_buffer"
        self.local_xml_dir.mkdir(exist_ok=True)

        self.ftp_host = "ftp.ncbi.nlm.nih.gov"
        self.baseline_ftp_path = "/pubmed/baseline/"
        self.updates_ftp_path = "/pubmed/updatefiles/"

        self.chunk_size = 25
        self.checkpoint_interval = 4  # Heavy DB checkpoint every 100 files (plus an early one after the first block)

        # Progress manifest written alongside the master DB. It records parse
        # progress independently of the database file, so a later quarantine of a
        # corrupt DB never erases the record of how far the build got.
        self.progress_path = Path(str(self.master_db_path) + ".progress.json")

    def _get_ftp_file_list(self, ftp_path: str) -> list:
        """Return the list of .xml.gz filenames available at the given NCBI FTP path."""
        try:
            ftp = ftplib.FTP(self.ftp_host)
            ftp.login()
            ftp.cwd(ftp_path)
            files = [f for f in ftp.nlst() if f.endswith(".xml.gz")]
            ftp.quit()
            return files
        except Exception as e:
            print(f"  [!] FTP Connection Error while listing files: {e}")
            return []

    def _download_files(self, ftp_path: str, local_dir: Path, file_list: list, limit: int = None):
        """Download the listed files to local_dir, skipping existing ones and reconnecting/retrying on failure."""
        files_to_download = file_list[:limit] if limit else file_list
        ftp = None

        def connect():
            nonlocal ftp
            if ftp:
                try:
                    ftp.quit()
                except:
                    pass
            ftp = ftplib.FTP(self.ftp_host)
            ftp.login()
            ftp.cwd(ftp_path)

        connect()
        total = len(files_to_download)
        print(f"  Starting download sequence for {total} files...")
        seq_start = time.time()

        for i, filename in enumerate(files_to_download, 1):
            local_filepath = local_dir / filename
            already_present = local_filepath.exists() and local_filepath.stat().st_size > 1024

            if not already_present:
                success = False
                for attempt in range(3):
                    try:
                        with open(local_filepath, 'wb') as f:
                            ftp.retrbinary(f"RETR {filename}", f.write)
                        success = True
                        break
                    except Exception as e:
                        print(f"    [!] Error on {filename} (Attempt {attempt+1}/3): {e}")
                        if local_filepath.exists():
                            local_filepath.unlink()
                        time.sleep(2)
                        connect()
                if not success:
                    print(f"    [!] FAILED to download {filename}. Skipping.")

            if i % 50 == 0 or i == total:
                print(f"    -> Downloaded {i}/{total} files... [{(time.time() - seq_start)/60:.1f} min]")

        if ftp:
            try:
                ftp.quit()
            except:
                pass

    def run_downloads(self, include_baseline=True, include_updates=False, limit=None):
        """Fetch the PubMed baseline and/or daily-update XML archives from NCBI over FTP."""
        print("\n" + "<"*30 + ">"*30)
        print("<<< NLM PubMed XML Data Download >>>")
        print("<"*30 + ">"*30)
        if include_baseline:
            print(f"  Fetching Baseline index from {self.ftp_host}...")
            baseline_files = self._get_ftp_file_list(self.baseline_ftp_path)
            if baseline_files:
                self._download_files(self.baseline_ftp_path, self.baseline_dir, baseline_files, limit)

        if include_updates:
            print(f"\n  Fetching Updates index from {self.ftp_host}...")
            update_files = self._get_ftp_file_list(self.updates_ftp_path)
            if update_files:
                self._download_files(self.updates_ftp_path, self.updates_dir, update_files, limit)

    def _save_progress_manifest(self, parsed_files, total_files, complete=False, record_count=None):
        """Persist a small JSON record of parse progress next to the master DB.

        Decoupled from the database file so progress stays visible and auditable
        even if the DB is later quarantined as corrupt. On completion the verified
        physical article count is stored so it can be cross-checked later.
        """
        try:
            payload = {
                "updated": time.strftime("%Y-%m-%dT%H:%M:%S"),
                "parsed_count": len(parsed_files),
                "total_files": total_files,
                "complete": complete,
                "record_count": record_count,
                "parsed_files": sorted(parsed_files),
            }
            tmp = Path(str(self.progress_path) + ".tmp")
            with open(tmp, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=1)
            os.replace(tmp, self.progress_path)
        except Exception as e:
            print(f"      [!] Could not write progress manifest: {e}")

    def _reset_local_db(self):
        """Delete a stale local-workspace DB (and its WAL/SHM sidecars) before a
        fresh build, with retries and a clear error if another pipeline instance
        (or an orphaned worker process) still holds the file open."""
        for suffix in ("", "-wal", "-shm"):
            p = Path(str(self.local_db_path) + suffix)
            for _ in range(3):
                try:
                    if p.exists():
                        p.unlink()
                    break
                except PermissionError:
                    time.sleep(1.0)
            else:
                raise RuntimeError(
                    "The ETL workspace is locked by another process - a previous pipeline "
                    f"run (or its orphaned Python workers) is still holding '{p.name}'. "
                    "Only one pipeline build can run at a time. Close any other running "
                    "pipeline, end stray python.exe processes, then delete this folder and "
                    f"re-run:\n    {self.local_workspace}"
                )

    def _checkpoint_to_target(self, conn):
        """Write the live DB to the master target as one self-contained file.

        Folds the WAL into the main database file and switches it to rollback
        (DELETE) journal mode in place, then verify-transfers that single file to
        the (possibly OneDrive-synced) target. Rollback mode matters because a
        WAL-mode file cannot be opened read-only on a OneDrive path (the -shm
        shared-memory file cannot be created there), which made the wizard and the
        relevance step misread a healthy DB as corrupt.

        Doing the consolidation in place avoids making a second full-size .snapshot
        copy in the workspace - halving the checkpoint footprint, which matters on
        a near-full disk. If the WAL cannot be fully consolidated (e.g. a reader
        holds it), it falls back to the SQLite backup API so a checkpoint is never
        written from an unflushed WAL.
        """
        conn.commit()
        conn.execute("PRAGMA wal_checkpoint(TRUNCATE);")
        mode = conn.execute("PRAGMA journal_mode=DELETE;").fetchone()[0]
        conn.commit()

        if mode == "delete":
            # Live file is now complete and self-contained: copy it directly.
            verified_safe_transfer(self.local_db_path, self.master_db_path)
        else:
            # Rare: WAL could not be consolidated in place. Fall back to the
            # backup API so we never transfer an incomplete main file.
            snap = Path(str(self.local_db_path) + ".snapshot")
            if snap.exists():
                snap.unlink()
            dst = sqlite3.connect(str(snap))
            try:
                conn.backup(dst)
                dst.execute("PRAGMA journal_mode=DELETE;")
                dst.commit()
            finally:
                dst.close()
            verified_safe_transfer(snap, self.master_db_path)
            if snap.exists():
                snap.unlink()

        # Restore WAL mode for fast continued bulk loading.
        conn.execute("PRAGMA journal_mode=WAL;")

    def _verify_built_database(self, expected_files):
        """Fast sanity check of the finished build before declaring success.

        The hard checks run against the LOCAL workspace database, not the OneDrive
        target. The target is a hash-verified byte copy of the local file, but
        OneDrive can dehydrate the freshly written target into a placeholder whose
        read-only view is transiently empty/partial - which previously false-failed
        a perfectly complete build. The local file is on a normal disk and is the
        true build output, so it is authoritative.

        Cheap by design: one open plus two counts (the article count walks the main
        table, so gross corruption surfaces here). Raises RuntimeError only when the
        build itself is genuinely empty/incomplete. Returns the article count.
        """
        print("\n" + "<"*30 + ">"*30)
        print("<<< Phase 5: Build Verification >>>")
        print("<"*30 + ">"*30)

        local = str(self.local_db_path)
        try:
            conn = sqlite3.connect(local, timeout=60.0)
            try:
                columns = [r[1] for r in conn.execute("PRAGMA table_info(master_mesh_annotations)").fetchall()]
                if not columns:
                    raise RuntimeError("[VERIFY FAILED] master_mesh_annotations table missing in the build output.")
                record_count = conn.execute("SELECT count(*) FROM master_mesh_annotations").fetchone()[0]
                parsed_count = conn.execute("SELECT count(*) FROM parsed_files").fetchone()[0]
            finally:
                conn.close()
        except sqlite3.DatabaseError as e:
            raise RuntimeError(f"[VERIFY FAILED] The build output is unreadable/corrupt: {e}. Rebuild.")

        size_gb = os.path.getsize(self.master_db_path) / (1024 ** 3) if os.path.exists(self.master_db_path) else 0.0
        print(f"  Build output (local) .. OK")
        print(f"  Target size ........... {size_gb:.2f} GB")
        print(f"  Files parsed .......... {parsed_count:,} / {expected_files:,} expected")
        print(f"  Articles in DB ........ {record_count:,}")

        problems = []
        if record_count <= 0:
            problems.append("the database contains zero articles")
        if expected_files and parsed_count < expected_files:
            problems.append(f"only {parsed_count}/{expected_files} input files are recorded as parsed")
        if problems:
            raise RuntimeError(
                "[VERIFY FAILED] " + "; ".join(problems) +
                ". The build is incomplete - do not proceed to scoring."
            )

        # Best-effort, NON-FATAL: confirm the OneDrive target is also readable
        # read-only (how Step 3 reads it), retrying to let OneDrive finish
        # hydrating/syncing the freshly written file. A failure here never fails
        # the build - the output is already verified above.
        _keep_on_device(self.master_db_path)
        target = str(self.master_db_path)
        readable = False
        for _ in range(6):
            try:
                c = sqlite3.connect(f"file:{target}?mode=ro", uri=True, timeout=30.0)
                if c.execute("PRAGMA table_info(master_mesh_annotations)").fetchall():
                    readable = True
                c.close()
            except sqlite3.DatabaseError:
                pass
            if readable:
                break
            time.sleep(5)

        if readable:
            print("  Target read-only ...... OK (this is how relevance scoring reads it)")
        else:
            print("  [!] NOTE: the OneDrive target is not yet readable read-only - OneDrive is")
            print("      still syncing/hydrating it. The build IS verified complete; if Step 3")
            print("      later reports it missing, let OneDrive finish (right-click > 'Always")
            print("      keep on this device') and re-run.")

        print("  [+] Verification PASSED - build output is complete and non-empty.")
        return record_count

    def _validate_synced_db(self):
        """Confirm the workspace copy synced from the OneDrive target is a readable
        database before resuming.

        OneDrive can serve a dehydrated/partial view of a large file (especially
        under disk pressure), producing a malformed local copy that would otherwise
        crash mid-parse with a cryptic "database disk image is malformed". Reading
        the annotations table walks its b-tree, so truncation/corruption surfaces
        here as a clear, actionable error instead.
        """
        try:
            conn = sqlite3.connect(str(self.local_db_path), timeout=60.0)
            try:
                cols = [r[1] for r in conn.execute("PRAGMA table_info(master_mesh_annotations)").fetchall()]
                if cols:
                    conn.execute("SELECT count(*) FROM master_mesh_annotations").fetchone()
            finally:
                conn.close()
        except sqlite3.DatabaseError as e:
            raise RuntimeError(
                f"[RESUME ABORTED] The existing master database could not be read after syncing "
                f"from OneDrive ({e}). OneDrive most likely served a dehydrated or partially-synced "
                "copy. Recovery: in Explorer, right-click the file in data/raw and choose 'Always "
                "keep on this device', wait for it to finish downloading, then re-run. If it stays "
                "unreadable the file is genuinely corrupt and must be rebuilt (delete it and re-run)."
            )
        if not cols:
            raise RuntimeError(
                "[RESUME ABORTED] The synced master database has no master_mesh_annotations table - "
                "OneDrive may have served a stale/empty placeholder. Ensure the file is fully "
                "downloaded ('Always keep on this device') and re-run."
            )

    def compile_database(self):
        """Parse all downloaded XML into the master SQLite database in parallel, with resumable checkpoints."""
        print("\n" + "<"*30 + ">"*30)
        print("<<< Phase 0: Environment Diagnostics & Bootstrapping >>>")
        print("<"*30 + ">"*30)

        os.makedirs(os.path.dirname(self.master_db_path), exist_ok=True)

        if self.master_db_path.exists():
            print("  Existing target database detected. Syncing to Local Workspace...")
            _keep_on_device(self.master_db_path)
            verified_safe_transfer(self.master_db_path, self.local_db_path)
            self._validate_synced_db()
        else:
            print("  No target database found. Initializing blank Local Workspace.")
            self._reset_local_db()

        conn = sqlite3.connect(self.local_db_path)
        conn.execute("PRAGMA journal_mode = WAL;")

        # --- SCHEMA SELF-HEALING CHECK ---
        cursor = conn.cursor()
        cursor.execute("PRAGMA table_info(master_mesh_annotations)")
        cols = cursor.fetchall()
        if cols and len(cols) < 4:
            print("  [!] Outdated 3-column schema detected. Purging and upgrading to 4-column schema...")
            conn.execute("DROP TABLE IF EXISTS master_mesh_annotations")
            conn.execute("DROP TABLE IF EXISTS parsed_files")

        # 4-column schema. pmid is the PRIMARY KEY so the INSERT OR REPLACE during
        # sharded load de-duplicates NLM baseline overlaps automatically and pmid
        # joins are indexed without a separate index build.
        conn.execute("CREATE TABLE IF NOT EXISTS master_mesh_annotations (pmid INTEGER PRIMARY KEY, pub_date TEXT, mesh_terms TEXT, source_file TEXT)")
        conn.execute("CREATE TABLE IF NOT EXISTS parsed_files (filename TEXT PRIMARY KEY)")

        cursor.execute("DROP INDEX IF EXISTS idx_pmid")
        cursor.execute("SELECT filename FROM parsed_files")
        completed_files = {row[0] for row in cursor.fetchall()}

        # If a progress manifest from a prior run exists but the database holds no
        # parsed files, the DB was reset/rebuilt - explain why a full re-parse follows
        # (the XML downloads are reused; only parsing repeats).
        if self.progress_path.exists() and not completed_files:
            try:
                prev = json.load(open(self.progress_path, encoding="utf-8"))
                if prev.get("parsed_count", 0) > 0 and not prev.get("complete", False):
                    print(f"  [i] Progress manifest from {prev.get('updated','?')} recorded "
                          f"{prev.get('parsed_count')} parsed files, but the database was reset "
                          f"(corrupt/rebuilt) so that parsed data is gone - re-parsing is required.")
                    print("      (Your downloaded XML files are reused; nothing is re-downloaded.)")
            except Exception:
                pass

        all_files = []
        if self.baseline_dir.exists():
            all_files.extend(list(self.baseline_dir.glob("*.xml.gz")))
        if self.updates_dir.exists():
            all_files.extend(list(self.updates_dir.glob("*.xml.gz")))

        pending_files = [f for f in all_files if f.name not in completed_files]

        if not pending_files:
            print("  All files parsed! Proceeding to Indexing Phase...")
        else:
            print(f"  Found {len(all_files)} total files. {len(completed_files)} already complete.")
            print(f"  Resuming with {len(pending_files)} files...")

            chunks = [pending_files[i:i + self.chunk_size] for i in range(0, len(pending_files), self.chunk_size)]
            cores = max(1, multiprocessing.cpu_count() - 1)

            print("\n" + "<"*30 + ">"*30)
            print(f"<<< Phase 1 & 2: Distributed Parsing & Aggregation ({cores} Cores) >>>")
            print("<"*30 + ">"*30)

            start_time = time.time()
            global_articles = 0

            # One worker pool for the whole run; re-creating it per chunk paid
            # process-spawn overhead on every block (dozens of times per ETL).
            pool = multiprocessing.Pool(cores)
            try:
                for chunk_idx, chunk in enumerate(chunks, 1):
                    local_chunk_paths = []
                    for raw_file in chunk:
                        local_dest = self.local_xml_dir / raw_file.name
                        shutil.copy2(raw_file, local_dest)
                        local_chunk_paths.append(local_dest)

                    sub_chunks = [local_chunk_paths[i::cores] for i in range(cores)]
                    sub_chunks = [(sc, str(self.local_workspace))
                                  for sc in sub_chunks if sc]

                    for shard_path_str, parsed_names, article_count in pool.imap_unordered(build_local_shard, sub_chunks):
                        cursor.execute("ATTACH DATABASE ? AS temp_shard", (shard_path_str,))
                        cursor.execute("BEGIN TRANSACTION;")

                        cursor.execute("INSERT OR REPLACE INTO master_mesh_annotations (pmid, pub_date, mesh_terms, source_file) SELECT pmid, pub_date, mesh_terms, source_file FROM temp_shard.shard_data")

                        for fname in parsed_names:
                            cursor.execute("INSERT OR REPLACE INTO parsed_files (filename) VALUES (?)", (fname,))

                        conn.commit()
                        cursor.execute("DETACH DATABASE temp_shard")
                        os.remove(shard_path_str)
                        global_articles += article_count

                    for local_file in local_chunk_paths:
                        if local_file.exists():
                            local_file.unlink()

                    cursor.execute("SELECT count(*) FROM parsed_files")
                    total_done = cursor.fetchone()[0]
                    elapsed = time.time() - start_time
                    print(f"  -> Processed Block {chunk_idx}/{len(chunks)}. (Total: {total_done}/{len(all_files)}) [+ {global_articles:,} articles] [{elapsed/60:.1f} min]")

                    # <<< ATOMIC CHECKPOINT >>>
                    # Checkpoint after the very first block (cheap - the DB is tiny)
                    # so a resumable copy exists within minutes, then on the normal
                    # cadence. This prevents an interrupt-before-first-checkpoint from
                    # leaving a tableless DB that gets flagged corrupt and wiped.
                    if chunk_idx == 1 or chunk_idx % self.checkpoint_interval == 0:
                        print(f"  <<< Executing Routine Checkpoint ({total_done} files) >>>")
                        conn.commit()
                        cursor.execute("SELECT filename FROM parsed_files")
                        parsed_now = {r[0] for r in cursor.fetchall()}
                        # Consistent snapshot via the backup API (keeps conn open).
                        self._checkpoint_to_target(conn)
                        self._save_progress_manifest(parsed_now, len(all_files))
            finally:
                # terminate() (not close()) so workers can't linger as orphaned
                # processes holding the workspace open if the run is interrupted.
                pool.terminate()
                pool.join()

        print("\n" + "<"*30 + ">"*30)
        print("<<< Phase 3: Post-Load Optimization & Indexing >>>")
        print("<"*30 + ">"*30)
        # De-duplication and pmid indexing are inherent to the schema: pmid is an
        # INTEGER PRIMARY KEY and rows are merged with INSERT OR REPLACE, so
        # duplicate pmids never exist and the primary key already provides the
        # lookup index. The previous explicit dedup-DELETE (a full-table scan) and
        # the separate UNIQUE INDEX duplicated that work for no benefit.
        print("  De-duplication and pmid indexing handled by the primary key (no extra pass needed).")

        print("\n" + "<"*30 + ">"*30)
        print("<<< Phase 4: Final Network Transfer >>>")
        print("<"*30 + ">"*30)
        conn.commit()
        conn.execute("PRAGMA wal_checkpoint(TRUNCATE);")
        # VACUUM removed: it rewrote the entire multi-GB database and required an
        # equal amount of free temp space, risking a disk-full failure at the very
        # end. With primary-key dedup there is negligible free space to reclaim.
        cursor.execute("SELECT filename FROM parsed_files")
        final_parsed = {r[0] for r in cursor.fetchall()}

        # Consistent snapshot of the finished DB to the target (backup API).
        self._checkpoint_to_target(conn)
        conn.close()

        # Verify the finished target the way downstream steps read it (read-only)
        # before declaring success; raises if the build is unusable/incomplete.
        record_count = self._verify_built_database(len(all_files))
        self._save_progress_manifest(final_parsed, len(all_files), complete=True, record_count=record_count)

        print("\n" + "<"*30 + ">"*30)
        print("MASTER DATABASE COMPILATION COMPLETE")
        print("<"*30 + ">"*30 + "\n")