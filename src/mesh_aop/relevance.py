# -*- coding: utf-8 -*-
"""
relevance.py - mean relevancy scoring / semantic re-ranking (pipeline Step 3).

Re-scores the network's nodes by how strongly they are supported across the
wider literature within a chosen time window, querying the local master database
so the calculation needs no API calls.

For each surviving node it aggregates per-article relevance scores into a
Mean Relevancy Score (MRS) - the plain mean article-relevance of a term:

    ARS(A) = D_a^-1/2 . X . D_t^-1/2 . w        (symmetrically normalised)
    MRS(i) = mean(ARS(A) for A in P_i)

then writes both the per-article scores and the relevance-annotated network to disk.

ARS uses the SYMMETRICALLY NORMALISED projection: each seed weight is damped by
the square root of that term's in-window document frequency, and each article
total by the square root of how many network terms it carries. A term therefore
cannot dominate merely by being common, nor an article merely by being broadly
indexed. This replaces the earlier plain sum over a constant denominator
(w-sum / total-w), which the projection comparison showed to be the weakest of
the projections tested; seeding MRS from it left MRS unable to beat the raw
centrality it was built from. Consequences of the change:
  - document frequency must be known before any article can be scored, so the
    corpus is now traversed TWICE (pass 1 counts, pass 2 scores). The scan
    dominates this step, so expect roughly double the previous runtime.
  - the resulting scores are ~1e-6 rather than ~1e-2, so each score column is
    max-scaled to [0, 1] after the pass. That is rank-preserving (the same
    convention already applied to MRS) and keeps ars_score on the range
    secondary_analysis assumes when blending it against a citation score.

Two multiplicative adjustments that earlier versions layered on top of this mean were
both dropped, leaving the mean itself (hence "non-adjusted for context"):
  - an evidence-volume / corpus-frequency factor, log10(|P_i| + 1); and
  - an information-content penalty, -log10(G_t / N_global), applied via a harmonic mean.
The mean-ARS term is already frequency-neutral by construction (averaging over P_i
dilutes common terms by their many low-scoring articles), so both factors double-corrected
for frequency: they inverted the association with corpus frequency and pushed discrimination
of externally attested terms down to chance. Frequency confounding is instead addressed by
conditional analysis at validation time, not by a multiplicative penalty here.
"""

import os
import math
import json
import sqlite3
from collections import defaultdict

from tqdm import tqdm

# Relative package imports
from .data_ops import parse_date_robust
from .baseline_manager import _keep_on_device


def _open_master_readonly(master_db_path):
    """Open the master DB for reading, resilient to OneDrive dehydration.

    Pins the file so OneDrive keeps it on-device, then opens read-only. If the
    read-only view shows no schema - a dehydrated placeholder reads back empty -
    it re-opens with a normal connection, which forces OneDrive to hydrate the
    file and exposes the real content. No writes are issued either way, and on a
    normal (already-hydrated) file the read-only path is used unchanged.
    """
    _keep_on_device(master_db_path)
    conn = sqlite3.connect(f'file:{master_db_path}?mode=ro', uri=True, timeout=120.0)
    if conn.execute("PRAGMA table_info(master_mesh_annotations)").fetchall():
        return conn
    conn.close()
    return sqlite3.connect(str(master_db_path), timeout=120.0)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# HELPER FUNCTIONS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def _extract_base_terms(mesh_str: str) -> set:
    """Safely extracts base MeSH terms by stripping asterisks and subheadings."""
    if not mesh_str or not isinstance(mesh_str, str):
        return set()
    bases = set()
    for t in mesh_str.split(';'):
        base = t.split('/')[0].lstrip('*').strip()
        if base:
            bases.add(base)
    return bases


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CORE PIPELINE OPERATIONS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def run_mean_relevancy_scoring(input_nodes_file: str, output_nodes_file: str,
                               master_db_path: str, relevance_db_path: str,
                               id_key: str,
                               start_date_param: str, end_date_param: str,
                               weightings: list,
                               entrez_email: str = "", entrez_api_key: str = "",
                               calculate_full_centrality: bool = True):
    """
    Scores each network node by how strongly the wider literature supports it within the
    chosen time window, querying the local Master DB so no API calls are made.

    `weightings` is a list of (node_attribute, output_attribute) pairs, one per node
    weighting to score. Each produces a per-article `score_<node_attribute>` column and a
    per-node `<output_attribute>` MRS (Mean Relevancy Score). Any number of weightings are
    scored from a single corpus pass - the scan dominates the cost, so extras are cheap.
    """
    if not weightings:
        raise ValueError("No node weightings supplied to relevance scoring.")

    if not calculate_full_centrality:
        print("\n[!] NOTICE: Betweenness centrality was estimated from a sampled subset of nodes "
              "rather than computed exactly. Eigenvector and PageRank centrality are unaffected, "
              "and MRS still reflects topological relevance.\n")

    print(f"\n<<< Loading Seed Terms from {os.path.basename(input_nodes_file)} >>>")
    n_w = len(weightings)
    seed_weights = [dict() for _ in range(n_w)]
    total_weight = [0.0] * n_w
    seed_terms = set()

    try:
        with open(input_nodes_file, 'r') as f:
            nodes = json.load(f).get('elements', {}).get('nodes', [])
        for node in nodes:
            data = node.get('data', {})
            term = data.get(id_key)
            if not term:
                continue
            seed_terms.add(term)
            for i, (w_key, _) in enumerate(weightings):
                w = data.get(w_key)
                if isinstance(w, (int, float)):
                    seed_weights[i][term] = float(w)
                    total_weight[i] += float(w)
        print(f"Loaded {len(seed_terms)} unique seed terms for relevance scoring.")
        print(f"  Weightings ({n_w}): " + ", ".join(w for w, _ in weightings))
        for i, (w_key, _) in enumerate(weightings):
            if total_weight[i] <= 0:
                print(f"  [!] WARNING: '{w_key}' is absent or zero on every node; "
                      f"its scores will be 0.")
    except Exception as e:
        raise RuntimeError(f"ERROR loading seed terms: {e}")

    if not seed_terms:
        return

    # <<< Enforce ISO Dates for High-Speed SQLite Query >>>
    start_iso = parse_date_robust(start_date_param).strftime('%Y-%m-%d')
    end_iso = parse_date_robust(end_date_param).strftime('%Y-%m-%d')

    print("\n" + "<"*30 + ">"*30)
    print("<<< Calculating Mean Relevancy (Local Data Lake) >>>")
    print("<"*30 + ">"*30)
    print(f"  Target Date Range: {start_iso} to {end_iso}")

    try:
        # Read-only to avoid lock issues, with a fallback that rehydrates the file
        # if OneDrive has dehydrated it into a placeholder (see _open_master_readonly).
        master_conn = _open_master_readonly(master_db_path)
        m_cursor = master_conn.cursor()

        # Schema compatibility check
        m_cursor.execute("PRAGMA table_info(master_mesh_annotations)")
        columns = [info[1] for info in m_cursor.fetchall()]
        if 'pub_date' not in columns:
            raise RuntimeError("CRITICAL ERROR: 'pub_date' column is missing from the Master Database. You MUST run Step 0 again to rebuild the database with the new 4-column schema before proceeding.")

        # In-window count drives the ARS pass (which articles get scored / form
        # each term's P_i) and the zero-check.
        m_cursor.execute("SELECT count(*) FROM master_mesh_annotations WHERE pub_date BETWEEN ? AND ?", (start_iso, end_iso))
        total_articles_in_range = m_cursor.fetchone()[0]

        if total_articles_in_range == 0:
            print("  [!] Warning: Zero articles found in Master Database for this date range.")
            master_conn.close()
            return

        print(f"  Found {total_articles_in_range:,} articles in window; scoring them now...")

        aggregators = [defaultdict(lambda: {'sum': 0.0, 'count': 0}) for _ in range(n_w)]
        # Peak score per weighting, so the normalised ARS can be max-scaled to
        # [0, 1] once the pass finishes (see the rescale step below).
        max_score = [0.0] * n_w

        # Article scores are streamed straight into SQLite in batches rather than
        # accumulated in a list and converted to a DataFrame at the end. On a corpus
        # this size the old approach held every row in RAM and then copied all of it
        # again for the DataFrame, which exhausted memory once several weightings were
        # requested; peak usage here is one batch regardless of corpus or weighting count.
        os.makedirs(os.path.dirname(relevance_db_path), exist_ok=True)
        score_cols = [f"score_{w_key}" for w_key, _ in weightings]
        out_conn = sqlite3.connect(relevance_db_path)
        out_conn.execute("PRAGMA journal_mode=WAL")
        out_conn.execute("PRAGMA synchronous=OFF")
        out_conn.execute("DROP TABLE IF EXISTS article_relevance_scores")
        out_conn.execute(
            "CREATE TABLE article_relevance_scores (pmid TEXT, "
            + ", ".join(f'"{c}" REAL' for c in score_cols)
            + ", contributing_seeds TEXT, pub_date TEXT)"
        )
        insert_sql = (
            "INSERT INTO article_relevance_scores VALUES ("
            + ", ".join("?" * (len(score_cols) + 3)) + ")"
        )
        batch = []
        rows_written = 0

        # <<< PASS 1 of 2: document frequency per seed term >>>
        # The symmetrically normalised projection divides each seed weight by the
        # square root of the number of in-window articles carrying that term, so
        # every df must be known before ANY article can be scored. That forces a
        # second corpus pass; the scan dominates Step 3, so expect roughly double
        # the previous runtime for this step.
        print("  Pass 1/2: counting document frequency per seed term ...")
        doc_freq = defaultdict(int)
        m_cursor.execute(
            "SELECT mesh_terms FROM master_mesh_annotations "
            "WHERE pub_date BETWEEN ? AND ?", (start_iso, end_iso))
        pbar_df = tqdm(total=total_articles_in_range, desc="Counting terms")
        while True:
            rows = m_cursor.fetchmany(100000)
            if not rows:
                break
            for (mesh_terms_str,) in rows:
                if not mesh_terms_str:
                    continue
                for term in _extract_base_terms(mesh_terms_str).intersection(seed_terms):
                    doc_freq[term] += 1
            pbar_df.update(len(rows))
        pbar_df.close()
        # A term absent from the window would divide by zero; 1 leaves its weight
        # untouched, which is the correct no-information default.
        inv_sqrt_df = {t: 1.0 / math.sqrt(doc_freq.get(t, 0) or 1) for t in seed_terms}
        print(f"    document frequencies collected for {len(doc_freq):,} seed terms")

        # <<< PASS 2 of 2: article scoring >>>
        # Only in-window articles are needed: they form each term's P_i, which drives
        # the mean ARS. Filtering in SQL rather than scanning all 33M records and
        # discarding in Python is what dropping the global-specificity term buys back.
        m_cursor.execute(
            "SELECT pmid, pub_date, mesh_terms FROM master_mesh_annotations "
            "WHERE pub_date BETWEEN ? AND ?", (start_iso, end_iso))
        chunk_size = 100000
        pbar = tqdm(total=total_articles_in_range, desc="Scoring articles")

        while True:
            rows = m_cursor.fetchmany(chunk_size)
            if not rows:
                break

            for row in rows:
                pmid, pub_date, mesh_terms_str = row
                if not mesh_terms_str:
                    continue

                article_terms = _extract_base_terms(mesh_terms_str)
                matching_seeds = article_terms.intersection(seed_terms)

                if not matching_seeds:
                    continue

                # Symmetric degree normalisation: 1/sqrt(number of network terms the
                # article carries), applied once per article.
                inv_sqrt_deg_a = 1.0 / math.sqrt(len(matching_seeds))

                record = [pmid]
                for i, (w_key, _) in enumerate(weightings):
                    score = 0.0
                    if total_weight[i] > 0:
                        # ARS = D_a^-1/2 . X . D_t^-1/2 . w
                        # Each seed weight is damped by the square root of its own
                        # document frequency, so a term cannot dominate simply by
                        # being common, and the article total is damped by how many
                        # network terms it carries, so long articles cannot dominate
                        # by breadth alone. This is the projection that performs best
                        # in the projection comparison; the previous plain sum over
                        # a constant denominator is its unnormalised counterpart.
                        score = inv_sqrt_deg_a * sum(
                            seed_weights[i].get(term, 0.0) * inv_sqrt_df[term]
                            for term in matching_seeds)
                        if score > max_score[i]:
                            max_score[i] = score
                        agg = aggregators[i]
                        for term in matching_seeds:
                            agg[term]['sum'] += score
                            agg[term]['count'] += 1
                    record.append(score)

                record.append(';'.join(sorted(matching_seeds)))
                record.append(pub_date)
                batch.append(record)

            if batch:
                out_conn.executemany(insert_sql, batch)
                # Commit per batch so the write-ahead log is checkpointed as we go.
                # Left to a single transaction the WAL grows to the size of the whole
                # table before any of it lands, and an interruption discards the lot.
                out_conn.commit()
                rows_written += len(batch)
                batch = []

            pbar.update(len(rows))

        if batch:
            out_conn.executemany(insert_sql, batch)
            rows_written += len(batch)

        pbar.close()
        master_conn.close()

        # <<< Max-scale the stored ARS to [0, 1] >>>
        # The normalised projection produces values around 1e-6, whereas the old
        # plain sum sat near 1e-2. Downstream consumers are not all scale-free:
        # secondary_analysis blends ars_score directly against a [0, 1] citation
        # score and clips at a dynamic floor, so an unscaled column would be
        # clipped flat and the relevance signal lost. Max-scaling is rank-
        # preserving (the same convention already used for MRS) and puts ARS on
        # the [0, 1] range those blends assume.
        for i, (w_key, _) in enumerate(weightings):
            if max_score[i] > 0:
                col = f"score_{w_key}"
                out_conn.execute(
                    f'UPDATE article_relevance_scores SET "{col}" = "{col}" / ?',
                    (max_score[i],))
        out_conn.commit()
        print(f"  Max-scaled {n_w} ARS column(s) to [0, 1].")

    except Exception as e:
        raise RuntimeError(f"ERROR calculating scores from Master DB: {e}")

    print("\n<<< Indexing Article Relevance Scores >>>")
    if rows_written:
        # The rows are already on disk (streamed during the scan); only the indexes
        # remain. A failure here is raised, not warned about: the previous code
        # swallowed save errors, so a run that wrote nothing still exited cleanly and
        # looked successful until the empty table was discovered downstream.
        try:
            cursor = out_conn.cursor()
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_pmid ON article_relevance_scores (pmid)")
            for i, (w_key, _) in enumerate(weightings, start=1):
                cursor.execute(
                    f"CREATE INDEX IF NOT EXISTS idx_score{i} "
                    f"ON article_relevance_scores ('score_{w_key}')"
                )
            out_conn.commit()
        finally:
            out_conn.close()
        _keep_on_device(relevance_db_path)
        print(f"  [+] Successfully saved {rows_written:,} contributing article scores to the database.")
    else:
        out_conn.close()
        raise RuntimeError(
            "No articles contained the target seed terms, so no relevance scores were "
            "written. Check the date window and that the network's node IDs match the "
            "MeSH terms in the master database."
        )

    # Calculate final Mean Relevancy Scores (MRS) with Max-Scaling
    print("\n<<< Calculating MRS (Mean Relevancy Score & Max-Scaling)... >>>")

    # MRS(i) = mean of S_article(A) over A in P_i  (the articles containing term i).
    # It is the plain mean article-relevance of a term, NOT adjusted for context: the
    # evidence-volume (CF) and information-content (IC) factors of the former CRS were
    # dropped because mean-ARS is already frequency-neutral by construction, so those
    # terms only distorted the node ranking rather than sharpening it.
    final_node_weights = []
    for i in range(len(weightings)):
        raw_node_weights = {}
        for term, data in aggregators[i].items():
            if data['count'] > 0:
                raw_node_weights[term] = data['sum'] / data['count']   # mean S_article over P_i

        # --- Relative max-scaling to [0.0, 1.0] (a display normalisation; rank-preserving) ---
        max_mrs = max(raw_node_weights.values()) if raw_node_weights else 1.0
        final_node_weights.append(
            {term: val / max_mrs for term, val in raw_node_weights.items()} if max_mrs else {}
        )

    print("\n<<< Generating Final Output JSON File >>>")
    try:
        with open(input_nodes_file, 'r') as f:
            network_data = json.load(f)

        for node in network_data.get('elements', {}).get('nodes', []):
            term = node.get('data', {}).get(id_key)
            for i, (_, f_key) in enumerate(weightings):
                node['data'][f_key] = final_node_weights[i].get(term, 0.0)

        os.makedirs(os.path.dirname(output_nodes_file), exist_ok=True)
        with open(output_nodes_file, 'w') as f:
            json.dump(network_data, f, indent=2)

        print(f"  [+] Success! Final relevance JSON created at: {os.path.basename(output_nodes_file)}")

    except Exception as e:
        raise RuntimeError(f"ERROR generating final relevance JSON: {e}")