# -*- coding: utf-8 -*-
"""
data_ops.py - NCBI data collection and project database management (pipeline Step 2).

Retrieves the article cohort for a MeSH search term from NCBI Entrez, expands it
across the requested number of citation generations, then stores and cleans the
result in SQLite for network construction.

Responsibilities include date-chunked PubMed searches with retry/back-off,
generational citation expansion via the iCite API, fetching MeSH annotations and
publication dates into the local master database, and propagating those
annotations back onto the project's PMID tables.
"""

import os
import sqlite3
import time
import requests
from collections import deque
from datetime import date, datetime, timedelta
from urllib.error import HTTPError

from Bio import Entrez, Medline
import pandas as pd
from tqdm import tqdm

from . import runcontrol as _runcontrol


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Helper and Utility Functions
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def get_generation_label(numeric_generation: int) -> str:
    """Converts numeric generation to a label ('P0', 'G1', 'G2'...)."""
    if numeric_generation == 0:
        return 'P0'
    elif numeric_generation > 0:
        return f'G{numeric_generation}'
    return 'Unknown'

def parse_date_robust(date_str: str) -> date:
    """Parses YYYY/MM/DD, YYYY/MM, or YYYY into a date object."""
    if not isinstance(date_str, str) or not date_str.strip():
        raise ValueError("Invalid date input: must be a non-empty string.")
    formats_to_try = ["%Y/%m/%d", "%Y-%m-%d", "%Y/%m", "%Y-%m", "%Y"]
    for fmt in formats_to_try:
        try:
            return datetime.strptime(date_str, fmt).date()
        except ValueError:
            continue
    raise ValueError(f"Invalid date format: '{date_str}'. Use YYYY/MM/DD, YYYY/MM, or YYYY.")

def format_date_for_entrez(date_obj: date) -> str:
    """Formats a date object into YYYY/MM/DD for Entrez."""
    return date_obj.strftime("%Y/%m/%d")


# NCBI's published E-utilities limits: 3 requests a second without an API key,
# 10 with one. The delay was a flat 0.11 s - about 9 a second - whichever it
# was, so an unregistered user, which is the default, ran at three times the
# permitted rate. NCBI blocks the IP for that rather than throttling it.
_DELAY_WITH_KEY = 0.11
_DELAY_NO_KEY = 0.34


def entrez_delay_for(api_key) -> float:
    """The polite interval between requests, given whether a key is set."""
    return _DELAY_WITH_KEY if str(api_key or "").strip() else _DELAY_NO_KEY


def set_entrez_credentials(email, api_key):
    """Hand the e-mail and key to Biopython, or None when they are blank.

    This is the fix for a fatal HTTP 400 that stopped data collection for every
    user who had not registered an API key - which is the default, and so was
    everyone on a fresh install.

    Assigning Entrez.api_key = "" is not the same as leaving it unset.
    Biopython includes any non-None value in the query string, so a blank key
    was sent as a real parameter with no value:

        ...&email=someone%40example.org&api_key=

    NCBI rejects that with "400 Bad Request" - not for the term, not for the
    dates, but for the empty parameter. The failure surfaced during the search,
    which is what made the query itself look like the culprit.

    Returns the delay that suits whichever credentials were actually set.
    """
    email = str(email or "").strip()
    api_key = str(api_key or "").strip()
    Entrez.email = email or None
    Entrez.api_key = api_key or None
    return entrez_delay_for(api_key)


def _announce_rate(api_key):
    """Say which limit is in force, once, where the user will see it."""
    if str(api_key or "").strip():
        print("  NCBI API key found - requesting at up to 10/second.")
    else:
        print("  No NCBI API key set - requesting at up to 3/second, which is "
              "the unregistered limit.\n"
              "  A free key roughly triples retrieval speed: "
              "https://account.ncbi.nlm.nih.gov/settings/")

def _entrez_search_with_retry(db: str, term: str, retmax: int, **kwargs) -> dict:
    """Wraps Entrez.esearch with a retry mechanism for common server errors."""
    max_retries = 5
    for attempt in range(max_retries):
        try:
            handle = Entrez.esearch(db=db, term=term, retmax=retmax, **kwargs)
            results = Entrez.read(handle)
            handle.close()
            return results
        except HTTPError as e:
            if e.code in [429, 500, 502, 503, 504]:
                wait_time = 2 ** attempt
                print(f"  - NCBI server returned error {e.code}. Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                # A 4xx is our fault, not the server's, and it used to arrive as
                # a bare "HTTP Error 400: Bad Request" with no hint of what was
                # sent - which left the query itself to be guessed at. Print the
                # term and the URL before re-raising: between them they name the
                # cause every time.
                print(f"\n  [!] NCBI rejected the request with HTTP {e.code}.")
                print(f"      query sent : {term!r}")
                print(f"      url        : {getattr(e, 'url', '(not reported)')}")
                if not str(term).strip():
                    print("      The query is empty. Set a search term on the "
                          "Search tab before running data collection.")
                raise
        except Exception as e:
            wait_time = 2 ** attempt
            print(f"  - An unexpected error occurred: {e}. Retrying in {wait_time} seconds...")
            time.sleep(wait_time)

    print(f"  - All {max_retries} retries failed for term '{term}'. Skipping.")
    return {"Count": "0", "IdList": []}

def parse_mesh_terms(mesh_terms_string: str) -> list:
    """Parses MeSH terms string into main terms and subheadings."""
    parsed_terms = []
    if not isinstance(mesh_terms_string, str) or not mesh_terms_string.strip():
        return parsed_terms

    mesh_terms = mesh_terms_string.split(';')
    for term_part in mesh_terms:
        term_part = term_part.strip()
        if not term_part: continue

        parts = term_part.split('/')
        main_term_component = parts[0].strip()
        if not main_term_component: continue

        is_main_major = main_term_component.startswith('*')
        main_term = main_term_component.lstrip('*').strip()

        if main_term:
            parsed_terms.append((main_term, 'mesh_term', is_main_major))

        if len(parts) > 1:
            for subheading_component in parts[1:]:
                subheading_component = subheading_component.strip()
                if not subheading_component: continue

                is_sub_major = subheading_component.startswith('*')
                subheading = subheading_component.lstrip('*').strip()

                if subheading:
                    parsed_terms.append((subheading, 'subheading', is_sub_major))
    return parsed_terms

def _fetch_ids_for_range(search_term: str, start_date_obj: date, end_date_obj: date, retmax_chunk: int, entrez_delay: float) -> set:
    """Performs a single esearch and returns PMIDs."""
    pmids = set()
    full_query = f"({search_term}) AND ({format_date_for_entrez(start_date_obj)}[EDAT] : {format_date_for_entrez(end_date_obj)}[EDAT])"

    search_results = _entrez_search_with_retry(db="pubmed", term=full_query, retmax=retmax_chunk, sort="relevance")
    id_list = search_results.get("IdList", [])
    pmids.update(id_list)
    time.sleep(entrez_delay)

    return pmids

def get_pmids_date_chunking(search_term: str, start_date_obj: date, end_date_obj: date, retmax_limit: int = 9999, entrez_delay: float = _DELAY_NO_KEY) -> list:
    """Fetches PMIDs by splitting large date ranges. Output is optimized for Colab rendering."""
    all_pmids = set()
    date_ranges_to_process = deque([(start_date_obj, end_date_obj)])
    chunks_processed = 0
    start_time = time.time()

    print(f"Processing PMIDs in date ranges for: '{search_term}'...")

    while date_ranges_to_process:
        current_start, current_end = date_ranges_to_process.popleft()

        date_query = f"({format_date_for_entrez(current_start)}[EDAT] : {format_date_for_entrez(current_end)}[EDAT])"
        full_query = f"({search_term}) AND {date_query}"

        search_results = _entrez_search_with_retry(db="pubmed", term=full_query, retmax="0")
        count = int(search_results["Count"])
        time.sleep(entrez_delay)

        if count <= retmax_limit:
            retrieved_ids = _fetch_ids_for_range(search_term, current_start, current_end, retmax_limit, entrez_delay)
            all_pmids.update(retrieved_ids)
            chunks_processed += 1
        else:
            if current_start == current_end:
                date_str = current_start.strftime('%Y-%m-%d')
                print(f"  Warning: Count for single day {date_str} ({count}) exceeds limit. Fetching first {retmax_limit}.")
                retrieved_ids = _fetch_ids_for_range(search_term, current_start, current_end, retmax_limit, entrez_delay)
                all_pmids.update(retrieved_ids)
                chunks_processed += 1
                continue

            mid_point_ord = (current_start.toordinal() + current_end.toordinal()) // 2
            mid_date = date.fromordinal(mid_point_ord)
            date_ranges_to_process.appendleft((mid_date + timedelta(days=1), current_end))
            date_ranges_to_process.appendleft((current_start, mid_date))

    total_time = time.time() - start_time
    time_per_chunk = total_time / chunks_processed if chunks_processed > 0 else 0
    print(f" -> Completed. Found {len(all_pmids):,} unique PMIDs across {chunks_processed} chunks in {total_time:.1f}s ({time_per_chunk:.2f}s/chunk).")

    return sorted(list(all_pmids))

def fetch_links_in_batches(parent_pmid_batch: list) -> dict:
    """Fetches citation data for a batch of PMIDs using the iCite API."""
    citation_results = {pmid: {'cited_by': "", 'cites': ""} for pmid in parent_pmid_batch}
    if not parent_pmid_batch:
        return citation_results

    api_url = "https://icite.od.nih.gov/api/pubs"
    params = {'pmids': ",".join(map(str, parent_pmid_batch))}

    try:
        response = requests.get(api_url, params=params)
        response.raise_for_status()
        data = response.json().get('data', [])

        for article_data in data:
            pmid = article_data.get('pmid')
            if pmid in citation_results:
                cited_by_data = article_data.get('cited_by') or []
                references_data = article_data.get('references') or []
                citation_results[pmid]['cited_by'] = ";".join([str(p) for p in cited_by_data])
                citation_results[pmid]['cites'] = ";".join([str(p) for p in references_data])

        return citation_results
    except requests.exceptions.RequestException as e:
        print(f"  ERROR: Failed to fetch data from iCite API: {e}")
        return citation_results

def _parse_medline_date(record: dict) -> str:
    """Safely extracts and formats the publication date from an Entrez Medline dictionary."""
    # 1. Try EDAT (Entrez Date) which is usually strictly formatted YYYY/MM/DD HH:MM
    edat = record.get("EDAT", "")
    if edat and len(edat) >= 10:
        date_part = edat[:10].replace("/", "-")
        if len(date_part.split("-")) == 3:
            return date_part

    # 2. Try DP (Date of Publication) which is often human-readable "2025 Jan 15"
    dp = record.get("DP", "")
    if dp:
        parts = dp.split()
        y = parts[0] if parts else "1900"
        m = '01'
        d = '01'

        if len(parts) > 1:
            m_str = parts[1][:3].lower()
            m_map = {
                'jan':'01', 'feb':'02', 'mar':'03', 'apr':'04',
                'may':'05', 'jun':'06', 'jul':'07', 'aug':'08',
                'sep':'09', 'oct':'10', 'nov':'11', 'dec':'12'
            }
            m = m_map.get(m_str, '01')

        if len(parts) > 2 and parts[2].isdigit():
            d = parts[2].zfill(2)

        return f"{y}-{m}-{d}"

    return "1900-01-01"


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Pipeline Operations
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def run_initial_data_collection(search_term_param: str, start_date_str: str, end_date_str: str,
                                generations_n_param: int, db_path: str, entrez_email: str,
                                entrez_api_key: str, batch_size: int = 100, entrez_delay: float = None):
    """Performs initial PMID search and generational expansion, storing results in SQLite."""
    resolved_delay = set_entrez_credentials(entrez_email, entrez_api_key)
    if entrez_delay is None:
        entrez_delay = resolved_delay
        _announce_rate(entrez_api_key)

    script_start_time = time.time()
    all_pmids_gen = {}
    conn = None
    generations_to_process = max(1, int(generations_n_param))

    try:
        # Checked here, before a single request goes out. An empty term is not
        # an empty result: "() AND (1950/01/01[EDAT] : ...)" is a legal PubMed
        # query that matches the entire database - forty million records, which
        # the date-chunking loop will then set about retrieving nine thousand at
        # a time. Failing in one line beats discovering that overnight.
        if not str(search_term_param or "").strip():
            raise ValueError(
                "No search term is set, so there is nothing to search for. "
                "Set one on the Search tab - for example "
                "'Dermatitis, Allergic Contact [Mesh]' - and run this again.")

        # A date word that never got resolved reaches PubMed as itself, and
        # PubMed does not reject it: it returns zero results and a warning
        # buried in the response, so the run ends with "No initial PMIDs found"
        # and no indication that the date was the reason.
        for label, value in (("Start date", start_date_str), ("End date", end_date_str)):
            if isinstance(value, str) and value.strip().lower() in ("today", "now"):
                raise ValueError(
                    f"{label} arrived as '{value.strip()}' rather than a date. "
                    "PubMed treats that as a search word, not a date, and "
                    "returns nothing. Enter a date as YYYY/MM/DD, or leave the "
                    "field empty to search every date PubMed has.")

        start_date_obj = parse_date_robust(start_date_str or "1900/01/01")
        end_date_obj = parse_date_robust(end_date_str or date.today().strftime("%Y/%m/%d"))
        if start_date_obj > end_date_obj:
            raise ValueError("Start date cannot be after end date.")

        print(f"  Search window: {format_date_for_entrez(start_date_obj)} "
              f"to {format_date_for_entrez(end_date_obj)} (EDAT - the date the "
              f"record entered PubMed)")

        initial_pmids_str_list = get_pmids_date_chunking(search_term_param, start_date_obj, end_date_obj, entrez_delay=entrez_delay)
        if not initial_pmids_str_list:
            print("No initial PMIDs found. Exiting.")
            return

        all_pmids_gen = {int(pmid): 0 for pmid in initial_pmids_str_list}
        print(f"Found {len(all_pmids_gen)} unique initial PMIDs (P0).")

        print(f"\n<<< Setting up Database: {db_path} >>>")
        os.makedirs(os.path.dirname(db_path), exist_ok=True)
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("PRAGMA journal_mode=WAL;")
        cursor.execute("PRAGMA synchronous=NORMAL;")
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS pmids_table (
                pmid INTEGER PRIMARY KEY,
                generation TEXT NOT NULL,
                cited_by TEXT,
                cites TEXT
            )
        ''')
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_generation ON pmids_table(generation)")
        conn.commit()

        p0_data = [(pmid, 'P0') for pmid in all_pmids_gen.keys()]
        cursor.executemany("INSERT OR IGNORE INTO pmids_table (pmid, generation) VALUES (?, ?)", p0_data)
        conn.commit()

        for current_gen_num in range(generations_to_process - 1):
            parent_gen_label = get_generation_label(current_gen_num)
            child_gen_label = get_generation_label(current_gen_num + 1)
            print(f"\n<<< Processing {parent_gen_label} PMIDs to find {child_gen_label} children >>>")

            cursor.execute("SELECT pmid FROM pmids_table WHERE generation = ?", (parent_gen_label,))
            parent_pmids_in_gen = [row[0] for row in cursor.fetchall()]

            if not parent_pmids_in_gen:
                print(f"No PMIDs for {parent_gen_label} in DB. Stopping expansion.")
                break

            for i in tqdm(range(0, len(parent_pmids_in_gen), batch_size), desc=f"Finding children of {parent_gen_label}"):
                # Between batches: the transaction for the previous one
                # is committed and the next has not begun.
                _runcontrol.checkpoint(f"retrieving {child_gen_label}")
                parent_batch = parent_pmids_in_gen[i:i + batch_size]
                citation_results = fetch_links_in_batches(parent_batch)
                time.sleep(entrez_delay)

                update_link_data, insert_new_pmid_data = [], []

                for pmid in parent_batch:
                    links = citation_results.get(pmid, {'cited_by': "", 'cites': ""})
                    cited_by_str, cites_str = links['cited_by'], links['cites']
                    update_link_data.append((cited_by_str, cites_str, pmid))

                    child_candidates = set((cited_by_str + ";" + cites_str).split(';'))
                    for child_pmid_str in filter(None, child_candidates):
                        try:
                            child_pmid = int(child_pmid_str)
                            if child_pmid not in all_pmids_gen:
                                all_pmids_gen[child_pmid] = current_gen_num + 1
                                insert_new_pmid_data.append((child_pmid, child_gen_label))
                        except ValueError:
                            continue

                if update_link_data:
                    cursor.executemany("UPDATE pmids_table SET cited_by = ?, cites = ? WHERE pmid = ?", update_link_data)
                if insert_new_pmid_data:
                    cursor.executemany("INSERT OR IGNORE INTO pmids_table (pmid, generation) VALUES (?, ?)", insert_new_pmid_data)
                conn.commit()

        if generations_to_process > 0:
            final_gen_label = get_generation_label(generations_to_process - 1)
            print(f"\n<<< Fetching links FOR final generation ({final_gen_label}) >>>")
            cursor.execute("SELECT pmid FROM pmids_table WHERE generation = ? AND (cited_by IS NULL OR cites IS NULL)", (final_gen_label,))
            final_gen_pmids = [row[0] for row in cursor.fetchall()]

            for i in tqdm(range(0, len(final_gen_pmids), batch_size), desc=f"Fetching Links for {final_gen_label}"):
                _runcontrol.checkpoint(f"fetching links for {final_gen_label}")
                final_batch = final_gen_pmids[i:i + batch_size]
                citation_results = fetch_links_in_batches(final_batch)
                time.sleep(entrez_delay)

                update_link_data = []
                for pmid in final_batch:
                    links = citation_results.get(pmid, {'cited_by': "", 'cites': ""})
                    update_link_data.append((links['cited_by'], links['cites'], pmid))

                if update_link_data:
                    cursor.executemany("UPDATE pmids_table SET cited_by = ?, cites = ? WHERE pmid = ?", update_link_data)
                conn.commit()

    except ValueError:
        # A settings problem the user can act on - an empty search term, a date
        # that is not a date. "An unexpected error occurred" is exactly the
        # wrong frame for those: it reads as a program fault and buries the one
        # sentence that says what to change. Let it through as it was written.
        raise
    except Exception as e:
        raise RuntimeError(f"An unexpected error occurred during data collection: {e}")
    finally:
        if conn:
            conn.close()
            print("\nDatabase connection closed.")

    script_end_time = time.time()
    print(f"\nTotal script execution time: {script_end_time - script_start_time:.2f} seconds")
    print(f"Total unique PMIDs found across {generations_to_process} generation(s): {len(all_pmids_gen)}")

def clean_database(original_db: str, cleaned_db: str, start_gen_to_remove_label: str):
    """Creates a cleaned copy of the database, removing specified generations and higher."""
    print(f"\n<<< Cleaning Database >>>")
    if not os.path.exists(original_db):
        print(f"Original database '{original_db}' not found. Skipping cleaning process.")
        return

    import shutil
    try:
        os.makedirs(os.path.dirname(cleaned_db), exist_ok=True)
        shutil.copyfile(original_db, cleaned_db)
    except IOError as e:
        raise RuntimeError(f"Error copying database file: {e}")

    conn = sqlite3.connect(cleaned_db)
    cursor = conn.cursor()

    try:
        start_num_to_remove = int(start_gen_to_remove_label.lstrip('G'))
        cursor.execute("SELECT DISTINCT generation FROM pmids_table")
        all_gens = [row[0] for row in cursor.fetchall()]

        gens_to_delete = []
        for gen_label in all_gens:
            if gen_label.startswith('G'):
                try:
                    if int(gen_label.lstrip('G')) >= start_num_to_remove:
                        gens_to_delete.append(gen_label)
                except ValueError:
                    continue

        if gens_to_delete:
            placeholders = ', '.join('?' for _ in gens_to_delete)
            sql_query = f"DELETE FROM pmids_table WHERE generation IN ({placeholders})"
            cursor.execute(sql_query, gens_to_delete)
            print(f"Successfully removed {cursor.rowcount} records from generations: {gens_to_delete}")
        else:
            print("No generations met the criteria for removal.")

    except (ValueError, IndexError):
        raise ValueError(f"Invalid format for `start_gen_to_remove`: '{start_gen_to_remove_label}'. Should be e.g., 'G2'.")
    except sqlite3.Error as e:
        raise RuntimeError(f"An error occurred during database cleaning: {e}")
    finally:
        conn.commit()
        conn.close()

def populate_master_mesh_database(source_pmids_input, master_db_path: str, entrez_email: str,
                                  entrez_api_key: str, fetch_batch_size: int = 9999,
                                  entrez_delay: float = None, failed_log: str = None, empty_log: str = None):
    """Checks source PMIDs against master DB, fetches missing MeSH terms, and updates master DB."""
    resolved_delay = set_entrez_credentials(entrez_email, entrez_api_key)
    if entrez_delay is None:
        entrez_delay = resolved_delay
        _announce_rate(entrez_api_key)

    if isinstance(source_pmids_input, (str, os.PathLike)) and os.path.exists(source_pmids_input):
        try:
            conn = sqlite3.connect(source_pmids_input)
            cursor = conn.cursor()
            cursor.execute("SELECT pmid FROM pmids_table")
            source_pmids = {row[0] for row in cursor.fetchall()}
            conn.close()
        except sqlite3.Error as e:
            raise RuntimeError(f"Database Error reading PMIDs: {e}")
    elif isinstance(source_pmids_input, (list, set)):
        source_pmids = set(source_pmids_input)
    else:
        raise ValueError(f"Invalid input. Got {type(source_pmids_input)}. Expecting DB path or list of PMIDs.")

    os.makedirs(os.path.dirname(master_db_path), exist_ok=True)
    master_conn = sqlite3.connect(master_db_path)
    m_cursor = master_conn.cursor()

    # 4-column schema. pmid is the PRIMARY KEY so INSERT OR IGNORE actually
    # de-duplicates (no unique constraint => OR IGNORE is a silent no-op) and so
    # pmid lookups/joins are indexed even when the baseline ETL never ran.
    m_cursor.execute("CREATE TABLE IF NOT EXISTS master_mesh_annotations (pmid INTEGER PRIMARY KEY, pub_date TEXT, mesh_terms TEXT, source_file TEXT)")

    m_cursor.execute("SELECT pmid FROM master_mesh_annotations")
    master_pmids = {row[0] for row in m_cursor.fetchall()}

    new_pmids_to_fetch = list(source_pmids - master_pmids)
    if not new_pmids_to_fetch:
        print("No new PMIDs from this analysis to add to the master database.")
        master_conn.close()
        return

    empty_pmids, failed_pmids = [], []

    for i in tqdm(range(0, len(new_pmids_to_fetch), fetch_batch_size), desc="Populating Master MeSH DB"):
        _runcontrol.checkpoint("fetching MeSH annotations")
        batch = new_pmids_to_fetch[i : i + fetch_batch_size]
        if not batch:
            continue

        try:
            epost_handle = Entrez.epost(db="pubmed", id=",".join(map(str, batch)))
            results = Entrez.read(epost_handle)
            epost_handle.close()
            webenv, query_key = results["WebEnv"], results["QueryKey"]
            efetch_handle = Entrez.efetch(db="pubmed", rettype="medline", retmode="text", webenv=webenv, query_key=query_key)
            records = Medline.parse(efetch_handle)

            all_batch_updates, processed_in_batch = [], set()

            for record in records:
                pmid_str = record.get("PMID")
                if pmid_str:
                    pmid_int = int(pmid_str)
                    processed_in_batch.add(pmid_int)
                    mesh_terms_list = record.get("MH", [])
                    pub_date = _parse_medline_date(record)

                    # UPDATED: Insert matches new 4-column format
                    all_batch_updates.append((pmid_int, pub_date, ";".join(mesh_terms_list), "Entrez_API"))

                    if not mesh_terms_list:
                        empty_pmids.append(pmid_int)

            missing_pmids = set(map(int, batch)) - processed_in_batch
            for pmid in missing_pmids:
                # UPDATED: Insert matches new 4-column format
                all_batch_updates.append((pmid, "1900-01-01", "", "Entrez_API"))
                failed_pmids.append(pmid)

            if all_batch_updates:
                m_cursor.executemany("INSERT OR IGNORE INTO master_mesh_annotations (pmid, pub_date, mesh_terms, source_file) VALUES (?, ?, ?, ?)", all_batch_updates)
                master_conn.commit()

            time.sleep(entrez_delay)

        except Exception as e:
            print(f"Error occurred during master DB population for batch starting with {batch[0]}: {e}")
            failed_pmids.extend(map(int, batch))

    master_conn.close()

    if failed_log and failed_pmids:
        try:
            os.makedirs(os.path.dirname(failed_log), exist_ok=True)
            with open(failed_log, 'w') as f:
                f.write("PMID\n" + "\n".join(map(str, failed_pmids)) + "\n")
        except Exception as e:
            print(f"Could not write failed PMIDs log: {e}")

    if empty_log and empty_pmids:
        try:
            os.makedirs(os.path.dirname(empty_log), exist_ok=True)
            with open(empty_log, 'w') as f:
                f.write("PMID\n" + "\n".join(map(str, empty_pmids)) + "\n")
        except Exception as e:
            print(f"Could not write empty PMIDs log: {e}")

def update_db_with_mesh_batch(local_db_path, master_db_path):
    """Updates the local project DB with MeSH terms from the master annotation DB."""
    # Convert Path objects to strings to prevent SQLite parameter binding errors
    local_db_path = str(local_db_path)
    master_db_path = str(master_db_path)

    if not os.path.exists(master_db_path):
        raise FileNotFoundError(f"Master database not found at {master_db_path}.")
    if not os.path.exists(local_db_path):
        raise FileNotFoundError(f"Local database not found at {local_db_path}.")

    import sqlite3
    conn = sqlite3.connect(local_db_path)
    cursor = conn.cursor()
    try:
        try:
            cursor.execute("ALTER TABLE pmids_table ADD COLUMN mesh_terms TEXT")
        except sqlite3.OperationalError:
            pass

        cursor.execute("ATTACH DATABASE ? AS master_db", (master_db_path,))
        update_query = """
        UPDATE pmids_table
        SET mesh_terms = (
            SELECT m.mesh_terms
            FROM master_db.master_mesh_annotations AS m
            WHERE m.pmid = pmids_table.pmid
        )
        WHERE
            EXISTS (
                SELECT 1
                FROM master_db.master_mesh_annotations AS m
                WHERE m.pmid = pmids_table.pmid
            )
            AND (pmids_table.mesh_terms IS NULL OR pmids_table.mesh_terms = '');
        """
        cursor.execute(update_query)
        conn.commit()
        print(f"Update complete. {cursor.rowcount:,} records were updated.")
        cursor.execute("DETACH DATABASE master_db")
        conn.commit()
    except sqlite3.Error as e:
        raise RuntimeError(f"An SQL error occurred during database annotation update: {e}")
    finally:
        conn.close()