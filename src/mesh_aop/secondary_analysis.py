# -*- coding: utf-8 -*-
"""
secondary_analysis.py - targeted literature queries and exports (pipeline Step 3.5).

Runs ad-hoc analyses against the finished network to surface the specific
publications behind a single node, an edge (relationship), or the network as a
whole, for manual review.

It ranks candidate articles with a dual-engine Article Impact Score (a
compensatory Linear blend, or a penalizing F1 harmonic mean, of relevance and
citation impact), hydrates the top results with live NCBI metadata (title, year,
publication type, with optional review filtering), and exports both the ranked
article tables and the full network as Excel/CSV files.
"""

import os
import json
import sqlite3
import time
import numpy as np
import pandas as pd
from Bio import Entrez, Medline

from .baseline_manager import _keep_on_device


def _open_readonly_resilient(db_path):
    """Open a database read-only, resilient to OneDrive dehydration.

    Pins the file and opens read-only; if the read-only view exposes no tables -
    a dehydrated OneDrive placeholder reads back empty - it reopens with a normal
    connection, which forces OneDrive to hydrate the real file. No writes are
    issued, and an already-hydrated file keeps the read-only path unchanged.
    """
    _keep_on_device(db_path)
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True, timeout=120.0)
    if conn.execute("SELECT name FROM sqlite_master WHERE type='table' LIMIT 1").fetchone():
        return conn
    conn.close()
    return sqlite3.connect(str(db_path), timeout=120.0)

def _calculate_impact_scores(raw_rows: list, sort_metric: str, linear_weight_ars: float, oversample_limit: int) -> list:
    """
    Pandas engine to apply dynamic scaling and dual-metric Article Impact Scoring.
    Replaces arbitrary floors with a dynamic baseline relative to the network's max citation count.
    """
    if not raw_rows:
        return []

    df = pd.DataFrame(raw_rows, columns=['pmid', 'ars_score', 'seeds', 'cited_by'])

    # 1. Extract raw citation counts (0, 1, 2, etc.)
    cited_by_str = df['cited_by'].fillna('').astype(str).str.strip()
    is_empty = (cited_by_str == '') | (cited_by_str.str.lower() == 'nan') | (cited_by_str.str.lower() == 'none')
    raw_counts = np.where(is_empty, 0, cited_by_str.str.count(';') + 1)

    # Bind the raw counts to the DataFrame immediately before any sorting can occur
    df['Incoming_Citations_Smoothed'] = raw_counts

    # 2. Log-Normalize (using standard +1 so 0 citations = log10(1) = 0.0 initially)
    df['log_cit'] = np.log10(raw_counts + 1)

    # 3. Dynamic Scaling based on the Network's Maximum
    max_log = df['log_cit'].max()

    if max_log > 0:
        one_citation_scale = np.log10(2) / max_log
        dynamic_floor = one_citation_scale / 2.0
        normalized = df['log_cit'] / max_log
        df['Normalized_Citation'] = np.where(raw_counts == 0, dynamic_floor, normalized)
    else:
        dynamic_floor = 1.0
        df['Normalized_Citation'] = 1.0

    df['ars_score'] = df['ars_score'].clip(lower=dynamic_floor)

    # 4. Dual Calculation Engine
    df['Linear_AIS'] = (df['ars_score'] * linear_weight_ars) + (df['Normalized_Citation'] * (1.0 - linear_weight_ars))
    df['F1_AIS'] = (2 * df['ars_score'] * df['Normalized_Citation']) / (df['ars_score'] + df['Normalized_Citation'])

    # 5. Sort by user's chosen metric and truncate
    sort_col = 'Linear_AIS' if sort_metric.lower() == 'linear' else 'F1_AIS'
    df = df.sort_values(by=sort_col, ascending=False)

    return df.head(oversample_limit).to_dict('records')


def _fetch_metadata_and_filter(pmid_score_dicts: list, exclude_reviews: bool, limit: int,
                               entrez_email: str, entrez_api_key: str, sort_metric: str) -> pd.DataFrame:
    """
    Hydrates PMIDs with title, year, and publication type via the NCBI API.
    Implements NCBI Official Pagination (retstart/retmax) to prevent payload timeouts.
    """
    if not pmid_score_dicts:
        return pd.DataFrame()

    Entrez.email = entrez_email
    Entrez.api_key = entrez_api_key

    pmid_str_list = [str(d['pmid']) for d in pmid_score_dicts]
    hydrated_results = []
    score_map = {str(d['pmid']): d for d in pmid_score_dicts}

    try:
        epost_handle = Entrez.epost(db="pubmed", id=",".join(pmid_str_list))
        post_results = Entrez.read(epost_handle)
        epost_handle.close()

        webenv = post_results["WebEnv"]
        query_key = post_results["QueryKey"]

        batch_size = 200
        total_candidates = len(pmid_str_list)

        for start in range(0, total_candidates, batch_size):
            if len(hydrated_results) >= limit:
                break

            fetch_handle = Entrez.efetch(db="pubmed", rettype="medline", retmode="text",
                                         retstart=start, retmax=batch_size,
                                         webenv=webenv, query_key=query_key)
            records = list(Medline.parse(fetch_handle))
            fetch_handle.close()

            for record in records:
                if len(hydrated_results) >= limit:
                    break

                pmid = record.get("PMID", "")
                if not pmid: continue

                if exclude_reviews:
                    pub_types = record.get("PT", [])
                    if any("Review" in pt for pt in pub_types):
                        continue

                title = record.get("TI", "No Title Available")
                dp = record.get("DP", "Unknown Year")
                year = dp.split()[0] if dp else "Unknown"

                base_data = score_map.get(pmid, {})

                hydrated_results.append({
                    "PMID": pmid,
                    "Year": year,
                    "Title": title,
                    "Raw_ARS": base_data.get('ars_score', 0.0),
                    "Incoming_Citations": base_data.get('Incoming_Citations_Smoothed', 0),
                    "Normalized_Citation": base_data.get('Normalized_Citation', 0.0),
                    "Linear_AIS": base_data.get('Linear_AIS', 0.0),
                    "F1_AIS": base_data.get('F1_AIS', 0.0),
                    "Contributing_Seeds": base_data.get('seeds', '')
                })

            time.sleep(0.35)

        sort_col = 'Linear_AIS' if sort_metric.lower() == 'linear' else 'F1_AIS'
        hydrated_results.sort(key=lambda x: x[sort_col], reverse=True)

        return pd.DataFrame(hydrated_results)

    except Exception as e:
        print(f"  [!] Error fetching metadata from NCBI: {e}")
        return pd.DataFrame()


def analyze_node_relevancy(node_name: str, db_path: str, cleaned_db_path: str, results_dir: str, file_prefix: str,
                           limit: int, exclude_reviews: bool, entrez_email: str, entrez_api_key: str,
                           sort_metric: str = 'Linear', linear_weight_ars: float = 0.5):
    """Pulls and scores all articles associated with a specific node."""
    if not os.path.exists(db_path) or not os.path.exists(cleaned_db_path):
        raise FileNotFoundError("Missing required databases.")

    print(f"\n<<< Searching for Node: '{node_name}' >>>")
    oversample_limit = limit * 5

    try:
        # Step 1: Query Relevance DB
        conn_rel = _open_readonly_resilient(db_path)
        cursor_rel = conn_rel.cursor()
        cursor_rel.execute("SELECT pmid, score_betweenness_centrality, contributing_seeds FROM article_relevance_scores WHERE contributing_seeds LIKE ?", (f'%{node_name}%',))
        rel_rows = cursor_rel.fetchall()
        conn_rel.close()

        if not rel_rows:
            print(f"  -> No articles found for node '{node_name}'.")
            return

        # Step 2: Clean PMIDs and map data
        clean_pmids = [str(r[0]).split('.')[0] for r in rel_rows]
        rel_map = {clean_pmids[i]: (rel_rows[i][1], rel_rows[i][2]) for i in range(len(clean_pmids))}

        # Step 3: Batch lookup in Master DB
        citation_map = {}
        conn_clean = _open_readonly_resilient(cleaned_db_path)
        cursor_clean = conn_clean.cursor()

        chunk_size = 900
        for i in range(0, len(clean_pmids), chunk_size):
            chunk = clean_pmids[i:i + chunk_size]
            placeholders = ','.join('?' for _ in chunk)
            cursor_clean.execute(f"SELECT pmid, cited_by FROM pmids_table WHERE pmid IN ({placeholders})", chunk)
            for row in cursor_clean.fetchall():
                citation_map[str(row[0])] = row[1]
        conn_clean.close()

        # Step 4: Reconstruct array
        raw_rows = [(pmid, rel_map[pmid][0], rel_map[pmid][1], citation_map.get(pmid, None)) for pmid in clean_pmids]

        pmid_dicts = _calculate_impact_scores(raw_rows, sort_metric, linear_weight_ars, oversample_limit)
        print(f"  -> Extracted & Scored top {len(pmid_dicts)} candidates. Hydrating metadata...")

        df = _fetch_metadata_and_filter(pmid_dicts, exclude_reviews, limit, entrez_email, entrez_api_key, sort_metric)

        if not df.empty:
            print(f"  Top {len(df)} results (Ranked by {sort_metric} Article Impact Score (AIS)):")
            sort_col = 'Linear_AIS' if sort_metric.lower() == 'linear' else 'F1_AIS'
            print(df[['PMID', 'Year', 'Raw_ARS', 'Incoming_Citations', sort_col]].head(5).to_string(index=False))

            os.makedirs(results_dir, exist_ok=True)
            safe_name = node_name.replace(" ", "_").replace("/", "-")
            out_file = os.path.join(results_dir, f"{file_prefix}_relevance_{safe_name}.csv")
            df.to_csv(out_file, index=False)
            print(f"  Saved full results to: {os.path.basename(out_file)}")
        else:
            print("  -> No valid articles remained after applying filters.")

    except Exception as e:
        raise RuntimeError(f"Error querying database for node relevancy: {e}")


def analyze_edge_relevancy(node1: str, node2: str, db_path: str, cleaned_db_path: str, results_dir: str, file_prefix: str,
                           limit: int, exclude_reviews: bool, entrez_email: str, entrez_api_key: str,
                           sort_metric: str = 'Linear', linear_weight_ars: float = 0.5):
    """Pulls and scores all articles associated with a specific EDGE (relationship)."""
    if not os.path.exists(db_path) or not os.path.exists(cleaned_db_path):
        raise FileNotFoundError("Missing required databases.")

    print(f"\n<<< Searching for Edge: '{node1}' <--> '{node2}' >>>")
    oversample_limit = limit * 5

    try:
        conn_rel = _open_readonly_resilient(db_path)
        cursor_rel = conn_rel.cursor()
        cursor_rel.execute("SELECT pmid, score_betweenness_centrality, contributing_seeds FROM article_relevance_scores WHERE contributing_seeds LIKE ? AND contributing_seeds LIKE ?", (f'%{node1}%', f'%{node2}%'))
        rel_rows = cursor_rel.fetchall()
        conn_rel.close()

        if not rel_rows:
            print(f"  -> No articles found linking '{node1}' and '{node2}'.")
            return

        clean_pmids = [str(r[0]).split('.')[0] for r in rel_rows]
        rel_map = {clean_pmids[i]: (rel_rows[i][1], rel_rows[i][2]) for i in range(len(clean_pmids))}

        citation_map = {}
        conn_clean = _open_readonly_resilient(cleaned_db_path)
        cursor_clean = conn_clean.cursor()

        chunk_size = 900
        for i in range(0, len(clean_pmids), chunk_size):
            chunk = clean_pmids[i:i + chunk_size]
            placeholders = ','.join('?' for _ in chunk)
            cursor_clean.execute(f"SELECT pmid, cited_by FROM pmids_table WHERE pmid IN ({placeholders})", chunk)
            for row in cursor_clean.fetchall():
                citation_map[str(row[0])] = row[1]
        conn_clean.close()

        raw_rows = [(pmid, rel_map[pmid][0], rel_map[pmid][1], citation_map.get(pmid, None)) for pmid in clean_pmids]

        pmid_dicts = _calculate_impact_scores(raw_rows, sort_metric, linear_weight_ars, oversample_limit)
        print(f"  -> Extracted & Scored top {len(pmid_dicts)} candidates. Hydrating metadata...")

        df = _fetch_metadata_and_filter(pmid_dicts, exclude_reviews, limit, entrez_email, entrez_api_key, sort_metric)

        if not df.empty:
            print(f"  Top {len(df)} results (Ranked by {sort_metric} Article Impact Score (AIS)):")
            sort_col = 'Linear_AIS' if sort_metric.lower() == 'linear' else 'F1_AIS'
            print(df[['PMID', 'Year', 'Raw_ARS', 'Incoming_Citations', sort_col]].head(5).to_string(index=False))

            os.makedirs(results_dir, exist_ok=True)
            safe_n1 = node1.replace(" ", "_").replace(",", "")
            safe_n2 = node2.replace(" ", "_").replace(",", "")
            out_file = os.path.join(results_dir, f"{file_prefix}_edge_relevance_{safe_n1}_{safe_n2}.csv")
            df.to_csv(out_file, index=False)
            print(f"  Saved full results to: {os.path.basename(out_file)}")
        else:
            print("  -> No valid articles remained after applying filters.")

    except Exception as e:
        raise RuntimeError(f"Error querying database for edge relevancy: {e}")


def get_top_network_articles(db_path: str, cleaned_db_path: str, results_dir: str, file_prefix: str, limit: int,
                             exclude_reviews: bool, entrez_email: str, entrez_api_key: str,
                             sort_metric: str = 'Linear', linear_weight_ars: float = 0.5):
    """Retrieves and scores articles across the entire refined network."""
    if not os.path.exists(db_path) or not os.path.exists(cleaned_db_path):
        raise FileNotFoundError("Missing required databases.")

    print(f"\n<<< Extracting Top {limit} Network-Wide Articles >>>")
    oversample_limit = limit * 5

    try:
        conn_rel = _open_readonly_resilient(db_path)
        cursor_rel = conn_rel.cursor()
        cursor_rel.execute("SELECT pmid, score_betweenness_centrality, contributing_seeds FROM article_relevance_scores ORDER BY score_betweenness_centrality DESC LIMIT 50000")
        rel_rows = cursor_rel.fetchall()
        conn_rel.close()

        if not rel_rows:
            print("  -> No articles found in the relevance database.")
            return

        clean_pmids = [str(r[0]).split('.')[0] for r in rel_rows]
        rel_map = {clean_pmids[i]: (rel_rows[i][1], rel_rows[i][2]) for i in range(len(clean_pmids))}

        citation_map = {}
        conn_clean = _open_readonly_resilient(cleaned_db_path)
        cursor_clean = conn_clean.cursor()

        chunk_size = 900
        for i in range(0, len(clean_pmids), chunk_size):
            chunk = clean_pmids[i:i + chunk_size]
            placeholders = ','.join('?' for _ in chunk)
            cursor_clean.execute(f"SELECT pmid, cited_by FROM pmids_table WHERE pmid IN ({placeholders})", chunk)
            for row in cursor_clean.fetchall():
                citation_map[str(row[0])] = row[1]
        conn_clean.close()

        raw_rows = [(pmid, rel_map[pmid][0], rel_map[pmid][1], citation_map.get(pmid, None)) for pmid in clean_pmids]

        pmid_dicts = _calculate_impact_scores(raw_rows, sort_metric, linear_weight_ars, oversample_limit)

        df = _fetch_metadata_and_filter(pmid_dicts, exclude_reviews, limit, entrez_email, entrez_api_key, sort_metric)

        if not df.empty:
            print(f"  Top {len(df)} results (Ranked by {sort_metric} Article Impact Score (AIS)):")
            sort_col = 'Linear_AIS' if sort_metric.lower() == 'linear' else 'F1_AIS'
            print(df[['PMID', 'Year', 'Raw_ARS', 'Incoming_Citations', sort_col]].head(5).to_string(index=False))

            os.makedirs(results_dir, exist_ok=True)
            out_file = os.path.join(results_dir, f"{file_prefix}_Top_Network_Articles.csv")
            df.to_csv(out_file, index=False)
            print(f"  Saved full results to: {os.path.basename(out_file)}")
        else:
            print("  -> No valid articles remained after applying filters.")

    except Exception as e:
        raise RuntimeError(f"Error retrieving top network articles: {e}")


def convert_network_json_to_excel(input_json_path: str, output_excel_path: str):
    """Converts a Cytoscape JSON network file into a structured Excel file."""
    print(f"\n<<< Converting JSON to Excel >>>")

    if not os.path.exists(input_json_path):
        raise FileNotFoundError(f"Input file does not exist: {input_json_path}")

    try:
        with open(input_json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)

        nodes = data.get('elements', {}).get('nodes', [])
        edges = data.get('elements', {}).get('edges', [])

        if not nodes:
            print("  Warning: No nodes found to export.")
            nodes_df = pd.DataFrame()
        else:
            nodes_df = pd.DataFrame([n.get('data', {}) for n in nodes])

        if not edges:
             print("  Warning: No edges found to export.")
             edges_df = pd.DataFrame()
        else:
            edges_df = pd.DataFrame([e.get('data', {}) for e in edges])

        os.makedirs(os.path.dirname(output_excel_path), exist_ok=True)
        with pd.ExcelWriter(output_excel_path, engine='openpyxl') as writer:
            nodes_df.to_excel(writer, sheet_name='Nodes', index=False)
            edges_df.to_excel(writer, sheet_name='Edges', index=False)

        print("  Conversion Successful.")

    except Exception as e:
        raise RuntimeError(f"Failed to convert JSON to Excel: {e}")