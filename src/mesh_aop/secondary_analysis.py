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

from . import paths as _paths


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


def _fetch_relevance_rows(conn, clause, params):
    """Fetch (pmid, ARS, contributing_seeds, pub_year) from article_relevance_scores.

    Includes the publication year when the relevance DB has a pub_date column; a
    legacy DB without it yields pub_year=None, so the impact engine falls back to a
    raw citation count instead of a citations-per-year rate. `clause` is the
    trailing SQL (WHERE/ORDER/LIMIT) appended to the SELECT.
    """
    cols = [r[1] for r in conn.execute("PRAGMA table_info(article_relevance_scores)").fetchall()]
    has_pubdate = 'pub_date' in cols
    select_cols = "pmid, score_pagerank_centrality, contributing_seeds" + (", pub_date" if has_pubdate else "")
    rows = conn.execute(f"SELECT {select_cols} FROM article_relevance_scores {clause}", params).fetchall()

    out = []
    for r in rows:
        year = None
        if has_pubdate and r[3]:
            try:
                year = int(str(r[3])[:4])
            except (ValueError, TypeError):
                year = None
        out.append((r[0], r[1], r[2], year))
    return out

def _calculate_impact_scores(raw_rows: list, sort_metric: str, linear_weight_ars: float, oversample_limit: int) -> list:
    """
    Pandas engine to apply dynamic scaling and dual-metric Article Impact Scoring.

    The citation axis is the citations-per-year RATE (total citations / years since
    publication, with a minimum-age floor) so that older papers are not favoured
    purely for having had longer to accumulate citations. This rate feeds both the
    Linear (weighted-average) and F1 (harmonic-mean) engines identically to the old
    raw-count signal - only the input changes. Rows whose publication year is
    unknown (e.g. a legacy relevance DB lacking a pub_date column) fall back to the
    raw citation count, so behaviour is unchanged for those.
    """
    if not raw_rows:
        return []

    df = pd.DataFrame(raw_rows, columns=['pmid', 'ars_score', 'seeds', 'cited_by', 'pub_year'])

    # 1. Raw incoming-citation counts (kept for display/transparency).
    cited_by_str = df['cited_by'].fillna('').astype(str).str.strip()
    is_empty = (cited_by_str == '') | (cited_by_str.str.lower() == 'nan') | (cited_by_str.str.lower() == 'none')
    raw_counts = np.where(is_empty, 0, cited_by_str.str.count(';') + 1)
    df['Incoming_Citations_Smoothed'] = raw_counts

    # 2. Citations-per-year rate with a minimum-age floor so very recent papers are
    #    not inflated. Unknown publication year -> fall back to the raw count.
    MIN_AGE_YEARS = 3
    current_year = time.localtime().tm_year
    pub_year = pd.to_numeric(df['pub_year'], errors='coerce')
    age = (current_year - pub_year).clip(lower=MIN_AGE_YEARS)
    citation_signal = np.where(pub_year.notna(), raw_counts / age, raw_counts)
    df['Citations_Per_Year'] = np.where(pub_year.notna(), citation_signal, np.nan)

    # 3. Log-Normalize the citation signal (rate where the year is known).
    df['log_cit'] = np.log10(citation_signal + 1)

    # 4. Dynamic Scaling based on the Network's Maximum
    max_log = df['log_cit'].max()

    if max_log > 0:
        one_citation_scale = np.log10(2) / max_log
        dynamic_floor = one_citation_scale / 2.0
        normalized = df['log_cit'] / max_log
        df['Normalized_Citation'] = np.where(citation_signal == 0, dynamic_floor, normalized)
    else:
        dynamic_floor = 1.0
        df['Normalized_Citation'] = 1.0

    df['ars_score'] = df['ars_score'].clip(lower=dynamic_floor)

    # 5. Dual Calculation Engine (unchanged structure; citation axis is now a rate)
    df['Linear_AIS'] = (df['ars_score'] * linear_weight_ars) + (df['Normalized_Citation'] * (1.0 - linear_weight_ars))
    df['F1_AIS'] = (2 * df['ars_score'] * df['Normalized_Citation']) / (df['ars_score'] + df['Normalized_Citation'])

    # 6. Sort by user's chosen metric and truncate
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

    # Blank credentials must reach Biopython as None, never as "" - an empty
    # api_key is sent as a real parameter and NCBI answers 400. See
    # data_ops.set_entrez_credentials.
    from .data_ops import set_entrez_credentials
    set_entrez_credentials(entrez_email, entrez_api_key)

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
                    "Citations_Per_Year": base_data.get('Citations_Per_Year', None),
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
                           sort_metric: str = 'F1', linear_weight_ars: float = 0.5):
    """Pulls and scores all articles associated with a specific node."""
    if not os.path.exists(db_path) or not os.path.exists(cleaned_db_path):
        raise FileNotFoundError("Missing required databases.")

    print(f"\n<<< Searching for Node: '{node_name}' >>>")
    oversample_limit = limit * 5

    try:
        # Step 1: Query Relevance DB
        conn_rel = _open_readonly_resilient(db_path)
        rel_rows = _fetch_relevance_rows(conn_rel, "WHERE contributing_seeds LIKE ?", (f'%{node_name}%',))
        conn_rel.close()

        if not rel_rows:
            print(f"  -> No articles found for node '{node_name}'.")
            return

        # Step 2: Clean PMIDs and map data
        clean_pmids = [str(r[0]).split('.')[0] for r in rel_rows]
        rel_map = {clean_pmids[i]: (rel_rows[i][1], rel_rows[i][2], rel_rows[i][3]) for i in range(len(clean_pmids))}

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
        raw_rows = [(pmid, rel_map[pmid][0], rel_map[pmid][1], citation_map.get(pmid, None), rel_map[pmid][2]) for pmid in clean_pmids]

        pmid_dicts = _calculate_impact_scores(raw_rows, sort_metric, linear_weight_ars, oversample_limit)
        print(f"  -> Extracted & Scored top {len(pmid_dicts)} candidates. Hydrating metadata...")

        df = _fetch_metadata_and_filter(pmid_dicts, exclude_reviews, limit, entrez_email, entrez_api_key, sort_metric)

        if not df.empty:
            print(f"  Top {len(df)} results (Ranked by {sort_metric} Article Impact Score (AIS)):")
            sort_col = 'Linear_AIS' if sort_metric.lower() == 'linear' else 'F1_AIS'
            print(df[['PMID', 'Year', 'Raw_ARS', 'Incoming_Citations', 'Citations_Per_Year', sort_col]].head(5).to_string(index=False))

            os.makedirs(results_dir, exist_ok=True)
            safe_name = node_name.replace(" ", "_").replace("/", "-")
            out_file = os.path.join(results_dir, f"{file_prefix}_relevance_{safe_name}.csv")
            df.to_csv(_paths.long_path(out_file), index=False)
            print(f"  Saved full results to: {os.path.basename(out_file)}")
        else:
            print("  -> No valid articles remained after applying filters.")

    except Exception as e:
        raise RuntimeError(f"Error querying database for node relevancy: {e}")


def analyze_edge_relevancy(node1: str, node2: str, db_path: str, cleaned_db_path: str, results_dir: str, file_prefix: str,
                           limit: int, exclude_reviews: bool, entrez_email: str, entrez_api_key: str,
                           sort_metric: str = 'F1', linear_weight_ars: float = 0.5):
    """Pulls and scores all articles associated with a specific EDGE (relationship)."""
    if not os.path.exists(db_path) or not os.path.exists(cleaned_db_path):
        raise FileNotFoundError("Missing required databases.")

    print(f"\n<<< Searching for Edge: '{node1}' <--> '{node2}' >>>")
    oversample_limit = limit * 5

    try:
        conn_rel = _open_readonly_resilient(db_path)
        rel_rows = _fetch_relevance_rows(conn_rel, "WHERE contributing_seeds LIKE ? AND contributing_seeds LIKE ?", (f'%{node1}%', f'%{node2}%'))
        conn_rel.close()

        if not rel_rows:
            print(f"  -> No articles found linking '{node1}' and '{node2}'.")
            return

        clean_pmids = [str(r[0]).split('.')[0] for r in rel_rows]
        rel_map = {clean_pmids[i]: (rel_rows[i][1], rel_rows[i][2], rel_rows[i][3]) for i in range(len(clean_pmids))}

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

        raw_rows = [(pmid, rel_map[pmid][0], rel_map[pmid][1], citation_map.get(pmid, None), rel_map[pmid][2]) for pmid in clean_pmids]

        pmid_dicts = _calculate_impact_scores(raw_rows, sort_metric, linear_weight_ars, oversample_limit)
        print(f"  -> Extracted & Scored top {len(pmid_dicts)} candidates. Hydrating metadata...")

        df = _fetch_metadata_and_filter(pmid_dicts, exclude_reviews, limit, entrez_email, entrez_api_key, sort_metric)

        if not df.empty:
            print(f"  Top {len(df)} results (Ranked by {sort_metric} Article Impact Score (AIS)):")
            sort_col = 'Linear_AIS' if sort_metric.lower() == 'linear' else 'F1_AIS'
            print(df[['PMID', 'Year', 'Raw_ARS', 'Incoming_Citations', 'Citations_Per_Year', sort_col]].head(5).to_string(index=False))

            os.makedirs(results_dir, exist_ok=True)
            safe_n1 = node1.replace(" ", "_").replace(",", "")
            safe_n2 = node2.replace(" ", "_").replace(",", "")
            out_file = os.path.join(results_dir, f"{file_prefix}_edge_relevance_{safe_n1}_{safe_n2}.csv")
            df.to_csv(_paths.long_path(out_file), index=False)
            print(f"  Saved full results to: {os.path.basename(out_file)}")
        else:
            print("  -> No valid articles remained after applying filters.")

    except Exception as e:
        raise RuntimeError(f"Error querying database for edge relevancy: {e}")


def get_top_network_articles(db_path: str, cleaned_db_path: str, results_dir: str, file_prefix: str, limit: int,
                             exclude_reviews: bool, entrez_email: str, entrez_api_key: str,
                             sort_metric: str = 'F1', linear_weight_ars: float = 0.5):
    """Retrieves and scores articles across the entire refined network."""
    if not os.path.exists(db_path) or not os.path.exists(cleaned_db_path):
        raise FileNotFoundError("Missing required databases.")

    print(f"\n<<< Extracting Top {limit} Network-Wide Articles >>>")
    oversample_limit = limit * 5

    try:
        conn_rel = _open_readonly_resilient(db_path)
        rel_rows = _fetch_relevance_rows(conn_rel, "ORDER BY score_pagerank_centrality DESC LIMIT 50000", ())
        conn_rel.close()

        if not rel_rows:
            print("  -> No articles found in the relevance database.")
            return

        clean_pmids = [str(r[0]).split('.')[0] for r in rel_rows]
        rel_map = {clean_pmids[i]: (rel_rows[i][1], rel_rows[i][2], rel_rows[i][3]) for i in range(len(clean_pmids))}

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

        raw_rows = [(pmid, rel_map[pmid][0], rel_map[pmid][1], citation_map.get(pmid, None), rel_map[pmid][2]) for pmid in clean_pmids]

        pmid_dicts = _calculate_impact_scores(raw_rows, sort_metric, linear_weight_ars, oversample_limit)

        df = _fetch_metadata_and_filter(pmid_dicts, exclude_reviews, limit, entrez_email, entrez_api_key, sort_metric)

        if not df.empty:
            print(f"  Top {len(df)} results (Ranked by {sort_metric} Article Impact Score (AIS)):")
            sort_col = 'Linear_AIS' if sort_metric.lower() == 'linear' else 'F1_AIS'
            print(df[['PMID', 'Year', 'Raw_ARS', 'Incoming_Citations', 'Citations_Per_Year', sort_col]].head(5).to_string(index=False))

            os.makedirs(results_dir, exist_ok=True)
            out_file = os.path.join(results_dir, f"{file_prefix}_Top_Network_Articles.csv")
            df.to_csv(_paths.long_path(out_file), index=False)
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
        with pd.ExcelWriter(_paths.long_path(output_excel_path), engine='openpyxl') as writer:
            nodes_df.to_excel(writer, sheet_name='Nodes', index=False)
            edges_df.to_excel(writer, sheet_name='Edges', index=False)

        print("  Conversion Successful.")

    except Exception as e:
        raise RuntimeError(f"Failed to convert JSON to Excel: {e}")


# Recognised graph formats besides the pipeline's own Cytoscape JSON, mapped to
# their networkx reader. Anything else is attempted as Cytoscape JSON.
_GRAPH_READERS = {
    ".graphml": "read_graphml",
    ".gml": "read_gml",
    ".gexf": "read_gexf",
    ".net": "read_pajek",
    ".pajek": "read_pajek",
    ".adjlist": "read_adjlist",
    ".edgelist": "read_edgelist",
}


def _parse_network_list(raw: str) -> list:
    """Split the wizard's comma-separated, quote-wrapped network list into names.

    Accepts `"a.json","b.graphml"` or a bare `a.json, b.graphml`; surrounding
    single/double quotes and whitespace are stripped and empties dropped.
    """
    if not raw:
        return []
    return [tok.strip().strip('"').strip("'").strip()
            for tok in raw.split(",") if tok.strip().strip('"').strip("'").strip()]


def _load_network_nodes(path: str) -> set:
    """Return the set of node ids in a network file (Cytoscape JSON or a networkx format)."""
    ext = os.path.splitext(path)[1].lower()
    if ext in _GRAPH_READERS:
        import networkx as nx
        graph = getattr(nx, _GRAPH_READERS[ext])(path)
        return set(str(n) for n in graph.nodes())
    # Default / fallback: the pipeline's Cytoscape JSON layout.
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    return {n.get("data", {}).get("id") for n in data.get("elements", {}).get("nodes", [])
            if n.get("data", {}).get("id") is not None}


def run_network_overlay_comparison(network_files, search_dir: str, results_dir: str,
                                   file_prefix: str, random_seed: int = 42):
    """Compare the node membership of an arbitrary set of networks.

    `network_files` is a list of file names (or a raw comma-separated string). Bare
    names are resolved against `search_dir` (data/processed, or data/reference_processed
    when reference data is in use); names carrying their own path are used as given.
    Networks may be the pipeline's Cytoscape JSON (default) or any networkx-readable
    format (.graphml/.gml/.gexf/...). Writes a membership matrix, a pairwise
    intersection/Jaccard table, and an overlay figure to `results_dir`. Non-JSON
    outputs (the PNG) stay in the results section, never in data/.
    """
    print("\n<<< Network Overlay Comparison >>>")

    names = _parse_network_list(network_files) if isinstance(network_files, str) else list(network_files)
    if not names:
        print("  [!] No networks were listed to compare; skipping.")
        return

    # Resolve each requested network and report the ones that cannot be found.
    resolved, missing = {}, []
    for name in names:
        path = name if os.path.dirname(name) else os.path.join(search_dir, name)
        if os.path.exists(path):
            resolved[name] = path
        else:
            missing.append(name)

    if missing:
        print(f"  [!] WARNING: {len(missing)} network file(s) could not be found:")
        for name in missing:
            print(f"        - {name}")
        print(f"      The system looked in: {os.path.abspath(search_dir)}")
        print("      Check the spelling and location of the file names above (bare names are")
        print("      resolved against that folder; include a path to load from elsewhere).")

    if len(resolved) < 2:
        print(f"  [!] Need at least two readable networks to compare; found {len(resolved)}. Skipping.")
        return

    # Load node sets, dropping any file that fails to parse (with a warning).
    sets = {}
    for name, path in resolved.items():
        try:
            sets[name] = _load_network_nodes(path)
            print(f"  Loaded '{name}': {len(sets[name]):,} nodes")
        except Exception as e:
            print(f"  [!] WARNING: could not read '{name}' ({e}); excluded from the comparison.")
    if len(sets) < 2:
        print(f"  [!] Fewer than two networks could be read; skipping.")
        return

    labels = list(sets.keys())
    os.makedirs(results_dir, exist_ok=True)

    # --- Membership matrix: every node x every network (Y/N + kept-in-N count) ---
    universe = sorted(set().union(*sets.values()))
    rows = []
    for node in universe:
        row = {"node": node}
        for lab in labels:
            row[lab] = "Y" if node in sets[lab] else "N"
        row["kept_in_N_networks"] = sum(node in sets[lab] for lab in labels)
        rows.append(row)
    membership_df = pd.DataFrame(rows, columns=["node"] + labels + ["kept_in_N_networks"])
    membership_df = membership_df.sort_values(["kept_in_N_networks", "node"],
                                              ascending=[False, True]).reset_index(drop=True)
    membership_path = os.path.join(results_dir, f"{file_prefix}_network_overlap_membership.csv")
    membership_df.to_csv(_paths.long_path(membership_path), index=False, encoding="utf-8-sig")

    # --- Pairwise intersection count and Jaccard similarity ---
    inter = pd.DataFrame(index=labels, columns=labels, dtype=float)
    jacc = pd.DataFrame(index=labels, columns=labels, dtype=float)
    for a in labels:
        for b in labels:
            common = len(sets[a] & sets[b])
            union = len(sets[a] | sets[b])
            inter.loc[a, b] = common
            jacc.loc[a, b] = (common / union) if union else 0.0
    matrix_path = os.path.join(results_dir, f"{file_prefix}_network_overlap_matrix.csv")
    with open(_paths.long_path(matrix_path), "w", encoding="utf-8-sig", newline="") as f:
        f.write("# Pairwise node-set intersection counts\n")
        inter.astype(int).to_csv(f)
        f.write("\n# Pairwise Jaccard similarity (|A n B| / |A u B|)\n")
        jacc.round(4).to_csv(f)

    # --- Figure: set sizes + pairwise Jaccard heatmap ---
    fig_path = os.path.join(results_dir, f"{file_prefix}_Network_Overlap.png")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        n = len(labels)
        fig, (ax_bar, ax_hm) = plt.subplots(
            1, 2, figsize=(max(9, 1.1 * n + 6), max(4.5, 0.5 * n + 3)),
            gridspec_kw={"width_ratios": [1, 1.3]})

        sizes = [len(sets[lab]) for lab in labels]
        ax_bar.barh(range(n), sizes, color="#4C72B0")
        ax_bar.set_yticks(range(n))
        ax_bar.set_yticklabels(labels, fontsize=9)
        ax_bar.invert_yaxis()
        ax_bar.set_xlabel("Node count")
        ax_bar.set_title("Network size")
        for i, s in enumerate(sizes):
            ax_bar.text(s, i, f" {s:,}", va="center", fontsize=8)

        jvals = jacc.values.astype(float)
        im = ax_hm.imshow(jvals, cmap="magma_r", vmin=0.0, vmax=1.0)
        ax_hm.set_xticks(range(n)); ax_hm.set_yticks(range(n))
        ax_hm.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
        ax_hm.set_yticklabels(labels, fontsize=8)
        ax_hm.set_title("Pairwise Jaccard overlap")
        for i in range(n):
            for j in range(n):
                ax_hm.text(j, i, f"{jvals[i, j]:.2f}", ha="center", va="center",
                           fontsize=7, color="white" if jvals[i, j] > 0.5 else "black")
        fig.colorbar(im, ax=ax_hm, fraction=0.046, pad=0.04)

        fig.suptitle("Network node-overlap comparison", fontsize=13)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        fig.savefig(fig_path, dpi=200)
        plt.close(fig)
    except Exception as e:
        print(f"  [!] WARNING: overlap figure could not be generated ({e}); tables were still written.")
        fig_path = None

    print(f"  [+] Compared {len(sets)} networks over {len(universe):,} distinct nodes.")
    print(f"      Membership matrix : {membership_path}")
    print(f"      Overlap matrix    : {matrix_path}")
    if fig_path:
        print(f"      Overlay figure    : {fig_path}")