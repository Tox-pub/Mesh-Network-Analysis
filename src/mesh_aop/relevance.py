# -*- coding: utf-8 -*-
"""
relevance.py - contextual relevance scoring / semantic re-ranking (pipeline Step 3).

Re-scores the network's nodes by how strongly they are supported across the
wider literature within a chosen time window, querying the local master database
so the calculation needs no API calls.

For each surviving node it aggregates per-article term-overlap scores into a
Contextual Relevance Score that blends topological weight, evidence volume, and
term specificity (an information-content penalty applied via a harmonic mean),
then writes both the per-article scores and the relevance-annotated network to disk.
"""

import os
import json
import sqlite3
from collections import defaultdict

import numpy as np
import pandas as pd
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

def run_contextual_relevance_scoring(input_nodes_file: str, output_nodes_file: str,
                                     master_db_path: str, relevance_db_path: str,
                                     id_key: str,
                                     start_date_param: str, end_date_param: str,
                                     weightings: list = None,
                                     weight_key_1: str = None, final_key_1: str = None,
                                     weight_key_2: str = None, final_key_2: str = None,
                                     weight_key_3: str = None, final_key_3: str = None,
                                     entrez_email: str = "", entrez_api_key: str = "",
                                     calculate_full_centrality: bool = True):
    """
    Calculates contextual relevance by querying the local Master Database
    over the specific time constraints, completely bypassing API limits.

    `weightings` is a list of (node_attribute, output_attribute) pairs, one per node
    weighting to score. Each produces a per-article `score_<node_attribute>` column
    and a per-node `<output_attribute>` CRS value, so any number of weightings can be
    compared from a single pass over the corpus - the corpus scan dominates the cost,
    so adding a weighting is far cheaper than re-running. The legacy
    weight_key_1/2/3 arguments are folded into this list when it is not supplied.
    """
    if weightings is None:
        weightings = [(w, f) for w, f in ((weight_key_1, final_key_1),
                                          (weight_key_2, final_key_2),
                                          (weight_key_3, final_key_3)) if w and f]
    if not weightings:
        raise ValueError("No node weightings supplied to relevance scoring.")

    if not calculate_full_centrality:
        print("\n[!] NOTICE: Betweenness centrality was estimated from a sampled subset of nodes "
              "rather than computed exactly. Eigenvector and PageRank centrality are unaffected, "
              "and CRS still reflects topological relevance.\n")

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
    print("<<< Calculating Contextual Relevance (Local Data Lake) >>>")
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

        # Global corpus size for Information Content (IC): the FULL master DB,
        # independent of the analysis window, so IC measures a term's specificity
        # across all of the literature rather than within the chosen date range.
        print(f"  Scanning 30-million record database...")
        m_cursor.execute("SELECT count(*) FROM master_mesh_annotations")
        N_global = float(m_cursor.fetchone()[0])

        # In-window count drives the ARS pass (which articles get scored / form
        # each term's P_i) and the zero-check.
        m_cursor.execute("SELECT count(*) FROM master_mesh_annotations WHERE pub_date BETWEEN ? AND ?", (start_iso, end_iso))
        total_articles_in_range = m_cursor.fetchone()[0]

        if total_articles_in_range == 0:
            print("  [!] Warning: Zero articles found in Master Database for this date range.")
            master_conn.close()
            return

        print(f"  Found {total_articles_in_range:,} articles in window; scanning the full {int(N_global):,}-article corpus for global term specificity...")

        aggregators = [defaultdict(lambda: {'sum': 0.0, 'count': 0}) for _ in range(n_w)]
        global_term_freq = defaultdict(int)
        article_scores_data = []

        # Single full-corpus pass: accumulate each seed term's GLOBAL document
        # frequency (G_t, for IC) over every article, and run ARS scoring only on
        # the in-window articles (which form each term's P_i for CF and mean ARS).
        m_cursor.execute("SELECT pmid, pub_date, mesh_terms FROM master_mesh_annotations")
        chunk_size = 100000
        pbar = tqdm(total=int(N_global), desc="Scanning corpus")

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

                # GLOBAL term frequency (all dates) -> G_t for the IC denominator.
                for term in matching_seeds:
                    global_term_freq[term] += 1

                # ARS scoring only for in-window articles.
                if not (start_iso <= pub_date <= end_iso):
                    continue

                article_row = {'pmid': pmid}
                for i, (w_key, _) in enumerate(weightings):
                    score = 0.0
                    if total_weight[i] > 0:
                        score = sum(seed_weights[i].get(term, 0)
                                    for term in matching_seeds) / total_weight[i]
                        agg = aggregators[i]
                        for term in matching_seeds:
                            agg[term]['sum'] += score
                            agg[term]['count'] += 1
                    article_row[f'score_{w_key}'] = score

                article_row['contributing_seeds'] = ';'.join(sorted(list(matching_seeds)))
                article_row['pub_date'] = pub_date
                article_scores_data.append(article_row)

            pbar.update(len(rows))

        pbar.close()
        master_conn.close()

    except Exception as e:
        raise RuntimeError(f"ERROR calculating scores from Master DB: {e}")

    print("\n<<< Saving Individual Article Relevance Scores to Database >>>")
    if article_scores_data:
        scores_df = pd.DataFrame(article_scores_data)
        try:
            os.makedirs(os.path.dirname(relevance_db_path), exist_ok=True)
            conn = sqlite3.connect(relevance_db_path)
            scores_df.to_sql('article_relevance_scores', conn, if_exists='replace', index=False)
            cursor = conn.cursor()
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_pmid ON article_relevance_scores (pmid)")
            for i, (w_key, _) in enumerate(weightings, start=1):
                cursor.execute(
                    f"CREATE INDEX IF NOT EXISTS idx_score{i} "
                    f"ON article_relevance_scores ('score_{w_key}')"
                )
            conn.commit()
            conn.close()
            _keep_on_device(relevance_db_path)
            print(f"  [+] Successfully saved {len(scores_df):,} contributing article scores to the database.")
        except Exception as e:
            print(f"  [!] WARNING: Could not save article scores to database: {e}")
    else:
        print("  [!] No articles contained the target seed terms. No scores generated.")

# Calculate final Contextual Relevance Scores (CRS) with Harmonic Mean & Max-Scaling
    print("\n<<< Calculating CRS (Harmonic Mean Penalty & Max-Scaling)... >>>")

    # CF (in-window evidence volume) and IC (global specificity) depend only on the
    # term, not on the weighting, so the harmonic multiplier is identical across
    # weightings and only the mean ARS differs.
    final_node_weights = []
    for i, (w_key, _) in enumerate(weightings):
        raw_node_weights = {}
        for term, data in aggregators[i].items():
            if data['count'] > 0:
                mean_ars = data['sum'] / data['count']
                cf_vol = np.log10(data['count'] + 1)             # CF: in-window evidence volume |P_i|
                g_t = global_term_freq.get(term, data['count'])  # G_t: global document frequency
                ic_spec = -np.log10(g_t / N_global)              # IC: specificity across the whole corpus

                # Harmonic Mean of Confidence Factor and Information Content
                denom = cf_vol + ic_spec
                harmonic_multiplier = (2 * cf_vol * ic_spec) / denom if denom > 0 else 0.0

                raw_node_weights[term] = mean_ars * harmonic_multiplier

        # --- Apply Relative Max-Scaling Normalization [0.0 to 1.0] ---
        max_crs = max(raw_node_weights.values()) if raw_node_weights else 1.0
        final_node_weights.append(
            {term: val / max_crs for term, val in raw_node_weights.items()} if max_crs else {}
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