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
                                     id_key: str, weight_key_1: str, final_key_1: str,
                                     weight_key_2: str, final_key_2: str,
                                     start_date_param: str, end_date_param: str,
                                     entrez_email: str = "", entrez_api_key: str = "",
                                     calculate_full_centrality: bool = True,
                                     weight_key_3: str = None, final_key_3: str = None):
    """
    Calculates contextual relevance by querying the local Master Database
    over the specific time constraints, completely bypassing API limits.
    """
    if not calculate_full_centrality:
        print("\n[!] NOTICE: User skipped centrality. CRS scores now reflect TOPIC DENSITY (term overlap), not topological relevance.\n")

    print(f"\n<<< Loading Seed Terms from {os.path.basename(input_nodes_file)} >>>")
    seed_weights_1, total_weight_1, seed_terms = {}, 0, set()
    seed_weights_2, total_weight_2 = {}, 0
    seed_weights_3, total_weight_3 = {}, 0

    try:
        with open(input_nodes_file, 'r') as f:
            nodes = json.load(f).get('elements', {}).get('nodes', [])
        for node in nodes:
            data = node.get('data', {})
            term, weight_1, weight_2 = data.get(id_key), data.get(weight_key_1), data.get(weight_key_2)
            weight_3 = data.get(weight_key_3) if weight_key_3 else None
            if term:
                seed_terms.add(term)
                if isinstance(weight_1, (int, float)):
                    seed_weights_1[term] = float(weight_1)
                    total_weight_1 += float(weight_1)
                if isinstance(weight_2, (int, float)):
                    seed_weights_2[term] = float(weight_2)
                    total_weight_2 += float(weight_2)
                if isinstance(weight_3, (int, float)):
                    seed_weights_3[term] = float(weight_3)
                    total_weight_3 += float(weight_3)
        print(f"Loaded {len(seed_terms)} unique seed terms for relevance scoring.")
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

        aggregator_1 = defaultdict(lambda: {'sum': 0.0, 'count': 0})
        aggregator_2 = defaultdict(lambda: {'sum': 0.0, 'count': 0})
        aggregator_3 = defaultdict(lambda: {'sum': 0.0, 'count': 0})
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

                score_1, score_2, score_3 = 0.0, 0.0, 0.0
                if total_weight_1 > 0:
                    score_1 = sum(seed_weights_1.get(term, 0) for term in matching_seeds) / total_weight_1
                    for term in matching_seeds:
                        aggregator_1[term]['sum'] += score_1
                        aggregator_1[term]['count'] += 1

                if total_weight_2 > 0:
                    score_2 = sum(seed_weights_2.get(term, 0) for term in matching_seeds) / total_weight_2
                    for term in matching_seeds:
                        aggregator_2[term]['sum'] += score_2
                        aggregator_2[term]['count'] += 1

                if total_weight_3 > 0:
                    score_3 = sum(seed_weights_3.get(term, 0) for term in matching_seeds) / total_weight_3
                    for term in matching_seeds:
                        aggregator_3[term]['sum'] += score_3
                        aggregator_3[term]['count'] += 1

                article_row = {
                    'pmid': pmid,
                    f'score_{weight_key_1}': score_1,
                    f'score_{weight_key_2}': score_2,
                    'contributing_seeds': ';'.join(sorted(list(matching_seeds))),
                    'pub_date': pub_date
                }
                if weight_key_3:
                    article_row[f'score_{weight_key_3}'] = score_3
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
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_score1 ON article_relevance_scores ('score_{weight_key_1}')")
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

    raw_node_weights_1 = {}
    for term, data in aggregator_1.items():
        if data['count'] > 0:
            mean_ars = data['sum'] / data['count']
            cf_vol = np.log10(data['count'] + 1)                 # CF: in-window evidence volume |P_i|
            g_t = global_term_freq.get(term, data['count'])      # G_t: global document frequency
            ic_spec = -np.log10(g_t / N_global)                  # IC: specificity across the whole corpus

            # Harmonic Mean of Confidence Factor and Information Content
            denom = cf_vol + ic_spec
            harmonic_multiplier = (2 * cf_vol * ic_spec) / denom if denom > 0 else 0.0

            raw_node_weights_1[term] = mean_ars * harmonic_multiplier

    raw_node_weights_2 = {}
    for term, data in aggregator_2.items():
        if data['count'] > 0:
            mean_ars = data['sum'] / data['count']
            cf_vol = np.log10(data['count'] + 1)
            g_t = global_term_freq.get(term, data['count'])
            ic_spec = -np.log10(g_t / N_global)

            denom = cf_vol + ic_spec
            harmonic_multiplier = (2 * cf_vol * ic_spec) / denom if denom > 0 else 0.0

            raw_node_weights_2[term] = mean_ars * harmonic_multiplier

    # Optional third weighting (e.g. PageRank). Same CF/IC multiplier - only the
    # mean-ARS weighting differs.
    raw_node_weights_3 = {}
    if weight_key_3:
        for term, data in aggregator_3.items():
            if data['count'] > 0:
                mean_ars = data['sum'] / data['count']
                cf_vol = np.log10(data['count'] + 1)
                g_t = global_term_freq.get(term, data['count'])
                ic_spec = -np.log10(g_t / N_global)
                denom = cf_vol + ic_spec
                harmonic_multiplier = (2 * cf_vol * ic_spec) / denom if denom > 0 else 0.0
                raw_node_weights_3[term] = mean_ars * harmonic_multiplier

    # --- Apply Relative Max-Scaling Normalization [0.0 to 1.0] ---
    final_node_weights_1 = {}
    max_crs_1 = max(raw_node_weights_1.values()) if raw_node_weights_1 else 1.0
    for term, val in raw_node_weights_1.items():
        final_node_weights_1[term] = val / max_crs_1

    final_node_weights_2 = {}
    max_crs_2 = max(raw_node_weights_2.values()) if raw_node_weights_2 else 1.0
    for term, val in raw_node_weights_2.items():
        final_node_weights_2[term] = val / max_crs_2

    final_node_weights_3 = {}
    max_crs_3 = max(raw_node_weights_3.values()) if raw_node_weights_3 else 1.0
    for term, val in raw_node_weights_3.items():
        final_node_weights_3[term] = val / max_crs_3

    print("\n<<< Generating Final Output JSON File >>>")
    try:
        with open(input_nodes_file, 'r') as f:
            network_data = json.load(f)

        for node in network_data.get('elements', {}).get('nodes', []):
            term = node.get('data', {}).get(id_key)
            node['data'][final_key_1] = final_node_weights_1.get(term, 0.0)
            node['data'][final_key_2] = final_node_weights_2.get(term, 0.0)
            if final_key_3:
                node['data'][final_key_3] = final_node_weights_3.get(term, 0.0)

        os.makedirs(os.path.dirname(output_nodes_file), exist_ok=True)
        with open(output_nodes_file, 'w') as f:
            json.dump(network_data, f, indent=2)

        print(f"  [+] Success! Final relevance JSON created at: {os.path.basename(output_nodes_file)}")

    except Exception as e:
        raise RuntimeError(f"ERROR generating final relevance JSON: {e}")