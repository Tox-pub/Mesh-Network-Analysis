"""
master_mesh_network.py

Master script used to execute the complete data collection and
processing pipeline for the MeSH AOP Network analysis.

This script requires access to the following files:
requirements.txt
scripts/python/config.py
scripts/python/mesh_stop_words.py

This script performs:
1. Initial Data Collection (PubMed Search & Citation Expansion)
2. Database Cleaning & Master MeSH Annotation
3. Full Network Construction (Co-occurrence & Weight Calculation)
4. Network Optimization Simulations (GLF & SA Algorithms)
5. Consensus Network Filtering & LCC Extraction
6. Community Detection (Louvain Algorithm)
7. Contextual Relevance Scoring (Semantic Re-ranking)

This saves processed networks and databases in:
data/processed/
results/logs/
figures/
"""
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 0. Preamble: Imports and Environment Check
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

import csv
import itertools
import json
import logging
import math
import os
import random
import re
import shutil
import sqlite3
import sys
import time
from collections import deque, defaultdict
from datetime import date, datetime, timedelta
from urllib.error import HTTPError

import community as community_louvain # You might need to uninstall and reinstall on VMs
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import requests
from Bio import Entrez, Medline
from scipy import stats
from scipy.stats import binom
from tqdm.notebook import tqdm # Notebook progress bars


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 1. Configuration Block
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# <<< Import Configuration (Step 1)>>>
if "__file__" in globals():
    current_dir = os.path.dirname(os.path.abspath(__file__))
else:
    current_dir = os.path.abspath(os.getcwd())

if current_dir not in sys.path:
    sys.path.append(current_dir)

import config # See config.py file for adjusting parameters
import mesh_stop_words # MeSH stop words list

# <<< Import MeSH Stop Words >>>
try:
    from mesh_stop_words import Mesh_stop_words # mesh_stop_words.py
    print(f"Loaded {len(Mesh_stop_words)} MeSH stop words.")
except ImportError:
    print("ERROR: Could not find 'mesh_stop_words.py'.")
    sys.exit(1)


# <<< Entrez Configuration (Step 2)>>>
Entrez.email = config.ENTREZ_EMAIL
Entrez.api_key = config.ENTREZ_API_KEY
if "YOUR_" in Entrez.email:
    print("WARNING: Entrez credentials not set in config.py!")

# <<< User Search Parameters (Step 3)>>>
search_term = config.SEARCH_TERM
start_date = config.START_DATE
end_date = config.END_DATE
generations_n = config.GENERATIONS_N
BATCH_SIZE = 100  # Internal batching
start_generation_to_remove = f"G{generations_n}" # Only sends designated generations for Step 2 onward.


# <<< Parameters for Consensus Network & Centrality (Step 5) >>>
BETWEENNESS_K_SAMPLES = config.BETWEENNESS_K_SAMPLES
TARGET_NUM_EDGES = config.TARGET_NUM_EDGES
GLF_ITERATIONS = config.GLF_ITERATIONS
SA_ITERATIONS = config.SA_ITERATIONS
SA_INITIAL_TEMPERATURE = config.SA_TEMP_START
SA_COOLING_RATE = config.SA_COOLING_RATE
RANDOM_SEED = config.RANDOM_SEED

# <<< Parameters for Contextual Relevance (Step 7) >>>
CONTEXT_START_DATE = config.CONTEXT_START_DATE
CONTEXT_END_DATE = config.CONTEXT_END_DATE
CONTEXT_INITIAL_WEIGHT_KEY_COL_1 = "betweenness_centrality" # CR Score parameter #1.
CONTEXT_INITIAL_WEIGHT_KEY_COL_2 = "eigenvector_centrality" # CR Score parameter #2.
CONTEXT_FINAL_WEIGHT_KEY_COL_1 = f"CRS_{CONTEXT_INITIAL_WEIGHT_KEY_COL_1}"
CONTEXT_FINAL_WEIGHT_KEY_COL_2 = f"CRS_{CONTEXT_INITIAL_WEIGHT_KEY_COL_2}"
MAX_PMIDS_PER_TERM = None # "None" = No PMID limit.

# <<< General NCBI Parameters (Step 8) >>>
ENTREZ_DELAY = 0.11           # Entrez limits to 10 calls a sec
RETMAX_CHUNK = 9999           # Chunks cannot exceed 9999 for Entrez calls
MESH_FETCH_BATCH_SIZE = 9999  # For contextual relevance step

# <<< Harmonized File Names (Step 9)>>>

FILE_PREFIX = config.FILE_PREFIX

original_db_path = config.FILES['pmids_db']
cleaned_db_path = config.FILES['cleaned_db']
MASTER_ANNOTATION_DB_PATH = config.FILES['master_db']
NETWORK_JSON_FULL = config.FILES['full_network']
GLF_OUTPUT_FILE = config.FILES['glf_subgraph']
SA_OUTPUT_FILE = config.FILES['sa_subgraph']
CONSENSUS_LCC_FILE = config.FILES['consensus_lcc']
CONTEXTUAL_RELEVANCE_DB_PATH = config.FILES['relevance_db']
FINAL_OUTPUT_WITH_RELEVANCE = config.FILES['final_network']
print(f"Project files will be stored in: {config.PROCESSED_DIR}")

# <<< Logging Setup (Step 10)>>>
LOG_FILE_PATH = config.FILES['log_file']
FAILED_MESH_FETCH_TSV = config.FILES['failed_mesh']
EMPTY_MESH_PMIDS_TSV = config.FILES['empty_mesh']

logging.basicConfig(filename=LOG_FILE_PATH, level=logging.WARNING,
                    format='%(asctime)s - %(levelname)s - %(message)s')
print(f"Logging output at: {LOG_FILE_PATH}")

# <<< Parameters for Network Building (Step 11) >>>
DB_PATH_FOR_NETWORK_BUILD = cleaned_db_path
LAMBDA_VALUE = 1.0        # Drives log value of generation_weight of nodes
NODE_WEIGHT_FACTORS = {   # Optional weights for adjusted_node_weight column
    'centrality': 0.45,
    'article_rank': 0.15,
    'rank_median_cit': 0.20,
    'rank_total_cit': 0.20
}

print(f"Loaded configuration for search: {search_term}")


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 2. Google Drive Mounting
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# This section was removed for the final published version

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 3. MeSH Stop Words List
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Combined MeSH Stop Word List
# Generated by parsing a MeSH ASCII descriptor file.
Mesh_stop_words = mesh_stop_words.Mesh_stop_words

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 4. Helper and Utility Functions
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def get_generation_label(numeric_generation):
    """Converts numeric generation (0, 1, 2...) to a label ('P0', 'G1', 'G2'...)."""
    if numeric_generation == 0:
        return 'P0'
    elif numeric_generation > 0:
        return f'G{numeric_generation}'
    return 'Unknown'

def parse_date_robust(date_str):
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

def format_date_for_entrez(date_obj):
    """Formats a date object into YYYY/MM/DD for Entrez."""
    return date_obj.strftime("%Y/%m/%d")

def _entrez_search_with_retry(db, term, retmax, **kwargs):
    """
    Wraps Entrez.esearch with a retry mechanism for common server errors.
    """
    max_retries = 5
    for attempt in range(max_retries):
        try:
            handle = Entrez.esearch(db=db, term=term, retmax=retmax, **kwargs)
            results = Entrez.read(handle)
            handle.close()
            return results
        except HTTPError as e:
            if e.code in [429, 500, 502, 503, 504]: # Common temp server errors
                wait_time = 2 ** attempt
                print(f"  - NCBI server returned error {e.code}. Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                raise
        except Exception as e:
            wait_time = 2 ** attempt
            print(f"  - An unexpected error occurred: {e}. Retrying in {wait_time} seconds...")
            time.sleep(wait_time)

    print(f"  - All {max_retries} retries failed for term '{term}'. Skipping.")
    return {"Count": "0", "IdList": []}

def parse_mesh_terms(mesh_terms_string):
    """
    Parses a string of MeSH terms, splitting them into main terms and subheadings,
    and identifies major topics. Handles '*' marker and '/'.
    """
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

def calculate_generation_weight(generation, lambda_val):
    """Calculates centrality weight based on generation string ('P0', 'G1', ...)."""
    try:
        if not isinstance(generation, str):
            if pd.notna(generation):
                distance = int(generation)
            else:
                distance = float('inf')
        elif generation.upper() == 'P0':
            distance = 0
        elif generation.upper().startswith('G'):
            distance = int(generation[1:])
        else:
            distance = float('inf')

        if distance == float('inf'):
            return 0.0
        else:
            distance = max(0, distance)
            return np.exp(-lambda_val * distance)
    except (ValueError, TypeError):
        print(f"Warning: Could not parse generation '{generation}'. Assigning 0 centrality.")
        return 0.0

def load_full_graph_data(filepath):
    """Loads a full graph from a Cytoscape.js JSON file, preserving all node and edge attributes."""
    all_nodes_data = {}
    all_edges_data = {}
    print(f"Loading full graph data from {os.path.basename(filepath)}...")
    try:
        with open(filepath, 'r') as f:
            cytoscape_data = json.load(f)
        nodes_list = cytoscape_data.get('elements', {}).get('nodes', [])
        for item in nodes_list:
            node_data = item.get('data', {})
            node_id = node_data.get('id')
            if node_id:
                all_nodes_data[node_id] = node_data
        edges_list = cytoscape_data.get('elements', {}).get('edges', [])
        for item in edges_list:
            edge_data = item.get('data', {})
            source = edge_data.get('source')
            target = edge_data.get('target')
            if source and target:
                edge_key = tuple(sorted((source, target)))
                all_edges_data[edge_key] = edge_data
        print(f" -> Loaded {len(all_nodes_data)} nodes and {len(all_edges_data)} edges.")
        return all_nodes_data, all_edges_data
    except FileNotFoundError:
        print(f" -> ERROR: Graph file not found at {filepath}")
        return {}, {}

def calculate_graph_stats(edges_data):
    """Calculates T and k_i for a graph from the full edge data dictionary."""
    node_strengths = defaultdict(int)
    total_strength = 0
    for edge_key, data in edges_data.items():
        weight = data.get('cooccurrence_count', 0)
        source, target = edge_key
        node_strengths[source] += weight
        node_strengths[target] += weight
        total_strength += weight
    return total_strength, node_strengths

def get_log_likelihood_term(edge_data, global_node_strengths, denominator):
    """Calculates the log-likelihood component for a single edge."""
    source = edge_data.get('source')
    target = edge_data.get('target')
    w_ij = edge_data.get('cooccurrence_count', 0)
    k_i = global_node_strengths.get(source, 0)
    k_j = global_node_strengths.get(target, 0)
    if k_i == 0 or k_j == 0: return 0
    p_ij = (k_i * k_j) / denominator
    if p_ij > 0 and w_ij > 0:
        return w_ij * np.log(p_ij) - math.lgamma(w_ij + 1)
    return 0

def run_glf_simulation(all_edges_data, global_node_strengths, global_T, target_num_edges, num_iterations, log_interval=10000):
    """Finds the minimum likelihood subgraph using the GLF Metropolis simulation."""
    print(f"\nStarting GLF simulation for {num_iterations:,} iterations...")
    history = []
    denominator = 2 * (global_T**2)
    all_edge_keys = sorted(list(all_edges_data.keys()))

    initial_selection = random.sample(all_edge_keys, target_num_edges)
    current_subgraph_keys = set(initial_selection)
    current_subgraph_list = sorted(initial_selection)

    current_T_sub = sum(all_edges_data[key].get('cooccurrence_count', 0) for key in current_subgraph_keys)
    current_sum_term = sum(get_log_likelihood_term(all_edges_data[key], global_node_strengths, denominator) for key in current_subgraph_keys)
    current_ll = math.lgamma(current_T_sub + 1) + current_sum_term

    best_subgraph_keys, min_ll = current_subgraph_keys.copy(), current_ll

    for i in tqdm(range(num_iterations), desc="GLF Simulation"):

        on_edge_index = random.randrange(len(current_subgraph_list))
        on_edge_key = current_subgraph_list[on_edge_index]

        off_edge_key = random.choice(all_edge_keys)
        while off_edge_key in current_subgraph_keys:
            off_edge_key = random.choice(all_edge_keys)

        w_on = all_edges_data[on_edge_key].get('cooccurrence_count', 0)
        w_off = all_edges_data[off_edge_key].get('cooccurrence_count', 0)

        proposed_T_sub = current_T_sub - w_on + w_off

        ll_term_on = get_log_likelihood_term(all_edges_data[on_edge_key], global_node_strengths, denominator)
        ll_term_off = get_log_likelihood_term(all_edges_data[off_edge_key], global_node_strengths, denominator)

        delta_ll = (math.lgamma(proposed_T_sub + 1) - math.lgamma(current_T_sub + 1)) + (ll_term_off - ll_term_on)

        if delta_ll < 0 or random.random() < math.exp(-delta_ll):
            current_subgraph_keys.remove(on_edge_key)
            current_subgraph_keys.add(off_edge_key)

            current_subgraph_list[on_edge_index] = off_edge_key

            current_ll += delta_ll
            current_T_sub = proposed_T_sub

            if current_ll < min_ll:
                min_ll = current_ll
                best_subgraph_keys = current_subgraph_keys.copy()

        if i % log_interval == 0:
            history.append((i, current_ll))

    return best_subgraph_keys, min_ll, history

def run_sa_simulation(all_edges_data, global_node_strengths, global_T, target_num_edges, num_iterations, initial_temp, cooling_rate, log_interval=10000):
    """Finds the minimum likelihood subgraph using Simulated Annealing."""
    print(f"\nStarting SA simulation for {num_iterations:,} iterations...")
    history, temperature = [], initial_temp
    denominator = 2 * (global_T**2)
    all_edge_keys = sorted(list(all_edges_data.keys()))

    initial_selection = random.sample(all_edge_keys, target_num_edges)
    current_subgraph_keys = set(initial_selection)
    current_subgraph_list = sorted(initial_selection)

    current_T_sub = sum(all_edges_data[key].get('cooccurrence_count', 0) for key in current_subgraph_keys)
    current_sum_term = sum(get_log_likelihood_term(all_edges_data[key], global_node_strengths, denominator) for key in current_subgraph_keys)
    current_ll = math.lgamma(current_T_sub + 1) + current_sum_term

    best_subgraph_keys, min_ll = current_subgraph_keys.copy(), current_ll

    for i in tqdm(range(num_iterations), desc="SA Simulation"):
        on_edge_index = random.randrange(len(current_subgraph_list))
        on_edge_key = current_subgraph_list[on_edge_index]

        off_edge_key = random.choice(all_edge_keys)
        while off_edge_key in current_subgraph_keys:
            off_edge_key = random.choice(all_edge_keys)

        w_on = all_edges_data[on_edge_key].get('cooccurrence_count', 0)
        w_off = all_edges_data[off_edge_key].get('cooccurrence_count', 0)

        proposed_T_sub = current_T_sub - w_on + w_off

        ll_term_on = get_log_likelihood_term(all_edges_data[on_edge_key], global_node_strengths, denominator)
        ll_term_off = get_log_likelihood_term(all_edges_data[off_edge_key], global_node_strengths, denominator)

        delta_ll = (math.lgamma(proposed_T_sub + 1) - math.lgamma(current_T_sub + 1)) + (ll_term_off - ll_term_on)

        if delta_ll < 0 or (temperature > 1e-9 and random.random() < math.exp(-delta_ll / temperature)):
            current_subgraph_keys.remove(on_edge_key)
            current_subgraph_keys.add(off_edge_key)

            current_subgraph_list[on_edge_index] = off_edge_key

            current_ll += delta_ll
            current_T_sub = proposed_T_sub

            if current_ll < min_ll:
                min_ll = current_ll
                best_subgraph_keys = current_subgraph_keys.copy()

        temperature *= cooling_rate
        if i % log_interval == 0:
            history.append((i, current_ll))

    return best_subgraph_keys, min_ll, history

def plot_trajectory(history, title, final_score):
    """Plots the log-likelihood trajectory and saves it to the figures folder."""
    if not history: return
    iterations, scores = zip(*history)

    plt.figure(figsize=(10, 6))
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.plot(iterations, scores, color='C0', alpha=0.7, label='Current Score')
    plt.axhline(y=final_score, color='C3', linestyle='--', label=f'Final Minimum: {final_score:.2f}')
    plt.title(f'{config.FILE_PREFIX}: {title} Optimization Trajectory', fontsize=14)
    plt.xlabel('Simulation Iteration', fontsize=12)
    plt.ylabel('Log-Likelihood Score', fontsize=12)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.legend()
    plt.tight_layout()

    filename = f"True_{config.FILE_PREFIX}_Optimization Trajectory.png"
    save_path = os.path.join(config.FIGURES_DIR, filename)

    plt.savefig(save_path, dpi=300)
    plt.close() # Closes the plot
    print(f"Saved trajectory plot to: {save_path}")

def save_subgraph_to_json(filepath, subgraph_keys, all_nodes_data, all_edges_data):
    """Saves a subgraph to a Cytoscape.js JSON, preserving all original attributes."""
    nodes_in_subgraph = {node for key in subgraph_keys for node in key}
    edges_list = [{'data': all_edges_data[key]} for key in subgraph_keys if key in all_edges_data]
    nodes_list = [{'data': all_nodes_data[node_id]} for node_id in sorted(list(nodes_in_subgraph)) if node_id in all_nodes_data]
    cytoscape_data = {'elements': {'nodes': nodes_list, 'edges': edges_list}}
    with open(filepath, 'w') as f:
        json.dump(cytoscape_data, f, indent=2)
    print(f"Successfully saved graph with {len(nodes_list)} nodes and {len(edges_list)} edges to {filepath}")

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 5. Script Core Functions (Organized by Step)
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# <<< STEP 1: Initial Data Collection >>>
def populate_master_mesh_database(source_pmids_input, master_db_path, failed_log=None, empty_log=None):
    """
    Checks a list of source PMIDs against the master MeSH DB, finds new PMIDs,
    fetches their MeSH terms, and populates the master DB.
    """
    if isinstance(source_pmids_input, str) and os.path.exists(source_pmids_input):
        conn = sqlite3.connect(source_pmids_input)
        cursor = conn.cursor()
        cursor.execute("SELECT pmid FROM pmids_table")
        source_pmids = {row[0] for row in cursor.fetchall()}
        conn.close()
    elif isinstance(source_pmids_input, list):
        source_pmids = set(source_pmids_input)
    else:
        print("Invalid input for populate_master_mesh_database. Expecting DB path or list of PMIDs.")
        return

    print(f"Connecting to master database at: {master_db_path}")
    master_conn = sqlite3.connect(master_db_path)
    m_cursor = master_conn.cursor()
    m_cursor.execute("CREATE TABLE IF NOT EXISTS master_mesh_annotations (pmid INTEGER PRIMARY KEY, mesh_terms TEXT)")
    m_cursor.execute("SELECT pmid FROM master_mesh_annotations")
    master_pmids = {row[0] for row in m_cursor.fetchall()}
    print(f"Master database currently contains {len(master_pmids):,} annotated PMIDs.")

    new_pmids_to_fetch = list(source_pmids - master_pmids)

    if not new_pmids_to_fetch:
        print("No new PMIDs from this analysis to add to the master database.")
        master_conn.close()
        return

    print(f"Found {len(new_pmids_to_fetch):,} new PMIDs to fetch and add to the master database.")

    empty_pmids, failed_pmids = [], []

    for i in tqdm(range(0, len(new_pmids_to_fetch), MESH_FETCH_BATCH_SIZE), desc="Populating Master MeSH DB"):
        batch = new_pmids_to_fetch[i : i + MESH_FETCH_BATCH_SIZE]
        if not batch: continue

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
                    mesh_terms_str = ";".join(mesh_terms_list)
                    all_batch_updates.append((pmid_int, mesh_terms_str))
                    if not mesh_terms_list: empty_pmids.append(pmid_int)

            missing_pmids = set(map(int, batch)) - processed_in_batch
            for pmid in missing_pmids:
                all_batch_updates.append((pmid, ""))
                failed_pmids.append(pmid)

            if all_batch_updates:
                m_cursor.executemany("INSERT OR IGNORE INTO master_mesh_annotations (pmid, mesh_terms) VALUES (?, ?)", all_batch_updates)
                master_conn.commit()
            time.sleep(ENTREZ_DELAY)
        except Exception as e:
            print(f"An error occurred during master DB population for batch starting with {batch[0]}: {e}")
            failed_pmids.extend(map(int, batch))

    master_conn.close()
    # Logging failed attempts
    if failed_log and failed_pmids:
        try:
            with open(failed_log, 'w') as f:
                f.write("PMID\n")
                for pmid in failed_pmids:
                    f.write(f"{pmid}\n")
            print(f"Saved {len(failed_pmids)} failed PMIDs to {failed_log}")
        except Exception as e:
            print(f"Could not write failed PMIDs log: {e}")

    if empty_log and empty_pmids:
        try:
            with open(empty_log, 'w') as f:
                f.write("PMID\n")
                for pmid in empty_pmids:
                    f.write(f"{pmid}\n")
            print(f"Saved {len(empty_pmids)} empty MeSH PMIDs to {empty_log}")
        except Exception as e:
            print(f"Could not write empty PMIDs log: {e}")
    print("Master database population complete.")

def _fetch_ids_for_range(search_term, start_date_obj, end_date_obj, retmax_chunk):
    """
    Performs a single esearch and returns PMIDs, using the retry mechanism.
    """
    pmids = set()
    full_query = f"({search_term}) AND ({format_date_for_entrez(start_date_obj)}[EDAT] : {format_date_for_entrez(end_date_obj)}[EDAT])"

    search_results = _entrez_search_with_retry(db="pubmed", term=full_query, retmax=retmax_chunk, sort="relevance")
    id_list = search_results.get("IdList", [])
    pmids.update(id_list)
    time.sleep(ENTREZ_DELAY)

    return pmids

def get_pmids_date_chunking(search_term, start_date_obj, end_date_obj, retmax_limit=9999):
    """
    Fetches PMIDs by splitting large date ranges and using a robust retry mechanism.
    """
    all_pmids = set()
    date_ranges_to_process = deque([(start_date_obj, end_date_obj)])
    pbar = tqdm(total=1, desc=f"Processing Date Ranges for '{search_term}'")

    while date_ranges_to_process:
        current_start, current_end = date_ranges_to_process.popleft()
        pbar.set_postfix_str(f"{current_start} to {current_end}", refresh=True)

        date_query = f"({format_date_for_entrez(current_start)}[EDAT] : {format_date_for_entrez(current_end)}[EDAT])"
        full_query = f"({search_term}) AND {date_query}"

        search_results = _entrez_search_with_retry(db="pubmed", term=full_query, retmax="0")
        count = int(search_results["Count"])
        time.sleep(ENTREZ_DELAY)

        if count <= retmax_limit:
            retrieved_ids = _fetch_ids_for_range(search_term, current_start, current_end, retmax_limit)
            all_pmids.update(retrieved_ids)
            pbar.update(1)
        else:
            if current_start == current_end:
                date_str = current_start.strftime('%Y-%m-%d')
                print(f"  Warning: Count for single day {date_str} ({count}) exceeds limit. Fetching only the first {retmax_limit}.")
                retrieved_ids = _fetch_ids_for_range(search_term, current_start, current_end, retmax_limit)
                all_pmids.update(retrieved_ids)
                pbar.update(1)
                continue

            mid_point_ord = (current_start.toordinal() + current_end.toordinal()) // 2
            mid_date = date.fromordinal(mid_point_ord)
            date_ranges_to_process.appendleft((mid_date + timedelta(days=1), current_end))
            date_ranges_to_process.appendleft((current_start, mid_date))
            pbar.total += 1

    pbar.close()
    print(f" -> Finished processing. Found {len(all_pmids):,} unique PMIDs for '{search_term}'.")
    return sorted(list(all_pmids))

def fetch_links_in_batches(parent_pmid_batch):
    """
    Fetches citation data for a batch of PMIDs using the iCite API
    with the correct GET method and robust handling for None values.
    """
    citation_results = {pmid: {'cited_by': "", 'cites': ""} for pmid in parent_pmid_batch}
    if not parent_pmid_batch:
        return citation_results

    api_url = "https://icite.od.nih.gov/api/pubs"

    params = {
        'pmids': ",".join(map(str, parent_pmid_batch))
    }

    try:
        response = requests.get(api_url, params=params)
        response.raise_for_status()

        data = response.json().get('data', [])

        for article_data in data:
            pmid = article_data.get('pmid')
            if pmid in citation_results:
                cited_by_data = article_data.get('cited_by') or []
                references_data = article_data.get('references') or []

                cited_by_list = [str(p) for p in cited_by_data]
                cites_list = [str(p) for p in references_data]

                citation_results[pmid]['cited_by'] = ";".join(cited_by_list)
                citation_results[pmid]['cites'] = ";".join(cites_list)

        return citation_results

    except requests.exceptions.RequestException as e:
        print(f"  ERROR: Failed to fetch data from iCite API: {e}")
        return citation_results

def run_initial_data_collection(search_term_param, start_date_str, end_date_str, generations_n_param, db_path):
    """
    Performs initial PMID search and generational expansion, storing results in an SQLite database.
    """
    script_start_time = time.time()
    all_pmids_gen = {}
    conn = None
    generations_to_process = max(1, int(generations_n_param))

    try:
        start_date_obj = parse_date_robust(start_date_str or "1900/01/01")
        end_date_obj = parse_date_robust(end_date_str or date.today().strftime("%Y/%m/%d"))
        if start_date_obj > end_date_obj: raise ValueError("Start date cannot be after end date.")

        initial_pmids_str_list = get_pmids_date_chunking(search_term_param, start_date_obj, end_date_obj, retmax_limit=RETMAX_CHUNK)
        if not initial_pmids_str_list:
            print("No initial PMIDs found. Exiting.")
            return

        all_pmids_gen = {int(pmid): 0 for pmid in initial_pmids_str_list}
        print(f"Found {len(all_pmids_gen)} unique initial PMIDs (P0).")

        print(f"\n<<< Setting up Database: {db_path} >>>")
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

        # <<< Generational Expansion Loop >>>
        for current_gen_num in range(generations_to_process - 1):
            parent_gen_label = get_generation_label(current_gen_num)
            child_gen_label = get_generation_label(current_gen_num + 1)
            print(f"\n<<< Processing {parent_gen_label} PMIDs to find {child_gen_label} children >>>")

            cursor.execute("SELECT pmid FROM pmids_table WHERE generation = ?", (parent_gen_label,))
            parent_pmids_in_gen = [row[0] for row in cursor.fetchall()]
            if not parent_pmids_in_gen:
                print(f"No PMIDs for {parent_gen_label} in DB. Stopping expansion.")
                break

            for i in tqdm(range(0, len(parent_pmids_in_gen), BATCH_SIZE), desc=f"Finding children of {parent_gen_label}"):
                parent_batch = parent_pmids_in_gen[i:i + BATCH_SIZE]

                # A single API call now gets all citation data for the batch. # modified
                citation_results = fetch_links_in_batches(parent_batch)
                time.sleep(ENTREZ_DELAY)

                update_link_data = []
                insert_new_pmid_data = []

                for pmid in parent_batch:
                    # Unpack the results for each pmid from the returned dictionary. # modified
                    links = citation_results.get(pmid, {'cited_by': "", 'cites': ""})
                    cited_by_str = links['cited_by']
                    cites_str = links['cites']
                    update_link_data.append((cited_by_str, cites_str, pmid))

                    child_candidates = set((cited_by_str + ";" + cites_str).split(';'))
                    for child_pmid_str in filter(None, child_candidates):
                        try:
                            child_pmid = int(child_pmid_str)
                            if child_pmid not in all_pmids_gen:
                                all_pmids_gen[child_pmid] = current_gen_num + 1
                                insert_new_pmid_data.append((child_pmid, child_gen_label))
                        except ValueError: continue

                if update_link_data:
                    cursor.executemany("UPDATE pmids_table SET cited_by = ?, cites = ? WHERE pmid = ?", update_link_data)
                if insert_new_pmid_data:
                    cursor.executemany("INSERT OR IGNORE INTO pmids_table (pmid, generation) VALUES (?, ?)", insert_new_pmid_data)
                conn.commit()

        # <<< Final Update Pass Loop >>>
        if generations_to_process > 0:
            final_gen_label = get_generation_label(generations_to_process - 1)
            print(f"\n<<< Fetching links FOR final generation ({final_gen_label}) >>>")
            cursor.execute("SELECT pmid FROM pmids_table WHERE generation = ? AND (cited_by IS NULL OR cites IS NULL)", (final_gen_label,))
            final_gen_pmids = [row[0] for row in cursor.fetchall()]

            for i in tqdm(range(0, len(final_gen_pmids), BATCH_SIZE), desc=f"Fetching Links for {final_gen_label}"):
                final_batch = final_gen_pmids[i:i + BATCH_SIZE]
                citation_results = fetch_links_in_batches(final_batch)
                time.sleep(ENTREZ_DELAY)

                update_link_data = []
                for pmid in final_batch:
                    links = citation_results.get(pmid, {'cited_by': "", 'cites': ""})
                    cited_by_str = links['cited_by']
                    cites_str = links['cites']
                    update_link_data.append((cited_by_str, cites_str, pmid))

                if update_link_data:
                    cursor.executemany("UPDATE pmids_table SET cited_by = ?, cites = ? WHERE pmid = ?", update_link_data)
                conn.commit()

    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")
        logging.error("An unexpected error occurred in main processing", exc_info=True)
    finally:
        if conn:
            conn.close()
            print("\nDatabase connection closed.")

    script_end_time = time.time()
    print(f"\nTotal script execution time: {script_end_time - script_start_time:.2f} seconds")
    print(f"Total unique PMIDs found across {generations_to_process} generation(s): {len(all_pmids_gen)}")

# <<< STEP 2: Database Cleaning >>>
def clean_database(original_db, cleaned_db, start_gen_to_remove_label):
    """
    Creates a cleaned copy of the database, removing specified generations and higher.
    """
    print(f"\n<<< Cleaning Database >>>")
    if not os.path.exists(original_db):
        print(f"Original database '{original_db}' not found. Skipping cleaning process.")
        return

    print(f"Original: {original_db}, Cleaned: {cleaned_db}")
    print(f"Removing generations at and after: {start_gen_to_remove_label}")

    shutil.copyfile(original_db, cleaned_db)
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
                    gen_num = int(gen_label.lstrip('G'))
                    if gen_num >= start_num_to_remove:
                        gens_to_delete.append(gen_label)
                except ValueError: continue

        if gens_to_delete:
            placeholders = ', '.join('?' for _ in gens_to_delete)
            sql_query = f"DELETE FROM pmids_table WHERE generation IN ({placeholders})"
            cursor.execute(sql_query, gens_to_delete)
            print(f"Successfully removed {cursor.rowcount} records from generations: {gens_to_delete}")
        else:
            print("No generations met the criteria for removal.")

    except (ValueError, IndexError):
        print(f"Error: Invalid format for `start_gen_to_remove`: '{start_gen_to_remove_label}'. Should be e.g., 'G2'.")
    except sqlite3.Error as e:
        print(f"An error occurred during database cleaning: {e}")
    finally:
        conn.commit()
        conn.close()
        print("Database cleaning complete.")

# <<< STEP 3: MeSH Term Annotation >>>
def update_db_with_mesh_batch(local_db_path, master_db_path):
    """
    Updates the local project DB with MeSH terms from the master annotation DB
    using attached SQL and updates accordingly.
    """
    print(f"Starting MeSH term update for local database: {local_db_path}")
    if not os.path.exists(master_db_path):
        print(f"ERROR: Master database not found at {master_db_path}. Cannot annotate local database.")
        return
    if not os.path.exists(local_db_path):
        print(f"ERROR: Local database not found at {local_db_path}.")
        return

    conn = None
    try:
        conn = sqlite3.connect(local_db_path)
        cursor = conn.cursor()
        try:
            cursor.execute("ALTER TABLE pmids_table ADD COLUMN mesh_terms TEXT")
            print("Added 'mesh_terms' column to local DB.")
        except sqlite3.OperationalError:
            pass

        print("Attaching master database for update...")
        cursor.execute("ATTACH DATABASE ? AS master_db", (master_db_path,))

        print("Performing update using SQL.")
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

        print("Detaching master database...")
        cursor.execute("DETACH DATABASE master_db")
        conn.commit()

    except sqlite3.Error as e:
        print(f"An SQL error occurred: {e}")
    finally:
        if conn:
            conn.close()

    print("Local DB annotation complete.")

# <<< STEP 4: Network Construction >>>
def normalize_cited_by_per_mesh(df, parse_mesh_terms_func, stop_words_set):
    """
    Calculates article citation ranks and MeSH term citation statistics using
    an optimized, vectorized approach.
    """
    print("  Calculating numeric citation counts per article...")
    # Convert 'cited_by' to a numeric count of citations.
    df['cited_by_count'] = df['cited_by'].apply(
        lambda x: len(x.split(';')) if isinstance(x, str) and x.strip() and x.lower() != 'nan' else 0
    ).astype(int)

    print("  Vectorizing MeSH terms for aggregation...")
    # Create a new DataFrame where each row corresponds to a single MeSH term from an article.
    # This is the key optimization step, replacing the slow iterrows() loop.

    # 1. Split the 'mesh_terms' string into a list of terms, handling potential None/NaN values.
    # We only care about main terms, so we'll filter out subheadings and clean the term name.
    def extract_main_mesh_terms(mesh_string):
        if not isinstance(mesh_string, str):
            return []
        return [
            term.strip().lstrip('*') for term in mesh_string.split(';')
            if '/' not in term and term.strip()
        ]

    df_exploded = df[['mesh_terms', 'cited_by_count']].copy()
    df_exploded['mesh_terms'] = df_exploded['mesh_terms'].apply(extract_main_mesh_terms)

    # 2. Explode the DataFrame so each MeSH term has its own row.
    df_exploded = df_exploded.explode('mesh_terms')
    df_exploded.rename(columns={'mesh_terms': 'mesh_term'}, inplace=True)

    # 3. Filter out empty terms and stop words in a single, vectorized step.
    df_exploded.dropna(subset=['mesh_term'], inplace=True)
    df_exploded = df_exploded[~df_exploded['mesh_term'].isin(stop_words_set)]

    print("  Aggregating citation counts per term using groupby...")
    if not df_exploded.empty:
        # 4. Use groupby().agg() for highly efficient, parallelizable calculation.
        mesh_stats_df = df_exploded.groupby('mesh_term')['cited_by_count'].agg(
            total_pmids='size',
            median_citations='median',
            mean_citations='mean',
            total_citations='sum'
        ).reset_index()
    else:
        print("  Warning: No valid MeSH terms found for statistics calculation.")
        mesh_stats_df = pd.DataFrame(columns=['mesh_term', 'total_pmids', 'median_citations', 'mean_citations', 'total_citations'])

    print("  Rank-normalizing MeSH term statistics...")
    num_mesh_terms = len(mesh_stats_df)
    stats_to_rank = ['total_pmids', 'median_citations', 'mean_citations', 'total_citations']
    if num_mesh_terms > 0:
        for col_name in stats_to_rank:
            rank_col_name = f"rank_{col_name.replace('citations','cit').replace('_','')}"
            values_to_rank = pd.to_numeric(mesh_stats_df[col_name], errors='coerce').fillna(0)
            ranks = stats.rankdata(values_to_rank.values, method='average')
            mesh_stats_df[rank_col_name] = ranks / num_mesh_terms
    else:
        for col_name in stats_to_rank:
            rank_col_name = f"rank_{col_name.replace('citations','cit').replace('_','')}"
            mesh_stats_df[rank_col_name] = 0.0
    print("  Rank-normalizing article-level citation counts...")
    if 'cited_by_count' in df.columns and pd.api.types.is_numeric_dtype(df['cited_by_count']):
        cited_counts_array = df['cited_by_count'].fillna(0).values
        num_articles = len(cited_counts_array)
        if num_articles > 1:
            ranks = stats.rankdata(cited_counts_array, method='average')
            df['normalized_cited_by'] = ranks / num_articles
        elif num_articles == 1:
            df['normalized_cited_by'] = 0.5 if cited_counts_array[0] > 1e-9 else 0.0
        else:
            df['normalized_cited_by'] = 0.0
    else:
        df['normalized_cited_by'] = 0.0

    print("Finished calculating statistics and weights.")
    return df, mesh_stats_df

def run_network_construction(db_path_param, output_json_path):
    print(f"Connecting to database: {db_path_param}")
    required_columns = ['mesh_terms', 'cited_by', 'generation']
    try:
        conn = sqlite3.connect(db_path_param)
        cursor = conn.cursor()
        cursor.execute("PRAGMA table_info(pmids_table)")
        available_columns = [info[1] for info in cursor.fetchall()]
        if not set(required_columns).issubset(available_columns):
            raise ValueError(f"Required columns missing from pmids_table: {set(required_columns) - set(available_columns)}")
        query = f"SELECT {', '.join(required_columns)} FROM pmids_table"
        df = pd.read_sql_query(query, conn)
        conn.close()
        print(f"Successfully loaded {len(df)} records.")
    except Exception as e:
        print(f"Data loading error: {e}")
        return

    if df.empty:
        print("Loaded DataFrame is empty. Exiting network construction.")
        return

    print("\n<<< Pre-computation >>>")
    df['generation_weight'] = df['generation'].apply(calculate_generation_weight, lambda_val=LAMBDA_VALUE)
    df, mesh_stats_df = normalize_cited_by_per_mesh(df, parse_mesh_terms, Mesh_stop_words)

    rank_norm_weight_dicts = {}
    if not mesh_stats_df.empty:
        rank_norm_weight_dicts['median_cit'] = pd.Series(mesh_stats_df.get('rank_mediancit', pd.Series(dtype='float64')).values, index=mesh_stats_df['mesh_term']).to_dict()
        rank_norm_weight_dicts['mean_cit'] = pd.Series(mesh_stats_df.get('rank_meancit', pd.Series(dtype='float64')).values, index=mesh_stats_df['mesh_term']).to_dict()
        rank_norm_weight_dicts['total_cit'] = pd.Series(mesh_stats_df.get('rank_totalcit', pd.Series(dtype='float64')).values, index=mesh_stats_df['mesh_term']).to_dict()

    print("\n<<< Initializing Network Data Structures >>>")
    node_data_agg = defaultdict(lambda: {'generation_weight': 0.0, 'term_type': 'Unknown', 'is_major_topic': False, 'generation': 'Unknown'})
    edges = defaultdict(lambda: {'edge_weight': 0.0, 'cooccurrence_count': 0})
    mesh_term_article_counts = defaultdict(int)
    articles_with_term_as_major_count = defaultdict(int)
    min_generation_for_term = defaultdict(lambda: float('inf'))

    print("\n<<< Processing Articles to Build Network Structure >>>")
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Processing Articles"):
        article_mesh_terms_str = row.get('mesh_terms', '')
        if not isinstance(article_mesh_terms_str, str):
            continue

        article_centrality = row.get('generation_weight', 0.0)
        article_norm_cited_by = row.get('normalized_cited_by', 0.0)
        article_generation_raw = row.get('generation', None)

        try:
            if isinstance(article_generation_raw, str):
                gen_str = article_generation_raw.upper()
                article_generation_num = 0 if gen_str == 'P0' else int(gen_str[1:]) if gen_str.startswith('G') else float('inf')
            else:
                article_generation_num = int(article_generation_raw) if pd.notna(article_generation_raw) else float('inf')
        except (ValueError, TypeError):
            article_generation_num = float('inf')

        parsed_terms = parse_mesh_terms(article_mesh_terms_str)

        main_terms_in_article = []
        for term_name, term_type, is_major_topic in parsed_terms:
            if term_type == 'mesh_term' and term_name and term_name not in Mesh_stop_words:
                main_terms_in_article.append({'name': term_name, 'is_major': is_major_topic})

        if not main_terms_in_article:
            continue

        unique_valid_terms_for_edges = {term['name'] for term in main_terms_in_article}
        major_terms_in_this_article = {term['name'] for term in main_terms_in_article if term['is_major']}
        for term_name in unique_valid_terms_for_edges:
            node_entry = node_data_agg[term_name]
            mesh_term_article_counts[term_name] += 1
            if term_name in major_terms_in_this_article:
                articles_with_term_as_major_count[term_name] += 1
                node_entry['is_major_topic'] = True
            if '_initialized' not in node_entry:
                node_entry['rank_norm_median_cit'] = rank_norm_weight_dicts.get('median_cit', {}).get(term_name, 0.0)
                node_entry['rank_norm_mean_cit'] = rank_norm_weight_dicts.get('mean_cit', {}).get(term_name, 0.0)
                node_entry['rank_norm_total_cit'] = rank_norm_weight_dicts.get('total_cit', {}).get(term_name, 0.0)
                node_entry['_initialized'] = True

            if pd.notna(article_centrality) and article_centrality > node_entry['generation_weight']:
                node_entry['generation_weight'] = article_centrality

            node_entry['term_type'] = 'mesh_term'

            if article_generation_num != float('inf'):
                min_generation_for_term[term_name] = min(min_generation_for_term[term_name], article_generation_num)

        if len(unique_valid_terms_for_edges) > 1:
            combined_article_weight = (article_norm_cited_by + article_centrality) / 2.0 if pd.notna(article_norm_cited_by) and pd.notna(article_centrality) else 0.0
            sorted_terms_list = sorted(list(unique_valid_terms_for_edges))
            for i in range(len(sorted_terms_list)):
                for j in range(i + 1, len(sorted_terms_list)):
                    edge_key = (sorted_terms_list[i], sorted_terms_list[j])
                    edges[edge_key]['edge_weight'] += combined_article_weight
                    edges[edge_key]['cooccurrence_count'] += 1

    for term_data in node_data_agg.values():
        term_data.pop('_initialized', None)

    print("\n<<< Calculating MLF P-Values for All Edges >>>")
    T = sum(data['cooccurrence_count'] for data in edges.values())
    node_strengths = defaultdict(int)
    for key, data in edges.items():
        source, target = key
        node_strengths[source] += data['cooccurrence_count']
        node_strengths[target] += data['cooccurrence_count']

    if T > 0:
        denominator = 2 * (T**2)
        for key, data in tqdm(edges.items(), desc="Calculating p-values"):
            source, target = key
            k_i = node_strengths[source]
            k_j = node_strengths[target]
            w_ij = data['cooccurrence_count']
            if w_ij > 0:
                p_ij = (k_i * k_j) / denominator
                p_value = binom.sf(w_ij - 1, n=T, p=p_ij)
                data['mlf_p_value'] = p_value
            else:
                data['mlf_p_value'] = 1.0

    print("\n<<< Post-computation and Finalization >>>")
    edge_keys = list(edges.keys())
    if edge_keys:
        for edge_key in edge_keys:
            source, target = edge_key
            cooccurrence = edges[edge_key]['cooccurrence_count']
            union = mesh_term_article_counts.get(source, 0) + mesh_term_article_counts.get(target, 0) - cooccurrence
            edges[edge_key]['normalized_cooccurrence'] = cooccurrence / union if union > 0 else 0.0

        cumulative_weights = np.array([edges[key]['edge_weight'] for key in edge_keys])

        num_edges = len(cumulative_weights)
        edge_ranks_norm = np.zeros(num_edges)
        if num_edges > 1:
            unique_weights = np.unique(cumulative_weights)
            if len(unique_weights) == 1:
                edge_ranks_norm.fill(0.5)
            else:
                ranks = stats.rankdata(cumulative_weights, method='average')
                edge_ranks_norm = ranks / num_edges
        elif num_edges == 1:
            edge_ranks_norm[0] = 0.5

        for i, key in enumerate(edge_keys):
            edges[key]['rank_normalized_weight'] = edge_ranks_norm[i]

        cooccurrence_counts = np.array([edges[key]['cooccurrence_count'] for key in edge_keys]).astype(float)
        log1p_counts = np.log1p(cooccurrence_counts)
        min_log, max_log = (np.min(log1p_counts), np.max(log1p_counts)) if len(log1p_counts) > 0 else (0,0)
        scaled_logs = (log1p_counts - min_log) / (max_log - min_log) if (max_log - min_log) > 1e-9 else np.full(log1p_counts.shape, 0.5)
        for i, key in enumerate(edge_keys):
            edges[key]['log1p_cooccurrence_minmax'] = scaled_logs[i]

    node_keys = list(node_data_agg.keys())
    if node_keys:
        for term_name in node_keys:
            node_data_agg[term_name]['generation'] = get_generation_label(min_generation_for_term.get(term_name, float('inf')))
        article_counts = np.array([mesh_term_article_counts.get(term, 0) for term in node_keys]).astype(float)
        if len(article_counts) > 1: ranks = stats.rankdata(article_counts, method='average') / len(article_counts)
        elif len(article_counts) == 1: ranks = [0.5]
        else: ranks = []
        for i, term_name in enumerate(node_keys):
            node_data_agg[term_name]['article_count_rank_normalized'] = ranks[i] if i < len(ranks) else 0.0
        total_factor = sum(NODE_WEIGHT_FACTORS.values())
        norm_factors = {k: v / total_factor if total_factor > 1e-9 else 0 for k, v in NODE_WEIGHT_FACTORS.items()}
        for term_name in node_keys:
            data = node_data_agg[term_name]
            article_count = mesh_term_article_counts.get(term_name, 0)
            data['article_count'] = article_count
            data['major_topic_proportion'] = articles_with_term_as_major_count.get(term_name, 0) / article_count if article_count > 0 else 0.0
            data['adjusted_node_weight'] = (data.get('generation_weight', 0.0) * norm_factors['centrality'] + data.get('article_count_rank_normalized', 0.0) * norm_factors['article_rank'] + data.get('rank_norm_median_cit', 0.0) * norm_factors['rank_median_cit'] + data.get('rank_norm_total_cit', 0.0) * norm_factors['rank_total_cit'])
    print("\n<<< Calculating Network Centrality Measures & Communities >>>")

    G = nx.Graph()
    if edges: G.add_edges_from(edges.keys())
    run_full_centrality = getattr(config, 'CALCULATE_FULL_CENTRALITY', True) # User choice in config file

    if run_full_centrality:
        # <<< OPTION TRUE: CALCULATIONS >>>
        print("\nThis may take awhile... Do NOT interrupt")
        print("This will not corrupt files, but calcs stored in RAM will disappear")

        G_analysis_conn = G
        if len(G) > 0 and not nx.is_connected(G):
            print("Graph is not connected. Using largest connected component for some metrics.")
            largest_cc = max(nx.connected_components(G), key=len)
            G_analysis_conn = G.subgraph(largest_cc).copy()

        if G_analysis_conn.number_of_nodes() > 0:
            degree_dict = dict(G.degree())
            betweenness_centrality = nx.betweenness_centrality(G_analysis_conn, k=BETWEENNESS_K_SAMPLES, normalized=True, seed=RANDOM_SEED)
            eigenvector_centrality = nx.eigenvector_centrality(G_analysis_conn, max_iter=500, tol=1.0e-06)
            clustering_coefficient = nx.clustering(G)
            edge_betweenness_centrality = nx.edge_betweenness_centrality(G_analysis_conn, k=BETWEENNESS_K_SAMPLES, normalized=True, seed=RANDOM_SEED)
            partition_map = community_louvain.best_partition(G_analysis_conn, random_state=RANDOM_SEED)
        else:
            degree_dict, betweenness_centrality, eigenvector_centrality, clustering_coefficient, edge_betweenness_centrality, partition_map = {}, {}, {}, {}, {}, {}

    else:
        # <<< OPTION FALSE: SKIPPED CALCULATIONS >>>
        print("  -> WARNING: Full centrality DISABLED.")
        print("  -> ACTION: Assigning 1.0 to LCC (Mainland), 0.0 to Disconnected (Islands).")

        all_nodes = list(G.nodes())
        all_edges = list(G.edges())

        # Initialize everyone to 0.0 (Island status)
        betweenness_centrality = {n: 0.0 for n in all_nodes}
        eigenvector_centrality = {n: 0.0 for n in all_nodes}
        clustering_coefficient = {n: 0.0 for n in all_nodes}
        partition_map = {n: -1 for n in all_nodes}
        edge_betweenness_centrality = {e: 0.0 for e in all_edges}

        degree_dict = dict(G.degree()) if len(G) > 0 else {}

        if len(G) > 0:
            # Identify LCC
            if nx.is_connected(G):
                lcc_nodes = all_nodes
            else:
                lcc_nodes = max(nx.connected_components(G), key=len)

            # Set Mainland to 1.0
            for n in lcc_nodes:
                betweenness_centrality[n] = 1.0
                eigenvector_centrality[n] = 1.0

            lcc_set = set(lcc_nodes)
            for u, v in all_edges:
                if u in lcc_set and v in lcc_set:
                    edge_betweenness_centrality[(u, v)] = 1.0

    # <<< MERGE LOGIC >>>
    # This matches your original code exactly and ensures data is saved for BOTH cases.
    for node_id, data in node_data_agg.items():
        data['degree'] = degree_dict.get(node_id, 0)
        data['betweenness_centrality'] = betweenness_centrality.get(node_id, 0.0)
        data['eigenvector_centrality'] = eigenvector_centrality.get(node_id, 0.0)
        data['clustering_coefficient'] = clustering_coefficient.get(node_id, 0.0)
        data['unfiltered_louvain_community_id'] = partition_map.get(node_id, -1)

    for edge_key, betweenness in edge_betweenness_centrality.items():
        sorted_key = tuple(sorted(edge_key))
        if sorted_key in edges:
            edges[sorted_key]['edge_betweenness'] = betweenness

    print("\n<<< Preparing Data for Cytoscape JSON Output >>>")
    nodes_list_for_cytoscape = []
    for term_name, data in node_data_agg.items():
        node_output_data = {
            'id': term_name, 'label': term_name, 'generation': data.get('generation', 'Unknown'),
            'term_type': data.get('term_type', 'Unknown'), 'is_major_topic': bool(data.get('is_major_topic', False)),
            'generation_weight': float(data.get('generation_weight', 0.0)), 'article_count': int(data.get('article_count', 0)),
            'article_count_rank_normalized': float(data.get('article_count_rank_normalized', 0.0)),
            'major_topic_proportion': float(data.get('major_topic_proportion', 0.0)),
            'rank_norm_median_cit': float(data.get('rank_norm_median_cit', 0.0)),
            'rank_norm_mean_cit': float(data.get('rank_norm_mean_cit', 0.0)),
            'rank_norm_total_cit': float(data.get('rank_norm_total_cit', 0.0)),
            'adjusted_node_weight': float(data.get('adjusted_node_weight', 0.0)), 'degree': int(data.get('degree', 0)),
            'betweenness_centrality': float(data.get('betweenness_centrality', 0.0)),
            'eigenvector_centrality': float(data.get('eigenvector_centrality', 0.0)),
            'clustering_coefficient': float(data.get('clustering_coefficient', 0.0)),
            'unfiltered_louvain_community_id': int(data.get('unfiltered_louvain_community_id', -1))
        }
        nodes_list_for_cytoscape.append({'data': node_output_data})
    edges_list_for_cytoscape = []
    for key, data_dict in edges.items():
        source, target = key
        edge_output_data = {
            'source': source, 'target': target, 'weight': float(data_dict.get('rank_normalized_weight', 0.0)),
            'cumulative_weight_original': float(data_dict.get('edge_weight', 0.0)),
            'cooccurrence_count': int(data_dict.get('cooccurrence_count', 0)),
            'normalized_cooccurrence': float(data_dict.get('normalized_cooccurrence', 0.0)),
            'log1p_cooccurrence_minmax': float(data_dict.get('log1p_cooccurrence_minmax', 0.0)),
            'edge_betweenness': float(data_dict.get('edge_betweenness', 0.0)),
            'mlf_p_value': float(data_dict.get('mlf_p_value', 1.0))
        }
        edges_list_for_cytoscape.append({'data': edge_output_data})

    cytoscape_data = {'elements': {'nodes': nodes_list_for_cytoscape, 'edges': edges_list_for_cytoscape}}

    print("\n<<< Full Network Summary >>>")
    print(f"Total Nodes in Full Graph: {len(nodes_list_for_cytoscape):,}")
    print(f"Total Edges in Full Graph: {len(edges_list_for_cytoscape):,}")

    temp_output_path = output_json_path + ".tmp"
    try:
        print(f"Writing network data to temporary file: {temp_output_path}...")
        with open(temp_output_path, 'w', encoding='utf-8') as f:
            json.dump(cytoscape_data, f, indent=2)
        os.rename(temp_output_path, output_json_path)
        print(f"Successfully saved network data to {output_json_path}")
    except Exception as e:
        print(f"ERROR: Failed to write final JSON file: {e}")
        if os.path.exists(temp_output_path):
            os.remove(temp_output_path)

    print("Network construction complete.")

# <<< STEP 5: Consensus Network Filtering and LCC Extraction >>>
def run_consensus_filtering_and_lcc(input_json_path, glf_output_path, sa_output_path, final_lcc_output_path):
    """
    Performs network filtering using a consensus of GLF and SA optimizers,
    then extracts and saves the Largest Connected Component (LCC).
    """
    print(f"\n<<< Running Consensus Network Filtering >>>")

    all_nodes_data, all_edges_data = load_full_graph_data(input_json_path)
    if not all_edges_data:
        print("No edges in the unfiltered graph. Aborting filtering.")
        return

    global_T, global_node_strengths = calculate_graph_stats(all_edges_data)

    # <<< Run Optimizations >>>
    glf_keys, glf_score, glf_history = run_glf_simulation(all_edges_data, global_node_strengths, global_T, TARGET_NUM_EDGES, GLF_ITERATIONS)
    sa_keys, sa_score, sa_history = run_sa_simulation(all_edges_data, global_node_strengths, global_T, TARGET_NUM_EDGES, SA_ITERATIONS, SA_INITIAL_TEMPERATURE, SA_COOLING_RATE)

    # <<< Display Plots and Save Individual Optimizer Graphs >>>
    plot_trajectory(glf_history, "GLF", glf_score)
    plot_trajectory(sa_history, "Simulated Annealing", sa_score)
    save_subgraph_to_json(glf_output_path, glf_keys, all_nodes_data, all_edges_data)
    save_subgraph_to_json(sa_output_path, sa_keys, all_nodes_data, all_edges_data)

    # <<< Build Temporary Consensus Network >>>
    print("\n<<< Building Consensus Network >>>")
    consensus_keys = glf_keys.intersection(sa_keys)
    print(f"Found {len(consensus_keys)} consensus edges (intersection of GLF and SA).")
    # FYI: The full consensus graph (with all subgraphs) is NOT saved to a separate file.

    # <<< Extract LCC from Consensus Network >>>
    print("\n<<< Extracting Largest Connected Component (LCC) from Consensus Network >>>")
    if not consensus_keys:
        print("Consensus network is empty. No LCC to extract.")
        with open(final_lcc_output_path, 'w') as f: json.dump({"elements": {"nodes": [], "edges": []}}, f)
        return

    G_consensus = nx.Graph()
    G_consensus.add_edges_from(list(consensus_keys))

    components = list(nx.connected_components(G_consensus))
    if not components:
        print("No connected components found in consensus network.")
        with open(final_lcc_output_path, 'w') as f: json.dump({"elements": {"nodes": [], "edges": []}}, f)
        return

    print(f"Found {len(components)} separate component(s).")
    sorted_components = sorted(components, key=len, reverse=True)
    lcc_nodes = sorted_components[0]

    smaller_components = sorted_components[1:]
    if smaller_components:
        print("\n<<< Top 5 Disconnected Networks (after LCC removal) >>>")
        for i, component_nodes in enumerate(smaller_components[:5]):
            if not component_nodes: continue
            central_node = max(component_nodes, key=lambda node_id: all_nodes_data.get(node_id, {}).get('betweenness_centrality', 0.0))
            centrality_val = all_nodes_data.get(central_node, {}).get('betweenness_centrality', 0.0)
            print(f"  #{i+1}: Size = {len(component_nodes)} nodes. Central Node = '{central_node}' (Betweenness: {centrality_val:.4f})")

    lcc_edge_keys = {key for key in consensus_keys if key[0] in lcc_nodes and key[1] in lcc_nodes}
    save_subgraph_to_json(final_lcc_output_path, lcc_edge_keys, all_nodes_data, all_edges_data)
    print("\nConsensus filtering and LCC extraction complete.")

# <<< STEP 6: Community Detection on Filtered Network >>>
def run_community_detection(network_file_path):
    """
    Loads a filtered network, runs Louvain, re-orders community IDs by size,
    adds the new IDs, and SAVES BACK to the same file.
    """
    print(f"Loading network for community detection from: {network_file_path}")
    if not os.path.exists(network_file_path):
        print(f"Error: Input file not found at '{network_file_path}'")
        return

    with open(network_file_path, 'r', encoding='utf-8') as f:
        network_data = json.load(f)

    nodes_json = network_data.get('elements', {}).get('nodes', [])
    edges_json = network_data.get('elements', {}).get('edges', [])

    G = nx.Graph()
    for node_obj in nodes_json:
        G.add_node(node_obj['data']['id'])
    for edge_obj in edges_json:
        G.add_edge(edge_obj['data']['source'], edge_obj['data']['target'])

    print(f"Built graph for Louvain with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")

    if G.number_of_nodes() > 0:
        print("Performing Louvain community detection...")
        partition_map = community_louvain.best_partition(G, random_state=RANDOM_SEED)
        print(f"Found {len(set(partition_map.values()))} communities.")

        # Reorder community IDs (largest to smallest)
        community_nodes = defaultdict(list)
        for node, community_id in partition_map.items():
            community_nodes[community_id].append(node)
        community_sizes = {cid: len(nodes) for cid, nodes in community_nodes.items()}
        sorted_communities = sorted(community_sizes.items(), key=lambda item: item[1], reverse=True)
        ranked_id_map = {old_id: new_rank for new_rank, (old_id, size) in enumerate(sorted_communities, 1)}
        for node_obj in nodes_json:
            node_id = node_obj['data']['id']
            original_community_id = partition_map.get(node_id, -1)

            if original_community_id != -1:
                ranked_community_id = ranked_id_map.get(original_community_id)
                node_obj['data']['filtered_louvain_community_id'] = ranked_community_id
            else:
                node_obj['data']['filtered_louvain_community_id'] = -1 # Assign -1 to nodes not in a community

    final_data = {"elements": {"nodes": nodes_json, "edges": edges_json}}

    print(f"Updating network with re-ordered community IDs and saving back to: {network_file_path}")
    with open(network_file_path, 'w', encoding='utf-8') as f:
        json.dump(final_data, f, indent=2)
    print("Community detection complete.")

# <<< STEP 7: Contextual Relevance Scoring >>>
def run_contextual_relevance_scoring(input_nodes_file, output_nodes_file, master_db_path, id_key,
                                     weight_key_1, final_key_1, weight_key_2, final_key_2,
                                     start_date_param, end_date_param, max_pmids):
    """
    Calculates contextual relevance by building a complete and resumable
    context pool of PMIDs before scoring.
    """
    if not getattr(config, 'CALCULATE_FULL_CENTRALITY', True):
        print("\n[!] NOTICE: User skipped centrality in config.py. CRS scores now reflect TOPIC DENSITY (term overlap), not topological relevance.\n")
    print(f"\n<<< Loading Seed Terms from {input_nodes_file} >>>")
    seed_weights_1, total_weight_1, seed_terms = {}, 0, set()
    seed_weights_2, total_weight_2 = {}, 0

    try:
        with open(input_nodes_file, 'r') as f:
            nodes = json.load(f).get('elements', {}).get('nodes', [])
        for node in nodes:
            data = node.get('data', {})
            term, weight_1, weight_2 = data.get(id_key), data.get(weight_key_1), data.get(weight_key_2)
            if term:
                seed_terms.add(term)
                if isinstance(weight_1, (int, float)):
                    seed_weights_1[term] = float(weight_1)
                    total_weight_1 += float(weight_1)
                if isinstance(weight_2, (int, float)):
                    seed_weights_2[term] = float(weight_2)
                    total_weight_2 += float(weight_2)
        print(f"Loaded {len(seed_terms)} unique seed terms for relevance scoring.")
    except Exception as e:
        print(f"ERROR loading seed terms: {e}"); return
    if not seed_terms: return

    # <<< Build Full PMID Pool >>>
    print("\n<<< Building Full PMID Pool via Entrez (with resumability) >>>")
    print("\nYou can interrupt this step and pick-up where you left off")
    print("\nRemeber to save your files if you want to resume this step")

    local_conn = sqlite3.connect(CONTEXTUAL_RELEVANCE_DB_PATH)
    l_cursor = local_conn.cursor()
    l_cursor.execute("CREATE TABLE IF NOT EXISTS relevance_pmid_pool (pmid INTEGER PRIMARY KEY)")
    l_cursor.execute("CREATE TABLE IF NOT EXISTS processed_pool_seeds (term TEXT PRIMARY KEY)")
    local_conn.commit()

    l_cursor.execute("SELECT term FROM processed_pool_seeds")
    processed_terms = {row[0] for row in l_cursor.fetchall()}
    terms_to_process = [term for term in seed_terms if term not in processed_terms]

    if not terms_to_process:
        print("All seed terms have already been processed. Loading PMID pool from cache.")
    else:
        print(f"Found {len(processed_terms)} previously processed terms. Resuming with {len(terms_to_process)} remaining terms.")
        for term in tqdm(terms_to_process, desc="Querying remaining Seed Terms"):
            sanitized_term = term.split('/')[0].strip('*')
            pmids_for_term = get_pmids_date_chunking(
                sanitized_term,
                parse_date_robust(start_date_param),
                parse_date_robust(end_date_param)
            )

            if pmids_for_term:
                pmid_tuples = [(int(p),) for p in pmids_for_term]
                l_cursor.executemany("INSERT OR IGNORE INTO relevance_pmid_pool (pmid) VALUES (?)", pmid_tuples)
            l_cursor.execute("INSERT OR IGNORE INTO processed_pool_seeds (term) VALUES (?)", (term,))
            local_conn.commit()

    l_cursor.execute("SELECT pmid FROM relevance_pmid_pool")
    comprehensive_pmid_pool = {row[0] for row in l_cursor.fetchall()}
    l_cursor.close()
    local_conn.close()
    print(f"Full relevance pool contains {len(comprehensive_pmid_pool):,} unique PMIDs.")

    print("\n<<< Checking for and Populating Missing Annotations in Master Database >>>")
    print("\nNOTICE [!] Best pracitce is to keep your master_mesh_database.db between independent runs")
    print("The first time you run this it might hurt your soul. FYI per 10^6MIDs \u2248 1.5hr")
    print("You can stop and resume anytime as long as you save all files associated")
    populate_master_mesh_database(list(comprehensive_pmid_pool), master_db_path)

    print("\n<<< Calculating Scores Directly from Master Database >>>")
    master_conn = sqlite3.connect(f'file:{master_db_path}?mode=ro', uri=True)
    m_cursor = master_conn.cursor()

    aggregator_1 = defaultdict(lambda: {'sum': 0.0, 'count': 0})
    aggregator_2 = defaultdict(lambda: {'sum': 0.0, 'count': 0})
    article_scores_data = []
    pmid_list = list(comprehensive_pmid_pool)
    chunk_size = 50000
    seed_term_set = set(seed_weights_1.keys())

    for i in tqdm(range(0, len(pmid_list), chunk_size), desc="Scoring Articles from Master Database"):
        batch_pmids = pmid_list[i : i + chunk_size]
        placeholders = ','.join('?' for _ in batch_pmids)
        query = f"SELECT pmid, mesh_terms FROM master_mesh_annotations WHERE pmid IN ({placeholders})"
        m_cursor.execute(query, batch_pmids)
        for row in m_cursor:
            pmid, mesh_terms_str = row
            if not isinstance(mesh_terms_str, str): continue
            article_terms = set(mesh_terms_str.split(';'))
            matching_seeds = article_terms.intersection(seed_term_set)
            if not matching_seeds: continue
            score_1, score_2 = 0.0, 0.0
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
            article_scores_data.append({
                'pmid': pmid,
                f'score_{weight_key_1}': score_1,
                f'score_{weight_key_2}': score_2,
                'contributing_seeds': ';'.join(sorted(list(matching_seeds)))
            })
    master_conn.close()

    print("\n<<< Saving Individual Article Relevance Scores to Database >>>")
    if article_scores_data:
        scores_df = pd.DataFrame(article_scores_data)
        try:
            conn = sqlite3.connect(CONTEXTUAL_RELEVANCE_DB_PATH)
            scores_df.to_sql('article_relevance_scores', conn, if_exists='replace', index=False)
            cursor = conn.cursor()
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_pmid ON article_relevance_scores (pmid)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_score1 ON article_relevance_scores ('score_{weight_key_1}')")
            conn.commit()
            conn.close()
            print(f"Successfully saved {len(scores_df):,} article scores to the database: {CONTEXTUAL_RELEVANCE_DB_PATH}")
        except Exception as e:
            print(f"ERROR: Could not save article scores to database: {e}")
    else:
        print("No article scores were generated to save.")

    final_node_weights_1 = {term: (data['sum'] / data['count']) * np.log10(data['count'] + 1) for term, data in aggregator_1.items() if data['count'] > 0}
    final_node_weights_2 = {term: (data['sum'] / data['count']) * np.log10(data['count'] + 1) for term, data in aggregator_2.items() if data['count'] > 0}

    print("\n<<< Generating Final Output JSON File >>>")
    with open(input_nodes_file, 'r') as f:
        network_data = json.load(f)
    for node in network_data.get('elements', {}).get('nodes', []):
        term = node.get('data', {}).get(id_key)
        node['data'][final_key_1] = final_node_weights_1.get(term, 0.0)
        node['data'][final_key_2] = final_node_weights_2.get(term, 0.0)

    with open(output_nodes_file, 'w') as f:
        json.dump(network_data, f, indent=2)
    print(f"Success! New file created at: {output_nodes_file}")

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 6. Main Execution Block
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if __name__ == "__main__":

    # <<< STEP 1: Initial Data Collection >>>
    print("\n\n<<< STEP 1: Initial Data Collection >>>")
    print("\nWARNING: Interruptions in STEP 1-3 may lead to incomplete " f"{FILE_PREFIX}_pmids.db")
    if os.path.exists(original_db_path):
        print(f"Output file '{original_db_path}' already exists. Skipping Step 1.")
    else:
        run_initial_data_collection(search_term, start_date, end_date, generations_n, original_db_path)

    # <<< STEP 2: Database Cleaning >>>
    print("\n\n<<< STEP 2: Database Cleaning >>>")
    if os.path.exists(cleaned_db_path):
        print(f"Output file '{cleaned_db_path}' already exists. Skipping Step 2.")
    else:
        clean_database(original_db_path, cleaned_db_path, start_generation_to_remove)

    # <<< 2.5: Populate Master MeSH Database >>>
    print("\n\n<<< STEP 2.5: Populate Master MeSH Database >>>")
    print("\nThis Step will take longer the bigger the Master db file is.")
    populate_master_mesh_database(
        cleaned_db_path,
        MASTER_ANNOTATION_DB_PATH,
        failed_log=FAILED_MESH_FETCH_TSV,
        empty_log=EMPTY_MESH_PMIDS_TSV
    )

    # <<< STEP 3: MeSH Term Annotation (from Master DB) >>>
    print("\n\n<<< STEP 3: Local MeSH Term Annotation >>>")
    update_db_with_mesh_batch(cleaned_db_path, MASTER_ANNOTATION_DB_PATH)

    # <<< STEP 4: Network Construction >>>
    print("\n\n<<< STEP 4: Network Construction >>>")
    if os.path.exists(NETWORK_JSON_FULL):
        print(f"Output file '{NETWORK_JSON_FULL}' already exists. Skipping Step 4.")
    else:
        run_network_construction(DB_PATH_FOR_NETWORK_BUILD, NETWORK_JSON_FULL)

    # <<< STEP 5: Network Filtering >>>
    print("\n\n<<< STEP 5: Consensus Network Filtering >>>")
    if os.path.exists(CONSENSUS_LCC_FILE):
        print(f"Output file '{os.path.basename(CONSENSUS_LCC_FILE)}' already exists. Skipping Step 5.")
    else:
        run_consensus_filtering_and_lcc(
            input_json_path=NETWORK_JSON_FULL,
            glf_output_path=GLF_OUTPUT_FILE,
            sa_output_path=SA_OUTPUT_FILE,
            final_lcc_output_path=CONSENSUS_LCC_FILE
        )

    # <<< STEP 6: Community Detection >>>
    print("\n\n<<< STEP 6: Community Detection >>>")
    run_step_6 = True
    if os.path.exists(CONSENSUS_LCC_FILE):
        try:
            with open(CONSENSUS_LCC_FILE, 'r') as f:
                temp_data = json.load(f)
                if temp_data.get('elements', {}).get('nodes', [{}])[0].get('data', {}).get('louvain_community') is not None:
                    print(f"Community IDs already present in '{CONSENSUS_LCC_FILE}'. Skipping Step 6.")
                    run_step_6 = False
        except (IOError, json.JSONDecodeError, IndexError):
            pass

    if run_step_6:
        run_community_detection(CONSENSUS_LCC_FILE)

    # <<< STEP 7: Contextual Relevance Scoring >>>
    print("\n\n<<< STEP 7: Contextual Relevance Scoring >>>")
    print("\nSave Master db between runs to significantly speed-up this step")
    if os.path.exists(FINAL_OUTPUT_WITH_RELEVANCE):
        print(f"Output file '{FINAL_OUTPUT_WITH_RELEVANCE}' already exists. Skipping Step 7.")
    else:
        run_contextual_relevance_scoring(
            input_nodes_file=CONSENSUS_LCC_FILE,
            output_nodes_file=FINAL_OUTPUT_WITH_RELEVANCE,
            master_db_path=MASTER_ANNOTATION_DB_PATH,
            id_key='id',
            weight_key_1=CONTEXT_INITIAL_WEIGHT_KEY_COL_1,
            final_key_1=CONTEXT_FINAL_WEIGHT_KEY_COL_1,
            weight_key_2=CONTEXT_INITIAL_WEIGHT_KEY_COL_2,
            final_key_2=CONTEXT_FINAL_WEIGHT_KEY_COL_2,
            start_date_param=CONTEXT_START_DATE,
            end_date_param=CONTEXT_END_DATE,
            max_pmids=MAX_PMIDS_PER_TERM
        )
    # <<< Cleanup Empty Log Files >>>
    logging.shutdown()
    if os.path.exists(LOG_FILE_PATH) and os.path.getsize(LOG_FILE_PATH) == 0:
        try:
            os.remove(LOG_FILE_PATH)
        except Exception:
            pass

    print("\n\nMaster script finished.")
