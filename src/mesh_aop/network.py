# -*- coding: utf-8 -*-
"""
network.py - co-occurrence network construction and filtering (pipeline Step 3).

Turns the annotated article database into a weighted MeSH co-occurrence network
and reduces it to an interpretable consensus subgraph.

The module builds the full NetworkX graph together with per-node and per-edge
statistics (generation-decayed weights, citation rank normalizations, MLF edge
p-values, degree, betweenness and eigenvector centrality), runs two independent
heuristic optimizers - Graph Likelihood Filtering and Simulated Annealing -
keeps the edges they agree on, extracts the largest connected component, and
assigns Louvain communities. All stochastic steps are seeded for reproducibility.
"""

import os
import json
import random
import sqlite3
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats
from scipy.stats import binom
from tqdm import tqdm

try:
    import community.community_louvain as community_louvain
except ImportError:
    raise ImportError("The 'python-louvain' library is missing. Install via pip.")

# Import utilities, stats, and stop words from our package.
# The generated stop-word file is a list; coerce to a frozenset once so the
# per-term membership tests in the article loop are O(1) instead of O(n).
from .mesh_stop_words import MESH_STOP_WORDS as _RAW_STOP_WORDS
MESH_STOP_WORDS = frozenset(_RAW_STOP_WORDS)
from .data_ops import parse_mesh_terms, get_generation_label
from .stats import calculate_graph_stats, run_simulation


class _ListSampleRandom(random.Random):
    """A seeded Random whose sample() tolerates set-like populations.

    Works around a networkx 2.8.8 bug: edge_betweenness_centrality passes
    G.nodes() (a set-like view) straight to random.sample(), which Python 3.11+
    rejects with "Population must be a sequence". Wrapping non-sequence
    populations in a list restores the pre-3.11 behaviour, leaving the already
    list-wrapped betweenness_centrality call unchanged. Seeding keeps the
    k-sample estimate reproducible.
    """
    def sample(self, population, k, *args, **kwargs):
        if not isinstance(population, (list, tuple)):
            population = list(population)
        return super().sample(population, k, *args, **kwargs)


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# I/O AND HELPER FUNCTIONS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def load_full_graph_data(filepath: str) -> tuple:
    """Loads a full graph from a Cytoscape.js JSON file, preserving all attributes."""
    all_nodes_data = {}
    all_edges_data = {}
    print(f"Loading full graph data from {os.path.basename(filepath)}...")
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Graph file not found at {filepath}")

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

        print(f"  -> Loaded {len(all_nodes_data)} nodes and {len(all_edges_data)} edges.")
        return all_nodes_data, all_edges_data
    except Exception as e:
        raise RuntimeError(f"Failed to parse graph file: {e}")

def save_subgraph_to_json(filepath: str, subgraph_keys: set, all_nodes_data: dict, all_edges_data: dict):
    """Saves a subgraph to a Cytoscape.js JSON, preserving all original attributes."""
    try:
        nodes_in_subgraph = {node for key in subgraph_keys for node in key}
        edges_list = [{'data': all_edges_data[key]} for key in subgraph_keys if key in all_edges_data]
        nodes_list = [{'data': all_nodes_data[node_id]} for node_id in sorted(list(nodes_in_subgraph)) if node_id in all_nodes_data]
        cytoscape_data = {'elements': {'nodes': nodes_list, 'edges': edges_list}}

        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, 'w') as f:
            json.dump(cytoscape_data, f, indent=2)
        print(f"Successfully saved graph with {len(nodes_list)} nodes and {len(edges_list)} edges to {filepath}")
    except Exception as e:
        raise RuntimeError(f"Failed to save subgraph JSON to {filepath}: {e}")

def calculate_generation_weight(generation, lambda_val: float) -> float:
    """Calculates centrality weight based on generation string ('P0', 'G1', ...)."""
    try:
        if pd.isna(generation):
            distance = float('inf')
        elif not isinstance(generation, str):
            distance = int(generation)
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
        return 0.0

def normalize_cited_by_per_mesh(df: pd.DataFrame, stop_words_set: set) -> tuple:
    """Calculates article citation ranks and MeSH term citation statistics."""
    print("  Calculating numeric citation counts per article...")
    df['cited_by_count'] = df['cited_by'].apply(
        lambda x: len(x.split(';')) if isinstance(x, str) and x.strip() and x.lower() != 'nan' else 0
    ).astype(int)

    print("  Vectorizing MeSH terms for aggregation...")
    def extract_main_mesh_terms(mesh_string):
        if not isinstance(mesh_string, str):
            return []
        return [
            term.strip().lstrip('*') for term in mesh_string.split(';')
            if '/' not in term and term.strip()
        ]

    df_exploded = df[['mesh_terms', 'cited_by_count']].copy()
    df_exploded['mesh_terms'] = df_exploded['mesh_terms'].apply(extract_main_mesh_terms)
    df_exploded = df_exploded.explode('mesh_terms')
    df_exploded.rename(columns={'mesh_terms': 'mesh_term'}, inplace=True)
    df_exploded.dropna(subset=['mesh_term'], inplace=True)

    stop_words_lower = {s.lower() for s in stop_words_set}
    df_exploded = df_exploded[~df_exploded['mesh_term'].str.lower().isin(stop_words_lower)]

    print("  Aggregating citation counts per term using groupby...")
    if not df_exploded.empty:
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

    return df, mesh_stats_df

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CORE PIPELINE OPERATIONS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def run_network_construction(db_path_param: str, output_json_path: str, lambda_val: float,
                             node_weight_factors: dict, run_full_centrality: bool,
                             betweenness_k_samples: int, random_seed: int,
                             eigenvector_max_iter: int, eigenvector_tol: float):
    """Constructs the co-occurrence network and calculates topological centralities."""
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
        raise RuntimeError(f"Data loading error: {e}")

    if df.empty:
        print("Loaded DataFrame is empty. Exiting network construction.")
        return

    print("\n<<< Pre-computation >>>")
    df['generation_weight'] = df['generation'].apply(calculate_generation_weight, lambda_val=lambda_val)
    df, mesh_stats_df = normalize_cited_by_per_mesh(df, MESH_STOP_WORDS)

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
    # itertuples is ~10-100x faster than iterrows for wide frames; access by attr.
    for row in tqdm(df.itertuples(index=False), total=len(df), desc="Processing Articles"):
        article_mesh_terms_str = getattr(row, 'mesh_terms', '')
        if not isinstance(article_mesh_terms_str, str):
            continue

        article_centrality = getattr(row, 'generation_weight', 0.0)
        article_norm_cited_by = getattr(row, 'normalized_cited_by', 0.0)
        article_generation_raw = getattr(row, 'generation', None)

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
            if term_type == 'mesh_term' and term_name and term_name not in MESH_STOP_WORDS:
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

        if len(article_counts) > 1:
            ranks = stats.rankdata(article_counts, method='average') / len(article_counts)
        elif len(article_counts) == 1:
            ranks = [0.5]
        else:
            ranks = []

        for i, term_name in enumerate(node_keys):
            node_data_agg[term_name]['article_count_rank_normalized'] = ranks[i] if i < len(ranks) else 0.0

        total_factor = sum(node_weight_factors.values())
        norm_factors = {k: v / total_factor if total_factor > 1e-9 else 0 for k, v in node_weight_factors.items()}

        for term_name in node_keys:
            data = node_data_agg[term_name]
            article_count = mesh_term_article_counts.get(term_name, 0)
            data['article_count'] = article_count
            data['major_topic_proportion'] = articles_with_term_as_major_count.get(term_name, 0) / article_count if article_count > 0 else 0.0
            data['adjusted_node_weight'] = (data.get('generation_weight', 0.0) * norm_factors['centrality'] +
                                            data.get('article_count_rank_normalized', 0.0) * norm_factors['article_rank'] +
                                            data.get('rank_norm_median_cit', 0.0) * norm_factors['rank_median_cit'] +
                                            data.get('rank_norm_total_cit', 0.0) * norm_factors['rank_total_cit'])

    print("\n<<< Calculating Network Centrality Measures & Communities >>>")
    G = nx.Graph()
    if edges:
        # NOTE: edges are added without a 'weight' attribute, so the centrality
        # calls below run UNWEIGHTED (this matches the prior behaviour: passing
        # weight='weight' with no such attribute silently defaulted to 1.0).
        # To weight the topology, populate per-metric attributes here, e.g.
        #   for (u, v), d in edges.items():
        #       G.add_edge(u, v, affinity=d['cooccurrence_count'],
        #                        distance=1.0 / d['cooccurrence_count'])
        # then use weight='affinity' for eigenvector and weight='distance' for
        # betweenness (betweenness treats weight as path length, so it must be
        # inverted). This is a modelling choice and is left disabled by default.
        G.add_edges_from(edges.keys())

    G_analysis_conn = G
    if len(G) > 0 and not nx.is_connected(G):
        print("Graph is not connected. Using largest connected component for some metrics.")
        largest_cc = max(nx.connected_components(G), key=len)
        G_analysis_conn = G.subgraph(largest_cc).copy()

    try:
        if G_analysis_conn.number_of_nodes() > 0:
            degree_dict = dict(G.degree())
            clustering_coefficient = nx.clustering(G)
            partition_map = community_louvain.best_partition(G_analysis_conn, random_state=random_seed)

            # Type-safety fallbacks for config omissions
            safe_max_iter = int(eigenvector_max_iter) if eigenvector_max_iter is not None else 1000
            safe_tol = float(eigenvector_tol) if eigenvector_tol is not None else 1.0e-6

            # weight=None => unweighted, stated explicitly (see note where G is built).
            print(f"  Calculating Eigenvector Centrality (Max Iterations: {safe_max_iter}, Tol: {safe_tol})...")
            try:
                eigenvector_centrality = nx.eigenvector_centrality_numpy(
                    G_analysis_conn, weight=None, max_iter=safe_max_iter, tol=safe_tol
                )
            except Exception:
                eigenvector_centrality = nx.eigenvector_centrality(
                    G_analysis_conn, weight=None, max_iter=safe_max_iter, tol=safe_tol
                )

            if run_full_centrality:
                print("  Calculating EXACT Betweenness Centrality (This may take a while)...")
                betweenness_centrality = nx.betweenness_centrality(G_analysis_conn, k=None, normalized=True, weight=None)
                edge_betweenness_centrality = nx.edge_betweenness_centrality(G_analysis_conn, k=None, normalized=True, weight=None)
            else:
                print(f"  Calculating ESTIMATED Betweenness Centrality (k={betweenness_k_samples})...")
                k_eff = min(betweenness_k_samples, G_analysis_conn.number_of_nodes())
                betweenness_centrality = nx.betweenness_centrality(G_analysis_conn, k=k_eff, normalized=True, weight=None, seed=_ListSampleRandom(random_seed))
                edge_betweenness_centrality = nx.edge_betweenness_centrality(G_analysis_conn, k=k_eff, normalized=True, weight=None, seed=_ListSampleRandom(random_seed))

        else:
            degree_dict, betweenness_centrality, eigenvector_centrality, clustering_coefficient, edge_betweenness_centrality, partition_map = {}, {}, {}, {}, {}, {}

    except MemoryError:
        raise RuntimeError("Memory Limit Exceeded during centrality calculations. Consider reducing generations.")

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

    final_output_path = Path(output_json_path)
    temp_output_path = final_output_path.with_suffix('.json.tmp')

    try:
        final_output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(temp_output_path, 'w', encoding='utf-8') as f:
            json.dump(cytoscape_data, f, indent=2)

        temp_output_path.replace(final_output_path)
        print(f"Successfully saved network data to {final_output_path}")

    except Exception as e:
        if temp_output_path.exists():
            temp_output_path.unlink()
        raise RuntimeError(f"Failed to write final JSON file: {e}")

def run_consensus_filtering_and_lcc(input_json_path: str, glf_output_path: str,
                                    sa_output_path: str, final_lcc_output_path: str,
                                    history_output_path: str,
                                    target_num_edges: int, glf_iterations: int,
                                    sa_iterations: int, sa_initial_temp: float,
                                    sa_cooling_rate: float, random_seed: int = 42):
    """
    Performs network filtering using a consensus of GLF and SA optimizers,
    saves the trajectory history to bypass redundancy,
    and extracts the Largest Connected Component (LCC).
    """
    print(f"\n" + "<"*30 + ">"*30)
    print(f"<<< Running Consensus Network Filtering >>>")
    print("<"*30 + ">"*30)

    all_nodes_data, all_edges_data = load_full_graph_data(input_json_path)
    if not all_edges_data:
        print("No edges in the unfiltered graph. Aborting filtering.")
        return

    global_T, global_node_strengths = calculate_graph_stats(all_edges_data)

    # Distinct (but reproducible) seeds so GLF and SA explore independently.
    glf_keys, glf_score, glf_history = run_simulation('GLF', all_edges_data, global_node_strengths, global_T, target_num_edges, glf_iterations, random_seed=random_seed)
    sa_keys, sa_score, sa_history = run_simulation('SA', all_edges_data, global_node_strengths, global_T, target_num_edges, sa_iterations, random_seed=random_seed + 1, initial_temp=sa_initial_temp, cooling_rate=sa_cooling_rate)

    try:
        os.makedirs(os.path.dirname(history_output_path), exist_ok=True)
        with open(history_output_path, 'w') as f:
            json.dump({"GLF": glf_history, "SA": sa_history}, f)
        print(f"  [+] Saved simulation trajectory data to bypass redundant calculations later.")
    except Exception as e:
        print(f"  [!] Warning: Could not save optimization history: {e}")

    save_subgraph_to_json(glf_output_path, glf_keys, all_nodes_data, all_edges_data)
    save_subgraph_to_json(sa_output_path, sa_keys, all_nodes_data, all_edges_data)

    print("\n<<< Building Consensus Network >>>")
    consensus_keys = glf_keys.intersection(sa_keys)
    print(f"Found {len(consensus_keys)} consensus edges (intersection of GLF and SA).")

    print("\n<<< Extracting Largest Connected Component (LCC) from Consensus Network >>>")
    if not consensus_keys:
        print("Consensus network is empty. No LCC to extract.")
        os.makedirs(os.path.dirname(final_lcc_output_path), exist_ok=True)
        with open(final_lcc_output_path, 'w') as f:
            json.dump({"elements": {"nodes": [], "edges": []}}, f)
        return

    G_consensus = nx.Graph()
    G_consensus.add_edges_from(list(consensus_keys))

    components = list(nx.connected_components(G_consensus))
    if not components:
        print("No connected components found in consensus network.")
        os.makedirs(os.path.dirname(final_lcc_output_path), exist_ok=True)
        with open(final_lcc_output_path, 'w') as f:
            json.dump({"elements": {"nodes": [], "edges": []}}, f)
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

def run_community_detection(network_file_path: str, random_seed: int):
    """Runs Louvain, re-orders community IDs by size, and saves back to JSON."""
    print(f"Loading network for community detection from: {network_file_path}")

    net_path = Path(network_file_path)

    if not net_path.exists():
        raise FileNotFoundError(f"Input file not found at '{net_path}'")

    try:
        with open(net_path, 'r', encoding='utf-8') as f:
            network_data = json.load(f)
    except Exception as e:
        raise RuntimeError(f"Error reading JSON for community detection: {e}")

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
        partition_map = community_louvain.best_partition(G, random_state=random_seed)
        print(f"Found {len(set(partition_map.values()))} communities.")

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
                node_obj['data']['filtered_louvain_community_id'] = -1

    final_data = {"elements": {"nodes": nodes_json, "edges": edges_json}}

    print(f"Updating network with re-ordered community IDs and saving back to: {net_path}")

    temp_file = net_path.with_suffix('.json.tmp')

    try:
        with open(temp_file, 'w', encoding='utf-8') as f:
            json.dump(final_data, f, indent=2)

        temp_file.replace(net_path)
        print("Community detection complete.")

    except Exception as e:
        if temp_file.exists():
            temp_file.unlink()
        raise RuntimeError(f"Error saving updated community data to JSON: {e}")