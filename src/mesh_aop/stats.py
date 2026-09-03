# -*- coding: utf-8 -*-
"""
stats.py - statistical models and graph-optimization routines.

Provides the mathematical core used during network filtering and reporting,
kept separate from the graph-assembly logic in network.py.

It contains the maximum-likelihood edge model and a unified optimizer that runs
either Global Likelihood Filter or Simulated Annealing to select a fixed-size,
high-likelihood subgraph, plus helpers that compute node strengths and summarize
the node and edge counts retained at each stage of the filtering cascade.
"""

import os
import json
import math
import random
from collections import defaultdict
import networkx as nx

def calculate_graph_stats(edges_data: dict) -> tuple:
    """Calculates total strength and individual node strengths from edge data."""
    node_strengths = defaultdict(int)
    total_strength = 0
    for edge_key, data in edges_data.items():
        w = data.get('cooccurrence_count', 0)
        node_strengths[edge_key[0]] += w
        node_strengths[edge_key[1]] += w
        total_strength += w
    return total_strength, node_strengths

def get_log_likelihood_term(edge_data: dict, node_strengths: dict, denominator: float) -> float:
    """Calculates the log-likelihood component for a specific edge."""
    w_ij = edge_data.get('cooccurrence_count', 0)
    k_i = node_strengths.get(edge_data.get('source'), 0)
    k_j = node_strengths.get(edge_data.get('target'), 0)
    if k_i == 0 or k_j == 0 or w_ij == 0:
        return 0.0
    p_ij = (k_i * k_j) / denominator
    return w_ij * math.log(p_ij) - math.lgamma(w_ij + 1) if p_ij > 0 else 0.0

def run_simulation(method: str, all_edges: dict, node_strengths: dict, total_T: float, target_edges: int, iterations: int, random_seed: int = None, **kwargs) -> tuple:
    """Unified runner for Global Likelihood Filter (GLF) or Simulated Annealing (SA).

    A local Random instance seeded with `random_seed` makes the stochastic
    search reproducible without mutating the global `random` state.
    """
    print(f"Starting {method} Simulation ({iterations:,} iters)...")
    try:
        rng = random.Random(random_seed)
        denominator = 2 * (total_T**2)
        keys = list(all_edges.keys())
        # Guard against requesting more edges than exist in the graph.
        target_edges = min(target_edges, len(keys))
        if target_edges <= 0:
            print(f"[!] Simulation {method}: no edges available to sample.")
            return set(), 0.0, []
        current_keys = set(rng.sample(keys, target_edges))

        curr_T = sum(all_edges[k].get('cooccurrence_count', 0) for k in current_keys)
        curr_term = sum(get_log_likelihood_term(all_edges[k], node_strengths, denominator) for k in current_keys)
        curr_ll = math.lgamma(curr_T + 1) + curr_term

        best_keys, min_ll = current_keys.copy(), curr_ll
        history = []

        temp = kwargs.get('initial_temp', 0.0)
        cooling = kwargs.get('cooling_rate', 0.0)

        interval = max(1, iterations // 5)
        for i in range(iterations):
            if i > 0 and i % interval == 0:
                print(f"    {method} Progress: {int((i / iterations) * 100)}% ({i:,} / {iterations:,})")
            on = rng.choice(list(current_keys))
            off = rng.choice(keys)
            while off in current_keys:
                off = rng.choice(keys)

            w_on = all_edges[on].get('cooccurrence_count', 0)
            w_off = all_edges[off].get('cooccurrence_count', 0)
            prop_T = curr_T - w_on + w_off

            ll_on = get_log_likelihood_term(all_edges[on], node_strengths, denominator)
            ll_off = get_log_likelihood_term(all_edges[off], node_strengths, denominator)

            delta = (math.lgamma(prop_T + 1) - math.lgamma(curr_T + 1)) + (ll_off - ll_on)

            accept = False
            if delta < 0:
                accept = True
            elif method == 'GLF' and rng.random() < math.exp(-delta):
                accept = True
            elif method == 'SA' and temp > 1e-9 and rng.random() < math.exp(-delta / temp):
                accept = True

            if accept:
                current_keys.remove(on)
                current_keys.add(off)
                curr_ll += delta
                curr_T = prop_T
                if curr_ll < min_ll:
                    min_ll = curr_ll
                    best_keys = current_keys.copy()

            if method == 'SA':
                temp *= cooling

            if i % 10000 == 0:
                history.append((i, curr_ll))

        return best_keys, min_ll, history
    except Exception as e:
        print(f"[!] Simulation {method} failed: {e}")
        return set(), 0.0, []

def load_basic_graph(filepath: str):
    """Loads a Cytoscape JSON into a basic NetworkX graph for structural analysis."""
    if not os.path.exists(filepath):
        print(f"[!] Warning: Missing intermediate file for Sankey: {filepath}")
        return None

    with open(filepath, 'r') as f:
        data = json.load(f)

    G_temp = nx.Graph()
    for n in data.get('elements', {}).get('nodes', []):
        G_temp.add_node(n['data']['id'])
    for e in data.get('elements', {}).get('edges', []):
        G_temp.add_edge(e['data']['source'], e['data']['target'])
    return G_temp

def calculate_dynamic_sankey_data(glf_path: str, sa_path: str, final_G: nx.Graph) -> dict:
    """Dynamically calculates node/edge counts across the filtering cascade."""
    print("\n<<< Calculating Network Filtering Statistics >>>")
    try:
        G_glf = load_basic_graph(glf_path)
        G_sa = load_basic_graph(sa_path)

        if G_glf is None or G_sa is None:
            return None

        glf_lcc_nodes = max(nx.connected_components(G_glf), key=len)
        G_glf_lcc = G_glf.subgraph(glf_lcc_nodes)

        sa_lcc_nodes = max(nx.connected_components(G_sa), key=len)
        G_sa_lcc = G_sa.subgraph(sa_lcc_nodes)

        final_nodes = final_G.number_of_nodes()
        final_edges = final_G.number_of_edges()

        return {
            'Initial Subgraph': {
                'GLF': {'Nodes': G_glf.number_of_nodes(), 'Edges': G_glf.number_of_edges()},
                'SA': {'Nodes': G_sa.number_of_nodes(), 'Edges': G_sa.number_of_edges()}
            },
            'LCC': {
                'GLF': {'Nodes': G_glf_lcc.number_of_nodes(), 'Edges': G_glf_lcc.number_of_edges()},
                'SA': {'Nodes': G_sa_lcc.number_of_nodes(), 'Edges': G_sa_lcc.number_of_edges()}
            },
            'Final Consensus LCC': {
                'GLF': {'Nodes': final_nodes, 'Edges': final_edges},
                'SA': {'Nodes': final_nodes, 'Edges': final_edges}
            }
        }
    except Exception as e:
        print(f"[!] Error calculating dynamic Sankey data: {e}")
        return None