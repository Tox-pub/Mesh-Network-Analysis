# -*- coding: utf-8 -*-
"""
gt_network_validation.py - node/edge convergent validation against a ground truth.

Where `benchmark.py` asks "does the per-article score rank the ground truth
highly?", this module asks the complementary structural question: "is the
network's own vocabulary and wiring independently reproduced by the ground-truth
literature?" A high article score could in principle come from a network that
merely names popular terms; recovering the network's *edges* from an independent
corpus is much harder to achieve by chance.

The ground-truth articles' MeSH terms are pulled from the master database and
filtered through the SAME node-eligibility stop-word list used to build the
network, so stop-listed headings (check tags, organisms, geographics) can never
appear as spurious "near misses" - they were never eligible to be nodes.

Reported at two levels, each calibrated against a permutation null of random
article draws of the same size (both corpora are "about the topic", so raw
overlap alone is not evidence):

  * node level - network nodes attested in the ground truth, their enrichment
    over corpus base rates, and the rank correlation between a node's
    ground-truth prominence and its network weight;
  * edge level - how many network edges reappear as ground-truth co-occurrences.

Also emits the "misses" (ground-truth-frequent eligible terms that are NOT nodes),
which localize where the network under-covers the domain.
"""

import os
import json
import sqlite3
import random
import datetime
from collections import Counter
from itertools import combinations

import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import spearmanr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D

from .mesh_stop_words import MESH_STOP_WORDS as _RAW_STOP_WORDS
from .benchmark import validate_ground_truth, _open_readonly_resilient

MESH_STOP_WORDS = frozenset(_RAW_STOP_WORDS)


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# HELPER FUNCTIONS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def _eligible_terms(mesh_str: str) -> set:
    """Base MeSH descriptors of one article, restricted to node-eligible terms.

    Strips the major-topic asterisk and any subheading, then removes stop-listed
    headings so the result is exactly the universe a network node could come from.
    """
    if not mesh_str:
        return set()
    out = set()
    for raw in mesh_str.split(';'):
        base = raw.split('/')[0].lstrip('*').strip()
        if base:
            out.add(base)
    return out - MESH_STOP_WORDS


def _raw_terms(mesh_str: str) -> set:
    """Base MeSH descriptors with no stop-word filtering (for the audit trail)."""
    if not mesh_str:
        return set()
    return {r.split('/')[0].lstrip('*').strip() for r in mesh_str.split(';') if r.strip()}


def _load_network(final_network_path: str, weight_key: str):
    """Load node ids, the chosen node weight, and the edge set from the final network."""
    with open(final_network_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    nodes = data.get('elements', {}).get('nodes', [])
    node_ids, weights, extras = set(), {}, {}
    for n in nodes:
        d = n.get('data', {})
        nid = d.get('id')
        if not nid:
            continue
        node_ids.add(nid)
        weights[nid] = float(d.get(weight_key, 0.0) or 0.0)
        extras[nid] = {
            'betweenness_centrality': float(d.get('betweenness_centrality', 0.0) or 0.0),
            'article_count': int(d.get('article_count', 0) or 0),
        }

    G = nx.Graph()
    for e in data.get('elements', {}).get('edges', []):
        d = e.get('data', {})
        s, t = d.get('source'), d.get('target')
        if s and t:
            G.add_edge(s, t, cooccurrence_count=int(d.get('cooccurrence_count', 0) or 0))

    return node_ids, weights, extras, G


def _sample_background_pool(master_db_path: str, pool_size: int, seed: int) -> list:
    """Draw a random article sample to estimate corpus base rates and the nulls.

    PMIDs are the table's INTEGER PRIMARY KEY, so random ids are cheap indexed
    lookups; misses (gaps in the id space) are simply skipped.
    """
    rng = random.Random(seed)
    conn = _open_readonly_resilient(master_db_path)
    cur = conn.cursor()
    max_pmid = cur.execute("SELECT MAX(pmid) FROM master_mesh_annotations").fetchone()[0] or 1

    pool = []
    while len(pool) < pool_size:
        ids = [rng.randint(1, max_pmid) for _ in range(int((pool_size - len(pool)) * 1.35) + 2000)]
        for i in range(0, len(ids), 1500):
            chunk = ids[i:i + 1500]
            q = ("SELECT mesh_terms FROM master_mesh_annotations WHERE pmid IN (%s)"
                 % ",".join("?" * len(chunk)))
            for (mt,) in cur.execute(q, chunk):
                pool.append(_eligible_terms(mt))
                if len(pool) >= pool_size:
                    break
            if len(pool) >= pool_size:
                break
    conn.close()
    return pool


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FIGURES
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def _plot_node_validation(shared, freq, weights, gt_terms_df, rho, n_gt, figures_dir, prefix):
    """Scatter of ground-truth prominence vs network weight, plus the top-term bar."""
    fig, ax = plt.subplots(1, 2, figsize=(15, 6.2))

    xs = [freq[t] for t in shared]
    ys = [weights.get(t, 0.0) for t in shared]
    ax[0].scatter(xs, ys, s=55, c=[freq[t] / n_gt for t in shared], cmap='viridis',
                  edgecolors='k', linewidths=.3, zorder=3)
    if any(y > 0 for y in ys):
        ax[0].set_yscale('log')
    ax[0].set_xlabel(f'Ground-truth article frequency (of {n_gt})')
    ax[0].set_ylabel('Network node weight (log)')
    ax[0].set_title(f'Shared nodes: GT prominence vs network weight\n'
                    f'Spearman rho = {rho:+.2f}  (n={len(shared)})')
    for t in sorted(shared, key=lambda x: -freq[x])[:10]:
        ax[0].annotate(t, (freq[t], weights.get(t, 0.0)), fontsize=7.5,
                       xytext=(4, 2), textcoords='offset points')
    ax[0].grid(True, ls='--', alpha=.35)

    top = gt_terms_df.head(22)[::-1]
    cols = ['#3b7d3b' if v == 'Y' else '#b0b0b0' for v in top['in_network']]
    ax[1].barh(range(len(top)), top['GT_articles'], color=cols)
    ax[1].set_yticks(range(len(top)))
    ax[1].set_yticklabels(top['term'], fontsize=8)
    ax[1].set_xlabel('Ground-truth article frequency')
    ax[1].set_title('Top eligible ground-truth terms\n(green = network node, grey = miss)')
    ax[1].grid(True, axis='x', ls='--', alpha=.35)

    plt.tight_layout()
    out = os.path.join(figures_dir, f"{prefix}_GT_Node_Validation.png")
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close('all')
    return out


def _plot_gt_network(GT, freq, node_ids, net_edges, n_recovered, figures_dir, prefix, min_articles):
    """Draw the ground-truth co-occurrence network, coloured by overlap with ours."""
    plt.figure(figsize=(13, 10.5))
    ax = plt.gca()
    pos = nx.spring_layout(GT, seed=42, k=0.55, iterations=200, weight='w')

    ncol = ['#2c7fb8' if n in node_ids else '#d95f0e' for n in GT.nodes()]
    nsize = [120 + freq[n] * 55 for n in GT.nodes()]
    ecol, ew = [], []
    for u, v in GT.edges():
        hit = frozenset((u, v)) in net_edges
        ecol.append('#2c7fb8' if hit else '#cccccc')
        ew.append(1.6 if hit else 0.6)

    nx.draw_networkx_edges(GT, pos, edge_color=ecol, width=ew, alpha=.6, ax=ax)
    nx.draw_networkx_nodes(GT, pos, node_color=ncol, node_size=nsize, alpha=.9,
                           linewidths=.5, edgecolors='#333', ax=ax)
    nx.draw_networkx_labels(GT, pos, labels={n: n for n in GT.nodes() if freq[n] >= 3},
                            font_size=7.5, ax=ax)

    ax.legend(handles=[
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#2c7fb8', markersize=11,
               label='term IS a network node'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#d95f0e', markersize=11,
               label='ground-truth-only term (network miss)'),
        Line2D([0], [0], color='#2c7fb8', lw=2.2, label='edge present in our network (recovered)'),
        Line2D([0], [0], color='#cccccc', lw=1.5, label='ground-truth-only co-occurrence'),
    ], loc='upper left', fontsize=9, framealpha=.95)

    shared_n = len(set(GT.nodes()) & node_ids)
    ax.set_title(f'Ground-truth MeSH co-occurrence network (terms in >={min_articles} GT articles)\n'
                 f'{GT.number_of_nodes()} nodes ({shared_n} are network nodes) - '
                 f'{n_recovered}/{len(net_edges)} network edges recovered', fontsize=12.5)
    ax.axis('off')
    plt.tight_layout()
    out = os.path.join(figures_dir, f"{prefix}_GT_Cooccurrence_Network.png")
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close('all')
    return out


def _plot_nulls(node_null, obs_nodes, z_n, p_n, edge_null, obs_edges, z_e, p_e,
                n_gt, figures_dir, prefix):
    """Observed overlap against the random-article-draw null, for nodes and edges."""
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    panels = [
        (node_null, obs_nodes, z_n, p_n, 'network nodes attested'),
        (edge_null, obs_edges, z_e, p_e, 'network edges recovered'),
    ]
    for a, (vals, obs, z, p, lab) in zip(ax, panels):
        a.hist(vals, bins=30, color='#bdbdbd', edgecolor='w')
        a.axvline(obs, color='#d7301f', lw=2.5)
        a.annotate(f'observed = {obs}\nnull {np.mean(vals):.0f}+/-{np.std(vals):.0f}\n'
                   f'z={z:.1f}, p={p:.3g}',
                   xy=(obs, a.get_ylim()[1] * .9), xytext=(-10, 0),
                   textcoords='offset points', ha='right', fontsize=10,
                   bbox=dict(boxstyle='round', fc='white', ec='#d7301f'))
        a.set_xlabel(lab)
        a.set_ylabel('permutations')
        a.set_title(f'{lab} vs random {n_gt}-article draws')
    plt.suptitle('Convergent validity against a permutation null', fontsize=12)
    plt.tight_layout()
    out = os.path.join(figures_dir, f"{prefix}_GT_Permutation_Nulls.png")
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close('all')
    return out


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ORCHESTRATION
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def run_gt_network_validation(ground_truth_path: str, master_db_path: str,
                              final_network_path: str, output_dir: str,
                              figures_dir: str, file_prefix: str = "network",
                              weight_key: str = "CRS_pagerank_centrality",
                              min_articles_per_node: int = 2,
                              pool_size: int = 50000,
                              n_perm_nodes: int = 1000,
                              n_perm_edges: int = 500,
                              random_seed: int = 42,
                              make_figures: bool = True) -> dict:
    """Validate the network's nodes and edges against a curated ground truth."""
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    rng = random.Random(random_seed)

    print("\n" + "<" * 30 + ">" * 30)
    print("<<< NODE/EDGE CONVERGENT VALIDATION vs GROUND TRUTH >>>")
    print("<" * 30 + ">" * 30)

    # <<< 1. Network >>>
    node_ids, weights, extras, G = _load_network(final_network_path, weight_key)
    net_edges = {frozenset((u, v)) for u, v in G.edges()}
    stopped_nodes = node_ids & MESH_STOP_WORDS
    print(f"  Network ....................... {len(node_ids)} nodes / {len(net_edges)} edges")
    print(f"  Node weight ................... {weight_key}")
    if stopped_nodes:
        print(f"  [!] WARNING: {len(stopped_nodes)} network nodes are stop-listed "
              f"(filter mismatch): {sorted(stopped_nodes)[:5]}")

    # <<< 2. Ground-truth articles -> eligible terms >>>
    gt = validate_ground_truth(ground_truth_path)
    gt_pmids = sorted(int(p) for p in gt["resolved"])
    conn = _open_readonly_resilient(master_db_path)
    cur = conn.cursor()
    gt_docs, gt_docs_raw = [], []
    for pm in gt_pmids:
        row = cur.execute("SELECT mesh_terms FROM master_mesh_annotations WHERE pmid=?",
                          (pm,)).fetchone()
        if row and row[0]:
            gt_docs.append(_eligible_terms(row[0]))
            gt_docs_raw.append(_raw_terms(row[0]))
    conn.close()

    n_gt = len(gt_docs)
    if n_gt == 0:
        raise RuntimeError("No ground-truth articles carried MeSH terms in the master database.")

    freq = Counter()
    for s in gt_docs:
        freq.update(s)
    eligible = set(freq)

    stop_removed = Counter()
    for raw in gt_docs_raw:
        for t in (raw & MESH_STOP_WORDS):
            stop_removed[t] += 1

    print(f"  Ground-truth articles ......... {n_gt} (of {len(gt_pmids)} resolved PMIDs)")
    print(f"  Eligible GT terms ............. {len(eligible)}  "
          f"(stop-listed removed: {len(stop_removed)} distinct)")
    print(f"  Stop-word leak into nodes/misses {len(set(stop_removed) & node_ids)} / "
          f"{len(set(stop_removed) & eligible)}   (both must be 0)")

    # <<< 3. Background pool for base rates + nulls >>>
    print(f"  Sampling {pool_size:,} random articles for base rates and nulls ...")
    pool = _sample_background_pool(master_db_path, pool_size, random_seed)
    pool_freq = Counter()
    for s in pool:
        pool_freq.update(s)
    n_pool = len(pool)

    def corpus_prev(t):
        return pool_freq.get(t, 0.5) / n_pool

    def enrichment(t):
        return round((freq.get(t, 0) / n_gt) / corpus_prev(t), 2) if t in eligible else 0.0

    # <<< 4. Node level >>>
    shared = sorted(eligible & node_ids)
    obs_nodes = len(shared)
    rho = float('nan')
    if len(shared) > 2:
        rho = spearmanr([freq[t] for t in shared], [weights.get(t, 0.0) for t in shared]).correlation

    node_null = []
    for _ in range(n_perm_nodes):
        drawn = set().union(*rng.sample(pool, n_gt)) if n_gt <= n_pool else set()
        node_null.append(len(drawn & node_ids))
    node_null = np.asarray(node_null, dtype=float)
    z_n = (obs_nodes - node_null.mean()) / (node_null.std() + 1e-9)
    p_n = (np.sum(node_null >= obs_nodes) + 1) / (n_perm_nodes + 1)

    print("\n<<< NODE LEVEL >>>")
    print(f"  Network nodes attested in GT .. {obs_nodes}/{len(node_ids)} "
          f"({100 * obs_nodes / max(len(node_ids), 1):.0f}%)")
    print(f"  Permutation null .............. {node_null.mean():.1f} +/- {node_null.std():.1f}"
          f"   z={z_n:.1f}  p={p_n:.4f}")
    print(f"  Spearman(GT freq, node weight)  {rho:+.3f}   (n={len(shared)})")

    # <<< 5. Edge level >>>
    gt_pair = Counter()
    for s in gt_docs:
        for a, b in combinations(sorted(s & node_ids), 2):
            gt_pair[frozenset((a, b))] += 1
    rec1 = sum(1 for e in net_edges if gt_pair.get(e, 0) >= 1)
    rec2 = sum(1 for e in net_edges if gt_pair.get(e, 0) >= 2)

    edge_null = []
    for _ in range(n_perm_edges):
        seen = set()
        for s in rng.sample(pool, n_gt):
            for a, b in combinations(sorted(s & node_ids), 2):
                seen.add(frozenset((a, b)))
        edge_null.append(len(seen & net_edges))
    edge_null = np.asarray(edge_null, dtype=float)
    z_e = (rec1 - edge_null.mean()) / (edge_null.std() + 1e-9)
    p_e = (np.sum(edge_null >= rec1) + 1) / (n_perm_edges + 1)

    print("\n<<< EDGE LEVEL >>>")
    print(f"  Network edges recovered ....... {rec1}/{len(net_edges)} "
          f"({100 * rec1 / max(len(net_edges), 1):.0f}%) at >=1 GT article; "
          f"{rec2} at >=2")
    print(f"  Permutation null .............. {edge_null.mean():.1f} +/- {edge_null.std():.1f}"
          f"   z={z_e:.1f}  p={p_e:.4f}")

    # <<< 6. Ground-truth co-occurrence network >>>
    gt_nodes = [t for t in eligible if freq[t] >= min_articles_per_node]
    co = Counter()
    gt_node_set = set(gt_nodes)
    for s in gt_docs:
        for a, b in combinations(sorted(s & gt_node_set), 2):
            co[frozenset((a, b))] += 1
    GT = nx.Graph()
    GT.add_nodes_from(gt_nodes)
    for e, c in co.items():
        u, v = tuple(e)
        GT.add_edge(u, v, w=c)
    n_comp = nx.number_connected_components(GT) if GT.number_of_nodes() else 0
    lcc = max((len(c) for c in nx.connected_components(GT)), default=0)
    print(f"\n  GT co-occurrence network ...... {GT.number_of_nodes()} nodes / "
          f"{GT.number_of_edges()} edges / {n_comp} component(s), LCC={lcc}")

    # <<< 7. Tables >>>
    gt_terms_df = pd.DataFrame([
        {'term': t, 'GT_articles': c, 'GT_prevalence': round(c / n_gt, 4),
         'in_network': ('Y' if t in node_ids else ''),
         'corpus_prevalence': round(corpus_prev(t), 6), 'enrichment': enrichment(t)}
        for t, c in freq.most_common()
    ])
    misses_df = gt_terms_df[gt_terms_df['in_network'] == ''].copy()
    nodes_df = pd.DataFrame([
        {'term': t, 'in_GT': ('Y' if t in eligible else ''), 'GT_articles': freq.get(t, 0),
         'node_weight': round(weights.get(t, 0.0), 6),
         'betweenness_centrality': round(extras[t]['betweenness_centrality'], 6),
         'network_article_count': extras[t]['article_count'],
         'corpus_prevalence': round(corpus_prev(t), 6), 'enrichment': enrichment(t)}
        for t in sorted(node_ids)
    ]).sort_values(['in_GT', 'GT_articles'], ascending=[False, False])
    edges_df = pd.DataFrame([
        {'source': u, 'target': v,
         'network_cooccurrence_count': G[u][v].get('cooccurrence_count', 0),
         'GT_cooccurrence_count': gt_pair.get(frozenset((u, v)), 0),
         'recovered_in_GT': ('Y' if gt_pair.get(frozenset((u, v)), 0) >= 1 else '')}
        for u, v in G.edges()
    ]).sort_values(['recovered_in_GT', 'GT_cooccurrence_count'], ascending=[False, False])
    stop_df = pd.DataFrame([
        {'stop_term': t, 'GT_occurrences': c, 'in_network': ('Y' if t in node_ids else '')}
        for t, c in stop_removed.most_common()
    ])
    summary_df = pd.DataFrame([
        ('ground_truth_articles', n_gt),
        ('eligible_gt_terms', len(eligible)),
        ('stop_listed_terms_removed', len(stop_removed)),
        ('stopword_leak_nodes_misses', f"{len(set(stop_removed) & node_ids)} / "
                                       f"{len(set(stop_removed) & eligible)}"),
        ('network_nodes', len(node_ids)),
        ('network_edges', len(net_edges)),
        ('nodes_attested_in_gt', obs_nodes),
        ('node_coverage_pct', round(100 * obs_nodes / max(len(node_ids), 1), 1)),
        ('node_null_mean_sd', f"{node_null.mean():.1f} +/- {node_null.std():.1f}"),
        ('node_overlap_z_p', f"{z_n:.1f} / {p_n:.4f}"),
        ('spearman_gtfreq_vs_weight', round(rho, 3) if rho == rho else 'nan'),
        ('edges_recovered_ge1', f"{rec1}/{len(net_edges)}"),
        ('edges_recovered_ge2', f"{rec2}/{len(net_edges)}"),
        ('edge_null_mean_sd', f"{edge_null.mean():.1f} +/- {edge_null.std():.1f}"),
        ('edge_overlap_z_p', f"{z_e:.1f} / {p_e:.4f}"),
        ('gt_network_nodes_edges_lcc', f"{GT.number_of_nodes()} / {GT.number_of_edges()} / {lcc}"),
        ('gt_frequent_misses', len(misses_df)),
    ], columns=['metric', 'value'])

    xlsx_path = os.path.join(output_dir, f"{file_prefix}_gt_network_validation.xlsx")
    with pd.ExcelWriter(xlsx_path, engine='openpyxl') as xw:
        summary_df.to_excel(xw, sheet_name='summary', index=False)
        stop_df.to_excel(xw, sheet_name='stopword_audit', index=False)
        nodes_df.to_excel(xw, sheet_name='network_nodes', index=False)
        gt_terms_df.to_excel(xw, sheet_name='GT_terms', index=False)
        misses_df.to_excel(xw, sheet_name='GT_misses', index=False)
        edges_df.to_excel(xw, sheet_name='network_edge_validation', index=False)
    print(f"\n  [+] Wrote {os.path.basename(xlsx_path)}")

    # <<< 8. Cytoscape.js export of the GT co-occurrence network >>>
    gt_deg = dict(GT.degree())
    cy = {
        'metadata': {
            'description': 'Ground-truth MeSH co-occurrence network annotated with '
                           'overlap against the concept network',
            'generated': datetime.date.today().isoformat(),
            'ground_truth': os.path.basename(ground_truth_path),
            'gt_articles': n_gt,
            'min_articles_per_node': min_articles_per_node,
            'stopword_filter': 'mesh_stop_words.MESH_STOP_WORDS (same list used to build the network)',
            'node_weight_key': weight_key,
            'node_overlap_z': round(float(z_n), 2),
            'edge_overlap_z': round(float(z_e), 2),
            'network_edges_recovered_ge1': rec1,
            'network_edges_total': len(net_edges),
        },
        'elements': {'nodes': [], 'edges': []},
    }
    for n in GT.nodes():
        inn = n in node_ids
        cy['elements']['nodes'].append({'data': {
            'id': n, 'label': n,
            'GT_articles': int(freq[n]), 'GT_prevalence': round(freq[n] / n_gt, 4),
            'in_network': bool(inn),
            'category': 'shared' if inn else 'GT_only_miss',
            'node_weight': round(weights.get(n, 0.0), 6) if inn else None,
            'corpus_prevalence': round(corpus_prev(n), 6),
            'enrichment': enrichment(n),
            'gt_degree': int(gt_deg.get(n, 0)),
        }})
    for u, v, d in GT.edges(data=True):
        inn = frozenset((u, v)) in net_edges
        cy['elements']['edges'].append({'data': {
            'source': u, 'target': v,
            'GT_cooccurrence_count': int(d.get('w', 0)),
            'in_network': bool(inn),
            'edge_type': 'recovered' if inn else 'GT_only',
            'network_cooccurrence_count': (int(G[u][v].get('cooccurrence_count', 0)) if inn else None),
        }})
    cy_path = os.path.join(output_dir, f"{file_prefix}_gt_cooccurrence_network.json")
    with open(cy_path, 'w', encoding='utf-8') as f:
        json.dump(cy, f, indent=1)
    print(f"  [+] Wrote {os.path.basename(cy_path)} (Cytoscape.js)")

    figures = []
    if make_figures:
        try:
            sns.set_context("paper", font_scale=1.15)
            sns.set_style("whitegrid")
            plt.rcParams['pdf.fonttype'] = 42
            figures.append(_plot_node_validation(shared, freq, weights, gt_terms_df,
                                                 rho, n_gt, figures_dir, file_prefix))
            if GT.number_of_nodes():
                figures.append(_plot_gt_network(GT, freq, node_ids, net_edges, rec1,
                                                figures_dir, file_prefix, min_articles_per_node))
            figures.append(_plot_nulls(node_null, obs_nodes, z_n, p_n,
                                       edge_null, rec1, z_e, p_e, n_gt, figures_dir, file_prefix))
            for p in figures:
                print(f"  [+] Wrote {os.path.basename(p)}")
        except Exception as e:
            print(f"  [!] WARNING: figure generation failed: {e}")

    if len(misses_df):
        top_misses = ", ".join(misses_df['term'].head(8))
        print(f"\n  Top ground-truth terms NOT in the network: {top_misses}")

    return {
        'gt_articles': n_gt,
        'nodes_attested': obs_nodes,
        'node_overlap_z': float(z_n),
        'node_overlap_p': float(p_n),
        'spearman_gtfreq_vs_weight': float(rho) if rho == rho else None,
        'edges_recovered_ge1': rec1,
        'edges_recovered_ge2': rec2,
        'edge_overlap_z': float(z_e),
        'edge_overlap_p': float(p_e),
        'workbook': xlsx_path,
        'cytoscape_json': cy_path,
        'figures': figures,
    }
