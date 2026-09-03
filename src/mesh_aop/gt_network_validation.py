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
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap

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
    """Load node ids, the chosen node weight, all node attributes, and the edge set."""
    with open(final_network_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    nodes = data.get('elements', {}).get('nodes', [])
    node_ids, weights, extras, node_attrs = set(), {}, {}, {}
    carried = 0
    for n in nodes:
        d = n.get('data', {})
        nid = d.get('id')
        if not nid:
            continue
        node_ids.add(nid)
        node_attrs[nid] = d
        if weight_key in d:
            carried += 1
        weights[nid] = float(d.get(weight_key, 0.0) or 0.0)

    # A weight key the network does not carry reads as 0.0 for every node, and
    # the run then completes with a full set of numbers scored against a
    # constant. That is worse than stopping, because nothing looks wrong. Say
    # so, loudly, and name what the network does have.
    if node_ids and carried == 0:
        available = sorted({k for d in node_attrs.values() for k, v in d.items()
                            if isinstance(v, (int, float)) and not isinstance(v, bool)})
        raise ValueError(
            f"The network carries no '{weight_key}'. Every node would score "
            f"zero and the validation below would be meaningless.\n"
            f"  This network offers: {', '.join(available) or '(no numeric attributes)'}\n"
            f"  Networks built before the subgraph centralities became "
            f"unconditional lack the six *_subgraph_centrality attributes; "
            f"rebuild the network step once to add them.")
    if node_ids and carried < len(node_ids):
        print(f"  [!] WARNING: only {carried} of {len(node_ids)} nodes carry "
              f"'{weight_key}'; the rest are scored as zero.")
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

    return node_ids, weights, extras, G, node_attrs


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# NODE-WEIGHTING COMPARISON
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Candidate node weightings, paired with the raw centrality each is derived from.
# An MRS is a transformation of its centrality, so the question is never "which
# correlates better" in isolation but "does the transformation add anything the
# centrality did not already carry" - hence the paired baseline.
_WEIGHT_CANDIDATES = [
    ("MRS_betweenness_centrality", "betweenness_centrality"),
    ("MRS_pagerank_centrality", "pagerank_centrality"),
    ("MRS_betweenness_subgraph_centrality", "betweenness_subgraph_centrality"),
    ("MRS_pagerank_subgraph_centrality", "pagerank_subgraph_centrality"),
    ("betweenness_centrality", None),
    ("pagerank_centrality", None),
    ("betweenness_subgraph_centrality", None),
    ("pagerank_subgraph_centrality", None),
    ("eigenvector_centrality", None),
    ("adjusted_node_weight", None),
]


def _auc_present_vs_absent(scores: np.ndarray, present: np.ndarray) -> float:
    """Probability a ground-truth-attested node outranks an unattested one (Mann-Whitney)."""
    n_pos = int(present.sum())
    n_neg = int(len(present) - n_pos)
    if n_pos == 0 or n_neg == 0:
        return float('nan')
    from scipy.stats import rankdata
    r = rankdata(scores)
    return float((r[present == 1].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def _partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    """Spearman correlation of x and y after linearly removing z, on ranks.

    Answers the incremental question: once the raw centrality (z) is accounted
    for, does the transformed weight (x) still track ground-truth prominence (y)?
    """
    from scipy.stats import rankdata
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    A = np.column_stack([np.ones_like(rz), rz])
    rx_res = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ry_res = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    if rx_res.std() == 0 or ry_res.std() == 0:
        return float('nan')
    return float(np.corrcoef(rx_res, ry_res)[0, 1])


def _compare_node_weightings(node_attrs: dict, freq: Counter, n_gt: int,
                             corpus_prev, n_boot: int, seed: int) -> pd.DataFrame:
    """Score every available node weighting against external ground-truth prominence.

    Three complementary criteria, because no single one is sufficient:
      * GT frequency  - how often ground-truth articles use the term (shared nodes);
      * GT enrichment - the same, divided by the corpus base rate, which removes any
        residual frequency component and so is the stricter test;
      * AUC           - whether the weighting ranks attested nodes above unattested
        ones across ALL nodes, using the terms the ground truth never mentions too.
    Paired weightings additionally get a partial correlation controlling for the raw
    centrality, and a bootstrap CI on the difference between the two.
    """
    rng = np.random.default_rng(seed)
    ids = sorted(node_attrs)
    present = np.array([1 if freq.get(t, 0) > 0 else 0 for t in ids])
    gt_freq = np.array([freq.get(t, 0) for t in ids], dtype=float)
    enrich = np.array([(freq.get(t, 0) / n_gt) / corpus_prev(t) for t in ids])
    shared = present == 1

    def col(key):
        """Column vector of one node attribute over the shared-node id list."""
        return np.array([float(node_attrs[t].get(key, 0.0) or 0.0) for t in ids])

    rows = []
    for key, base_key in _WEIGHT_CANDIDATES:
        if not any(key in node_attrs[t] for t in ids):
            continue
        w = col(key)
        if w.std() == 0:
            continue

        rho_f = spearmanr(w[shared], gt_freq[shared]).correlation if shared.sum() > 2 else float('nan')
        rho_e = spearmanr(w[shared], enrich[shared]).correlation if shared.sum() > 2 else float('nan')
        auc = _auc_present_vs_absent(w, present)

        partial = diff_lo = diff_hi = float('nan')
        if base_key and any(base_key in node_attrs[t] for t in ids):
            b = col(base_key)
            if b.std() > 0:
                partial = _partial_spearman(w[shared], gt_freq[shared], b[shared])
                # Paired bootstrap of the correlation difference: the two weightings
                # are measured on the same nodes, so they must be resampled together.
                idx_pool = np.flatnonzero(shared)
                diffs = []
                for _ in range(n_boot):
                    s = rng.choice(idx_pool, size=len(idx_pool), replace=True)
                    if gt_freq[s].std() == 0 or w[s].std() == 0 or b[s].std() == 0:
                        continue
                    d = (spearmanr(w[s], gt_freq[s]).correlation
                         - spearmanr(b[s], gt_freq[s]).correlation)
                    if d == d:
                        diffs.append(d)
                if diffs:
                    diff_lo, diff_hi = np.percentile(diffs, [2.5, 97.5])

        rows.append({
            'weighting': key,
            'vs_raw_centrality': base_key or '',
            'spearman_gt_freq': round(rho_f, 4) if rho_f == rho_f else None,
            'spearman_gt_enrichment': round(rho_e, 4) if rho_e == rho_e else None,
            'auc_attested_vs_not': round(auc, 4) if auc == auc else None,
            'partial_spearman_given_raw': round(partial, 4) if partial == partial else None,
            'delta_vs_raw_ci95_lo': round(diff_lo, 4) if diff_lo == diff_lo else None,
            'delta_vs_raw_ci95_hi': round(diff_hi, 4) if diff_hi == diff_hi else None,
        })

    return pd.DataFrame(rows)


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

# Full paper colour palette (GIMP Fig1.gpl export), kept here so any figure can
# reference the exact publication colours later. Order follows the .gpl file.
FIG1_PALETTE = [
    '#000000', '#003A46', '#00819C', '#1A1A1A', '#30526C', '#353564', '#39D0EE',
    '#3C476C', '#4D4D4D', '#4D4D9F', '#55306C', '#666666', '#6CE9FF', '#6DDBC2',
    '#85DB6D', '#8686BF', '#999999', '#AFAFDE', '#BB6DDB', '#C863B1', '#D7D7FF',
    '#D8DA6C', '#D9F4F9', '#DAF4F9', '#DB6D78', '#DB6DC0', '#DBAC6D', '#E9E9FF',
    '#FFFFFF',
]

# Semantic role -> colour for the three ground-truth validation figures, all taken
# from FIG1_PALETTE. Change a role here to recolour every GT figure at once.
GT_C_ATTESTED = '#00819C'   # teal      - in-network node / recovered edge
GT_C_MISS     = '#DB6D78'   # coral     - ground-truth-only node / network miss
GT_C_GT_EDGE  = '#DB6D78'   # coral     - ground-truth-only edge (coordinated with the miss nodes)
GT_C_NULL     = '#8686BF'   # purple    - permutation-null histogram (lighter)
GT_C_OBSERVED = '#D62728'   # red       - observed-statistic line (deliberately NOT from
                            #             Fig1.gpl: a strong red the eye can't mistake for
                            #             the right axis spine, since the observed value sits
                            #             far to the right of the null)
# Sequential ramp for the GT-frequency scatter: slate -> teal -> gold (low -> high).
GT_FREQ_CMAP  = LinearSegmentedColormap.from_list('gt_freq', ['#3C476C', '#00819C', '#DBAC6D'])


def _adjust_labels(texts, ax):
    """De-overlap a set of matplotlib text labels, if adjustText is installed.

    adjustText is already a visualisation dependency; if it is somehow missing the
    labels are simply left in place rather than failing the whole figure.
    """
    if not texts:
        return
    try:
        import io
        from contextlib import redirect_stdout
        from adjustText import adjust_text
        with redirect_stdout(io.StringIO()):
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#999', lw=.5))
    except Exception:
        pass


def _plot_node_validation(shared, freq, weights, gt_terms_df, rho, n_gt, figures_dir, prefix):
    """Scatter of ground-truth prominence vs network weight, plus the top-term bar."""
    fig, ax = plt.subplots(1, 2, figsize=(15.5, 6.6))

    xs = [freq[t] for t in shared]
    ys = [weights.get(t, 0.0) for t in shared]
    sc = ax[0].scatter(xs, ys, s=60, c=[freq[t] / n_gt for t in shared], cmap=GT_FREQ_CMAP,
                       edgecolors='k', linewidths=.3, zorder=3)
    if any(y > 0 for y in ys):
        ax[0].set_yscale('log')
    ax[0].set_xlabel(f'Ground-truth article frequency (of {n_gt})')
    ax[0].set_ylabel('Network node weight (log)')
    ax[0].set_title(f'Shared nodes: GT prominence vs network weight\n'
                    f'Spearman rho = {rho:+.2f}  (n={len(shared)})')
    cbar = fig.colorbar(sc, ax=ax[0], fraction=.046, pad=.02)
    cbar.set_label('GT frequency (fraction of articles)', fontsize=9)
    # Label only the ten most prominent shared terms, de-overlapped.
    texts = [ax[0].text(freq[t], weights.get(t, 0.0), t, fontsize=7.5)
             for t in sorted(shared, key=lambda x: -freq[x])[:10]]
    _adjust_labels(texts, ax[0])
    ax[0].grid(True, ls='--', alpha=.35)

    top = gt_terms_df.head(22)[::-1]
    cols = [GT_C_ATTESTED if v == 'Y' else GT_C_MISS for v in top['in_network']]
    ax[1].barh(range(len(top)), top['GT_articles'], color=cols)
    ax[1].set_yticks(range(len(top)))
    ax[1].set_yticklabels(top['term'], fontsize=8)
    ax[1].set_ylim(-0.6, len(top) - 0.4)
    ax[1].set_xlabel('Ground-truth article frequency')
    ax[1].set_title('Top eligible ground-truth terms')
    ax[1].legend(handles=[Patch(facecolor=GT_C_ATTESTED, label='network node'),
                          Patch(facecolor=GT_C_MISS, label='network miss')],
                 loc='lower right', fontsize=9, framealpha=.95)
    ax[1].grid(True, axis='x', ls='--', alpha=.35)

    plt.tight_layout()
    out = os.path.join(figures_dir, f"{prefix}_GT_Node_Validation.png")
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close('all')
    return out


def _plot_gt_network(GT, freq, node_ids, net_edges, n_recovered, figures_dir, prefix, min_articles):
    """Draw the ground-truth co-occurrence network, coloured by overlap with ours."""
    plt.figure(figsize=(14, 11.5))
    ax = plt.gca()
    # A slightly larger k spreads a dense (~100-node) graph so labels/edges breathe.
    pos = nx.spring_layout(GT, seed=42, k=0.75, iterations=300, weight='w')

    ncol = [GT_C_ATTESTED if n in node_ids else GT_C_MISS for n in GT.nodes()]
    nsize = [120 + freq[n] * 55 for n in GT.nodes()]
    ecol, ew = [], []
    for u, v in GT.edges():
        hit = frozenset((u, v)) in net_edges
        ecol.append(GT_C_ATTESTED if hit else GT_C_GT_EDGE)
        ew.append(1.7 if hit else 0.7)

    nx.draw_networkx_edges(GT, pos, edge_color=ecol, width=ew, alpha=.55, ax=ax)
    nx.draw_networkx_nodes(GT, pos, node_color=ncol, node_size=nsize, alpha=.92,
                           linewidths=.5, edgecolors='#333', ax=ax)
    # Label only the most prominent terms (top 22, freq>=3) with a white halo so a
    # dense graph stays readable instead of a wall of overlapping text.
    to_label = {n for n in sorted(GT.nodes(), key=lambda x: -freq[x])[:22] if freq[n] >= 3}
    lbls = nx.draw_networkx_labels(GT, pos, labels={n: n for n in to_label},
                                   font_size=7.5, ax=ax)
    for t in lbls.values():
        t.set_bbox(dict(facecolor='white', alpha=.55, edgecolor='none', pad=.4))

    ax.legend(handles=[
        Line2D([0], [0], marker='o', color='w', markerfacecolor=GT_C_ATTESTED, markersize=11,
               label='term IS a network node'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=GT_C_MISS, markersize=11,
               label='ground-truth-only term (network miss)'),
        Line2D([0], [0], color=GT_C_ATTESTED, lw=2.2, label='edge present in our network (recovered)'),
        Line2D([0], [0], color=GT_C_GT_EDGE, lw=1.5, label='ground-truth-only co-occurrence'),
    ], loc='upper left', fontsize=9, framealpha=.95, edgecolor='#cccccc').set_zorder(5)

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
        a.hist(vals, bins=30, color=GT_C_NULL, edgecolor='w', label='permutation null')
        a.axvline(obs, color=GT_C_OBSERVED, lw=2.5, label='observed')
        a.annotate(f'observed = {obs}\nnull = {np.mean(vals):.0f} +/- {np.std(vals):.0f} '
                   f'(mean +/- SD)\nz={z:.1f}, p={p:.3g}',
                   xy=(obs, a.get_ylim()[1] * .78), xytext=(-12, 0),
                   textcoords='offset points', ha='right', fontsize=9.5,
                   bbox=dict(boxstyle='round', fc='white', ec=GT_C_OBSERVED, alpha=.95))
        a.set_xlabel(lab)
        a.set_ylabel('permutations')
        a.set_title(f'{lab} vs random {n_gt}-article draws')
        a.legend(loc='upper left', fontsize=9, framealpha=.95)
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
                              weight_key: str = "MRS_pagerank_centrality",
                              min_articles_per_node: int = 2,
                              pool_size: int = 50000,
                              n_perm_nodes: int = 1000,
                              n_perm_edges: int = 500,
                              n_boot_weighting: int = 2000,
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
    node_ids, weights, extras, G, node_attrs = _load_network(final_network_path, weight_key)
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
        """A term's prevalence in the random background pool (fraction of its articles)."""
        return pool_freq.get(t, 0.5) / n_pool

    def enrichment(t):
        """GT frequency of a term over its background prevalence (base-rate-corrected)."""
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

    # <<< 5b. Which node weighting best tracks ground-truth prominence? >>>
    # An MRS is derived from its centrality, so the two are collinear and comparing
    # their raw correlations alone would be misleading. The partial correlation and
    # the paired bootstrap difference are the parts that actually answer whether the
    # transformation contributes anything the centrality did not already carry.
    weight_cmp = _compare_node_weightings(node_attrs, freq, n_gt, corpus_prev,
                                          n_boot_weighting, random_seed)
    if len(weight_cmp):
        print("\n<<< NODE WEIGHTING vs GROUND-TRUTH PROMINENCE >>>")
        print("  (rho_freq = vs GT article frequency; rho_enrich = base-rate corrected;")
        print("   AUC = attested vs unattested nodes; partial = controlling for raw centrality)")
        for _, r in weight_cmp.iterrows():
            base = f" | vs {r['vs_raw_centrality']}" if r['vs_raw_centrality'] else ""
            line = (f"  {r['weighting']:<34} rho_freq {str(r['spearman_gt_freq']):>7}"
                    f"  rho_enrich {str(r['spearman_gt_enrichment']):>7}"
                    f"  AUC {str(r['auc_attested_vs_not']):>6}")
            if r['partial_spearman_given_raw'] is not None:
                line += (f"  partial {r['partial_spearman_given_raw']:+.3f}"
                         f"  delta95 [{r['delta_vs_raw_ci95_lo']:+.3f}, {r['delta_vs_raw_ci95_hi']:+.3f}]")
            print(line + base)
        print("  A delta interval spanning 0 means the transformation is not")
        print("  distinguishable from the centrality it was built from.")

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

    # <<< 6b. Set-overlap between the two networks >>>
    # Precision / recall / F1 with the ground-truth network as the reference and the
    # consensus network as the prediction. All three are reported rather than F1
    # alone: the two networks are built from corpora differing by orders of
    # magnitude, so a consensus term the ground-truth articles never mention is not
    # necessarily an error. Precision and recall keep that asymmetry visible where a
    # single harmonic mean would hide it. Jaccard is included because it is the
    # measure already used elsewhere for the GLF/SA intersection.
    def _prf(n_shared: int, n_pred: int, n_ref: int):
        """Precision (of the prediction), recall (of the reference), and their F1."""
        p = n_shared / n_pred if n_pred else 0.0
        r = n_shared / n_ref if n_ref else 0.0
        f = (2 * p * r / (p + r)) if (p + r) else 0.0
        return round(p, 4), round(r, 4), round(f, 4)

    gt_edge_set = {frozenset(e) for e in GT.edges()}
    shared_nodes_net = node_ids & gt_node_set
    shared_edges_net = net_edges & gt_edge_set
    p_nn, r_nn, f_nn = _prf(len(shared_nodes_net), len(node_ids), len(gt_node_set))
    p_ne, r_ne, f_ne = _prf(obs_nodes, len(node_ids), len(eligible))
    p_ee, r_ee, f_ee = _prf(len(shared_edges_net), len(net_edges), len(gt_edge_set))
    union_n = node_ids | gt_node_set
    union_e = net_edges | gt_edge_set
    jacc_n = round(len(shared_nodes_net) / len(union_n), 4) if union_n else 0.0
    jacc_e = round(len(shared_edges_net) / len(union_e), 4) if union_e else 0.0

    print("\n  <<< NETWORK OVERLAP >>>  (reference = GT network, prediction = consensus)")
    print(f"    nodes vs GT network ({len(gt_node_set)})".ljust(38)
          + f"shared {len(shared_nodes_net):4d}   P {p_nn:.3f}  R {r_nn:.3f}  "
            f"F1 {f_nn:.3f}  Jaccard {jacc_n:.3f}")
    print(f"    nodes vs eligible GT terms ({len(eligible)})".ljust(38)
          + f"shared {obs_nodes:4d}   P {p_ne:.3f}  R {r_ne:.3f}  F1 {f_ne:.3f}")
    print(f"    edges vs GT network ({len(gt_edge_set)})".ljust(38)
          + f"shared {len(shared_edges_net):4d}   P {p_ee:.3f}  R {r_ee:.3f}  "
            f"F1 {f_ee:.3f}  Jaccard {jacc_e:.3f}")

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
        # Overlap of the two networks, one scalar per row so each value can be
        # referenced directly. Precision is over the consensus network, recall
        # over the ground-truth reference; both are reported because the corpora
        # differ in size by orders of magnitude (see section 6b).
        ('node_overlap_shared', len(shared_nodes_net)),
        ('node_overlap_precision', p_nn),
        ('node_overlap_recall', r_nn),
        ('node_overlap_f1', f_nn),
        ('node_overlap_jaccard', jacc_n),
        ('node_overlap_precision_vs_eligible_terms', p_ne),
        ('node_overlap_recall_vs_eligible_terms', r_ne),
        ('node_overlap_f1_vs_eligible_terms', f_ne),
        ('edge_overlap_shared', len(shared_edges_net)),
        ('edge_overlap_precision', p_ee),
        ('edge_overlap_recall', r_ee),
        ('edge_overlap_f1', f_ee),
        ('edge_overlap_jaccard', jacc_e),
    ], columns=['metric', 'value'])

    xlsx_path = os.path.join(output_dir, f"{file_prefix}_gt_network_validation.xlsx")
    with pd.ExcelWriter(xlsx_path, engine='openpyxl') as xw:
        summary_df.to_excel(xw, sheet_name='summary', index=False)
        stop_df.to_excel(xw, sheet_name='stopword_audit', index=False)
        nodes_df.to_excel(xw, sheet_name='network_nodes', index=False)
        gt_terms_df.to_excel(xw, sheet_name='GT_terms', index=False)
        misses_df.to_excel(xw, sheet_name='GT_misses', index=False)
        edges_df.to_excel(xw, sheet_name='network_edge_validation', index=False)
        if len(weight_cmp):
            weight_cmp.to_excel(xw, sheet_name='node_weighting_comparison', index=False)
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
        'node_weighting_comparison': weight_cmp.to_dict('records') if len(weight_cmp) else [],
    }
