"""
figures.py

Script used to generate all visualizations for the
MeSH AOP Network analysis.

This script requires access to the following files:
requirements.txt
scripts/python/config.py
data/raw/aop_annotations_master.csv
data/raw/FILE_PREFIX_full_network_data.json
data/processed/FILE_PREFIX_final_network_with_relevance.json
data/processed/FILE_PREFIX_glf_optimal_subgraph.json
data/processed/FILE_PREFIX_sa_optimal_subgraph.json

This script performs:
1. Data Loading & Annotations
2. Statistical Analysis (Dispersion, GLF/SA Simulation Comparisons)
3. Visualization Generation:
   - Figure 1: Edge Weight Distribution
   - Figure 2: Optimization Trajectory Plot (GLF vs SA)
   - Figure 3: Community Composition
   - Figure 4: CRS Correlation Joint Plot
   - Figure 5: t-SNE with Communities
   - Figure 6: AOP Alluvial Flow (Labeled & Unlabeled)
   - Figure 7: Dumbbell Plot
   - Figure 8: Centrality Scatter Plots (CRS vs Raw)
   - Node2Vec Dendrogram (with adjusted font sizes)
   - Sankey Diagrams (Node & Edge Filtering) + Summary Table

This saves finalized figures in:
results/figures/
"""

import os
import sys
import json
import math
import random
from collections import defaultdict
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import seaborn as sns
from scipy.stats import nbinom
from scipy.spatial import ConvexHull
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score, davies_bouldin_score
import statsmodels.api as sm
import plotly.graph_objects as go
from adjustText import adjust_text
from tqdm import tqdm

# Handle optional Node2Vec dependency (this needs to be reinstalled in Colab to work)
try:
    from node2vec import Node2Vec
except ImportError:
    print("Node2Vec not installed. Dendrogram generation will be skipped.")
    Node2Vec = None

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 1. GLOBAL CONFIGURATION & STYLE
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# <<< Import Configuration >>>
# Get the directory where this script/notebook is running
if "__file__" in globals():
    current_dir = os.path.dirname(os.path.abspath(__file__))
else:
    current_dir = os.path.abspath(os.getcwd())

# Add this folder to the system path so Python finds 'config.py' here
if current_dir not in sys.path:
    sys.path.append(current_dir)

import config

# <<< Map Configuration to Local Variables >>>
OUTPUT_DIR = config.FIGURES_DIR
INPUT_JSON_FILE = config.FILES['final_network']
FULL_NETWORK_JSON = config.FILES['full_network']
AOP_ANNOTATION_FILE = config.FILES['annotations']
GLF_SUBGRAPH_PATH = config.FILES['glf_subgraph']
SA_SUBGRAPH_PATH = config.FILES['sa_subgraph']

print(f"Loading Main Network: {INPUT_JSON_FILE}")
print(f"Saving Figures To:    {OUTPUT_DIR}")

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FILE PATHS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
FILE_PREFIX = config.FILE_PREFIX

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CONSTANT VARIABLES
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Plotting Constants
AOP_ORDER = ['Stressor', 'Molecular', 'Cellular', 'Tissue', 'Organ', 'Adverse Outcome', 'Uncategorized']

# Visual Style
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
sns.set_context("paper", font_scale=1.2)
sns.set_style("whitegrid")

# Palettes
MAGMA_PALETTE = sns.color_palette("magma_r", n_colors=len(AOP_ORDER))
AOP_COLOR_MAP = dict(zip(AOP_ORDER, MAGMA_PALETTE))

# Markers (Ensure enough for all levels)
MARKERS_LIST = ['o', 's', '^', 'D', 'v', 'X', 'P']

# <<< Simulation Parameters >>>
TARGET_NUM_EDGES = config.TARGET_NUM_EDGES
RANDOM_SEED      = config.RANDOM_SEED
GLF_ITERATIONS   = config.GLF_ITERATIONS
SA_ITERATIONS    = config.SA_ITERATIONS
SA_INITIAL_TEMP  = config.SA_TEMP_START
SA_COOLING_RATE  = config.SA_COOLING_RATE

# Set Seeds
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

# Ensure Output Directory Exists
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Output directory set to: {OUTPUT_DIR}")

# <<< Hardcoded Data for Sankey (Filtering Summary) >>>
SANKEY_DATA = {
    'Initial Subgraph': {'GLF': {'Nodes': 295, 'Edges': 500}, 'SA': {'Nodes': 292, 'Edges': 500}},
    'LCC': {'GLF': {'Nodes': 211, 'Edges': 451}, 'SA': {'Nodes': 217, 'Edges': 456}},
    'Final Consensus LCC': {'GLF': {'Nodes': 175, 'Edges': 336}, 'SA': {'Nodes': 175, 'Edges': 336}}
}


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 2. DATA LOADING HELPER FUNCTIONS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def load_and_prepare_data(json_path, annotation_path):
    """Loads final network data and merges AOP annotations."""
    print(f"\n<<< Loading Final Network Data: {os.path.basename(json_path)} >>>")
    if not os.path.exists(json_path):
        print(f"ERROR: File not found: {json_path}")
        return None, None, None

    with open(json_path, 'r') as f:
        data = json.load(f)

    # Process Nodes
    nodes_data = [n['data'] for n in data['elements']['nodes']]
    node_df = pd.DataFrame(nodes_data)
    node_df.rename(columns={
        'id': 'mesh_term',
        'adjusted_node_weight': 'ARS',
        'CRS_betweenness_centrality': 'CRS_betweenness',
        'CRS_eigenvector_centrality': 'CRS_eigenvector'
    }, inplace=True)
    node_df.set_index('mesh_term', inplace=True)

    # Process Edges & Graph
    edges_data = [e['data'] for e in data['elements']['edges']]
    edge_df = pd.DataFrame(edges_data)
    G = nx.Graph()
    for node_id, attrs in node_df.iterrows():
        G.add_node(node_id, **attrs)
    for idx, edge in edge_df.iterrows():
        G.add_edge(edge['source'], edge['target'], **edge)

    # Merge Annotations
    if os.path.exists(annotation_path):
        print(f"Merging annotations from {os.path.basename(annotation_path)}...")
        try:
            aop_df = pd.read_csv(annotation_path)

            # CHECK: If using new data, ensure annotations actually exist
            if not getattr(config, 'USE_REFERENCE_DATA', True):
                if aop_df.empty or 'aop_level' not in aop_df.columns or aop_df['aop_level'].dropna().empty:
                    print("\n" + "!*"*30)
                    print("[!] CRITICAL WARNING: Biological Strata Annotations Missing/Empty")
                    print(f"    File: {annotation_path}")
                    print("    You are running a NEW analysis, but 'aop_level' assignments seem to be missing.")
                    print("    Figures 3, 4, 6, 7, & 8 will be incorrect or uncolored.")
                    print("    Please manually assign levels (Molecular, Cellular, etc.) in the CSV.")
                    print("!*"*30 + "\n")

            node_df = node_df.reset_index().merge(aop_df, on='mesh_term', how='left').set_index('mesh_term')

        except Exception as e:
            print(f"WARNING: Could not read annotation file: {e}")
            node_df['aop_level'] = 'Uncategorized'

        node_df['aop_level'] = pd.Categorical(node_df['aop_level'], categories=AOP_ORDER, ordered=True)
        if 'Uncategorized' in AOP_ORDER:
             node_df['aop_level'] = node_df['aop_level'].fillna('Uncategorized')

        nx.set_node_attributes(G, node_df['aop_level'].to_dict(), 'aop_level')
    else:
        print("WARNING: Annotation file not found.")
        node_df['aop_level'] = 'Uncategorized'

    print(f"Graph Loaded: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges.")
    return node_df, edge_df, G

def load_full_raw_data(filepath):
    """Loads raw Cytoscape JSON for statistical comparison (Pre-filtering)."""
    print(f"Loading raw network data: {os.path.basename(filepath)}...")
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)

        edges = [e['data'] for e in data['elements']['edges']]
        edge_df = pd.DataFrame(edges)

        all_edges_data = {}
        for item in data.get('elements', {}).get('edges', []):
            ed = item.get('data', {})
            if ed.get('source') and ed.get('target'):
                key = tuple(sorted((ed['source'], ed['target'])))
                all_edges_data[key] = ed

        return edge_df, all_edges_data
    except FileNotFoundError:
        print(f"Error: {filepath} not found.")
        return None, None

def save_high_res(filename_base):
    """Helper to save JPEG preview and TIFF high-res."""
    plt.savefig(os.path.join(OUTPUT_DIR, f"{config.FILE_PREFIX}_{filename_base}.jpeg"),
                dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, f"{config.FILE_PREFIX}_{filename_base}.tif"),
                dpi=600, pil_kwargs={'compression': 'tiff_lzw'}, bbox_inches='tight')
    print(f"Saved figures: {filename_base}")

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 3. ANALYSIS FUNCTIONS (Dispersion, GLF/SA)
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# This code is repeated again from the Master_MeSH file to generate a plot

def analyze_dispersion(edge_df, output_dir):
    """Performs Negative Binomial Regression to assess over-dispersion."""
    print("\n<<< Running Dispersion Analysis (HNB) >>>")
    counts = edge_df['cooccurrence_count'].dropna()

    # Plot Distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(counts, bins=100, kde=False, color='#2c3e50')
    plt.yscale('log')
    plt.title('Edge Weight Distribution (Log Scale)', fontweight='bold')
    plt.xlabel('Co-occurrence Count')
    plt.ylabel('Frequency')
    save_high_res("distribution_plot")
    plt.close()

    # Regression
    df_reg = pd.DataFrame({'counts': counts, 'intercept': 1.0})
    try:
        model = sm.NegativeBinomial(df_reg['counts'], df_reg[['intercept']]).fit(disp=0)
        print(model.summary())
        alpha = model.params.get('alpha', model.params.iloc[1] if len(model.params)>1 else 0)
        print(f"Dispersion Parameter (Alpha): {alpha:.4f}")
    except Exception as e:
        print(f"Regression failed: {e}")

# <<< Simulation Helpers >>>
def calculate_graph_stats(edges_data):
    node_strengths = defaultdict(int)
    total_strength = 0
    for edge_key, data in edges_data.items():
        w = data.get('cooccurrence_count', 0)
        node_strengths[edge_key[0]] += w
        node_strengths[edge_key[1]] += w
        total_strength += w
    return total_strength, node_strengths

def get_log_likelihood_term(edge_data, node_strengths, denominator):
    w_ij = edge_data.get('cooccurrence_count', 0)
    k_i = node_strengths.get(edge_data.get('source'), 0)
    k_j = node_strengths.get(edge_data.get('target'), 0)
    if k_i == 0 or k_j == 0 or w_ij == 0: return 0
    p_ij = (k_i * k_j) / denominator
    return w_ij * np.log(p_ij) - math.lgamma(w_ij + 1) if p_ij > 0 else 0

def run_simulation(method, all_edges, node_strengths, total_T, target_edges, iterations, **kwargs):
    """Unified runner for GLF or SA simulations."""
    print(f"Starting {method} Simulation ({iterations:,} iters)...")

    denominator = 2 * (total_T**2)
    keys = list(all_edges.keys())
    current_keys = set(random.sample(keys, target_edges))

    curr_T = sum(all_edges[k].get('cooccurrence_count', 0) for k in current_keys)
    curr_term = sum(get_log_likelihood_term(all_edges[k], node_strengths, denominator) for k in current_keys)
    curr_ll = math.lgamma(curr_T + 1) + curr_term

    best_keys, min_ll = current_keys.copy(), curr_ll
    history = []

    temp = kwargs.get('initial_temp', 0)
    cooling = kwargs.get('cooling_rate', 0)

    for i in tqdm(range(iterations), desc=f"{method}"):
        on = random.choice(list(current_keys))
        off = random.choice(keys)
        while off in current_keys: off = random.choice(keys)

        w_on = all_edges[on].get('cooccurrence_count', 0)
        w_off = all_edges[off].get('cooccurrence_count', 0)
        prop_T = curr_T - w_on + w_off

        ll_on = get_log_likelihood_term(all_edges[on], node_strengths, denominator)
        ll_off = get_log_likelihood_term(all_edges[off], node_strengths, denominator)

        delta = (math.lgamma(prop_T + 1) - math.lgamma(curr_T + 1)) + (ll_off - ll_on)

        # fixed GLF accept probability
        accept = False
        if delta < 0:

            accept = True
        elif method == 'GLF' and random.random() < math.exp(-delta):
            accept = True
        elif method == 'SA' and temp > 1e-9 and random.random() < math.exp(-delta / temp):
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


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 4. PLOTTING FUNCTIONS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def plot_cooccurrance_distribution(full_edge_df, filtered_edge_df):
    print("\n<<< Generating Figure 1: Edge Weight Distribution >>>")
    plt.figure(figsize=(12, 8))

    global_max = full_edge_df['cooccurrence_count'].max()
    bins = np.arange(0, global_max + 25, 25)

    ax = sns.histplot(data=full_edge_df, x='cooccurrence_count', bins=bins,
                     color='gray', label=f'Before Filtering ({len(full_edge_df):,} edges)')
    sns.histplot(data=filtered_edge_df, x='cooccurrence_count', bins=bins,
                color=MAGMA_PALETTE[0], alpha=0.7,
                label=f'After Filtering ({len(filtered_edge_df):,} edges)', ax=ax)

    ax.set_yscale('log')
    plt.title('Effect of Filtering on Edge Co-occurrence Count Distribution', fontsize=16)
    plt.xlabel('Co-occurrence Count (Linear)', fontsize=12)
    plt.ylabel('Frequency (Log)', fontsize=12)
    plt.xlim(0, 2000)
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    save_high_res("EdgeWeightDistribution_Comparison")
    plt.close()

def run_optimization_comparison(raw_all_edges):
    """
    Runs GLF and SA simulations, plots the trajectory comparison,
    and prints the metrics table with F1 and Jaccard scores.
    """
    print("\n<<< Running GLF vs SA Simulation Comparison >>>")
    tot_str, node_strs = calculate_graph_stats(raw_all_edges)

    # 1. Run Simulations
    glf_keys, glf_sc, glf_hist = run_simulation('GLF', raw_all_edges, node_strs, tot_str, TARGET_NUM_EDGES, GLF_ITERATIONS)
    sa_keys, sa_sc, sa_hist = run_simulation('SA', raw_all_edges, node_strs, tot_str, TARGET_NUM_EDGES, SA_ITERATIONS,
                                             initial_temp=SA_INITIAL_TEMP, cooling_rate=SA_COOLING_RATE)

    # 2. Plot Trajectory
    plt.figure(figsize=(10, 6))
    plt.plot(*zip(*glf_hist), label='GLF')
    plt.plot(*zip(*sa_hist), label='SA')

    # Add Final Score Lines
    plt.axhline(y=glf_sc, color='C0', linestyle='--', label=f'GLF Final: {glf_sc:,.2f}')
    plt.axhline(y=sa_sc, color='C1', linestyle='--', label=f'SA Final: {sa_sc:,.2f}')

    plt.title('Optimization Trajectory: GLF vs SA (Likelihood Score)')
    plt.xlabel('Iterations')
    plt.ylabel('Log-Likelihood Score')
    plt.legend()

    print("\n<<< Generating Optimization Trajectory Plot (GLF vs SA) >>>")
    traj_path = os.path.join(OUTPUT_DIR, "GLF_SA_Trajectory.png")
    plt.savefig(traj_path, dpi=600)
    plt.show()
    plt.close()

    # 3. Generate Metrics Table (with F1 Score)
    print("\n<<< Generating Comparison Metrics Table >>>")

    glf_set = set(glf_keys)
    sa_set = set(sa_keys)
    consensus_set = glf_set.intersection(sa_set)
    union_set = glf_set.union(sa_set)

    # Calculations
    jaccard = len(consensus_set) / len(union_set) if len(union_set) > 0 else 0
    total_size = len(glf_set) + len(sa_set)
    f1_score = (2 * len(consensus_set)) / total_size if total_size > 0 else 0

    G_glf = nx.Graph(list(glf_set))
    G_sa = nx.Graph(list(sa_set))

    table_data = {
        "Metric": [
            "Final Log-Likelihood",
            "Number of Nodes",
            "Number of Edges",
            "Network Density",
            "Avg. Clustering Coeff.",
            "Jaccard Similarity",
            "F1 Score (Overlap)"
        ],
        "GLF Result": [
            f"{glf_sc:,.2f}",
            f"{G_glf.number_of_nodes():,}",
            f"{G_glf.number_of_edges():,}",
            f"{nx.density(G_glf):.4f}",
            f"{nx.average_clustering(G_glf):.4f}",
            f"{jaccard:.4f}",
            f"{f1_score:.4f}"
        ],
        "SA Result": [
            f"{sa_sc:,.2f}",
            f"{G_sa.number_of_nodes():,}",
            f"{G_sa.number_of_edges():,}",
            f"{nx.density(G_sa):.4f}",
            f"{nx.average_clustering(G_sa):.4f}",
            f"{jaccard:.4f}",
            f"{f1_score:.4f}"
        ]
    }

    print("\n" + "="*60)
    print("      Comparison of Resulting Subgraph Metrics")
    print("="*60)
    print(pd.DataFrame(table_data).to_string(index=False))
    print("="*60 + "\n")

def plot_louvain_community_bars(node_df):
    print("\n<<< Generating Figure 3: Community Composition >>>")
    if 'filtered_louvain_community_id' not in node_df.columns: return

    cross_tab = pd.crosstab(node_df['filtered_louvain_community_id'], node_df['aop_level'])
    cross_tab = cross_tab.reindex(columns=AOP_ORDER, fill_value=0)

    # Sort
    sorted_idx = cross_tab.sum(axis=1).sort_values(ascending=False).index
    cross_tab = cross_tab.loc[sorted_idx]

    # Plot
    ax = cross_tab.plot(kind='bar', stacked=True, figsize=(14, 8),
                        color=[AOP_COLOR_MAP[c] for c in cross_tab.columns],
                        alpha=0.7, edgecolor='black', linewidth=0.5, width=0.8)

    ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.grid(False)
    ax.tick_params(axis='x', length=5)
    ax.tick_params(axis='y', length=5)
    sns.despine(ax=ax, left=False, bottom=False)

    plt.xlabel('Louvain Community ID', fontsize=12)
    plt.ylabel('Number of Nodes', fontsize=12)
    plt.legend(title=r"$\bf{AOP\ Level}$", bbox_to_anchor=(0.81, 1), loc='upper left')
    plt.subplots_adjust(right=0.75)
    plt.tight_layout(rect=[0, 0, 0.85, 1])

    save_high_res("Louvain_Community_Composition")
    plt.close()

def plot_joint_plot(node_df):
    print("\n<<< Generating Figure 4: CRS Correlation Joint Plot >>>")
    fig = plt.figure(figsize=(12, 12))
    gs = fig.add_gridspec(2, 2, width_ratios=(4, 1), height_ratios=(1, 4),
                          left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.05, hspace=0.05)

    ax_scat = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_scat)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax_scat)

    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # Scatter
    sns.scatterplot(data=node_df, x='CRS_betweenness', y='CRS_eigenvector',
                    hue='aop_level', style='aop_level', markers=MARKERS_LIST,
                    palette='magma_r', hue_order=AOP_ORDER, s=60, alpha=0.7, ax=ax_scat)
    ax_scat.grid(True)
    ax_scat.legend(title=r"$\bf{AOP\ Level}$", bbox_to_anchor=(1.02, 1), loc='lower left')
    ax_scat.set_xlabel('CRS (Betweenness)', fontsize=12)
    ax_scat.set_ylabel('CRS (Eigenvector)', fontsize=12)

    # Hists
    sns.histplot(data=node_df, x='CRS_betweenness', hue='aop_level', multiple="stack",
                 palette='magma_r', hue_order=AOP_ORDER, legend=False, ax=ax_histx, bins=10)
    sns.histplot(data=node_df, y='CRS_eigenvector', hue='aop_level', multiple="stack",
                 palette='magma_r', hue_order=AOP_ORDER, legend=False, ax=ax_histy, bins=10)

    # Label top nodes
    top = pd.concat([node_df.nlargest(5, 'CRS_betweenness'), node_df.nlargest(5, 'CRS_eigenvector')])
    texts = [ax_scat.text(row['CRS_betweenness'], row['CRS_eigenvector'], idx, fontsize=9)
             for idx, row in top.iterrows()]
    adjust_text(texts, ax=ax_scat, arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

    save_high_res("Joint_Plot_Final")
    plt.close()

def plot_tsne_louvain_overlap(node_df, G):
    print("\n<<< Generating Figure 5: t-SNE with Communities >>>")
    # Path matrix
    try:
        # Try using weighted degree if available, else unweighted
        degree_dict = dict(G.degree(weight='cooccurrence_count'))
    except Exception:
        degree_dict = dict(G.degree())
        print("\n Weighted degree not available, proceeding with unweighted degree")

    path_len = dict(nx.all_pairs_shortest_path_length(G))
    nodes = list(G.nodes())
    matrix = pd.DataFrame(index=nodes, columns=nodes, dtype=float)
    for u in nodes:
        for v in nodes: matrix.loc[u, v] = path_len.get(u, {}).get(v, np.nan)
    matrix.fillna(matrix.max().max() * 2, inplace=True)

    # TSNE
    tsne = TSNE(n_components=2, perplexity=30, metric='precomputed', init='random', random_state=42)
    emb = tsne.fit_transform(matrix.values)
    df_tsne = pd.DataFrame(emb, index=nodes, columns=['X', 'Y'])
    df_tsne = df_tsne.join(node_df[['filtered_louvain_community_id', 'CRS_betweenness']])

    # Plot
    plt.figure(figsize=(12, 10))
    ax = plt.gca()

    unique_comms = sorted(df_tsne['filtered_louvain_community_id'].unique())
    comm_palette = sns.color_palette('viridis', n_colors=len(unique_comms))
    comm_map = dict(zip(unique_comms, comm_palette))

    # Community with top node legend
    community_labels = {}
    for comm_id in unique_comms:
        subset = df_tsne[df_tsne['filtered_louvain_community_id'] == comm_id]
        if not subset.empty:
            # Find representative node (highest degree in graph)
            nodes_in_comm = subset.index.tolist()
            top_node = max(nodes_in_comm, key=lambda x: degree_dict.get(x, 0))
            if pd.notna(node_df.loc[top_node, 'aop_level']):
                aop_tag = node_df.loc[top_node, 'aop_level']
            else:
                aop_tag = "?"

            label = f"{comm_id}: {top_node} ({aop_tag})"
        else:
            label = f"{comm_id}: (Empty)"
        community_labels[comm_id] = label

    # Hulls
    for cid in unique_comms:
        pts = df_tsne[df_tsne['filtered_louvain_community_id'] == cid][['X', 'Y']].values
        if len(pts) >= 3:
            try:
                hull = ConvexHull(pts)
                poly = plt.Polygon(pts[hull.vertices], closed=True, facecolor=comm_map[cid],
                                  alpha=0.1, edgecolor=comm_map[cid], linestyle='--')
                ax.add_patch(poly)
            except:
                  #Straight line = Zero area
                  continue

    sns.scatterplot(data=df_tsne, x='X', y='Y', hue='filtered_louvain_community_id',
                    palette=comm_map, size='CRS_betweenness', sizes=(20, 300), alpha=0.9, legend=False, ax=ax)

    # Custom Legend
    legend_patches = [mpatches.Patch(color=comm_map[cid], label=community_labels[cid], alpha=0.9)
                      for cid in unique_comms]
    ax.legend(handles=legend_patches,
                  title='Louvain Communities & Topics',
                  loc='upper right',
                  frameon=True,
                  facecolor='white',
                  edgecolor='black',
                  framealpha=0.95,
                  fontsize=10,
                  title_fontsize=11)
    ax.set_xlim(right=15)
    plt.xlabel("t-SNE Dimension 1")
    plt.ylabel("t-SNE Dimension 2")

    save_high_res("T-SNE_Final")
    plt.close()

def plot_sankey_alluvial(G, node_df):
    print("\n<<< Figure 6: AOP Alluvial Flow (Labeled & Unlabeled) >>>")

    # Filter only relevant levels
    levels = [l for l in AOP_ORDER if l != 'Uncategorized']
    lvl_map = {l: i for i, l in enumerate(levels)}

    flows = defaultdict(int)
    for u, v in G.edges():
        lu, lv = node_df.loc[u, 'aop_level'], node_df.loc[v, 'aop_level']
        if lu in lvl_map and lv in lvl_map:
            u_i, v_i = lvl_map[lu], lvl_map[lv]
            if u_i <= v_i: flows[(u_i, v_i)] += 1
            else: flows[(v_i, u_i)] += 1

    # <<< ADDED: Export CSV of relationships >>>
    csv_data = []
    for (src_idx, tgt_idx), count in flows.items():
        csv_data.append({
            "Source": levels[src_idx],
            "Target": levels[tgt_idx],
            "Weight": count
        })

    csv_path = os.path.join(OUTPUT_DIR, "Alluvial_Connections_Table_F6.csv")
    pd.DataFrame(csv_data).sort_values(['Source', 'Target']).to_csv(csv_path, index=False)
    print(f"Exported connection data to: {csv_path}")

    # Define the config for SVG export
    config = {
        'toImageButtonOptions': {
            'format': 'svg', # one of png, svg, jpeg, webp
            'filename': 'custom_sankey_layout',
            'height': 700,
            'width': 1200,
            'scale': 1
        }
    }

    # <<< PLOT 1: Labeled (For manual reference) >>>
    fig_labels = go.Figure(go.Sankey(
        node=dict(
            pad=25, thickness=20, line=dict(color="black", width=0.5),
            label=levels,  # SHOW LABELS
            color=sns.color_palette("magma_r", len(levels)).as_hex()
        ),
        link=dict(
            source=[k[0] for k in flows], target=[k[1] for k in flows],
            value=list(flows.values()),
            label=[str(v) for v in flows.values()] # SHOW COUNTS
        )
    ))
    fig_labels.update_layout() #LABELED - Reference Only
    fig_labels.show(config=config)

    # <<< PLOT 2: Clean (Unlabeled - For publication) >>>
    fig_clean = go.Figure(go.Sankey(
        node=dict(
            pad=25, thickness=20, line=dict(color="black", width=0.5),
            label=[""] * len(levels),  # HIDE LABELS
            color=sns.color_palette("magma_r", len(levels)).as_hex()
        ),
        link=dict(
            source=[k[0] for k in flows], target=[k[1] for k in flows],
            value=list(flows.values()),
            label=[""] * len(flows)    # HIDE COUNTS
        )
    ))
    fig_clean.update_layout() # Clean NO LABELS for figure and ps
    fig_clean.show(config=config)

def plot_dumbell_plot(node_df):
    print("\n<<< Figure 7: Dumbbell Plot >>>")
    df_plot = node_df.copy()
    df_plot['diff'] = (df_plot['CRS_betweenness'] - df_plot['CRS_eigenvector']).abs()
    top_diff = df_plot.nlargest(20, 'diff').sort_values('CRS_betweenness')

    plt.figure(figsize=(10, 8))
    plt.hlines(y=top_diff.index, xmin=top_diff['CRS_betweenness'], xmax=top_diff['CRS_eigenvector'],
               color='gray', alpha=0.5, linewidth=1.5, zorder=1)

    # Neutral colors for method dots
    colors = ["#2C3E50", "#B0B0B0"]
    plt.scatter(top_diff['CRS_betweenness'], top_diff.index, color=colors[0], s=100, zorder=2, label='Betweenness')
    plt.scatter(top_diff['CRS_eigenvector'], top_diff.index, color=colors[1], s=100, zorder=2, label='Eigenvector')

    # Color Y-labels
    ax = plt.gca()
    for lbl in ax.get_yticklabels():
        if lbl.get_text() in top_diff.index:
            lvl = top_diff.loc[lbl.get_text(), 'aop_level']
            if pd.notna(lvl) and lvl in AOP_COLOR_MAP:
                lbl.set_color(AOP_COLOR_MAP[lvl])
                lbl.set_fontweight('bold')

    # Legend
    leg_h = ([mpatches.Patch(color='none', label=r"$\bf{AOP\ Level}$")] +
             [mpatches.Patch(color=c, label=l) for l, c in AOP_COLOR_MAP.items()] +
             [mpatches.Patch(color='none', label="")] +
             [mpatches.Patch(color='none', label=r"$\bf{Method}$")] +
             [mpatches.Patch(color=colors[0], label='Betweenness'), mpatches.Patch(color=colors[1], label='Eigenvector')])

    plt.legend(handles=leg_h, loc='lower right', bbox_to_anchor=(1.0, 0.0))
    plt.xlabel("Contextual Relevance Score (CRS)")
    plt.tight_layout()
    save_high_res("CRS_Dumbbell_Plot_Top20")
    plt.close()

def plot_scatter_panels(node_df):
    print("\n<<< Figure 8: Centrality Scatter Plots (CRS vs Raw) >>>")

    # PANEL C: Betweenness
    print("  - Generating Panel C: Raw Betweenness vs CRS (Betweenness)...")
    plot_df_b = node_df[(node_df['betweenness_centrality'] > 0) & (node_df['CRS_betweenness'] > 0)]

    g_b = sns.lmplot(data=plot_df_b, x='betweenness_centrality', y='CRS_betweenness',
                     hue='aop_level', palette='magma_r', hue_order=AOP_ORDER,
                     height=8, aspect=1.1, ci=None, legend=False,
                     markers=MARKERS_LIST, scatter_kws={'s': 50, 'alpha': 0.7})

    ax_b = g_b.ax
    ax_b.set_xscale('log')
    ax_b.set_yscale('log')
    ax_b.set_xlabel('Raw Betweenness Centrality (Log Scale)', fontsize=12)
    ax_b.set_ylabel('CRS (Betweenness) (Log Scale)', fontsize=12)
    ax_b.grid(True, which="both", ls="--", alpha=0.5)

    output_file_b = os.path.join(OUTPUT_DIR, "Panel_C_Scatter_Betweenness.tif")
    plt.savefig(output_file_b, dpi=1200, pil_kwargs={'compression': 'tiff_lzw'}, bbox_inches='tight')
    print(f"    Saved: {output_file_b}")
    plt.close()

    # PANEL D: Eigenvector
    print("  - Generating Panel D: Raw Eigenvector vs CRS (Eigenvector)...")
    plot_df_e = node_df[(node_df['eigenvector_centrality'] > 0) & (node_df['CRS_eigenvector'] > 0)]

    g_e = sns.lmplot(data=plot_df_e, x='eigenvector_centrality', y='CRS_eigenvector',
                     hue='aop_level', palette='magma_r', hue_order=AOP_ORDER,
                     height=8, aspect=1.1, ci=None, legend=False,
                     markers=MARKERS_LIST, scatter_kws={'s': 50, 'alpha': 0.7})

    ax_e = g_e.ax
    ax_e.set_xscale('log')
    ax_e.set_yscale('log')
    ax_e.set_xlabel('Raw Eigenvector Centrality (Log Scale)', fontsize=12)
    ax_e.set_ylabel('CRS (Eigenvector) (Log Scale)', fontsize=12)
    ax_e.grid(True, which="both", ls="--", alpha=0.5)

    output_file_e = os.path.join(OUTPUT_DIR, "Panel_D_Scatter_Eigenvector.tif")
    plt.savefig(output_file_e, dpi=600, pil_kwargs={'compression': 'tiff_lzw'}, bbox_inches='tight')
    print(f"    Saved: {output_file_e}")
    plt.close()

def plot_dendrogram(G, node_df):
    if Node2Vec is None: return
    print("\n<<< Generating Node2Vec Dendrogram >>>")

    n2v = Node2Vec(G, dimensions=64, walk_length=30, num_walks=200, workers=4, quiet=True)
    model = n2v.fit(window=10, min_count=1)

    labels = model.wv.index_to_key
    vectors = model.wv.vectors
    Z = linkage(vectors, method='ward', metric='euclidean')

    plt.figure(figsize=(16, max(15, len(labels)*0.25)))

    dendrogram(Z, labels=labels, orientation='right',
               color_threshold=0, above_threshold_color='black',
               leaf_font_size=12)

    ax = plt.gca()
    for lbl in ax.get_yticklabels():
        term = lbl.get_text()
        if term in node_df.index:
            lvl = node_df.loc[term, 'aop_level']
            if lvl in AOP_COLOR_MAP: lbl.set_color(AOP_COLOR_MAP[lvl])

    aop_legend_patches = [mpatches.Patch(color=color, label=level) for level, color in AOP_COLOR_MAP.items()]
    ax.legend(
        handles=aop_legend_patches,
        title='AOP Level (Label Color)',
        loc='upper left',
        bbox_to_anchor=(0.82, 1.0)
    )
    plt.xlabel("Cluster Distance", fontsize=12) #Ward method
    plt.ylabel("MeSH Terms", fontsize=12)
    plt.grid(axis='x', linestyle='--', alpha=0.5)

    save_high_res("Dendrogram_Node2Vec_AOP")
    plt.close()

def generate_filtering_summary_sankeys(data):
    """Generates Filtering Sankeys (Nodes & Edges) and summary table."""
    print("\n<<< Generating Filtering Cascade Sankeys >>>")

    def make_sankey(item_type, output_name):
        labels = [
            "GLF Initial", "SA Initial", "GLF LCC", "SA LCC",
            "Final Consensus LCC", "Discarded (Not in LCC)", "Discarded (Not in Consensus)"
        ]
        sources = [0, 0, 1, 1, 2, 2, 3, 3]
        targets = [2, 5, 3, 5, 4, 6, 4, 6]

        glf_i = data['Initial Subgraph']['GLF'][item_type]
        sa_i = data['Initial Subgraph']['SA'][item_type]
        glf_lcc = data['LCC']['GLF'][item_type]
        sa_lcc = data['LCC']['SA'][item_type]
        final = data['Final Consensus LCC']['GLF'][item_type]

        values = [
            glf_lcc, glf_i - glf_lcc,
            sa_lcc, sa_i - sa_lcc,
            final, glf_lcc - final,
            final, sa_lcc - final
        ]

        # Colors: Blue (GLF), Orange (SA), Purple (Final), Grey (Discarded)
        link_colors = [
            'rgba(55, 126, 184, 0.6)', 'rgba(200, 200, 200, 0.4)',
            'rgba(255, 127, 0, 0.6)', 'rgba(200, 200, 200, 0.4)',
            'rgba(55, 126, 184, 0.6)', 'rgba(200, 200, 200, 0.4)',
            'rgba(255, 127, 0, 0.6)', 'rgba(200, 200, 200, 0.4)'
        ]
        node_colors = ["#377eb8", "#ff7f00", "#377eb8", "#ff7f00", "#984ea3", "grey", "grey"]

        fig = go.Figure(data=[go.Sankey(
            node=dict(pad=25, thickness=20, line=dict(color="black", width=0.5),
                      label=[f"{l}\n({v:,})" for l, v in zip(labels,
                             [glf_i, sa_i, glf_lcc, sa_lcc, final, (glf_i-glf_lcc)+(sa_i-sa_lcc), (glf_lcc-final)+(sa_lcc-final)])],
                      color=node_colors),
            link=dict(source=sources, target=targets, value=values, color=link_colors)
        )])
        fig.update_layout(title_text=f"Filtering Cascade: {item_type}", font_size=12)
        fig.show()

    make_sankey('Nodes', 'Node_Sankey')
    make_sankey('Edges', 'Edge_Sankey')

    # Print Summary Table
    print("\n" + "="*80)
    print("                Summary of Network Filtering Cascade")
    print("="*80)

    # Map raw data to table rows
    rows = []
    stages = [('Initial Subgraph', 'Initial Subgraph'),
              ('LCC', 'Largest Connected Component'),
              ('Final Consensus LCC', 'Final Consensus Largest Connected Component')]

    for key, name in stages:
        r = {'Stage': name}
        for method in ['GLF', 'SA']:
            for item in ['Nodes', 'Edges']:
                r[f"{method} {item}"] = data[key][method][item]
        rows.append(r)

    df = pd.DataFrame(rows).set_index('Stage')

    # Retention Calcs
    initial = 'Initial Subgraph'
    lcc = 'Largest Connected Component'
    final = 'Final Consensus Largest Connected Component'

    df['GLF Node Retention'] = ''
    df['SA Node Retention'] = ''

    # LCC vs Initial
    df.loc[lcc, 'GLF Node Retention'] = f"{df.loc[lcc, 'GLF Nodes']/df.loc[initial, 'GLF Nodes']:.1%}"
    df.loc[lcc, 'SA Node Retention'] = f"{df.loc[lcc, 'SA Nodes']/df.loc[initial, 'SA Nodes']:.1%}"

    # Final vs LCC
    df.loc[final, 'GLF Node Retention'] = f"{df.loc[final, 'GLF Nodes']/df.loc[lcc, 'GLF Nodes']:.1%}"
    df.loc[final, 'SA Node Retention'] = f"{df.loc[final, 'SA Nodes']/df.loc[lcc, 'SA Nodes']:.1%}"

    print(df[['GLF Nodes', 'GLF Edges', 'GLF Node Retention',
              'SA Nodes', 'SA Edges', 'SA Node Retention']].to_string())
    print("-" * 80 + "\n")


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 5. MAIN EXECUTION BLOCK
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if __name__ == "__main__":
    print(f"<<< Starting Analysis Pipeline >>>\n")

    # 1. LOAD MAIN DATA
    node_df, edge_df, G = load_and_prepare_data(INPUT_JSON_FILE, AOP_ANNOTATION_FILE)
    if node_df is None:
        print("Critical Error: Main data not found. Exiting.")
        sys.exit()

    # 2. LOAD RAW DATA (For Pre-filtering Stats)
    raw_edge_df, raw_all_edges = load_full_raw_data(FULL_NETWORK_JSON)

    # 3. DISPERSION & FILT DIST (Figure 1A)
    if raw_edge_df is not None:
        analyze_dispersion(raw_edge_df, OUTPUT_DIR)
        plot_cooccurrance_distribution(raw_edge_df, edge_df)

    # 4. GLF vs SA COMPARISON (Trajectory only - Venn Removed)
    if raw_all_edges:
        run_optimization_comparison(raw_all_edges)

    # 5. GENERATE FILTERING CASCADE SUMMARY (Sankeys + Table)
    generate_filtering_summary_sankeys(SANKEY_DATA)

    # 6. GENERATE ALL OTHER VISUALIZATIONS
    plot_louvain_community_bars(node_df)
    plot_joint_plot(node_df)
    plot_tsne_louvain_overlap(node_df, G)
    plot_sankey_alluvial(G, node_df)
    plot_dumbell_plot(node_df)
    plot_scatter_panels(node_df)
    plot_dendrogram(G, node_df)

    print("\n<<< Pipeline Completed Successfully >>>")

