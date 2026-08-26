# -*- coding: utf-8 -*-
"""
viz.py - figure generation and visualization (pipeline Step 4).

Produces the publication-quality figures and interactive exports from the
finished, AOP-annotated network.

Outputs include edge-weight and dispersion (zero-truncated negative binomial)
distributions, the GLF-versus-SA optimization trajectory, Louvain community
composition bars, centrality correlation and dumbbell plots, a t-SNE community
projection, a Node2Vec dendrogram, and the Sankey/alluvial flows that connect
stressors to adverse outcomes. Static images are written as high-resolution
JPEG/TIFF and interactive figures as self-contained HTML.
"""

import os
import json
import hashlib
import warnings
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import seaborn as sns
import scipy.stats as stats
from scipy.spatial import ConvexHull
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.manifold import TSNE
import statsmodels.api as sm
from statsmodels.discrete.truncated_model import TruncatedLFNegativeBinomialP
import plotly.graph_objects as go
from adjustText import adjust_text
from collections import defaultdict
import io
from contextlib import redirect_stdout

# Relative imports from our package
from .stats import calculate_graph_stats
from . import paths as _paths

# Node2Vec embedding for the dendrogram figure. We use a vendored, dependency-light
# implementation (node2vec_embedding) instead of the `node2vec` package, so the
# pipeline carries none of that package's import constraints - its edges.py imports
# pkg_resources (removed in setuptools 81). The vendored module needs only gensim,
# networkx and numpy, and works with modern versions of each. Capture any failure
# so the dendrogram step can explain a skip instead of silently omitting the figure.
try:
    from .node2vec_embedding import Node2Vec
    _NODE2VEC_IMPORT_ERROR = None
except Exception as e:
    Node2Vec = None
    _NODE2VEC_IMPORT_ERROR = e

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PLOTTING CONSTANTS & STYLE
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
AOP_ORDER = ['Stressor', 'Molecular', 'Cellular', 'Tissue', 'Organ', 'Adverse Outcome', 'Uncategorized']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
sns.set_context("paper", font_scale=1.2)
sns.set_style("whitegrid")
MAGMA_PALETTE = sns.color_palette("magma_r", n_colors=len(AOP_ORDER))
AOP_COLOR_MAP = dict(zip(AOP_ORDER, MAGMA_PALETTE))
MARKERS_LIST = ['o', 's', '^', 'D', 'v', 'X', 'P']

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DATA IO & HELPER FUNCTIONS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Figure output, overridable from the settings form. Set by configure_output()
# before any figure is drawn; the defaults stand if nothing calls it, so the
# module still behaves when driven directly from a notebook.
_FIGURE_DPI = 300
_FIGURE_FORMATS = ('jpeg', 'tif')

# TIFF is written with LZW because journals ask for lossless at submission and
# an uncompressed 600 dpi panel runs to tens of megabytes.
_FORMAT_KWARGS = {'tif': {'pil_kwargs': {'compression': 'tiff_lzw'}},
                  'tiff': {'pil_kwargs': {'compression': 'tiff_lzw'}}}


def configure_output(dpi=None, formats=None):
    """Set the resolution and file types every figure is written at.

    Called once from the CLI with the user's settings. Kept as module state
    rather than threaded through twenty plotting functions, each of which would
    have to pass it on unchanged.
    """
    global _FIGURE_DPI, _FIGURE_FORMATS
    if dpi:
        try:
            _FIGURE_DPI = max(50, min(1200, int(dpi)))
        except (TypeError, ValueError):
            pass
    if formats:
        if isinstance(formats, str):
            formats = [f.strip().lower().lstrip('.') for f in formats.split(',')]
        keep = [f for f in formats if f in ('jpeg', 'jpg', 'png', 'tif', 'tiff', 'pdf', 'svg')]
        if keep:
            _FIGURE_FORMATS = tuple(keep)
    return _FIGURE_DPI, _FIGURE_FORMATS


def save_high_res(filename_base: str, output_dir: str, file_prefix: str):
    """Write the current figure in every configured format."""
    os.makedirs(_paths.long_path(output_dir), exist_ok=True)
    written = []
    for ext in _FIGURE_FORMATS:
        path = os.path.join(output_dir, f"{file_prefix}_{filename_base}.{ext}")
        try:
            # long_path: a redirected Documents folder plus a descriptive figure
            # name clears 260 characters easily, and Windows then refuses the
            # write with a "file not found" naming a path that exists.
            plt.savefig(_paths.long_path(path), dpi=_FIGURE_DPI, bbox_inches='tight',
                        **_FORMAT_KWARGS.get(ext, {}))
            written.append(ext.upper())
        except Exception as e:
            print(f"    [!] Error saving {filename_base}.{ext}: {e}")
    if written:
        print(f"    Saved figures: {filename_base} ({', '.join(written)} at {_FIGURE_DPI} dpi)")

# <<< Download resolution for the interactive HTML figures >>>
# Plotly's camera button captures the plot at its ON-SCREEN pixel size unless it
# is told otherwise, so a figure that looks sharp in the browser downloads at
# roughly screen resolution (~96 ppi) and is unusable in print. The export is
# pinned here instead: an explicit width/height in CSS pixels makes the exported
# canvas deterministic regardless of the viewer's window size, and `scale`
# multiplies the raster on top of it.
#
# Sankey node coordinates are normalised, so a reader who drags nodes into the
# arrangement they want before pressing download keeps that arrangement - it is
# re-rendered at the larger size rather than recomputed.
HTML_EXPORT_DPI = 600
HTML_EXPORT_SIZE_IN = (11.0, 7.5)   # intended print size of the downloaded image
_CSS_PPI = 96                        # the reference pixels-per-inch Plotly assumes


def export_plotly_figure(fig, filename_base: str, output_dir: str, file_prefix: str,
                         export_dpi: int = HTML_EXPORT_DPI,
                         export_size_in: tuple = HTML_EXPORT_SIZE_IN):
    """
    Bypasses the broken Kaleido static exporter entirely to avoid Python 3.12
    compatibility hangs. Outputs a highly stable, interactive HTML dashboard instead.

    The embedded download button is configured to emit `export_size_in` inches at
    `export_dpi`; at the defaults that is 11.0 x 7.5 in at 600 dpi (6600 x 4500 px).
    """
    os.makedirs(output_dir, exist_ok=True)
    html_path = os.path.join(output_dir, f"{file_prefix}_{filename_base}.html")

    w_px = int(round(export_size_in[0] * _CSS_PPI))
    h_px = int(round(export_size_in[1] * _CSS_PPI))
    scale = export_dpi / _CSS_PPI
    config = {
        "toImageButtonOptions": {
            "format": "png",
            "filename": f"{file_prefix}_{filename_base}",
            "width": w_px,
            "height": h_px,
            "scale": scale,
        },
        "displaylogo": False,
        "responsive": True,
    }

    try:
        # HTML generation relies on pure Python string formatting (100% stable on Py 3.12)
        fig.write_html(_paths.long_path(html_path), config=config)
        print(f"    [+] Saved interactive HTML figure: {os.path.basename(html_path)}"
              f"  (downloads at {export_size_in[0]:g}x{export_size_in[1]:g} in / "
              f"{export_dpi} dpi = {int(w_px * scale):,}x{int(h_px * scale):,} px)")
    except Exception as e:
        print(f"    [!] HTML export failed: {e}")

def load_and_prepare_data(json_path: str, annotation_path: str):
    """Loads final network data and merges AOP annotations via Semicolon CSV."""
    print(f"\n<<< Loading Final Network Data: {os.path.basename(json_path)} >>>")
    if not os.path.exists(json_path):
        raise FileNotFoundError(f"Final network file not found: {json_path}")

    with open(json_path, 'r') as f:
        data = json.load(f)

    nodes_data = [n['data'] for n in data['elements']['nodes']]
    node_df = pd.DataFrame(nodes_data)
    node_df.rename(columns={
        'id': 'mesh_term',
        'adjusted_node_weight': 'ARS',
        'MRS_betweenness_centrality': 'MRS_betweenness',
        'MRS_pagerank_centrality': 'MRS_pagerank'
    }, inplace=True)
    node_df.set_index('mesh_term', inplace=True)

    edges_data = [e['data'] for e in data['elements']['edges']]
    edge_df = pd.DataFrame(edges_data)
    G = nx.Graph()
    for node_id, attrs in node_df.iterrows():
        G.add_node(node_id, **attrs)
    for idx, edge in edge_df.iterrows():
        G.add_edge(edge['source'], edge['target'], **edge)

    if os.path.exists(annotation_path):
        print(f"Merging annotations from {os.path.basename(annotation_path)}...")
        try:
            # Read explicitly with semicolon
            try:
                aop_df = pd.read_csv(annotation_path, sep=';')
                if len(aop_df.columns) == 1 and ',' in aop_df.columns[0]:
                    aop_df = pd.read_csv(annotation_path, sep=',')
            except Exception:
                aop_df = pd.read_csv(annotation_path, sep=',')

            # Standardize column casing
            col_map = {c: c.lower() for c in aop_df.columns}
            aop_df.rename(columns=col_map, inplace=True)

            node_df = node_df.reset_index().merge(aop_df, on='mesh_term', how='left').set_index('mesh_term')
        except Exception as e:
            print(f"[!] WARNING: Could not read annotation file: {e}")
            node_df['aop_level'] = 'Uncategorized'

        node_df['aop_level'] = pd.Categorical(node_df['aop_level'], categories=AOP_ORDER, ordered=True)
        if 'Uncategorized' in AOP_ORDER:
            node_df['aop_level'] = node_df['aop_level'].fillna('Uncategorized')
        nx.set_node_attributes(G, node_df['aop_level'].to_dict(), 'aop_level')
    else:
        print(f"[!] WARNING: Annotation file not found at {annotation_path}. Proceeding 'Uncategorized'.")
        node_df['aop_level'] = 'Uncategorized'

    print(f"Graph Loaded: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges.")
    return node_df, edge_df, G

def load_full_raw_data(filepath: str):
    """Loads raw Cytoscape JSON for statistical comparison."""
    print(f"Loading raw network data: {os.path.basename(filepath)}...")
    if not os.path.exists(filepath):
        print(f"[!] WARNING: {filepath} not found. Pre-filtering statistics will be skipped.")
        return None, None
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
    except Exception as e:
        print(f"[!] ERROR loading raw data: {e}")
        return None, None


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ANALYSIS & PLOTTING FUNCTIONS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def analyze_dispersion(edge_df: pd.DataFrame, output_dir: str, file_prefix: str):
    """Performs Zero-Truncated Negative Binomial Regression and plots distribution."""
    print("\n<<< Running Dispersion Analysis (ZTNB) >>>")
    try:
        counts = edge_df['cooccurrence_count'].dropna()
        plt.figure(figsize=(10, 6))
        sns.histplot(counts, bins=100, kde=False, color='#2c3e50')
        plt.yscale('log')
        plt.title('Edge Weight Distribution (Log Scale)', fontweight='bold')
        plt.xlabel('Co-occurrence Count')
        plt.ylabel('Frequency')
        save_high_res("distribution_plot", output_dir, file_prefix)
        plt.close('all')

        N_POP = len(counts)
        N_SAMPLE = 10000

        if N_POP > N_SAMPLE:
            print(f"    [!] Notice: Computing ZTNB parameters using a random sample of {N_SAMPLE:,} edges (from {N_POP:,}).")
            counts_fit = counts.sample(n=N_SAMPLE, random_state=42)
        else:
            counts_fit = counts

        df_reg = pd.DataFrame({'counts': counts_fit, 'intercept': 1.0})
        intercept_start = np.log(counts_fit.mean())
        alpha_start = 0.1

        original_filters = warnings.filters[:]
        warnings.simplefilter("ignore")
        os.environ['PYTHONWARNINGS'] = 'ignore'

        try:
            with np.errstate(all='ignore'):
                model = TruncatedLFNegativeBinomialP(
                    df_reg['counts'],
                    df_reg[['intercept']],
                    truncation=0
                ).fit(
                    start_params=[intercept_start, alpha_start],
                    method='lbfgs', maxiter=200, pgtol=1e-3, disp=0
                )
            print(model.summary())
            log_intercept = model.params['intercept']
            latent_mean = np.exp(log_intercept)
            alpha = model.params.get('alpha', model.params.iloc[1] if len(model.params)>1 else 0)

            p_zero = (1 + alpha * latent_mean) ** (-1 / alpha) if alpha > 0 else np.exp(-latent_mean)
            real_intercept = latent_mean / (1 - p_zero)

            print(f"\n--- REAL-WORLD METRICS ---")
            print(f"Calculated Average Edge Weight: {real_intercept:.4f} connections")
            print(f"Dispersion Parameter (Alpha):   {alpha:.4f}")

        finally:
            warnings.filters = original_filters
            if 'PYTHONWARNINGS' in os.environ:
                del os.environ['PYTHONWARNINGS']

    except Exception as e:
        print(f"[!] ZTNB regression failed: {e}")
        plt.close('all')

def plot_cooccurrance_distribution(full_edge_df: pd.DataFrame, filtered_edge_df: pd.DataFrame, output_dir: str, file_prefix: str):
    """Figure 1: overlaid co-occurrence count distributions before vs after filtering."""
    print("\n<<< Generating Figure 1: Edge Weight Distribution >>>")
    try:
        plt.figure(figsize=(12, 8))
        global_max = full_edge_df['cooccurrence_count'].max()
        bins = np.arange(0, global_max + 25, 25)

        ax = sns.histplot(data=full_edge_df, x='cooccurrence_count', bins=bins, color='gray', label=f'Before Filtering ({len(full_edge_df):,} edges)')
        sns.histplot(data=filtered_edge_df, x='cooccurrence_count', bins=bins, color=MAGMA_PALETTE[0], alpha=0.7, label=f'After Filtering ({len(filtered_edge_df):,} edges)', ax=ax)

        ax.set_yscale('log')
        plt.title('Effect of Filtering on Edge Co-occurrence Count Distribution', fontsize=16)
        plt.xlabel('Co-occurrence Count (Linear)', fontsize=12)
        plt.ylabel('Frequency (Log)', fontsize=12)
        plt.xlim(0, 2000)
        plt.legend()
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        save_high_res("EdgeWeightDistribution_Comparison", output_dir, file_prefix)
    except Exception as e:
        print(f"[!] Error generating Figure 1: {e}")
    finally:
        plt.close('all')

def run_optimization_comparison(history_json_path: str, output_dir: str, file_prefix: str):
    """Parses precomputed simulation history and generates comparative trajectory."""
    print("\n<<< Generating Optimization Trajectory Plot (GLF vs SA) >>>")
    try:
        if not os.path.exists(history_json_path):
            print(f"    [!] History JSON not found: {os.path.basename(history_json_path)}. Optimization trajectory generation bypassed.")
            return

        with open(history_json_path, 'r') as f:
            hist_data = json.load(f)

        glf_hist = hist_data.get('GLF', [])
        sa_hist = hist_data.get('SA', [])

        if not glf_hist or not sa_hist:
            print("    [!] Missing arrays in optimization history structure.")
            return

        glf_sc = glf_hist[-1][1]
        sa_sc = sa_hist[-1][1]

        plt.figure(figsize=(10, 6))
        plt.plot(*zip(*glf_hist), label='GLF')
        plt.plot(*zip(*sa_hist), label='SA')
        plt.axhline(y=glf_sc, color='C0', linestyle='--', label=f'GLF Final: {glf_sc:,.2f}')
        plt.axhline(y=sa_sc, color='C1', linestyle='--', label=f'SA Final: {sa_sc:,.2f}')
        plt.title('Optimization Trajectory: GLF vs SA (Likelihood Score)')
        plt.xlabel('Iterations')
        plt.ylabel('Log-Likelihood Score')
        plt.legend()

        save_high_res("GLF_SA_Trajectory", output_dir, file_prefix)

    except Exception as e:
        print(f"[!] Error in Simulation Comparison processing: {e}")
    finally:
        plt.close('all')

def plot_louvain_community_bars(node_df: pd.DataFrame, output_dir: str, file_prefix: str):
    """Figure 3: stacked bars of AOP-level composition within each Louvain community."""
    print("\n<<< Generating Figure 3: Community Composition >>>")
    try:
        if 'filtered_louvain_community_id' not in node_df.columns:
            print("    [!] Missing 'filtered_louvain_community_id' column. Skipping Figure 3.")
            return

        cross_tab = pd.crosstab(node_df['filtered_louvain_community_id'], node_df['aop_level'])
        cross_tab = cross_tab.reindex(columns=AOP_ORDER, fill_value=0)
        sorted_idx = cross_tab.sum(axis=1).sort_values(ascending=False).index
        cross_tab = cross_tab.loc[sorted_idx]

        ax = cross_tab.plot(kind='bar', stacked=True, figsize=(14, 8), color=[AOP_COLOR_MAP[c] for c in cross_tab.columns], alpha=0.7, edgecolor='black', linewidth=0.5, width=0.8)
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
        save_high_res("Louvain_Community_Composition", output_dir, file_prefix)
    except Exception as e:
        print(f"[!] Error generating Figure 3: {e}")
    finally:
        plt.close('all')

def plot_joint_plot(node_df: pd.DataFrame, output_dir: str, file_prefix: str):
    """Figure 4: joint scatter + marginal histograms of betweenness vs PageRank MRS."""
    print("\n<<< Generating Figure 4: MRS Correlation Joint Plot >>>")
    try:
        plot_df = node_df[(node_df['MRS_betweenness'] > 0) & (node_df['MRS_pagerank'] > 0)].copy()

        if plot_df.empty:
            print("[!] Error: No strictly positive MRS values available.")
            return

        fig = plt.figure(figsize=(12, 12))
        gs = fig.add_gridspec(2, 2, width_ratios=(4, 1), height_ratios=(1, 4), left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.05, hspace=0.05)
        ax_scat = fig.add_subplot(gs[1, 0])
        ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_scat)
        ax_histy = fig.add_subplot(gs[1, 1], sharey=ax_scat)

        ax_scat.set_xscale('log')
        ax_scat.set_yscale('log')

        ax_histx.tick_params(axis="x", which="both", bottom=True, labelbottom=False)
        ax_histy.tick_params(axis="y", which="both", left=True, labelleft=False)

        sns.scatterplot(data=plot_df, x='MRS_betweenness', y='MRS_pagerank', hue='aop_level', style='aop_level', markers=MARKERS_LIST, palette='magma_r', hue_order=AOP_ORDER, s=60, alpha=0.7, ax=ax_scat)
        
        ax_scat.grid(True, which="both", ls="--", alpha=0.5)
        ax_histx.grid(True, axis="x", which="both", ls="--", alpha=0.5)
        ax_histy.grid(True, axis="y", which="both", ls="--", alpha=0.5)

        ax_scat.legend(title=r"$\bf{AOP\ Level}$", bbox_to_anchor=(1.02, 1), loc='lower left')
        ax_scat.set_xlabel('MRS (Betweenness)', fontsize=12)
        ax_scat.set_ylabel('MRS (PageRank)', fontsize=12)

        sns.histplot(data=plot_df, x='MRS_betweenness', hue='aop_level', multiple="stack", palette='magma_r', hue_order=AOP_ORDER, legend=False, ax=ax_histx, log_scale=True, bins=15)
        sns.histplot(data=plot_df, y='MRS_pagerank', hue='aop_level', multiple="stack", palette='magma_r', hue_order=AOP_ORDER, legend=False, ax=ax_histy, log_scale=True, bins=15)

        top = pd.concat([plot_df.nlargest(5, 'MRS_betweenness'), plot_df.nlargest(5, 'MRS_pagerank')]).drop_duplicates()
        texts = [ax_scat.text(row['MRS_betweenness'], row['MRS_pagerank'], idx, fontsize=9) for idx, row in top.iterrows()]

        with redirect_stdout(io.StringIO()):
            adjust_text(texts, ax=ax_scat, arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

        save_high_res("Joint_Plot_Final", output_dir, file_prefix)
    except Exception as e:
        print(f"[!] Error generating Figure 4: {e}")
    finally:
        plt.close('all')

def plot_tsne_louvain_overlap(node_df: pd.DataFrame, G: nx.Graph, output_dir: str, file_prefix: str):
    """Figure 5: t-SNE projection of the graph distance matrix, coloured by Louvain community."""
    print("\n<<< Generating Figure 5: t-SNE with Communities >>>")
    try:
        if 'filtered_louvain_community_id' not in node_df.columns:
            print("    [!] Missing 'filtered_louvain_community_id'. Skipping Figure 5.")
            return

        degree_dict = dict(G.degree(weight='cooccurrence_count')) if nx.is_weighted(G, weight='cooccurrence_count') else dict(G.degree())
        nodes = list(G.nodes())
        N = len(nodes)
        if N < 3:
            print("    [!] Too few nodes for a t-SNE projection. Skipping Figure 5.")
            return
        node_idx = {n: i for i, n in enumerate(nodes)}

        # Build the shortest-path distance matrix directly in NumPy. The prior
        # implementation filled an NxN DataFrame with scalar .loc assignments,
        # which is dramatically slower than indexing a preallocated array.
        dist = np.full((N, N), np.nan, dtype=float)
        for u, dmap in nx.all_pairs_shortest_path_length(G):
            iu = node_idx[u]
            for v, d in dmap.items():
                dist[iu, node_idx[v]] = d
        finite_vals = dist[np.isfinite(dist)]
        fill_value = (finite_vals.max() * 2) if finite_vals.size else 1.0
        dist = np.where(np.isnan(dist), fill_value, dist)

        # perplexity must be strictly less than n_samples; clamp for small graphs.
        perplexity = min(30, max(2, N - 1))
        tsne = TSNE(n_components=2, perplexity=perplexity, metric='precomputed', init='random', random_state=42)
        emb = tsne.fit_transform(dist)
        df_tsne = pd.DataFrame(emb, index=nodes, columns=['X', 'Y']).join(node_df[['filtered_louvain_community_id', 'MRS_betweenness']])

        plt.figure(figsize=(12, 10))
        ax = plt.gca()
        unique_comms = sorted(df_tsne['filtered_louvain_community_id'].unique())
        comm_map = dict(zip(unique_comms, sns.color_palette('viridis', n_colors=len(unique_comms))))

        community_labels = {}
        for comm_id in unique_comms:
            subset = df_tsne[df_tsne['filtered_louvain_community_id'] == comm_id]
            if not subset.empty:
                top_node = max(subset.index.tolist(), key=lambda x: degree_dict.get(x, 0))
                aop_tag = node_df.loc[top_node, 'aop_level'] if pd.notna(node_df.loc[top_node, 'aop_level']) else "?"
                community_labels[comm_id] = f"{comm_id}: {top_node} ({aop_tag})"
            else:
                community_labels[comm_id] = f"{comm_id}: (Empty)"

        for cid in unique_comms:
            pts = df_tsne[df_tsne['filtered_louvain_community_id'] == cid][['X', 'Y']].values
            if len(pts) >= 3:
                try:
                    hull = ConvexHull(pts)
                    poly = plt.Polygon(pts[hull.vertices], closed=True, facecolor=comm_map[cid], alpha=0.1, edgecolor=comm_map[cid], linestyle='--')
                    ax.add_patch(poly)
                except Exception:
                    continue

        sns.scatterplot(data=df_tsne, x='X', y='Y', hue='filtered_louvain_community_id', palette=comm_map, size='MRS_betweenness', sizes=(20, 300), alpha=0.9, legend=False, ax=ax)
        legend_patches = [mpatches.Patch(color=comm_map[cid], label=community_labels[cid], alpha=0.9) for cid in unique_comms]
        ax.legend(handles=legend_patches, title='Louvain Communities & Topics', loc='upper right', frameon=True, facecolor='white', edgecolor='black', framealpha=0.95, fontsize=10, title_fontsize=11)
        ax.set_xlim(right=15)
        plt.xlabel("t-SNE Dimension 1")
        plt.ylabel("t-SNE Dimension 2")

        save_high_res("T-SNE_Final", output_dir, file_prefix)
    except Exception as e:
        print(f"[!] Error generating Figure 5: {e}")
    finally:
        plt.close('all')

def plot_sankey_alluvial(G: nx.Graph, node_df: pd.DataFrame, output_dir: str, file_prefix: str):
    """Figure 6: interactive AOP Sankey/alluvial flow from stressors to adverse outcomes (labeled and blank)."""
    print("\n<<< Figure 6: AOP Alluvial Flow (Labeled & Unlabeled) >>>")
    try:
        levels = [l for l in AOP_ORDER if l != 'Uncategorized']
        lvl_map = {l: i for i, l in enumerate(levels)}
        flows = defaultdict(int)

        unassigned_edges = 0

        for u, v in G.edges():
            lu, lv = node_df.loc[u, 'aop_level'], node_df.loc[v, 'aop_level']
            if lu in lvl_map and lv in lvl_map:
                u_i, v_i = lvl_map[lu], lvl_map[lv]
                flows[(u_i, v_i) if u_i <= v_i else (v_i, u_i)] += 1
            else:
                unassigned_edges += 1

        if unassigned_edges > 0:
            print(f"    [!] Note: {unassigned_edges:,} edge connections involved 'Uncategorized' nodes and were not plotted.")

        # Ensure we always pass valid arrays to Plotly, even if flows is completely empty
        sources = [k[0] for k in flows] if flows else []
        targets = [k[1] for k in flows] if flows else []
        values = list(flows.values()) if flows else []
        labels = [str(v) for v in flows.values()] if flows else []

        # If empty, this will create an empty DataFrame and still output the CSV structurally intact
        csv_data = [{"Source": levels[s], "Target": levels[t], "Weight": w} for (s, t), w in flows.items()]
        csv_path = os.path.join(output_dir, f"{file_prefix}_Alluvial_Connections_Table_F6.csv")
        pd.DataFrame(csv_data, columns=["Source", "Target", "Weight"]).sort_values(['Source', 'Target']).to_csv(_paths.long_path(csv_path), index=False)
        print(f"    Exported connection data to: {os.path.basename(csv_path)}")

        # Labeled Plot
        fig_labels = go.Figure(go.Sankey(
            node=dict(pad=25, thickness=20, line=dict(color="black", width=0.5), label=levels, color=sns.color_palette("magma_r", len(levels)).as_hex()),
            link=dict(source=sources, target=targets, value=values, label=labels)
        ))
        fig_labels.update_layout(title="AOP Labeled Alluvial Flow")
        export_plotly_figure(fig_labels, "Alluvial_Labeled", output_dir, file_prefix)

        # Unlabeled Plot
        blank_labels = [""] * len(flows) if flows else []
        fig_clean = go.Figure(go.Sankey(
            node=dict(pad=25, thickness=20, line=dict(color="black", width=0.5), label=[""] * len(levels), color=sns.color_palette("magma_r", len(levels)).as_hex()),
            link=dict(source=sources, target=targets, value=values, label=blank_labels)
        ))
        fig_clean.update_layout(title="AOP Unlabeled Alluvial Flow")
        export_plotly_figure(fig_clean, "Alluvial_Clean", output_dir, file_prefix)

    except Exception as e:
        print(f"[!] Error generating Figure 6 (Alluvial Flow): {e}")

def plot_dumbell_plot(node_df: pd.DataFrame, output_dir: str, file_prefix: str):
    """Figure 7: dumbbell plot of the nodes with the largest shift between betweenness and PageRank MRS."""
    print("\n<<< Figure 7: Dumbbell Plot >>>")
    try:
        df_plot = node_df.copy()
        df_plot['diff'] = (df_plot['MRS_betweenness'] - df_plot['MRS_pagerank']).abs()
        top_diff = df_plot.nlargest(20, 'diff').sort_values('MRS_betweenness')

        plt.figure(figsize=(10, 8))
        plt.hlines(y=top_diff.index, xmin=top_diff['MRS_betweenness'], xmax=top_diff['MRS_pagerank'], color='gray', alpha=0.5, linewidth=1.5, zorder=1)
        colors = ["#2C3E50", "#B0B0B0"]
        plt.scatter(top_diff['MRS_betweenness'], top_diff.index, color=colors[0], s=100, zorder=2, label='Betweenness')
        plt.scatter(top_diff['MRS_pagerank'], top_diff.index, color=colors[1], s=100, zorder=2, label='PageRank')

        ax = plt.gca()
        for lbl in ax.get_yticklabels():
            if lbl.get_text() in top_diff.index:
                lvl = top_diff.loc[lbl.get_text(), 'aop_level']
                if pd.notna(lvl) and lvl in AOP_COLOR_MAP:
                    lbl.set_color(AOP_COLOR_MAP[lvl])
                    lbl.set_fontweight('bold')

        leg_h = ([mpatches.Patch(color='none', label=r"$\bf{AOP\ Level}$")] +
                 [mpatches.Patch(color=c, label=l) for l, c in AOP_COLOR_MAP.items()] +
                 [mpatches.Patch(color='none', label="")] +
                 [mpatches.Patch(color='none', label=r"$\bf{Method}$")] +
                 [mpatches.Patch(color=colors[0], label='Betweenness'), mpatches.Patch(color=colors[1], label='PageRank')])

        plt.legend(handles=leg_h, loc='lower right', bbox_to_anchor=(1.0, 0.0))
        plt.xlabel("Mean Relevancy Score (MRS)")
        plt.tight_layout()
        save_high_res("MRS_Dumbbell_Plot_Top20", output_dir, file_prefix)
    except Exception as e:
        print(f"[!] Error generating Figure 7 (Dumbbell Plot): {e}")
    finally:
        plt.close('all')

def _overlay_loglog_fits(ax, plot_df: pd.DataFrame, x_col: str, y_col: str, min_points: int = 3):
    """Draw per-stratum power-law trend lines on an axis that is already log-log.

    Seaborn fits in LINEAR space, and its `logx` option fits only y ~ log(x), so
    letting lmplot draw the line here produced a straight line in linear units
    rendered on log axes: it looked like a power law, was dominated by the few
    largest-x points, and its slope was not the exponent quoted in the caption.
    Fitting log10(y) on log10(x) makes the drawn line and the reported exponent
    the same object. Colours come from AOP_COLOR_MAP, which is built from the same
    magma_r ramp lmplot uses for the points, so line and markers stay matched.
    """
    fits = []
    for level in AOP_ORDER:
        grp = plot_df[plot_df['aop_level'] == level]
        if len(grp) < min_points:
            if len(grp) > 0:
                print(f"      {level:<18} n={len(grp):<4} (fewer than {min_points} points, not fitted)")
            continue
        lx = np.log10(grp[x_col].to_numpy(dtype=float))
        ly = np.log10(grp[y_col].to_numpy(dtype=float))
        slope, intercept = np.polyfit(lx, ly, 1)
        ss_tot = float(((ly - ly.mean()) ** 2).sum())
        r2 = 1.0 - float(((ly - (intercept + slope * lx)) ** 2).sum()) / ss_tot if ss_tot > 0 else float('nan')
        xs = np.linspace(lx.min(), lx.max(), 50)
        ax.plot(10 ** xs, 10 ** (intercept + slope * xs), color=AOP_COLOR_MAP.get(level),
                lw=2.0, alpha=0.9, zorder=4)
        fits.append((level, len(grp), slope, r2))
        print(f"      {level:<18} n={len(grp):<4} slope {slope:+.3f}   R2 {r2:.3f}")
    return fits


def plot_scatter_panels(node_df: pd.DataFrame, output_dir: str, file_prefix: str):
    """Figure 8: log-log scatter panels of each raw centrality against its mean relevancy score."""
    print("\n<<< Figure 8: Centrality Scatter Plots (MRS vs Raw) >>>")
    try:
        print("  - Generating Panel C: Raw Betweenness vs MRS (Betweenness)...")
        plot_df_b = node_df[(node_df['betweenness_centrality'] > 0) & (node_df['MRS_betweenness'] > 0)]
        g_b = sns.lmplot(data=plot_df_b, x='betweenness_centrality', y='MRS_betweenness', hue='aop_level', palette='magma_r', hue_order=AOP_ORDER, height=8, aspect=1.1, fit_reg=False, legend=False, markers=MARKERS_LIST, scatter_kws={'s': 50, 'alpha': 0.7})
        ax_b = g_b.ax
        ax_b.set_xscale('log'); ax_b.set_yscale('log')
        _overlay_loglog_fits(ax_b, plot_df_b, 'betweenness_centrality', 'MRS_betweenness')
        ax_b.set_xlabel('Raw Betweenness Centrality (Log Scale)', fontsize=12)
        ax_b.set_ylabel('MRS (Betweenness) (Log Scale)', fontsize=12)
        ax_b.grid(True, which="both", ls="--", alpha=0.5)
        save_high_res("Panel_C_Scatter_Betweenness", output_dir, file_prefix)
        plt.close()

        print("  - Generating Panel D: Raw PageRank vs MRS (PageRank)...")
        plot_df_p = node_df[(node_df['pagerank_centrality'] > 0) & (node_df['MRS_pagerank'] > 0)]
        g_p = sns.lmplot(data=plot_df_p, x='pagerank_centrality', y='MRS_pagerank', hue='aop_level', palette='magma_r', hue_order=AOP_ORDER, height=8, aspect=1.1, fit_reg=False, legend=False, markers=MARKERS_LIST, scatter_kws={'s': 50, 'alpha': 0.7})
        ax_p = g_p.ax
        ax_p.set_xscale('log'); ax_p.set_yscale('log')
        _overlay_loglog_fits(ax_p, plot_df_p, 'pagerank_centrality', 'MRS_pagerank')
        ax_p.set_xlabel('Raw PageRank Centrality (Log Scale)', fontsize=12)
        ax_p.set_ylabel('MRS (PageRank) (Log Scale)', fontsize=12)
        ax_p.grid(True, which="both", ls="--", alpha=0.5)
        save_high_res("Panel_D_Scatter_PageRank", output_dir, file_prefix)
        plt.close()
    except Exception as e:
        print(f"[!] Error generating Figure 8 (Scatter Panels): {e}")
    finally:
        plt.close('all')

def _reproducible_hash(value):
    """Deterministic word hash so Word2Vec initialization is reproducible regardless
    of PYTHONHASHSEED (Python's built-in str hash is randomized per process)."""
    return int(hashlib.md5(str(value).encode('utf-8')).hexdigest(), 16)


def plot_dendrogram(G: nx.Graph, node_df: pd.DataFrame, output_dir: str, file_prefix: str, random_seed: int = None):
    """Ward dendrogram of Node2Vec embeddings, with leaf labels coloured by AOP level (warns clearly if node2vec is unavailable).

    When `random_seed` is provided the embedding is fully reproducible run-to-run:
    the biased walks are seeded, gensim trains single-threaded, and a deterministic
    word-hash removes the PYTHONHASHSEED dependence. Without a seed it stays
    stochastic (as the original library was)."""
    if Node2Vec is None:
        print("\n<<< Node2Vec Dendrogram >>>")
        print("  [!] SKIPPED - the embedding backend could not be loaded, so this figure")
        print("      was not generated. Every other figure is unaffected.")
        print(f"      Reason: {type(_NODE2VEC_IMPORT_ERROR).__name__}: {_NODE2VEC_IMPORT_ERROR}")
        if 'gensim' in str(_NODE2VEC_IMPORT_ERROR).lower():
            print("      -> This figure needs gensim.  Fix:  python -m pip install gensim")
        else:
            print("      -> Ensure gensim, networkx and numpy are installed and importable.")
        return
    print("\n<<< Generating Node2Vec Dendrogram >>>")
    try:
        # workers=1 + seed + deterministic hashfxn make the embedding reproducible
        # (gensim training is non-deterministic across workers / hash seeds).
        n2v = Node2Vec(G, dimensions=64, walk_length=30, num_walks=200, workers=1, quiet=True, seed=random_seed)
        w2v_params = {'window': 10, 'min_count': 1}
        if random_seed is not None:
            w2v_params.update(seed=random_seed, hashfxn=_reproducible_hash)
        model = n2v.fit(**w2v_params)
        labels = model.wv.index_to_key
        vectors = model.wv.vectors
        Z = linkage(vectors, method='ward', metric='euclidean')

        plt.figure(figsize=(16, max(15, len(labels)*0.25)))
        dendrogram(Z, labels=labels, orientation='right', color_threshold=0, above_threshold_color='black', leaf_font_size=12)

        ax = plt.gca()
        for lbl in ax.get_yticklabels():
            term = lbl.get_text()
            if term in node_df.index:
                lvl = node_df.loc[term, 'aop_level']
                if lvl in AOP_COLOR_MAP: lbl.set_color(AOP_COLOR_MAP[lvl])

        aop_legend_patches = [mpatches.Patch(color=color, label=level) for level, color in AOP_COLOR_MAP.items()]
        ax.legend(handles=aop_legend_patches, title='AOP Level (Label Color)', loc='upper left', bbox_to_anchor=(0.82, 1.0))
        plt.xlabel("Cluster Distance", fontsize=12)
        plt.ylabel("MeSH Terms", fontsize=12)
        plt.grid(axis='x', linestyle='--', alpha=0.5)

        save_high_res("Dendrogram_Node2Vec_AOP", output_dir, file_prefix)
    except Exception as e:
        print(f"[!] Error generating Dendrogram: {e}")
    finally:
        plt.close('all')

def generate_filtering_summary_sankeys(data: dict, output_dir: str, file_prefix: str):
    """Generates Filtering Sankeys (Nodes & Edges) and summary table."""
    print("\n<<< Generating Filtering Cascade Sankeys >>>")
    try:
        def make_sankey(item_type, output_name):
            labels = ["GLF Initial", "SA Initial", "GLF LCC", "SA LCC", "Final Consensus LCC", "Discarded (Not in LCC)", "Discarded (Not in Consensus)"]
            sources = [0, 0, 1, 1, 2, 2, 3, 3]
            targets = [2, 5, 3, 5, 4, 6, 4, 6]

            glf_i = data['Initial Subgraph']['GLF'][item_type]
            sa_i = data['Initial Subgraph']['SA'][item_type]
            glf_lcc = data['LCC']['GLF'][item_type]
            sa_lcc = data['LCC']['SA'][item_type]
            final = data['Final Consensus LCC']['GLF'][item_type]

            values = [glf_lcc, glf_i - glf_lcc, sa_lcc, sa_i - sa_lcc, final, glf_lcc - final, final, sa_lcc - final]
            link_colors = ['rgba(55, 126, 184, 0.6)', 'rgba(200, 200, 200, 0.4)', 'rgba(255, 127, 0, 0.6)', 'rgba(200, 200, 200, 0.4)', 'rgba(55, 126, 184, 0.6)', 'rgba(200, 200, 200, 0.4)', 'rgba(255, 127, 0, 0.6)', 'rgba(200, 200, 200, 0.4)']
            node_colors = ["#377eb8", "#ff7f00", "#377eb8", "#ff7f00", "#984ea3", "grey", "grey"]

            fig = go.Figure(data=[go.Sankey(
                node=dict(pad=25, thickness=20, line=dict(color="black", width=0.5), label=[f"{l}\n({v:,})" for l, v in zip(labels, [glf_i, sa_i, glf_lcc, sa_lcc, final, (glf_i-glf_lcc)+(sa_i-sa_lcc), (glf_lcc-final)+(sa_lcc-final)])], color=node_colors),
                link=dict(source=sources, target=targets, value=values, color=link_colors)
            )])
            fig.update_layout(title_text=f"Filtering Cascade: {item_type}", font_size=12)
            export_plotly_figure(fig, output_name, output_dir, file_prefix)

        make_sankey('Nodes', 'Node_Sankey')
        make_sankey('Edges', 'Edge_Sankey')

        print("\n" + "="*80)
        print("                Summary of Network Filtering Cascade")
        print("="*80)
        rows = []
        stages = [('Initial Subgraph', 'Initial Subgraph'), ('LCC', 'Largest Connected Component'), ('Final Consensus LCC', 'Final Consensus Largest Connected Component')]
        for key, name in stages:
            r = {'Stage': name}
            for method in ['GLF', 'SA']:
                for item in ['Nodes', 'Edges']:
                    r[f"{method} {item}"] = data[key][method][item]
            rows.append(r)

        df = pd.DataFrame(rows).set_index('Stage')
        initial, lcc, final = 'Initial Subgraph', 'Largest Connected Component', 'Final Consensus Largest Connected Component'
        df['GLF Node Retention'] = ''
        df['SA Node Retention'] = ''
        df.loc[lcc, 'GLF Node Retention'] = f"{df.loc[lcc, 'GLF Nodes']/df.loc[initial, 'GLF Nodes']:.1%}"
        df.loc[lcc, 'SA Node Retention'] = f"{df.loc[lcc, 'SA Nodes']/df.loc[initial, 'SA Nodes']:.1%}"
        df.loc[final, 'GLF Node Retention'] = f"{df.loc[final, 'GLF Nodes']/df.loc[lcc, 'GLF Nodes']:.1%}"
        df.loc[final, 'SA Node Retention'] = f"{df.loc[final, 'SA Nodes']/df.loc[lcc, 'SA Nodes']:.1%}"

        print(df[['GLF Nodes', 'GLF Edges', 'GLF Node Retention', 'SA Nodes', 'SA Edges', 'SA Node Retention']].to_string())
        print("-" * 80 + "\n")
    except Exception as e:
        print(f"[!] Error generating Filtering Cascade Sankeys: {e}")

def plot_mrs_ubiquity_bias(node_df: pd.DataFrame, output_dir: str, file_prefix: str):
    """Calculates Spearman rank correlation to evaluate ubiquity bias and generates scatter plots."""
    print("\n<<< Generating Figure: MRS Ubiquity Bias Evaluation >>>")
    try:
        # Map to the renamed columns in load_and_prepare_data
        mrs_cols = ["MRS_pagerank", "MRS_betweenness"]
        count_col = "article_count"

        missing_cols = [col for col in [count_col] + mrs_cols if col not in node_df.columns]
        if missing_cols:
            print(f"    [!] Missing columns for bias evaluation: {missing_cols}")
            return

        df_clean = node_df.dropna(subset=[count_col] + mrs_cols).copy()

        df_clean[count_col] = pd.to_numeric(df_clean[count_col], errors='coerce')
        for col in mrs_cols:
            df_clean[col] = pd.to_numeric(df_clean[col], errors='coerce')

        df_clean = df_clean.dropna(subset=[count_col] + mrs_cols)

        print(f"    Extracted {len(df_clean)} nodes for statistical evaluation.\n")

        fig, axes = plt.subplots(1, len(mrs_cols), figsize=(6 * len(mrs_cols), 5))
        if len(mrs_cols) == 1: axes = [axes]

        for i, mrs_col in enumerate(mrs_cols):
            rho, p_val = stats.spearmanr(df_clean[count_col], df_clean[mrs_col])

            print(f"--- Analysis for: {mrs_col} ---")
            print(f"  Spearman Rho (ρ): {rho:.4f}")
            print(f"  P-value:          {p_val:.4e}")

            if rho >= 0.7:
                print("  Status: HIGH Probability of Bias (Ubiquity bias is still present).")
            elif rho <= 0.4:
                print("  Status: DECOUPLED (Neutralized volume bias).")
            else:
                print("  Status: MODERATE correlation (Bias confounders may still exist).")
            print("-" * 40)

            ax = axes[i]
            sns.regplot(
                x=df_clean[count_col],
                y=df_clean[mrs_col],
                ax=ax,
                scatter_kws={'alpha': 0.5, 'color': '#2c3e50'},
                line_kws={'color': '#e74c3c'},
                logx=True
            )

            ax.set_xscale('log')
            ax.set_title(f"{mrs_col} vs Publication Count\nSpearman ρ = {rho:.3f}", fontsize=12, fontweight='bold')
            ax.set_xlabel(f"Publication Count (Log Scale)", fontsize=10)
            ax.set_ylabel(f"Score [{mrs_col}]", fontsize=10)

        plt.tight_layout()
        save_high_res("MRS_Ubiquity_Bias_Evaluation", output_dir, file_prefix)
    except Exception as e:
        print(f"[!] Error generating MRS Bias Evaluation: {e}")
    finally:
        plt.close('all')