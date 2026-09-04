# -*- coding: utf-8 -*-
"""
viz.py - figure generation and visualization (pipeline Step 4).

Produces the publication-quality figures and interactive exports from the
finished, annotated network.

Outputs include edge-weight and dispersion (zero-truncated negative binomial)
distributions, the GLF-versus-SA optimization trajectory, Louvain community
composition bars, centrality correlation and dumbbell plots, a t-SNE community
projection, a Node2Vec dendrogram, and the Sankey/alluvial flows between the
strata. Static images are written as high-resolution JPEG/TIFF and interactive
figures as self-contained HTML.

Three figures read the annotation: the community composition bars, the alluvial
flow, and the dendrogram's label colours. None of them assumes what the strata
are. The names come from the annotation file, their order from the settings,
and the palette is sized to however many the run turns out to have - so a
network divided into four organ systems plots exactly as well as one divided
into the seven levels of an adverse outcome pathway.
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
import matplotlib.patheffects as path_effects
import seaborn as sns
import scipy.stats as stats
from scipy.spatial import ConvexHull
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.manifold import TSNE
import statsmodels.api as sm
from statsmodels.discrete.truncated_model import TruncatedLFNegativeBinomialP
import plotly.graph_objects as go
from collections import defaultdict

# Relative imports from our package
from .stats import calculate_graph_stats
from . import paths as _paths
from . import strata as _strata

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
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
sns.set_context("paper", font_scale=1.2)
sns.set_style("whitegrid")
MARKERS_LIST = ['o', 's', '^', 'D', 'v', 'X', 'P']
# A fixed sample of the same palette, for the figures that colour by something
# other than stratum and so must not move when the strata count changes.
MAGMA_PALETTE = sns.color_palette("magma_r", n_colors=7)

# The strata this run divides the network into, and a colour for each. Neither
# can be a constant any more: the names come from the annotation file, which
# may say anything, the order comes from the settings, and the palette has to
# be sized to however many there turn out to be. configure_strata() takes the
# configured order before any figure is drawn; _adopt_strata() reconciles it
# with what the file actually contained. The defaults stand if nothing calls
# them, so the module still behaves when driven directly from a notebook.
_CONFIGURED_STRATA = ''
STRATA_ORDER = list(_strata.DEFAULT_ORDER)
STRATA_COLOR_MAP = dict(zip(_strata.DEFAULT_ORDER,
                            sns.color_palette("magma_r",
                                              n_colors=len(_strata.DEFAULT_ORDER))))


def configure_strata(order_text):
    """The strata order from the settings. Call before drawing anything."""
    global _CONFIGURED_STRATA
    _CONFIGURED_STRATA = order_text or ''


def _adopt_strata(present):
    """Fix the order and the palette for the strata this run actually has.

    Anything in the file but not in the setting is appended and reported. A
    stratum the user typed is one they meant, and a figure that silently
    dropped it would be a figure that lied about the annotation.
    """
    global STRATA_ORDER, STRATA_COLOR_MAP
    order, appended = _strata.resolve_order(
        _CONFIGURED_STRATA or list(_strata.DEFAULT_ORDER), present)
    if not order:
        order = [_strata.UNASSIGNED]
    STRATA_ORDER = order
    STRATA_COLOR_MAP = dict(
        zip(order, sns.color_palette("magma_r", n_colors=len(order))))
    if appended:
        print(f"    [i] Not in the Strata order setting, so placed at the end: "
              f"{', '.join(appended)}")
    return order

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
    """Loads the final network and merges the strata annotations (semicolon CSV)."""
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

    col = _strata.COLUMN
    if os.path.exists(annotation_path):
        print(f"Merging annotations from {os.path.basename(annotation_path)}...")
        try:
            # Read explicitly with semicolon
            try:
                anno_df = pd.read_csv(annotation_path, sep=';')
                if len(anno_df.columns) == 1 and ',' in anno_df.columns[0]:
                    anno_df = pd.read_csv(annotation_path, sep=',')
            except Exception:
                anno_df = pd.read_csv(annotation_path, sep=',')

            # Standardize column casing, then accept either spelling of the
            # stratum column: a file annotated before the rename says
            # 'aop_level' and must keep working without being re-annotated.
            anno_df.rename(columns={c: str(c).lower() for c in anno_df.columns},
                           inplace=True)
            if col not in anno_df.columns and _strata.LEGACY_COLUMN in anno_df.columns:
                anno_df.rename(columns={_strata.LEGACY_COLUMN: col}, inplace=True)

            node_df = node_df.reset_index().merge(anno_df, on='mesh_term', how='left').set_index('mesh_term')
        except Exception as e:
            print(f"[!] WARNING: Could not read annotation file: {e}")
            node_df[col] = _strata.UNASSIGNED

        if col not in node_df.columns:
            print(f"[!] WARNING: the annotation file has no '{col}' column. "
                  f"Proceeding with every term {_strata.UNASSIGNED}.")
            node_df[col] = _strata.UNASSIGNED

        # Whatever the user wrote is a stratum. Blanks, the template's
        # placeholder and the word the figures use all mean the same thing.
        node_df[col] = node_df[col].map(_strata.normalise)
        order = _adopt_strata(node_df[col].tolist())
        node_df[col] = pd.Categorical(node_df[col], categories=order, ordered=True)
        nx.set_node_attributes(G, node_df[col].to_dict(), col)
    else:
        print(f"[!] WARNING: Annotation file not found at {annotation_path}. "
              f"Proceeding {_strata.UNASSIGNED}.")
        node_df[col] = _strata.UNASSIGNED
        _adopt_strata([_strata.UNASSIGNED])

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
    """Figure 1: how often terms co-occur, before and after consensus filtering."""
    print("\n<<< Figure 1 of 7: Edge weight distribution - what filtering removed >>>")
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
    print("\n<<< Figure 2 of 7: Optimisation trajectory - GLF against simulated annealing >>>")
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
    """Figure 3: which strata make up each Louvain community."""
    print("\n<<< Figure 3 of 7: Community composition - strata per community >>>")
    try:
        if 'filtered_louvain_community_id' not in node_df.columns:
            print("    [!] Missing 'filtered_louvain_community_id' column. Skipping Figure 3.")
            return

        cross_tab = pd.crosstab(node_df['filtered_louvain_community_id'],
                                node_df[_strata.COLUMN])
        cross_tab = cross_tab.reindex(columns=STRATA_ORDER, fill_value=0)
        sorted_idx = cross_tab.sum(axis=1).sort_values(ascending=False).index
        cross_tab = cross_tab.loc[sorted_idx]

        ax = cross_tab.plot(kind='bar', stacked=True, figsize=(14, 8), color=[STRATA_COLOR_MAP[c] for c in cross_tab.columns], alpha=0.7, edgecolor='black', linewidth=0.5, width=0.8)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
        ax.grid(False)
        ax.tick_params(axis='x', length=5)
        ax.tick_params(axis='y', length=5)
        sns.despine(ax=ax, left=False, bottom=False)

        plt.xlabel('Louvain Community ID', fontsize=12)
        plt.ylabel('Number of Nodes', fontsize=12)
        plt.legend(title=r"$\bf{Stratum}$", bbox_to_anchor=(0.81, 1), loc='upper left')
        plt.subplots_adjust(right=0.75)
        plt.tight_layout(rect=[0, 0, 0.85, 1])
        save_high_res("Louvain_Community_Composition", output_dir, file_prefix)
    except Exception as e:
        print(f"[!] Error generating Figure 3: {e}")
    finally:
        plt.close('all')


def plot_tsne_louvain_overlap(node_df: pd.DataFrame, G: nx.Graph, output_dir: str, file_prefix: str):
    """Figure 4: t-SNE of the graph distance matrix, coloured by Louvain community."""
    print("\n<<< Figure 4 of 7: t-SNE projection - do the communities separate? >>>")
    try:
        if 'filtered_louvain_community_id' not in node_df.columns:
            print("    [!] Missing 'filtered_louvain_community_id'. Skipping Figure 4.")
            return

        degree_dict = dict(G.degree(weight='cooccurrence_count')) if nx.is_weighted(G, weight='cooccurrence_count') else dict(G.degree())
        nodes = list(G.nodes())
        N = len(nodes)
        if N < 3:
            print("    [!] Too few nodes for a t-SNE projection. Skipping Figure 4.")
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
                tag = node_df.loc[top_node, _strata.COLUMN] \
                    if pd.notna(node_df.loc[top_node, _strata.COLUMN]) else "?"
                community_labels[comm_id] = f"{comm_id}: {top_node} ({tag})"
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
        print(f"[!] Error generating Figure 4: {e}")
    finally:
        plt.close('all')

def plot_sankey_alluvial(G: nx.Graph, node_df: pd.DataFrame, output_dir: str, file_prefix: str):
    """Figure 5: interactive alluvial flow between the strata.

    The order of the strata is the flow: it is what makes the left-to-right
    reading mean anything. It comes from the Strata order setting, because no
    order can be inferred from the names alone - only the person who chose them
    knows which comes first.
    """
    print("\n<<< Figure 5 of 7: alluvial flow between strata >>>")
    try:
        levels = [l for l in STRATA_ORDER if l != _strata.UNASSIGNED]
        if len(levels) < 2:
            print(f"    [!] Only {len(levels)} stratum in this network. An "
                  f"alluvial flow needs at least two; skipping Figure 5.")
            return
        lvl_map = {l: i for i, l in enumerate(levels)}
        flows = defaultdict(int)

        unassigned_edges = 0

        for u, v in G.edges():
            lu, lv = node_df.loc[u, _strata.COLUMN], node_df.loc[v, _strata.COLUMN]
            if lu in lvl_map and lv in lvl_map:
                u_i, v_i = lvl_map[lu], lvl_map[lv]
                flows[(u_i, v_i) if u_i <= v_i else (v_i, u_i)] += 1
            else:
                unassigned_edges += 1

        if unassigned_edges > 0:
            print(f"    [!] Note: {unassigned_edges:,} edge connections involved "
                  f"'{_strata.UNASSIGNED}' nodes and were not plotted.")

        # Ensure we always pass valid arrays to Plotly, even if flows is completely empty
        sources = [k[0] for k in flows] if flows else []
        targets = [k[1] for k in flows] if flows else []
        values = list(flows.values()) if flows else []
        labels = [str(v) for v in flows.values()] if flows else []

        # If empty, this will create an empty DataFrame and still output the CSV structurally intact
        csv_data = [{"Source": levels[s], "Target": levels[t], "Weight": w} for (s, t), w in flows.items()]
        csv_path = os.path.join(output_dir, f"{file_prefix}_Alluvial_Connections_Table_F5.csv")
        pd.DataFrame(csv_data, columns=["Source", "Target", "Weight"]).sort_values(['Source', 'Target']).to_csv(_paths.long_path(csv_path), index=False)
        print(f"    Exported connection data to: {os.path.basename(csv_path)}")

        # Labeled Plot
        fig_labels = go.Figure(go.Sankey(
            node=dict(pad=25, thickness=20, line=dict(color="black", width=0.5), label=levels, color=sns.color_palette("magma_r", len(levels)).as_hex()),
            link=dict(source=sources, target=targets, value=values, label=labels)
        ))
        fig_labels.update_layout(title="Labelled Alluvial Flow Between Strata")
        export_plotly_figure(fig_labels, "Alluvial_Labeled", output_dir, file_prefix)

        # Unlabeled Plot
        blank_labels = [""] * len(flows) if flows else []
        fig_clean = go.Figure(go.Sankey(
            node=dict(pad=25, thickness=20, line=dict(color="black", width=0.5), label=[""] * len(levels), color=sns.color_palette("magma_r", len(levels)).as_hex()),
            link=dict(source=sources, target=targets, value=values, label=blank_labels)
        ))
        fig_clean.update_layout(title="Unlabelled Alluvial Flow Between Strata")
        export_plotly_figure(fig_clean, "Alluvial_Clean", output_dir, file_prefix)

    except Exception as e:
        print(f"[!] Error generating Figure 5 (alluvial flow): {e}")


def _reproducible_hash(value):
    """Deterministic word hash so Word2Vec initialization is reproducible regardless
    of PYTHONHASHSEED (Python's built-in str hash is randomized per process)."""
    return int(hashlib.md5(str(value).encode('utf-8')).hexdigest(), 16)


def plot_dendrogram(G: nx.Graph, node_df: pd.DataFrame, output_dir: str, file_prefix: str, random_seed: int = None):
    """Ward dendrogram of Node2Vec embeddings, with leaf labels coloured by stratum (warns clearly if node2vec is unavailable).

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
    print("\n<<< Figure 6 of 7: Node2Vec dendrogram - clustering by learned position >>>")
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
                lvl = node_df.loc[term, _strata.COLUMN]
                if lvl in STRATA_COLOR_MAP: lbl.set_color(STRATA_COLOR_MAP[lvl])

        legend_patches = [mpatches.Patch(color=color, label=level)
                          for level, color in STRATA_COLOR_MAP.items()]
        ax.legend(handles=legend_patches, title='Stratum (Label Colour)', loc='upper left', bbox_to_anchor=(0.82, 1.0))
        plt.xlabel("Cluster Distance", fontsize=12)
        plt.ylabel("MeSH Terms", fontsize=12)
        plt.grid(axis='x', linestyle='--', alpha=0.5)

        save_high_res("Dendrogram_Node2Vec_Strata", output_dir, file_prefix)
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
        print("<" * 30 + ">" * 30)
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


# =============================================================================
# FIGURE 7 - THE NETWORK ITSELF
# =============================================================================

# Metrics offered for colouring, in the order they are preferred when the one
# asked for is not in the network. Every one of these is a per-node float
# written by the network or relevance step.
NETWORK_COLOUR_METRICS = [
    'MRS_pagerank_centrality',
    'MRS_betweenness_centrality',
    'MRS_eigenvector_centrality',
    'MRS_pagerank_subgraph_centrality',
    'MRS_betweenness_subgraph_centrality',
    'MRS_eigenvector_subgraph_centrality',
    'pagerank_centrality',
    'betweenness_centrality',
    'eigenvector_centrality',
    'pagerank_subgraph_centrality',
    'betweenness_subgraph_centrality',
    'eigenvector_subgraph_centrality',
    'adjusted_node_weight',
    'article_count',
    'degree',
    'clustering_coefficient',
    'major_topic_proportion',
    'article_count_rank_normalized',
    'rank_norm_mean_cit',
    'rank_norm_median_cit',
    'rank_norm_total_cit',
]

# load_and_prepare_data renames three columns on the way in. A user picking a
# metric sees the name as it appears in the JSON, so both spellings resolve.
_METRIC_ALIASES = {
    'adjusted_node_weight': 'ARS',
    'MRS_betweenness_centrality': 'MRS_betweenness',
    'MRS_pagerank_centrality': 'MRS_pagerank',
}


def _network_scale(n_nodes):
    """Canvas, node size, label size and edge width for a graph of this size.

    Interpolated rather than stepped, so a 90-node network does not look like a
    50-node one stretched. The canvas is capped: past a few hundred terms this
    is a picture of the shape of the thing, and reading individual labels is a
    job for Cytoscape - which is what the caption says.
    """
    n = max(n_nodes, 1)
    side = float(np.clip(7.0 + 1.35 * np.sqrt(n), 9.0, 40.0))
    # Markers big enough that the colour is legible at a glance - the colour is
    # the whole point of the figure, and a 3-pixel dot carries none of it.
    node_size = float(np.clip(52000.0 / n, 90.0, 2200.0))
    font_size = float(np.clip(135.0 / np.sqrt(n), 3.0, 12.0))
    edge_width = float(np.clip(20.0 / np.sqrt(n), 0.15, 1.4))
    edge_alpha = float(np.clip(26.0 / np.sqrt(n), 0.10, 0.45))
    # Spring constant: the bigger the graph, the more room each node needs
    # relative to the unit square, or the middle collapses into a hairball.
    k = float(np.clip(4.2 / np.sqrt(n), 0.06, 0.85))
    iterations = int(np.clip(12000 / np.sqrt(n), 150, 900))
    # How far under its node a label sits. Derived from the marker itself, not
    # picked: matplotlib sizes markers by AREA in points squared, so the radius
    # is sqrt(area/pi), and the canvas after _spread covers 2.0 data units
    # across side*72 points. A constant drop looked right at one network size
    # and printed the label through the marker at every other.
    radius_pts = float(np.sqrt(node_size / np.pi))
    label_drop = radius_pts / (side * 72.0) * 2.0 + 0.006
    return dict(side=side, node_size=node_size, font_size=font_size,
                edge_width=edge_width, edge_alpha=edge_alpha, k=k,
                iterations=iterations, label_drop=label_drop)


def _spread(pos):
    """Push a layout out to fill its canvas in both directions.

    spring_layout returns something roughly circular, so a square canvas is
    left with empty corners while the middle is crowded. Rescaling each axis
    to its own full range buys back real separation for nothing.
    """
    if not pos:
        return pos
    xs = np.array([p[0] for p in pos.values()], dtype=float)
    ys = np.array([p[1] for p in pos.values()], dtype=float)
    span_x = xs.max() - xs.min() or 1.0
    span_y = ys.max() - ys.min() or 1.0
    return {n: np.array([2.0 * (p[0] - xs.min()) / span_x - 1.0,
                         2.0 * (p[1] - ys.min()) / span_y - 1.0])
            for n, p in pos.items()}


def resolve_colour_metric(node_df, requested=''):
    """The metric actually used to colour, given what the network carries.

    Returns (column in node_df, name to show). A metric that is absent - the
    subgraph centralities when the ground-truth analysis was off, say - falls
    back rather than failing the figure, and says which it used.
    """
    requested = (requested or '').strip()
    candidates = ([requested] if requested else []) + NETWORK_COLOUR_METRICS
    for name in candidates:
        for column in (_METRIC_ALIASES.get(name, name), name):
            if column in node_df.columns and pd.api.types.is_numeric_dtype(node_df[column]):
                if node_df[column].notna().any():
                    return column, name
    return None, requested


def plot_network_graph(G, node_df, output_dir: str, file_prefix: str,
                       metric: str = '', random_seed: int = 42):
    """Figure 7: the consensus network drawn, nodes coloured by one metric.

    Deliberately a still image. It exists so a reader can see the shape of the
    network - where the hubs are, how the communities sit against each other,
    whether the metric they chose picks out anything structural - without
    loading anything into another program. Anything closer than that is a job
    for Cytoscape, and the caption says so.

    Colour is viridis, min-max scaled across the nodes actually drawn, so the
    full range of the colour map is always used however narrow the spread is.
    """
    print("\n<<< Figure 7 of 7: The network, coloured by a chosen metric >>>")
    try:
        if G is None or G.number_of_nodes() == 0:
            print("    [!] The network has no nodes. Skipping Figure 7.")
            return

        column, shown = resolve_colour_metric(node_df, metric)
        if column is None:
            print(f"    [!] No numeric metric available to colour by "
                  f"(asked for {metric!r}). Skipping Figure 7.")
            return
        if metric and shown != metric:
            print(f"    [!] '{metric}' is not in this network - colouring by "
                  f"'{shown}' instead.")

        # Draw the largest connected component. The final network is already an
        # LCC, but a filtered or hand-edited one may not be, and a scatter of
        # detached singletons around the edge wastes most of the canvas.
        if not nx.is_connected(G):
            largest = max(nx.connected_components(G), key=len)
            dropped = G.number_of_nodes() - len(largest)
            G = G.subgraph(largest).copy()
            print(f"    Drawing the largest connected component "
                  f"({G.number_of_nodes():,} of {G.number_of_nodes() + dropped:,} terms).")

        nodes = list(G.nodes())
        values = pd.to_numeric(node_df.reindex(nodes)[column], errors='coerce')
        if values.notna().sum() == 0:
            print(f"    [!] '{shown}' is empty for every node. Skipping Figure 7.")
            return
        values = values.fillna(values.min())

        lo, hi = float(values.min()), float(values.max())
        if hi <= lo:
            # Every node identical: a min-max scale is undefined, so colour them
            # all mid-viridis rather than divide by zero.
            print(f"    [!] Every node has the same '{shown}' ({lo:g}); "
                  f"colouring uniformly.")
            scaled = np.full(len(nodes), 0.5)
        else:
            scaled = (values.to_numpy(dtype=float) - lo) / (hi - lo)

        s = _network_scale(len(nodes))
        print(f"    {len(nodes):,} terms, {G.number_of_edges():,} relations; "
              f"colouring by {shown}")

        # Kamada-Kawai separates a few hundred nodes far better than a spring
        # layout does, because it works from graph distances rather than
        # repulsion alone. It is O(n^2) though, so past a few hundred terms the
        # spring layout is the only affordable choice.
        layout_used = 'Kamada-Kawai'
        pos = None
        if len(nodes) <= 400:
            try:
                pos = nx.kamada_kawai_layout(G, weight=None)
            except Exception:
                pos = None
        if pos is None:
            layout_used = 'spring'
            pos = nx.spring_layout(G, k=s['k'], iterations=s['iterations'],
                                   seed=random_seed, weight=None)
        else:
            # Relax it briefly with the spring model to open out any pairs the
            # distance solution left sitting on top of each other.
            pos = nx.spring_layout(G, pos=pos, k=s['k'], iterations=60,
                                   seed=random_seed, weight=None)
        pos = _spread(pos)
        print(f"    layout: {layout_used}")

        fig, ax = plt.subplots(figsize=(s['side'], s['side']))
        cmap = plt.get_cmap('viridis')

        nx.draw_networkx_edges(G, pos, ax=ax, width=s['edge_width'],
                               edge_color='#9aa0a6', alpha=s['edge_alpha'])
        drawn = nx.draw_networkx_nodes(
            G, pos, ax=ax, nodelist=nodes, node_size=s['node_size'],
            node_color=scaled, cmap=cmap, vmin=0.0, vmax=1.0,
            linewidths=0.5, edgecolors='#2b2e31')
        drawn.set_zorder(2)

        # Labels sit UNDER their node rather than on it. Printed over the
        # markers they hid the very colour the figure exists to show.
        label_pos = {n: (p[0], p[1] - s['label_drop']) for n, p in pos.items()}
        texts = nx.draw_networkx_labels(
            G, label_pos, ax=ax, font_size=s['font_size'], font_color='#16181a',
            font_family='DejaVu Sans',
            verticalalignment='top')
        for t in texts.values():
            t.set_zorder(3)
            # A white halo, so a label crossing an edge or a neighbouring node
            # stays readable without hiding what is underneath.
            t.set_path_effects([
                path_effects.Stroke(linewidth=max(1.4, s['font_size'] / 3.0),
                                    foreground='white', alpha=0.85),
                path_effects.Normal()])

        # The colour bar carries the real values, not the 0-1 the colours were
        # scaled to - the scaling is a drawing device, not the measurement.
        mappable = plt.cm.ScalarMappable(cmap=cmap,
                                         norm=plt.Normalize(vmin=lo, vmax=hi))
        cbar = fig.colorbar(mappable, ax=ax, fraction=0.035, pad=0.02)
        cbar.set_label(shown, fontsize=11)

        ax.set_title(
            f"Consensus network - {len(nodes):,} MeSH terms, "
            f"{G.number_of_edges():,} co-occurrence relations\n"
            f"colour: {shown} (min-max scaled, viridis)",
            fontsize=13, fontweight='bold')
        ax.text(0.5, -0.015,
                "An overview, not a working view - open the network JSON in "
                "Cytoscape to inspect it properly.",
                transform=ax.transAxes, ha='center', va='top',
                fontsize=9, color='#5f6368')
        ax.set_axis_off()
        fig.tight_layout()

        save_high_res("Network_Graph", output_dir, file_prefix)
    except Exception as e:
        print(f"[!] Error generating Figure 7 (network graph): {e}")
    finally:
        plt.close('all')
