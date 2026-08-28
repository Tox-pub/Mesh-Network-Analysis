"""Draw Figure 7 at several sizes and metrics, and check what comes out."""
import os
import sys
import tempfile

sys.path.insert(0, 'src')
import networkx as nx
import numpy as np

from mesh_aop.viz import (plot_network_graph, resolve_colour_metric,
                          _network_scale, NETWORK_COLOUR_METRICS,
                          configure_output, load_and_prepare_data)

FAILS = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + (f'   -- {extra}' if extra and not ok else ''))
    if not ok:
        FAILS.append(msg)


configure_output(110, 'png')
out = tempfile.mkdtemp()

NET = 'data/processed/DAC_Mesh_final_network_with_relevance.json'
ANN = 'data/raw/aop_annotations_master.csv'
node_df, edge_df, G = load_and_prepare_data(NET, ANN)
print(f'\n  real network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges')

print('\n=== the scaling adapts ===')
prev = None
for n in (20, 60, 173, 500, 2000):
    s = _network_scale(n)
    print(f'   n={n:<5} canvas {s["side"]:.1f}in  node {s["node_size"]:.0f}  '
          f'label {s["font_size"]:.1f}pt  k {s["k"]:.3f}  iters {s["iterations"]}')
    if prev:
        ck(s['node_size'] <= prev['node_size'], f'n={n}: nodes no larger than at the smaller size')
        ck(s['font_size'] <= prev['font_size'], f'n={n}: labels no larger')
        ck(s['side'] >= prev['side'], f'n={n}: canvas no smaller')
    prev = s
ck(_network_scale(20)['font_size'] > _network_scale(2000)['font_size'],
   'a small network gets much bigger labels than a huge one')
ck(8.0 <= _network_scale(100000)['side'] <= 40.0,
   'the canvas is capped for absurd sizes',
   f"{_network_scale(100000)['side']:.1f} in")

print('\n=== metric resolution ===')
col, shown = resolve_colour_metric(node_df, 'MRS_pagerank_centrality')
ck(col is not None, 'the default metric resolves', str(col))
col2, shown2 = resolve_colour_metric(node_df, 'not_a_real_metric')
ck(col2 is not None and shown2 != 'not_a_real_metric',
   'an unknown metric falls back rather than failing', f'{col2} / {shown2}')
col3, shown3 = resolve_colour_metric(node_df, 'adjusted_node_weight')
ck(col3 == 'ARS', 'a renamed column resolves through its alias', str(col3))
for m in NETWORK_COLOUR_METRICS:
    c, _ = resolve_colour_metric(node_df, m)
    assert c is not None, m
ck(True, f'all {len(NETWORK_COLOUR_METRICS)} offered metrics resolve to something')

print('\n=== draw the real network, three different metrics ===')
for metric in ('MRS_pagerank_centrality', 'degree', 'article_count'):
    d = os.path.join(out, metric)
    plot_network_graph(G, node_df, d, 'DAC_Mesh', metric=metric, random_seed=42)
    f = os.path.join(d, 'DAC_Mesh_Network_Graph.png')
    ok = os.path.exists(f) and os.path.getsize(f) > 50_000
    ck(ok, f'drawn for {metric}', f'{f} {os.path.getsize(f) if os.path.exists(f) else 0} bytes')

print('\n=== it is deterministic for a given seed ===')
a = os.path.join(out, 'seed_a'); b = os.path.join(out, 'seed_b')
plot_network_graph(G, node_df, a, 'X', metric='degree', random_seed=7)
plot_network_graph(G, node_df, b, 'X', metric='degree', random_seed=7)
ck(open(os.path.join(a, 'X_Network_Graph.png'), 'rb').read()
   == open(os.path.join(b, 'X_Network_Graph.png'), 'rb').read(),
   'the same seed draws the identical image')

print('\n=== small and large synthetic networks ===')
import pandas as pd
for n, name in ((18, 'tiny'), (420, 'large')):
    H = nx.barabasi_albert_graph(n, 2, seed=1)
    H = nx.relabel_nodes(H, {i: f'MeSH term {i}' for i in H.nodes()})
    df = pd.DataFrame({'degree': [d for _, d in H.degree()],
                       'MRS_pagerank_centrality': list(nx.pagerank(H).values())},
                      index=list(H.nodes()))
    d = os.path.join(out, name)
    plot_network_graph(H, df, d, name, metric='MRS_pagerank_centrality', random_seed=42)
    f = os.path.join(d, f'{name}_Network_Graph.png')
    ck(os.path.exists(f) and os.path.getsize(f) > 20_000,
       f'{name} network ({n} nodes) draws',
       f'{os.path.getsize(f) if os.path.exists(f) else 0} bytes')

print('\n=== the awkward cases do not crash ===')
H = nx.path_graph(5)
H = nx.relabel_nodes(H, {i: f'T{i}' for i in H.nodes()})
flat = pd.DataFrame({'m': [1.0] * 5}, index=list(H.nodes()))
plot_network_graph(H, flat, os.path.join(out, 'flat'), 'F', metric='m')
ck(os.path.exists(os.path.join(out, 'flat', 'F_Network_Graph.png')),
   'every node the same value (no min-max range) still draws')

D = nx.Graph(); D.add_nodes_from(['a', 'b']); D.add_edge('a', 'b')
D.add_nodes_from(['x', 'y']); D.add_edge('x', 'y')
dfd = pd.DataFrame({'m': [1.0, 2.0, 3.0, 4.0]}, index=['a', 'b', 'x', 'y'])
plot_network_graph(D, dfd, os.path.join(out, 'disc'), 'D', metric='m')
ck(os.path.exists(os.path.join(out, 'disc', 'D_Network_Graph.png')),
   'a disconnected graph draws its largest component')

E = nx.Graph()
plot_network_graph(E, pd.DataFrame(), os.path.join(out, 'empty'), 'E', metric='m')
ck(not os.path.exists(os.path.join(out, 'empty', 'E_Network_Graph.png')),
   'an empty network is skipped, not crashed')

print(f'\n  images written under: {out}')
print()
print('FAILED' if FAILS else 'FIGURE 7 WORKS')
for f in FAILS:
    print('  -', f)
sys.exit(1 if FAILS else 0)
