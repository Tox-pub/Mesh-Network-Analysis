"""Every network JSON the pipeline produces: are the attributes there, and the right type?

Checks each file against what its stage is supposed to write, and checks that
each stage's output carries everything the next stage reads from it.
"""
import json
import os
import sys

sys.path.insert(0, 'src')
FAILS = []
NUM = {'int', 'float'}


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + (f'   -- {extra}' if extra and not ok else ''))
    if not ok:
        FAILS.append(msg)


def profile(path):
    """attribute -> (count, {types}) for nodes and edges."""
    with open(path, encoding='utf-8') as fh:
        d = json.load(fh)
    el = d.get('elements', {})
    out = {}
    for kind in ('nodes', 'edges'):
        items = el.get(kind, [])
        prof = {}
        for it in items:
            for k, v in it.get('data', {}).items():
                e = prof.setdefault(k, [0, set()])
                e[0] += 1
                e[1].add(type(v).__name__)
        out[kind] = (len(items), prof)
    return out


# What each stage is contracted to write, and the type each value must be.
# 'num' means int or float - counts come back as either depending on the value.
CONTRACTS = {
    'full_network': {
        'nodes': {'id': 'str', 'label': 'str', 'article_count': 'num',
                  'degree': 'num', 'term_type': 'str', 'is_major_topic': 'bool'},
        'edges': {'source': 'str', 'target': 'str', 'cooccurrence_count': 'num',
                  'weight': 'num'},
    },
    'consensus_lcc': {
        'nodes': {'id': 'str', 'label': 'str', 'degree': 'num',
                  'filtered_louvain_community_id': 'num',
                  'unfiltered_louvain_community_id': 'num',
                  'betweenness_centrality': 'num', 'pagerank_centrality': 'num',
                  'eigenvector_centrality': 'num', 'clustering_coefficient': 'num'},
        'edges': {'source': 'str', 'target': 'str', 'weight': 'num',
                  'edge_betweenness': 'num', 'mlf_p_value': 'num'},
    },
    'final_network': {
        # 'stratum' is what is written now; 'aop_level' is what files annotated
        # before the rename carry, and both are still read.
        'nodes': {'id': 'str', 'label': 'str', 'stratum': None,     # optional
                  'MRS_betweenness_centrality': 'num',
                  'MRS_pagerank_centrality': 'num',
                  'MRS_eigenvector_centrality': 'num',
                  'adjusted_node_weight': 'num',
                  'generation_weight': 'num',
                  'rank_norm_mean_cit': 'num', 'rank_norm_median_cit': 'num',
                  'rank_norm_total_cit': 'num',
                  'article_count_rank_normalized': 'num',
                  'major_topic_proportion': 'num'},
        'edges': {'source': 'str', 'target': 'str', 'weight': 'num',
                  'normalized_cooccurrence': 'num',
                  'log1p_cooccurrence_minmax': 'num',
                  'cumulative_weight_original': 'num'},
    },
    'glf_subgraph': {'nodes': {'id': 'str'}, 'edges': {'source': 'str', 'target': 'str'}},
    'sa_subgraph': {'nodes': {'id': 'str'}, 'edges': {'source': 'str', 'target': 'str'}},
}

FILES = {
    'consensus_lcc': 'data/processed/DAC_Mesh_consensus_lcc_network.json',
    'final_network': 'data/processed/DAC_Mesh_final_network_with_relevance.json',
    'glf_subgraph': 'data/processed/DAC_Mesh_glf_optimal_subgraph.json',
    'sa_subgraph': 'data/processed/DAC_Mesh_sa_optimal_subgraph.json',
}

for key, path in FILES.items():
    print(f'\n=== {key}   {os.path.basename(path)} ===')
    if not os.path.exists(path):
        ck(False, f'{key} is present to check')
        continue
    prof = profile(path)
    contract = CONTRACTS[key]
    for kind in ('nodes', 'edges'):
        count, attrs = prof[kind]
        print(f'  {count:,} {kind}, {len(attrs)} attributes')
        for name, want in contract.get(kind, {}).items():
            if want is None:
                continue
            if name not in attrs:
                ck(False, f'{kind}.{name} present')
                continue
            n, types = attrs[name]
            ok_count = (n == count)
            ok_type = (types <= NUM) if want == 'num' else (types == {want})
            ck(ok_count and ok_type,
               f'{kind}.{name}  ({want})',
               f'{n}/{count} present, types {sorted(types)}')

print('\n=== what the next stage reads must be written by the one before ===')
lcc = profile(FILES['consensus_lcc'])
fin = profile(FILES['final_network'])
carried = ['id', 'label', 'degree', 'filtered_louvain_community_id',
           'betweenness_centrality', 'pagerank_centrality', 'eigenvector_centrality']
for a in carried:
    ck(a in lcc['nodes'][1] and a in fin['nodes'][1],
       f'{a} survives from the consensus network into the final one',
       f"lcc={a in lcc['nodes'][1]} final={a in fin['nodes'][1]}")

print('\n=== the figures only read attributes that exist ===')
src = open('src/mesh_aop/viz.py', encoding='utf-8').read()
final_attrs = (set(fin['nodes'][1]) | set(fin['edges'][1])
               | {'stratum', 'aop_level'})
import re
referenced = set(re.findall(r"(?:node_df|edge_df|df|plot_df|top_diff)\[['\"]([A-Za-z_]+)['\"]\]", src))
unknown = sorted(a for a in referenced
                 if a not in final_attrs and not a.startswith(('diff', 'count', 'x', 'y')))
ck(not unknown, 'every node/edge attribute the figures index is in the network',
   f'not found: {unknown}')

print('\n=== the benchmark weight key names a real attribute ===')
from mesh_aop.config_parser import MeshConfig
key = MeshConfig.__init__  # just to import cleanly
import json as _j
default_key = 'MRS_pagerank_centrality'
ck(default_key in fin['nodes'][1],
   f'the default validation weight key ({default_key}) exists on every node',
   str(fin['nodes'][1].get(default_key)))

print()
print('FAILED' if FAILS else 'EVERY OUTPUT MATCHES ITS CONTRACT')
for f in FAILS:
    print('  -', f)
sys.exit(1 if FAILS else 0)
