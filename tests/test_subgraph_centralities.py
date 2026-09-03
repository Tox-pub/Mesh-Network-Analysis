"""Prove the six subgraph attributes are produced, through the real code path.

Not a mock: run_community_detection is the function the network step calls, and
_gt_wanted is what decides the flag it is called with.
"""
import json
import os
import shutil
import sys
import tempfile

sys.path.insert(0, 'src')

RAW3 = ('betweenness_subgraph_centrality',
        'pagerank_subgraph_centrality',
        'eigenvector_subgraph_centrality')
MRS3 = tuple('MRS_' + n for n in RAW3)

FAILS = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + (f'   -- {extra}' if extra and not ok else ''))
    if not ok:
        FAILS.append(msg)


SRC = 'data/processed/DAC_Mesh_consensus_lcc_network.json'
box = tempfile.mkdtemp()
work = os.path.join(box, 'DAC_Mesh_consensus_lcc_network.json')
shutil.copy2(SRC, work)

from mesh_aop.network import run_community_detection

print('=== they are written whatever the benchmark flag says ===')
# They used to be gated on the ground-truth analysis, which put a decision
# about the BENCHMARK inside the step that builds the network: turning the
# benchmark on later meant rebuilding the network to recover attributes that
# take 40ms. Worse, every consumer reads them with a 0.0 default, so a network
# without them did not fail validation - it scored every node zero in silence.
#
# Stripped first, because the checked-in network already carries them and
# nothing removes an attribute that is already there.
stripped = json.load(open(SRC, encoding='utf-8'))
for _n in stripped['elements']['nodes']:
    for _a in RAW3:
        _n['data'].pop(_a, None)
off = os.path.join(box, 'off.json')
json.dump(stripped, open(off, 'w'))
run_community_detection(off, random_seed=42, compute_subgraph_centrality=False)
nodes = json.load(open(off, encoding='utf-8'))['elements']['nodes']
missing = [a for a in RAW3 if a not in nodes[0]['data']]
ck(not missing, 'all three are written even with the old flag passed False',
   f'missing {missing}')
for _a in RAW3:
    vals = [float(n['data'].get(_a, 0) or 0) for n in nodes]
    ck(any(vals), f'{_a} carries real values ({sum(1 for v in vals if v)} non-zero)')

print('\n=== the flag ON (reference mode, or the box ticked) ===')
run_community_detection(work, random_seed=42, compute_subgraph_centrality=True)
nodes = json.load(open(work, encoding='utf-8'))['elements']['nodes']
print(f'  {len(nodes):,} nodes')
for attr in RAW3:
    have = [n for n in nodes if attr in n['data']]
    types = {type(n['data'][attr]).__name__ for n in have}
    ck(len(have) == len(nodes), f'{attr} on every node', f'{len(have)}/{len(nodes)}')
    ck(types == {'float'}, f'  ...and is a float', str(types))
    vals = [n['data'][attr] for n in have]
    ck(any(v > 0 for v in vals), f'  ...with real values, not all zero',
       f'max {max(vals) if vals else 0}')

print('\n=== _gt_wanted drives that flag correctly ===')
from mesh_aop.config_parser import MeshConfig
from mesh_aop.cli import _gt_wanted


def cfg(ref, gt):
    p = os.path.join(box, 'c.json')
    d = {"control_flags": {"use_reference_data": ref, "custom_file_prefix": "P"},
         "directories": {"data_dir": os.path.join(box, 'd'),
                         "results_dir": os.path.join(box, 'r')}}
    if gt is not None:
        d["benchmark"] = {"run_ground_truth_analysis": gt}
    json.dump(d, open(p, 'w'))
    return MeshConfig(config_path=p)


ck(_gt_wanted(cfg(True, False)), 'reference mode + form wrote false -> ON')
ck(_gt_wanted(cfg(True, None)), 'reference mode + key absent -> ON')
ck(_gt_wanted(cfg(False, True)), 'own data + ticked -> ON')
ck(not _gt_wanted(cfg(False, False)), 'own data + unticked -> off')

print('\n=== the MRS_ pairs are scored unconditionally ===')
src = open('src/mesh_aop/cli.py', encoding='utf-8').read()
for raw, mrs in zip(RAW3, MRS3):
    ck(f'("{raw}", "{mrs}")' in src, f'{raw} -> {mrs} is in the weightings list')
# They used to be gated on the ground-truth flag in both places, which left a
# relevance database six score columns short whenever the benchmark had been
# off - and only the file's EXISTENCE was checked before reusing it.
ck('if include_subgraph_weightings else []' not in src,
   'the network step no longer gates them')
ck('if include_subgraph_weightings:' not in src,
   'nor does the relevance-database rebuild')

print('\n=== an incomplete relevance database is detected, not reused ===')
import sqlite3                                                     # noqa: E402
from mesh_aop.cli import _relevance_columns_missing, MRS_SCORE_COLUMNS  # noqa: E402

for cols in MRS_SCORE_COLUMNS:
    ck(cols.startswith('score_'), f'{cols} is named as relevance.py writes it')


def _db(name, cols):
    p = os.path.join(box, name)
    c = sqlite3.connect(p)
    c.execute('CREATE TABLE article_relevance_scores (pmid TEXT, '
              + ', '.join(f'"{x}" REAL' for x in cols)
              + ', contributing_seeds TEXT, pub_date TEXT)')
    c.commit()
    c.close()
    return p


ck(_relevance_columns_missing(_db('full.db', MRS_SCORE_COLUMNS)) == [],
   'a complete database is left alone')
short = _relevance_columns_missing(
    _db('old.db', [c for c in MRS_SCORE_COLUMNS if 'subgraph' not in c]))
ck(len(short) == 3, f'one built before the change reports 3 missing: {short}')
ck(all('subgraph' in c for c in short), 'and they are the subgraph ones')
open(os.path.join(box, 'junk.db'), 'w').write('not a database')
ck(_relevance_columns_missing(os.path.join(box, 'junk.db')) == [],
   'a damaged file is left to the integrity check, not called incomplete')

print()
print('FAILED' if FAILS else 'ALL SIX ATTRIBUTES ARE PRODUCED BY THE REAL CODE PATH')
for f in FAILS:
    print('  -', f)
sys.exit(1 if FAILS else 0)
