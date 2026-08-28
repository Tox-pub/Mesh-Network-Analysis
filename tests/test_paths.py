"""Every path the pipeline resolves, and whether it lands where it should.

Static reasoning is not enough here: the layout changed, several paths are
derived from other paths, and reference mode resolves differently from own-data
mode. So this builds real configurations and interrogates them.
"""
import json
import os
import sys
import tempfile

sys.path.insert(0, 'src')

FAILS = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + (f'   -- {extra}' if extra and not ok else ''))
    if not ok:
        FAILS.append(msg)


def make(box, ref=False, prefix='P'):
    from mesh_aop.config_parser import MeshConfig
    cfgp = os.path.join(box, 'c.json')
    json.dump({"control_flags": {"use_reference_data": ref, "custom_file_prefix": prefix},
               "directories": {"data_dir": os.path.join(box, 'data'),
                               "results_dir": os.path.join(box, 'results')}},
              open(cfgp, 'w'))
    return MeshConfig(config_path=cfgp)


box = tempfile.mkdtemp()
c = make(box)

print('=== 1. every file the pipeline names resolves under a folder it owns ===')
OWNED = {
    'master_db':        'databases',
    'pmids_db':         'databases',
    'mesh_marc':        'raw',
    'mesh_ascii':       'raw',
    'full_network':     'networks',
    'glf_subgraph':     'networks',
    'sa_subgraph':      'networks',
    'consensus_lcc':    'networks',
    'final_network':    'networks',
    'cleaned_db':       'databases',
    'relevance_db':     'databases',
    'mesh_terms_csv':   'processed',
    'failed_mesh':      'logs',
    'empty_mesh':       'logs',
    'log_file':         'logs',
}
where = {
    'raw':        c.active_raw_dir,
    'processed':  c.active_source_dir,
    'networks':   c.networks_dir,
    'databases':  c.databases_dir,
    'logs':       c.log_dir,
}
for key, folder in OWNED.items():
    got = c.files[key].parent
    ck(got == where[folder], f'{key:<16} -> {folder}', f'landed in {got}')

# The two that are deliberately elsewhere.
import mesh_aop
pkg = os.path.dirname(os.path.abspath(mesh_aop.__file__))
ck(str(c.files['mesh_stopwords_py']) == os.path.join(pkg, 'mesh_stop_words.py'),
   'mesh_stopwords_py -> the package module the pipeline actually imports',
   str(c.files['mesh_stopwords_py']))
ck(c.files['mesh_stopwords_py'].exists(), '  ...and it exists, so Step 1 no longer rebuilds every run')
ck(c.files['annotations'].name == 'aop_annotations_master.csv',
   'annotations -> the AOP stratum dictionary', str(c.files['annotations']))

print('\n=== 2. no two artefacts collide on one path ===')
paths = {}
for key, p in c.files.items():
    paths.setdefault(str(p), []).append(key)
dupes = {p: k for p, k in paths.items() if len(k) > 1}
ck(not dupes, 'every artefact has its own path', str(dupes))

print('\n=== 3. derived paths agree with the file they derive from ===')
opt_hist = os.path.join(os.path.dirname(c.files['full_network']),
                        f'{c.prefix}_optimization_history.json')
ck(os.path.dirname(opt_hist) == str(c.networks_dir),
   'the optimisation history sits beside the network it describes', opt_hist)

print('\n=== 4. the results tree ===')
ck(c.figures_dir.parent == c.results_dir, 'figures/ is under results')
ck(c.secondary_dir.parent == c.results_dir, 'secondary_analysis/ is under results')
ck(c.benchmark_dir.parent == c.results_dir, 'benchmark/ is under results')
ck(c.networks_dir.parent == c.active_source_dir, 'networks/ is under processed, NOT results')
ck(c.databases_dir.parent == c.active_source_dir, 'databases/ is under processed')
ck(not str(c.results_dir).startswith(str(c.active_raw_dir)),
   'results are not inside the data folder')

print('\n=== 4b. the databases never move with the reference checkbox ===')
import json as _json
_p = os.path.join(box, 'c.json')
_d = _json.load(open(_p)); _d['control_flags']['use_reference_data'] = True
_json.dump(_d, open(_p, 'w'))
from mesh_aop.config_parser import MeshConfig as _MC
_ref = _MC(config_path=_p)
# The FOLDER, not the filename: the checkbox also forces the project prefix to
# DAC_Mesh, so the per-project database names legitimately differ. What must not
# change is that they stay in the user's databases folder rather than following
# the corpus into the read-only program directory.
for _k in ('master_db', 'pmids_db', 'cleaned_db', 'relevance_db'):
    ck(_ref.files[_k].parent == c.databases_dir,
       f'{_k} still lives in the user databases folder with the checkbox on',
       str(_ref.files[_k]))
ck(_ref.files['master_db'] == c.files['master_db'],
   'and the master database is the very same file (it carries no prefix)')
_d['control_flags']['use_reference_data'] = False
_json.dump(_d, open(_p, 'w'))

print('\n=== 5. writability: everything the pipeline writes to exists and is writable ===')
for name, d in (('raw', c.active_raw_dir), ('processed', c.active_source_dir),
                ('networks', c.networks_dir), ('databases', c.databases_dir),
                ('figures', c.figures_dir), ('secondary_analysis', c.secondary_dir),
                ('benchmark', c.benchmark_dir), ('logs', c.log_dir)):
    ok = d.is_dir() and os.access(d, os.W_OK)
    ck(ok, f'{name}/ exists and is writable', str(d))

print('\n=== 6. reference mode resolves its own corpus, and writes nowhere in it ===')
box2 = tempfile.mkdtemp()
cr = make(box2, ref=True, prefix='DAC_Mesh')
for key in ('final_network', 'consensus_lcc', 'glf_subgraph', 'sa_subgraph'):
    p = cr.files[key]
    ck(p.exists(), f'reference {key} is found', str(p))
ck(cr.results_dir != cr.active_source_dir,
   'reference results still go to the user results folder', str(cr.results_dir))
ck(not (cr.active_source_dir / 'networks').exists(),
   'no empty subfolders created inside the shipped corpus')

print('\n=== 7. a legacy project keeps working ===')
box3 = tempfile.mkdtemp()
c3 = make(box3)
legacy_names = {
    'final_network':  f'{c3.prefix}_final_network_with_relevance.json',
    'consensus_lcc':  f'{c3.prefix}_consensus_lcc_network.json',
    'cleaned_db':     f'{c3.prefix}_cleaned_pmids.db',
    'relevance_db':   f'{c3.prefix}_mean_relevancy.db',
}
for name in legacy_names.values():
    (c3.active_source_dir / name).write_text('x')
c3b = make(box3)
for key, name in legacy_names.items():
    ck(c3b.files[key] == c3b.active_source_dir / name,
       f'legacy {key} used where it already is', str(c3b.files[key]))
ck(c3b.files['glf_subgraph'].parent == c3b.networks_dir,
   'anything not yet built still goes to the new folder')

print('\n=== 8. two projects never collide ===')
box4 = tempfile.mkdtemp()
a = make(box4, prefix='ProjA')
b = make(box4, prefix='ProjB')
shared = {str(v) for v in a.files.values()} & {str(v) for v in b.files.values()}
expected_shared = {str(a.files['master_db']), str(a.files['mesh_terms_csv']),
                   str(a.files['mesh_stopwords_py']), str(a.files['annotations']),
                   str(a.files['mesh_marc']), str(a.files['mesh_ascii'])}
unexpected = shared - expected_shared
ck(not unexpected, 'only the shared reference assets are shared between prefixes',
   str(unexpected))

print()
print('FAILED' if FAILS else 'EVERY PATH RESOLVES WHERE IT SHOULD')
for f in FAILS:
    print('  -', f)
sys.exit(1 if FAILS else 0)
