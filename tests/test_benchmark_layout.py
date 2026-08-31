"""The benchmark folder groups what went in and what each step asked.

Benchmarking runs four things against one network - it keeps the ground truth,
ranks articles, validates that ranking, and validates the network's structure -
and they all used to land in one folder where the only clue to which output came
from which was the filename. Each now has its own folder, and the ground truth
the run actually used is kept beside them: a benchmark number cannot be read
without knowing what it was scored against, and the original lives wherever the
user put it, editable and not named after this project.

Checked here: the folders resolve where they say, they nest under benchmark/
rather than scattering through results/, the kept input carries the project
prefix, an already-kept file is not recopied onto itself, and the Results screen
recognises the new names as this project's.
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


def make(box, prefix='DAC_Mesh'):
    from mesh_aop.config_parser import MeshConfig
    cfgp = os.path.join(box, 'c.json')
    json.dump({"control_flags": {"custom_file_prefix": prefix},
               "directories": {"data_dir": os.path.join(box, 'data'),
                               "results_dir": os.path.join(box, 'results')}},
              open(cfgp, 'w'))
    return MeshConfig(config_path=cfgp)


box = tempfile.mkdtemp()
cfg = make(box)

print('=== 1. every benchmark folder is inside results/benchmark ===')
SUBS = {
    'benchmark_inputs_dir': 'inputs',
    'benchmark_ranking_dir': 'ranking',
    'benchmark_validation_dir': 'ranking_validation',
    'benchmark_network_dir': 'network_validation',
}
for attr, leaf in SUBS.items():
    p = getattr(cfg, attr, None)
    ck(p is not None, f'config exposes {attr}')
    if p is None:
        continue
    ck(p.name == leaf, f'{attr} is named {leaf!r}', f'got {p.name!r}')
    ck(p.parent == cfg.benchmark_dir,
       f'{attr} sits directly under benchmark/', f'parent is {p.parent}')

ck(cfg.benchmark_dir.parent == cfg.results_dir,
   'benchmark/ sits under results/', f'parent is {cfg.benchmark_dir.parent}')
ck(len({str(getattr(cfg, a)) for a in SUBS}) == 4,
   'the four folders are four distinct paths')

print('\n=== 2. the ground truth used is kept, under the project prefix ===')
from mesh_aop.cli import _keep_benchmark_input                      # noqa: E402

for d in (cfg.benchmark_dir, cfg.benchmark_inputs_dir):
    os.makedirs(d, exist_ok=True)

src_dir = os.path.join(box, 'somewhere else')
os.makedirs(src_dir, exist_ok=True)
gt = os.path.join(src_dir, 'oecd_ground_truth_curated.xlsx')
with open(gt, 'wb') as fh:
    fh.write(b'PMID\n123\n456\n')

kept = _keep_benchmark_input(gt, cfg, 'ground_truth')
ck(os.path.exists(kept), 'the kept file exists')
ck(os.path.dirname(kept) == str(cfg.benchmark_inputs_dir),
   'it landed in benchmark/inputs', f'landed in {os.path.dirname(kept)}')
ck(os.path.basename(kept) == 'DAC_Mesh_ground_truth.xlsx',
   'it carries the prefix and keeps the extension',
   f'named {os.path.basename(kept)}')
ck(open(kept, 'rb').read() == b'PMID\n123\n456\n', 'the contents are unchanged')
ck(os.path.exists(gt), 'the original is left where the user put it')

print('\n=== 3. keeping it twice is harmless ===')
again = _keep_benchmark_input(kept, cfg, 'ground_truth')
ck(again == kept, 'a file already in inputs/ is not copied onto itself')
ck(open(kept, 'rb').read() == b'PMID\n123\n456\n', 'and it is not truncated')
again2 = _keep_benchmark_input(gt, cfg, 'ground_truth')
ck(again2 == kept, 'rerunning refreshes the same path rather than making a second')

print('\n=== 4. a missing or unset input is not an error ===')
ck(_keep_benchmark_input(None, cfg, 'negative_control') is None,
   'no negative control configured: returns None, does not raise')
absent = os.path.join(src_dir, 'not_there.csv')
ck(_keep_benchmark_input(absent, cfg, 'negative_control') == absent,
   'a path that does not exist is handed back unchanged')

print('\n=== 5. the results screen claims the kept files for this project ===')
from mesh_workbench.app import Workbench                            # noqa: E402

for name, want in (('DAC_Mesh_ground_truth.xlsx', True),
                   ('DAC_Mesh_negative_control.csv', True),
                   ('DAC_Mesh_benchmark_results.json', True),
                   ('DAC_Mesh_1_ground_truth.xlsx', False),
                   ('Other_ground_truth.xlsx', False)):
    got = Workbench._belongs_to(name, 'DAC_Mesh')
    ck(got is want, f'{name!r} belongs to DAC_Mesh: {got}', f'expected {want}')

print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('all benchmark-layout checks passed')
