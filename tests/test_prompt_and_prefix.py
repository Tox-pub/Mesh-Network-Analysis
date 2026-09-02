"""Three ways the pipeline failed a user who had done nothing wrong.

1. A brand-new project prefix was told it had "already been used", naming
   benchmark files that belonged to a different project. The benchmark glob was
   the only one in STEP_RESULT_GLOBS with no {p} in it, so it matched every file
   any run had ever left in that folder.

2. Pause for annotation killed the run. The pipeline asked whether to sync AOP
   levels back to the master library with input(), and the Workbench runs it as
   a subprocess with no console, so input() got EOF and the step died - with the
   question printed just above the traceback, which read as though answering it
   was what crashed.

3. Target nodes typed as the help showed them - "Skin; Haptens" - matched
   nothing. The quotes were stripped per entry AFTER splitting on ';', so the
   two halves were '"Skin' and 'Haptens"', neither a matched pair.
"""
import io
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


print('=== 1. a fresh prefix is not warned about another project\'s files ===')
from mesh_aop import integrity                                     # noqa: E402
from mesh_aop.config_parser import MeshConfig                      # noqa: E402

box = tempfile.mkdtemp()


def cfg(prefix):
    p = os.path.join(box, f'{prefix}.json')
    json.dump({"control_flags": {"custom_file_prefix": prefix},
               "directories": {"data_dir": os.path.join(box, 'data'),
                               "results_dir": os.path.join(box, 'results')}},
              open(p, 'w'))
    c = MeshConfig(config_path=p)
    c.refresh_paths()
    return c


old = cfg('DAC_Mesh')
# what the earlier project left behind, in the real folder layout
for sub, name in (('inputs', 'DAC_Mesh_ground_truth.xlsx'),
                  ('ranking', 'DAC_Mesh_benchmark_results.json'),
                  ('ranking_validation', 'DAC_Mesh_validation_report.xlsx'),
                  ('network_validation', 'DAC_Mesh_gt_network_validation.xlsx')):
    d = os.path.join(str(old.benchmark_dir), sub)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name), 'w') as fh:
        fh.write('x' * 1000)

fresh = cfg('Derm_Mesh_run')
at_risk = integrity.would_overwrite(fresh, 'benchmark')
ck(at_risk == [], f'a never-used prefix risks nothing: {[a[0] for a in at_risk]}')

same = integrity.would_overwrite(old, 'benchmark')
ck(len(same) >= 1, f'the prefix that DID write them is still warned: {[a[0] for a in same]}')
found = sum(1 for a in same for _ in [a])
ck(any('Benchmark' in a[0] for a in same), 'and the warning names the benchmark outputs')

# a prefix that is a strict prefix of another must not inherit its files
sub_pref = cfg('DAC')
ck(integrity.would_overwrite(sub_pref, 'benchmark') == [],
   "'DAC' is not warned about 'DAC_Mesh' files")

print('\n=== 2. a question nobody can answer does not kill the run ===')
from mesh_aop.cli import _console_answer                           # noqa: E402


class NoConsole:
    def isatty(self):
        return False

    def readline(self):
        raise EOFError


real_stdin = sys.stdin
try:
    sys.stdin = NoConsole()
    ans = _console_answer('  [?] sync? [y/n]: ')
    ck(ans == '', f'returns the default instead of raising: {ans!r}')
    ans2 = _console_answer('  [?] go? [y/n]: ', default='n')
    ck(ans2 == 'n', f'an explicit default is honoured: {ans2!r}')

    sys.stdin = None
    ck(_console_answer('  [?] anything? ') == '', 'stdin=None is survivable too')
finally:
    sys.stdin = real_stdin

# and the sync itself must not fire when nobody said yes
import inspect                                                     # noqa: E402
from mesh_aop import cli                                           # noqa: E402
src = inspect.getsource(cli._sync_run_to_master)
ck('input(' not in src, 'the sync no longer calls input() directly')
ck('_console_answer' in src, 'it goes through the guarded helper')
ck("ans in ['y', 'yes']" in src, 'and still only syncs on an explicit yes')

# no other prompt is left unguarded in the pipeline module
cli_src = io.open('src/mesh_aop/cli.py', encoding='utf-8').read()
import re                                                          # noqa: E402
bare = [ln.strip() for ln in cli_src.splitlines()
        if re.search(r'(?<![\w_])input\s*\(', ln)
        and not ln.strip().startswith('#')
        and 'return input(' not in ln          # the one inside the helper
        and not ln.strip().startswith('"""')
        and 'input() ' not in ln]              # prose in docstrings
ck(not bare, f'no bare input() left in cli.py: {bare}')

print('\n=== 3. quotes typed round the list are understood ===')


def targets(value, key='target_nodes'):
    c = cfg('Q')
    c.update('secondary_analysis', key, value)
    return c.get('secondary_analysis', key)


for typed, want in (
    ('"Skin; Haptens"',           'Skin;Haptens'),   # the case that failed
    ("'Skin; Haptens'",           'Skin;Haptens'),
    ('"Skin"; "Haptens"',         'Skin;Haptens'),   # quoted individually
    ('Skin; Haptens',             'Skin;Haptens'),   # unquoted
    ('  Skin ;  Haptens  ',       'Skin;Haptens'),
    ('"Dermatitis, Allergic Contact"', 'Dermatitis, Allergic Contact'),
    ('Skin',                      'Skin'),
    ('',                          ''),
):
    got = targets(typed)
    ck(got == want, f'{typed!r:34} -> {got!r}', f'expected {want!r}')

for typed, want in (
    ('"Skin - Dermatitis, Allergic Contact"', 'Skin - Dermatitis, Allergic Contact'),
    ('"Skin - Haptens; Skin - Allergens"',    'Skin - Haptens;Skin - Allergens'),
    ('Skin - Haptens',                        'Skin - Haptens'),
):
    got = targets(typed, 'target_edges')
    ck(got == want, f'edges {typed!r:40} -> {got!r}', f'expected {want!r}')

# the help must not show an example that breaks when typed verbatim
from mesh_workbench.settings_schema import SETTINGS                # noqa: E402
tn = [f for _t, f in SETTINGS if f.key.endswith('target_nodes')][0]
ck('Skin; Haptens' in (tn.note or ''), 'the field shows the example')
# It used to tell people NOT to use quotes, which was backwards once quotes
# started working. The help now says they are accepted either way round.
ck('Quotes are accepted' in (tn.note or ''),
   'and says quotes are accepted rather than forbidden')

print('\n=== 4. the sync answer comes from the application, not a dead prompt ===')
import csv                                                         # noqa: E402
import subprocess                                                  # noqa: E402
from mesh_workbench.app import Workbench                           # noqa: E402

# the flag exists and only accepts the three answers
help_out = subprocess.run(
    [sys.executable, '-m', 'mesh_aop.cli', '--help'],
    capture_output=True, text=True,
    env={**os.environ, 'PYTHONPATH': 'src'}).stdout
ck('--sync-annotations {ask,yes,no}' in help_out,
   'the pipeline takes --sync-annotations {ask,yes,no}')

# a decided answer must not consult the console at all
sync_src = inspect.getsource(cli._sync_run_to_master)
ck("decided in ('yes', 'no')" in sync_src,
   'a decided answer short-circuits the prompt')
ck(sync_src.index("decided in ('yes', 'no')") < sync_src.index('_console_answer'),
   'and it is checked BEFORE the console is asked')

# the GUI asks on both routes into step 4
start_src = inspect.getsource(Workbench.start_run)
# Only a run that STARTS at the figures asks up front. Under "all" the merge
# happens after the annotation pause, an hour or more in, so asking at the
# start puts the question before the user has seen a single term to assign.
ck("step == 'viz'" in start_src,
   'a run starting at the figures asks up front')
ck("step in ('viz', 'all')" not in start_src,
   'a full run does NOT ask at the very beginning')
ck("'--sync-annotations' not in extra" in start_src,
   'and it is not asked twice when the answer is already known')
ck('_ask_annotation_sync' in start_src, 'via the dialog helper')

ask_src = inspect.getsource(Workbench._ask_annotation_sync)
ck('askyesnocancel' in ask_src, 'the dialog offers yes, no and cancel')
ck("return 'no'" in ask_src, 'AFK runs answer no without a dialog')
ck('return None' in ask_src, 'cancel means nothing runs')

# the unassigned count the dialog shows
anno = os.path.join(box, 'run_annotations.csv')
with open(anno, 'w', encoding='utf-8', newline='') as fh:
    w = csv.writer(fh, delimiter=';')
    w.writerow(['mesh_term', 'aop_level'])
    w.writerow(['Skin', 'Key Event'])
    w.writerow(['Haptens', 'Unassigned'])
    w.writerow(['Allergens', 'unassigned'])
ck(Workbench._unassigned_count(anno) == 2,
   f'counts Unassigned case-insensitively: {Workbench._unassigned_count(anno)}')
ck(Workbench._unassigned_count(os.path.join(box, 'nope.csv')) == 0,
   'a missing file counts zero rather than raising')

print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('all prompt and prefix checks passed')
