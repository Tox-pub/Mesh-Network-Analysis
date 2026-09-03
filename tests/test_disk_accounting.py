"""The Database screen has to account for what is actually using the disk.

Every search builds its own PMID, citation and relevance databases, named for
its project prefix, and nothing removes them when the prefix changes. Three
projects leave nine files behind and the screen said nothing about any of them,
so a user wondering where forty gigabytes went had no way to find out from the
program that put them there.

Checked here: databases belonging to OTHER prefixes are found, grouped and
sized; the current project's are not double-counted; the two folder totals are
reported; and the prefix is parsed correctly, which is the part that is easy to
get wrong - '_cleaned_pmids.db' also ends with '_pmids.db'.
"""
import json
import os
import sys
import tempfile

sys.path.insert(0, 'src')

import tkinter as tk                                               # noqa: E402
from mesh_aop.config_parser import MeshConfig                       # noqa: E402
from mesh_workbench.app import Workbench, _size_gb                  # noqa: E402

FAILS = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + (f'   -- {extra}' if extra else ''))
    if not ok:
        FAILS.append(msg)


MARKERS = ('_pmids.db', '_cleaned_pmids.db', '_mean_relevancy.db')

box = tempfile.mkdtemp()
cfgp = os.path.join(box, 'mesh_config.json')
json.dump({'control_flags': {'custom_file_prefix': 'Current'},
           'directories': {'data_dir': os.path.join(box, 'data'),
                           'results_dir': os.path.join(box, 'res')}},
          open(cfgp, 'w'))
cfg = MeshConfig(config_path=cfgp)
cfg.refresh_paths()

db_dir = str(cfg.databases_dir)
os.makedirs(db_dir, exist_ok=True)
OTHERS = ('DAC_Mesh', 'DAC_Mesh_1', 'Pilot')
for pref in OTHERS + ('Current',):
    for m in MARKERS:
        with open(os.path.join(db_dir, pref + m), 'wb') as fh:
            fh.write(b'x' * 500_000)
with open(os.path.join(db_dir, 'master_mesh_database.db'), 'wb') as fh:
    fh.write(b'x' * 2_000_000)
os.makedirs(str(cfg.active_raw_dir), exist_ok=True)
with open(os.path.join(str(cfg.active_raw_dir), 'desc2025.xml'), 'wb') as fh:
    fh.write(b'x' * 1_000_000)

app = Workbench(os.getcwd(), sys.executable)
app.withdraw()
app.show('setup')
app.update_idletasks()
for w in app.setup_rows.winfo_children():
    w.destroy()

print('=== 1. other projects are found, grouped and sized ===')
app._setup_other_projects(cfg, 'Current')
app.update_idletasks()
rows = [r for r in app.setup_rows.winfo_children() if r.winfo_children()]
ck(len(rows) == 1, f'one row is added: {len(rows)}')

labels = [w for w in rows[0].winfo_children() if isinstance(w, tk.Label)]
head = next(str(w.cget('text')) for w in labels
            if int(w.grid_info().get('column', -1)) == 0
            and int(w.grid_info().get('row', -1)) == 0)
note = next(str(w.cget('text')) for w in labels
            if int(w.grid_info().get('row', -1)) == 1)
size = next(str(w.cget('text')) for w in labels
            if int(w.grid_info().get('column', -1)) == 2)
print(f'    {head}  |  {note}  |  {size}')

ck('other projects' in head.lower(), f'the row says what it is: {head!r}')
ck('9 file(s)' in note, f'nine files counted (3 prefixes x 3): {note!r}')
ck('3 prefix' in note, f'three prefixes: {note!r}')
for p in OTHERS:
    ck(p in note, f'{p} is named')
ck('Current' not in note, "the CURRENT project is not counted as another")

print('\n=== 2. the prefix is parsed off the longest matching suffix ===')
# '_cleaned_pmids.db' also ends with '_pmids.db'. Matching the short one turns
# 'DAC_Mesh_cleaned_pmids.db' into a project called 'DAC_Mesh_cleaned'.
for bogus in ('DAC_Mesh_cleaned', 'DAC_Mesh_1_cleaned', 'Pilot_cleaned'):
    ck(bogus not in note, f'{bogus!r} is not mistaken for a project')
ck('1.5' in size or '0.0' in size or 'GB' in size, f'a size is shown: {size!r}')

print('\n=== 3. nothing to report means no row ===')
for w in app.setup_rows.winfo_children():
    w.destroy()
solo = tempfile.mkdtemp()
json.dump({'control_flags': {'custom_file_prefix': 'Only'},
           'directories': {'data_dir': os.path.join(solo, 'data'),
                           'results_dir': os.path.join(solo, 'res')}},
          open(os.path.join(solo, 'c.json'), 'w'))
lone = MeshConfig(config_path=os.path.join(solo, 'c.json'))
lone.refresh_paths()
os.makedirs(str(lone.databases_dir), exist_ok=True)
for m in MARKERS:
    with open(os.path.join(str(lone.databases_dir), 'Only' + m), 'wb') as fh:
        fh.write(b'x' * 1000)
app._setup_other_projects(lone, 'Only')
app.update_idletasks()
ck(not [r for r in app.setup_rows.winfo_children() if r.winfo_children()],
   'a single project adds no row at all')

print('\n=== 4. both folder totals are reported ===')
for w in app.setup_rows.winfo_children():
    w.destroy()
app._setup_folder_totals(cfg)
app.update_idletasks()
rows = [r for r in app.setup_rows.winfo_children() if r.winfo_children()]
ck(len(rows) == 2, f'two rows: {len(rows)}')
heads = []
for r in rows:
    labs = [w for w in r.winfo_children() if isinstance(w, tk.Label)]
    h = next(str(w.cget('text')) for w in labs
             if int(w.grid_info().get('column', -1)) == 0
             and int(w.grid_info().get('row', -1)) == 0)
    v = next(str(w.cget('text')) for w in labs
             if int(w.grid_info().get('column', -1)) == 2)
    heads.append(h)
    print(f'    {h:<24}{v}')
    ck('GB' in v, f'{h}: a size is shown', v)
ck(any('Raw' in h for h in heads), 'the raw folder is one of them')
ck(any('Working' in h for h in heads), 'the working folder is the other')

print('\n=== 5. the totals are real, not zero ===')
gb = _size_gb(str(cfg.active_source_dir))
ck(gb is not None and gb > 0, f'the working folder measures above zero: {gb:.4f} GB')

print('\n=== 6. the archive figure is one number everywhere ===')
import re                                                          # noqa: E402
stale = []
for root, _d, files in os.walk('src'):
    if 'notebooks' in root or '__pycache__' in root:
        continue
    for f in files:
        if not f.endswith('.py'):
            continue
        body = open(os.path.join(root, f), encoding='utf-8').read()
        if re.search(r'\b44\s*GB\b', body):
            stale.append(f)
ck(not stale, 'no source file still says 44 GB', f'found in {stale}')

app.destroy()
print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('all disk-accounting checks passed')
