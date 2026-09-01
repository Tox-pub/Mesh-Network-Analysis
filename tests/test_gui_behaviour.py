"""Two new behaviours, driven for real rather than inspected.

The annotation-sync dialog and the post-abort damage check. Both were
written, reviewed and committed with a NameError in them that a broad
except swallowed - the abort check reported "could not check" and did
nothing at all, quietly, on every abort. Inspecting the source would not
have caught it; running it did.
"""
import json
import os
import sqlite3
import sys
import tempfile

sys.path.insert(0, 'src')

import tkinter as tk                                               # noqa: E402
from mesh_workbench import app as appmod                            # noqa: E402
from mesh_workbench.app import Workbench                            # noqa: E402

FAILS = []
ASKED = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + (f'   -- {extra}' if extra else ''))
    if not ok:
        FAILS.append(msg)


box = tempfile.mkdtemp()
cfgp = os.path.join(box, 'mesh_config.json')
json.dump({"control_flags": {"custom_file_prefix": "AUDIT",
                             "pause_for_annotation": True},
           "directories": {"data_dir": os.path.join(box, 'data'),
                           "results_dir": os.path.join(box, 'results')}},
          open(cfgp, 'w'))

w = Workbench(os.getcwd(), sys.executable)
w.withdraw()
w.cfg_path = cfgp
w.update_idletasks()

# ------------------------------------------------------------------ C
print('=== C. the annotation-sync dialog ===')
from mesh_aop.config_parser import MeshConfig                       # noqa: E402
cfg = MeshConfig(config_path=cfgp)
cfg.refresh_paths()
anno = str(cfg.results_dir / f'{cfg.prefix}_run_annotations.csv')
os.makedirs(os.path.dirname(anno), exist_ok=True)
with open(anno, 'w', encoding='utf-8', newline='') as fh:
    fh.write('mesh_term;aop_level\nSkin;Key Event\nHaptens;Unassigned\n')

w.vars['control_flags.pause_for_annotation'].set(True)
w._paused_annotation = anno

for answer, want in ((True, 'yes'), (False, 'no'), (None, None)):
    ASKED.clear()

    def fake(*a, **k):
        ASKED.append(a[0] if a else '')
        return answer
    appmod.messagebox.askyesnocancel = fake
    got = w._ask_annotation_sync('viz')
    ck(got == want, f'dialog answer {answer!r} -> {got!r}', f'wanted {want!r}')
    ck(len(ASKED) == 1, 'the dialog was shown exactly once')

# the unassigned count reaches the dialog text
shown = []
appmod.messagebox.askyesnocancel = lambda *a, **k: (shown.append(a[1]), False)[1]
w._ask_annotation_sync('viz')
ck('1 term(s) are still Unassigned' in shown[0],
   'the dialog names the real Unassigned count', shown[0][-120:])

# AFK: no dialog at all
w.vars['control_flags.pause_for_annotation'].set(False)
ASKED.clear()
appmod.messagebox.askyesnocancel = lambda *a, **k: ASKED.append(1) or True
got = w._ask_annotation_sync('viz')
ck(got == 'no' and not ASKED, f'AFK answers no with no dialog: {got!r}, shown={len(ASKED)}')

# no annotations file: also no dialog
w.vars['control_flags.pause_for_annotation'].set(True)
w._paused_annotation = None
os.remove(anno)
ASKED.clear()
got = w._ask_annotation_sync('viz')
ck(got == 'no' and not ASKED,
   f'no annotations file -> no dialog: {got!r}, shown={len(ASKED)}')

# ------------------------------------------------------------------ E
print('\n=== E. the post-abort check, against real damage ===')
cfg2 = MeshConfig(config_path=cfgp)
cfg2.refresh_paths()

# a truncated network JSON - what a kill mid-write actually leaves
bad_json = str(cfg2.files['final_network'])
os.makedirs(os.path.dirname(bad_json), exist_ok=True)
with open(bad_json, 'w', encoding='utf-8') as fh:
    fh.write('{"elements": {"nodes": [{"data": {"id": "Sk')
# a zero-byte database
bad_db = str(cfg2.files['cleaned_db'])
os.makedirs(os.path.dirname(bad_db), exist_ok=True)
open(bad_db, 'w').close()
# and one that is perfectly fine, which must survive
good_db = str(cfg2.files['relevance_db'])
c = sqlite3.connect(good_db)
c.execute('CREATE TABLE article_relevance_scores (pmid INTEGER PRIMARY KEY)')
c.execute('INSERT INTO article_relevance_scores VALUES (1)')
c.commit()
c.close()

from mesh_aop import integrity                                      # noqa: E402
found = integrity.scan(cfg2, deep=False)
damaged = [a for a in found
           if a.status in (integrity.EMPTY, integrity.CORRUPT, integrity.ORPHAN)]
names = sorted(os.path.basename(a.path) for a in damaged)
ck(os.path.basename(bad_json) in names, f'the truncated JSON is found: {names}')
ck(os.path.basename(bad_db) in names, f'the zero-byte database is found: {names}')
ck(os.path.basename(good_db) not in names, 'the healthy database is NOT flagged')

# drive the real handler, answering yes to removal
logged = []
w._log = lambda msg, tag=None: logged.append(msg)
w.run_status.config = lambda **k: None
appmod.messagebox.askyesno = lambda *a, **k: True
w._check_after_abort()
ck(not os.path.exists(bad_json), 'the truncated JSON was removed')
ck(not os.path.exists(bad_db), 'the zero-byte database was removed')
ck(os.path.exists(good_db), 'the healthy database was left alone')
ck(any('half-written' in m for m in logged),
   f'and it said so in the log', str(logged[:4]))

# declining leaves everything in place
with open(bad_json, 'w', encoding='utf-8') as fh:
    fh.write('{"elements": {"nod')
logged.clear()
appmod.messagebox.askyesno = lambda *a, **k: False
w._check_after_abort()
ck(os.path.exists(bad_json), 'declining leaves the file in place')
ck(any('left in place' in m for m in logged), 'and says where to find it later')

# nothing damaged -> no dialog at all
os.remove(bad_json)
shown2 = []
appmod.messagebox.askyesno = lambda *a, **k: shown2.append(1) or True
logged.clear()
w._check_after_abort()
ck(not shown2, f'a clean abort shows no dialog: {len(shown2)}')
ck(any('nothing was left half-written' in m for m in logged),
   f'and says so', str(logged))

w.destroy()
print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('every behavioural checklist item passed')
