"""The parts of the window a user actually looks at.

Logic has unit tests. These do not: whether a label is really bold, whether
the description pane is really short, whether two buttons really read the
same thing. Each of these was a defect someone had to notice by eye. So the
window is built and the widget tree is walked.
"""
import os
import sys
import tempfile

sys.path.insert(0, 'src')
os.environ.setdefault('MESH_WORKBENCH_TEST', '1')

import tkinter as tk                                               # noqa: E402

FAILS = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + (f'   -- {extra}' if extra else ''))
    if not ok:
        FAILS.append(msg)


def walk(w, out=None):
    out = [] if out is None else out
    out.append(w)
    for c in w.winfo_children():
        walk(c, out)
    return out


from mesh_workbench.app import Workbench                            # noqa: E402

box = tempfile.mkdtemp()
cfgp = os.path.join(box, 'mesh_config.json')

app = Workbench(os.getcwd(), sys.executable)
app.withdraw()
app.update_idletasks()

widgets = walk(app)
print(f'  built the window: {len(widgets)} widgets\n')

# ---------------------------------------------------------------- A1
print('=== A1. Pipeline step is bold ===')
labels = [w for w in widgets if isinstance(w, tk.Label)
          and 'Pipeline step' in str(w.cget('text'))]
ck(bool(labels), 'the "Pipeline step:" label exists')
if labels:
    f = str(labels[0].cget('font'))
    ck(f == str(app.f_bold), f'its font is the bold one', f'font={f}, bold={app.f_bold}')
oms = [w for w in widgets if isinstance(w, tk.OptionMenu)]
if oms:
    f = str(oms[0].cget('font'))
    ck(f == str(app.f_bold), 'the step dropdown is bold too', f'font={f}')

# ---------------------------------------------------------------- A2/A3
print('\n=== A2/A3. the description pane is not half the window ===')
h = app._help_height()
H = app.winfo_height() or 720
try:
    import tkinter.font as tkfont
    cap = int(H * 0.46)
except Exception:
    cap = 331
ck(h < cap, f'pane height {h}px is below the {cap}px ceiling', f'h={h} cap={cap}')
ck(h < 220, f'and is modest in absolute terms ({h}px)', f'h={h}')

# which field drives it now
from mesh_workbench import settings_schema as schema                # noqa: E402
longest = max(((len(f.note or ''), f.key) for _t, fs in schema.TABS for f in fs))
ck('use_reference_data' not in longest[1],
   f'the reference checkbox no longer drives the height (longest is {longest[1]})',
   f'longest={longest}')

# ---------------------------------------------------------------- A7
print('\n=== A7. the tab note boxes scroll ===')
import inspect                                                      # noqa: E402
src = inspect.getsource(Workbench)
tab_note = src[src.index('note = schema.TAB_NOTES.get'):]
tab_note = tab_note[:tab_note.index('def _browse')]
ck('tk.Scrollbar' in tab_note, 'a scrollbar is created for the note box')
ck('yscrollcommand' in tab_note, 'and wired to the text widget')

# ---------------------------------------------------------------- E1/E2
print('\n=== E1/E2. one View results button, not two ===')
btns = [w for w in widgets if isinstance(w, tk.Button)]
texts = [str(b.cget('text')) for b in btns]
vr = [t for t in texts if t == 'View results']
ck(len(vr) <= 1, f'at most one "View results" button exists: {len(vr)}', f'texts={vr}')

# now simulate a finished run landing on results, which is when it duplicated
app._run_step = 'viz'
app.current_screen = 'running'
app._offer_next_screen()
app.update_idletasks()
back_text = str(app.btn_back.cget('text'))
res_text = str(app.btn_results.cget('text'))
ck(back_text != res_text,
   f'after a figures run the two buttons differ: back={back_text!r} results={res_text!r}')
ck(back_text == 'Back to settings',
   f'the back button keeps its own job', f'text={back_text!r}')

# and the database case still relabels
app._run_step = 'baseline'
app._offer_next_screen()
app.update_idletasks()
ck(str(app.btn_back.cget('text')) == 'Back to the database',
   'a database step still relabels it', f'text={app.btn_back.cget("text")!r}')

# ---------------------------------------------------------------- E3
print('\n=== E3/E4. abort checks for half-written files ===')
ck(hasattr(app, '_check_after_abort'), 'the post-abort check exists')
fin = inspect.getsource(Workbench._run_finished) if hasattr(Workbench, '_run_finished') else src
seg = fin[fin.index('elif rc == -1'):fin.index('elif rc == 130')]
ck('_check_after_abort()' in seg, 'the hard abort calls it')
seg130 = fin[fin.index('elif rc == 130'):]
seg130 = seg130[:seg130.index('else:')]
ck('_check_after_abort' not in seg130, 'the clean checkpoint stop does not')

# ---------------------------------------------------------------- A5/A6
print('\n=== A5/A6. the two tab boxes say the right things ===')
sec = schema.TAB_NOTES.get('Secondary', ('', ''))[1]
ben = schema.TAB_NOTES.get('Benchmark', ('', ''))[1]
ck('Cytoscape' in sec, 'Secondary names Cytoscape')
ck('Excel' in sec, 'and the Excel export')
ck('not' in sec.lower() and 'full' in sec.lower(), 'and says it is not part of a full run')
ck('Cytoscape' in ben, 'Benchmark names Cytoscape too')
from mesh_aop.benchmark import KNOWN_GT_FILENAMES                   # noqa: E402
missing = [n for n in KNOWN_GT_FILENAMES if n not in ben]
ck(not missing, 'and lists every auto-detected ground-truth filename', f'missing {missing}')
ck('RAW DATA' in ben or 'raw data' in ben.lower(), 'and says which folder to put it in')

app.destroy()
print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('every visual checklist item passed')
