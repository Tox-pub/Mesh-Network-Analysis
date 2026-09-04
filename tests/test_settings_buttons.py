"""The action buttons have to be ON SCREEN on every tab, not merely created.

They were created once for the whole Settings screen, so every check that
looked for them found them - and on Folders, Secondary, Benchmark and Figures
they were nowhere to be seen. The sheet was packed with expand=True BEFORE the
description pane and the button bar, so on any tab whose content is tall it
took the cavity and pushed both past the bottom edge of a fixed-size window.

Existence is therefore not the test. Position is: the bar has to sit inside the
window, below the sheet, with a real height, on every tab.
"""
import os
import sys

sys.path.insert(0, 'src')

import tkinter as tk                                               # noqa: E402
from mesh_workbench import settings_schema as schema                # noqa: E402
from mesh_workbench.app import Workbench                            # noqa: E402

FAILS = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + (f'   -- {extra}' if extra else ''))
    if not ok:
        FAILS.append(msg)


app = Workbench(os.getcwd(), sys.executable)
app.update_idletasks()
app.show('settings')
app.update_idletasks()

# The window that matters is a SHORT one. _geometry_for_font scales the window
# by the font it got, then clamps to the screen - so on a large font or a small
# display the content is laid out for more height than the window ends up with,
# and whatever was packed last falls off the bottom. At the reference font on a
# tall screen nothing overflows and this check passes on broken code, which is
# exactly what it did before this window was shrunk deliberately.
app.minsize(1, 1)
app.geometry('980x430')
app.update()
app.update_idletasks()
H = app.winfo_height()
print(f'  window shrunk to {app.winfo_width()}x{H}, as a small screen would\n')

print('=== 1. the bar is visible on every tab, even when space is short ===')
for name, _fields in schema.TABS:
    app._tab(name)
    app.update()
    app.update_idletasks()
    bar = app.btn_run.master
    y, h = bar.winfo_y(), bar.winfo_height()
    bottom = y + h
    on_screen = bar.winfo_ismapped() and h > 1 and bottom <= H
    ck(on_screen, f'{name:<12} bar at y={y}, height={h}, bottom={bottom} (window {H})',
       'off the bottom' if bottom > H else ('not mapped' if not bar.winfo_ismapped()
                                            else f'height {h}'))

print('\n=== 2. all four buttons, on every tab ===')
WANT = {'Save', 'Defaults', 'Data setup…'}
for name, _fields in schema.TABS:
    app._tab(name)
    app.update_idletasks()
    labels = {str(w.cget('text')) for w in app.btn_run.master.winfo_children()
              if isinstance(w, tk.Button)}
    missing = WANT - labels
    ck(not missing, f'{name:<12} has Save/Defaults/Data setup', f'missing {missing}')
    ck(str(app.btn_run.cget('text')).startswith('Run'),
       f'{name:<12} run button reads {str(app.btn_run.cget("text"))!r}')

print('\n=== 3. the step tabs run their own step ===')
for tab, (step, label) in Workbench.TAB_STEP.items():
    app._tab(tab)
    app.update_idletasks()
    ck(str(app.btn_run.cget('text')) == label,
       f'{tab:<12} button reads {label!r}', f'got {app.btn_run.cget("text")!r}')
    ck(app.step_var.get() == step,
       f'{tab:<12} moved the step list to {step!r}', f'got {app.step_var.get()!r}')

print('\n=== 4. an ordinary tab leaves the step list alone ===')
app._tab('Benchmark')
app.update_idletasks()
app.step_var.set('all')
app._tab('Network')
app.update_idletasks()
ck(app.step_var.get() == 'all', f'Network did not change it: {app.step_var.get()!r}')
ck(str(app.btn_run.cget('text')) == 'Run',
   f'and the button is plain Run: {app.btn_run.cget("text")!r}')

print('\n=== 5. Run sits beside the step list as well ===')
top_buttons = []


def walk(w):
    for c in w.winfo_children():
        if isinstance(c, tk.Button) and str(c.cget('text')) == 'Run':
            top_buttons.append(c)
        walk(c)


walk(app)
ck(len(top_buttons) >= 1, f'a Run button exists near the top: {len(top_buttons)}')
near_top = [b for b in top_buttons if b.winfo_rooty() - app.winfo_rooty() < H // 3]
ck(bool(near_top), f'and one of them is in the top third of the window')

print('\n=== 6. the step label no longer carries the benchmark caveat ===')
labels = dict(schema.STEPS)
ck('not the benchmark' not in labels['all'], f"'all' reads {labels['all']!r}")
ck(labels['all'] == 'all - full pipeline', f"exactly: {labels['all']!r}")
app._step_changed()
ck(str(app.step_desc.cget('text')) == 'full pipeline',
   f'and the description beside it: {app.step_desc.cget("text")!r}')

print('\n=== 7. the step menu shows labels, not argument names ===')
labels = [lab for _k, lab in schema.STEPS]
ck(all('_' not in lab.split(' - ')[0] or lab.startswith('data_ops')
       for lab in labels), 'labels are readable')
ck(not any(lab.strip() == 'viz' for lab in labels), "'viz' is never shown alone")
ck(any(lab.startswith('figures') for lab in labels),
   f'the figures step is called figures: {[l for l in labels if l.startswith("figures")]}')

# the menu widget itself must carry the labels, not the keys
menu = None
# The STEP menu, not the first OptionMenu found: every 'choice' field builds
# one too, and the last one in the tree is a weight-key list.
_want_var = str(app._step_label)


def find_menu(w):
    global menu
    for c in w.winfo_children():
        if isinstance(c, tk.OptionMenu) and str(c.cget('textvariable')) == _want_var:
            menu = c
        find_menu(c)


find_menu(app)
ck(menu is not None, 'the step menu was found')
if menu is not None:
    entries = [menu['menu'].entrycget(i, 'label')
               for i in range(menu['menu'].index('end') + 1)]
    ck('viz' not in entries, f'the menu lists no bare keys: {entries[:3]}...')
    ck(any(e.startswith('figures') for e in entries),
       'and does list the figures step by name')
    ck(entries == labels, 'in the order the schema declares')

# selecting by label must set the KEY the pipeline takes
app._step_picked('figures - plots and network images')
ck(app.step_var.get() == 'viz',
   f'picking the figures label sets the viz key: {app.step_var.get()!r}')
app.set_step('benchmark')
ck(app._step_label.get().startswith('benchmark'),
   f'and setting a key updates the label: {app._step_label.get()!r}')

print('\n=== 7b. no tab hides a control by running off the bottom ===')
# A field clipped by an ancestor still reports winfo_ismapped() == 1, so the
# checks above cannot see this. Figures lost its figure-formats box by 23px and
# Stop words lost three controls, and nothing said so: the fields were laid
# out, reported as mapped, and impossible to see or click. Each tab scrolls
# now, so the invariants are that the scrollbar appears exactly when the
# content overflows, and that the scroll region covers all of it.
for wide, tall in ((980, 430), (1100, 700)):
    app.geometry(f'{wide}x{tall}')
    app.update()
    app.update_idletasks()
    for name, fields in schema.TABS:
        app._tab(name)
        app.update()
        app.update_idletasks()
        anchor = next((app.widgets[f.key] for f in fields
                       if f.key in app.widgets), None)
        if anchor is None:
            continue
        inner, canvas = anchor, anchor.master
        while canvas is not None and not isinstance(canvas, tk.Canvas):
            inner, canvas = canvas, canvas.master
        ck(isinstance(canvas, tk.Canvas), f'{name:<12} sits in a scrollable pane')
        if not isinstance(canvas, tk.Canvas):
            continue
        overflows = inner.winfo_reqheight() > canvas.winfo_height()
        bars = [w for w in canvas.master.winfo_children()
                if isinstance(w, tk.Scrollbar)]
        shown = bool(bars) and bars[0].winfo_ismapped()
        ck(shown == overflows,
           f'{name:<12} {wide}x{tall}: bar shown={shown} for overflow={overflows}',
           f'content {inner.winfo_reqheight()}px in {canvas.winfo_height()}px')
        region = canvas.cget('scrollregion').split()
        ck(len(region) == 4 and int(float(region[3])) >= inner.winfo_reqheight(),
           f'{name:<12} {wide}x{tall}: the scroll region covers the content',
           f'region={canvas.cget("scrollregion")!r} content={inner.winfo_reqheight()}')

app.geometry('980x430')
app.update_idletasks()

print('\n=== 8. step order follows the tabs ===')
keys = [k for k, _lab in schema.STEPS]
ck(keys.index('secondary') < keys.index('benchmark') < keys.index('viz'),
   f'secondary, then benchmark, then figures: {keys}')

app.destroy()
print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('all settings-button checks passed')
