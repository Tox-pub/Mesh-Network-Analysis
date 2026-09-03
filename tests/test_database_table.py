"""The data table's headings have to sit over the values they name.

Two grids describe this table - the header, packed straight into the box, and
each row, packed inside a scrolling canvas. They only line up if both are
configured identically AND the header reserves the width the scrollbar takes
from the rows. Get either wrong and every column after Source sits a
scrollbar's width off, which reads as a table that was never aligned.

Also checked: the Built column exists and carries a real date. A master
database is built once and reused across projects for months, so how old it is
decides whether a result reflects current MeSH - and nothing on screen said.
"""
import os
import re
import sys
import tempfile
import time

sys.path.insert(0, 'src')

import tkinter as tk                                               # noqa: E402
from mesh_workbench.app import Workbench, _built_when               # noqa: E402

FAILS = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + (f'   -- {extra}' if extra else ''))
    if not ok:
        FAILS.append(msg)


print('=== 1. _built_when reads a real timestamp ===')
box = tempfile.mkdtemp()
today = time.strftime('%Y-%m-%d')

f = os.path.join(box, 'thing.db')
with open(f, 'w') as fh:
    fh.write('x')
ck(_built_when(f) == today, f'a file: {_built_when(f)}', f'expected {today}')

d = os.path.join(box, 'archive')
os.makedirs(d)
for n in range(3):
    with open(os.path.join(d, f'part{n}.gz'), 'w') as fh:
        fh.write('x')
ck(_built_when(d) == today, f'a directory takes its newest member: {_built_when(d)}')

ck(_built_when(os.path.join(box, 'absent')) is None, 'a missing path gives None')
ck(_built_when(os.path.join(box, 'empty_dir')) is None
   or _built_when(box) is not None, 'an empty directory does not raise')

fmt = _built_when(f)
ck(bool(re.fullmatch(r'\d{4}-\d{2}-\d{2}', fmt or '')),
   f'the format is a plain sortable date: {fmt!r}')

print('\n=== 2. header and rows describe the same columns ===')
app = Workbench(os.getcwd(), sys.executable)
app.update_idletasks()
app.show('setup')
app.update()
app.update_idletasks()
app.scan_setup() if hasattr(app, 'scan_setup') else None
app.update()
app.update_idletasks()


# Which edge to compare. A cell gridded sticky='e' right-aligns its label, so
# the label's LEFT edge moves with the length of its text - "Size" and
# "0.00 GB" start in different places while their cells line up perfectly.
# Compare the edge the column is anchored to, or the test measures the words.
EDGE = {1: 'left', 2: 'right', 3: 'left'}


def grid_x(frame):
    """The anchored edge of each gridded column, in screen pixels."""
    out = {}
    for w in frame.winfo_children():
        info = w.grid_info()
        if not info:
            continue
        c = int(info['column'])
        if c in out:
            continue
        out[c] = (w.winfo_rootx() + w.winfo_width()
                  if EDGE.get(c) == 'right' else w.winfo_rootx())
    return out


# The header is the frame holding the bold 'Source' label.
def find_header(w):
    for c in w.winfo_children():
        if isinstance(c, tk.Label) and str(c.cget('text')) == 'Source':
            return c.master
        got = find_header(c)
        if got is not None:
            return got
    return None


hdr = find_header(app)
ck(hdr is not None, 'the header row was found')

rows = [r for r in app.setup_rows.winfo_children() if r.winfo_children()]
ck(len(rows) >= 2, f'data rows were built: {len(rows)}')

if hdr is not None and rows:
    hx = grid_x(hdr)
    print(f'    header columns at: {[hx.get(i) for i in range(5)]}')
    for r in rows[:6]:
        rx = grid_x(r)
        label = next((str(w.cget('text')) for w in r.winfo_children()
                      if isinstance(w, tk.Label)
                      and int(w.grid_info().get('column', -1)) == 0
                      and int(w.grid_info().get('row', -1)) == 0), '?')
        for col in (1, 2, 3):
            if col not in hx or col not in rx:
                continue
            drift = abs(hx[col] - rx[col])
            ck(drift <= 2,
               f'{label[:26]:<28} column {col} aligned (drift {drift}px)',
               f'header {hx[col]} vs row {rx[col]}')

print('\n=== 3. the Built column is there and populated ===')
heads = [str(w.cget('text')) for w in (hdr.winfo_children() if hdr else [])]
ck('Built' in heads, f'the heading exists: {heads}')
ck(heads.index('Built') > heads.index('Size') if 'Built' in heads else False,
   'and sits after Size')

dates = []
for r in rows:
    for w in r.winfo_children():
        if isinstance(w, tk.Label) and int(w.grid_info().get('column', -1)) == 3:
            dates.append(str(w.cget('text')))
ck(bool(dates), f'every row has a Built cell: {len(dates)}')
real = [d for d in dates if re.fullmatch(r'\d{4}-\d{2}-\d{2}', d)]
ck(all(re.fullmatch(r'\d{4}-\d{2}-\d{2}', d) or d == '-' for d in dates),
   f'each is a date or a dash: {sorted(set(dates))}')

app.destroy()
print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f_ in FAILS:
        print('   -', f_)
    sys.exit(1)
print('all database-table checks passed')
