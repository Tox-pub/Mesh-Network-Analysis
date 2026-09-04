"""Any annotation scheme is valid, and the MeSH trees are a choice.

The pipeline used to hard-code the seven levels of an adverse outcome pathway
and the four MeSH trees an AOP is built from. Both are now settings, and both
have a failure mode that is silent rather than loud:

  * a stratum the user typed that the settings did not list could be dropped
    from the figures, which would make the figure disagree with the file it
    was drawn from;
  * a choice of trees that never reaches the network would leave the run
    looking correct while analysing the old vocabulary.

So neither is tested by asking whether the code runs. Both are tested by
changing the setting and checking the OUTPUT changed with it.
"""
import csv
import os
import sys
import tempfile

sys.path.insert(0, 'src')

from mesh_aop import strata, vocabulary, guides                     # noqa: E402

FAILS = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + ('' if ok or not extra
                                                    else f'   -- {extra}'))
    if not ok:
        FAILS.append(msg)


print('=== 1. the strata are whatever the file says ===')
order, appended = strata.resolve_order(
    'Stressor; Molecular; Organ',
    ['Molecular', 'Rogue Group', 'Uncategorized', 'Stressor'])
ck(order == ['Stressor', 'Molecular', 'Rogue Group', 'Uncategorized'],
   f'configured order kept, unlisted appended: {order}')
ck(appended == ['Rogue Group'], f'and reported: {appended}')
ck('Organ' not in order, 'a configured stratum nobody used is not plotted')

order, _ = strata.resolve_order('A; B', ['Uncategorized', 'B', 'A'])
ck(order[-1] == strata.UNASSIGNED, f'unplaced always sits last: {order}')

order, appended = strata.resolve_order('', ['Only One'])
ck(order == ['Only One'], f'no setting at all still plots the file: {order}')

ck(strata.parse_order('A;;B; A ;C') == ['A', 'B', 'C'],
   'blanks and repeats are dropped')
ck(strata.parse_order('A, B, C') == ['A', 'B', 'C'],
   'a comma list is accepted too')

print('\n=== 2. unplaced means unplaced, however it is spelled ===')
for value in ('', None, 'Unassigned', 'unassigned', 'Uncategorized', 'nan'):
    ck(strata.is_unplaced(value), f'{value!r} counts as unplaced')
ck(not strata.is_unplaced('Molecular'), 'a real stratum does not')
ck(strata.normalise(' Molecular ') == 'Molecular', 'space is trimmed')
ck(strata.normalise('molecular') == 'molecular',
   'but the spelling is left alone - it is the user\'s name, not ours')
ck(strata.normalise('') == strata.UNASSIGNED, 'blank becomes the unplaced name')

print('\n=== 3. a file annotated before the rename still opens ===')
ck(strata.column_of(['mesh_term', 'aop_level']) == 'aop_level',
   'the old column is found')
ck(strata.column_of(['mesh_term', 'stratum']) == 'stratum',
   'so is the new one')
ck(strata.column_of(['stratum', 'aop_level']) == 'stratum',
   'and the new one wins if a file somehow has both')
ck(strata.column_of(['mesh_term']) is None, 'a file with neither says so')

print('\n=== 4. changing the trees changes what is excluded ===')
tmp = tempfile.mkdtemp()
terms_csv = os.path.join(tmp, 'mesh_terms.csv')
ROWS = [('Skin', 'A'), ('Dermatitis', 'C'), ('Nickel', 'D'),
        ('Inflammation', 'CG'), ('Escherichia coli', 'B'), ('Sweden', 'Z'),
        ('Review', 'V'), ('Nurses', 'MN'), ('Orphan Term', '')]
with open(terms_csv, 'w', encoding='utf-8', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['DescriptorName', 'DescriptorUI', 'MeSH_stop_term',
                vocabulary.TREE_COLUMN])
    for name, trees in ROWS:
        w.writerow([name, 'D000000', 'False', trees])

default, prov = vocabulary.resolve(mesh_csv_path=terms_csv)
ck('Sweden' in default, 'by default a Geographical is excluded')
ck('Skin' not in default, 'and an Anatomy term is not')
ck('Orphan Term' in default, 'a term in no tree at all is excluded')
ck(terms_csv.endswith('.csv') and 'mesh_terms.csv' in prov['source'],
   f'the source is named in the provenance: {prov["source"]}')

kept_z = vocabulary.resolve(excluded_trees='B;V', mesh_csv_path=terms_csv)[0]
ck('Sweden' not in kept_z, 'un-ticking Geographicals brings Sweden back')
ck('Nurses' not in kept_z, 'and Named Groups brings Nurses back')
ck('Escherichia coli' in kept_z, 'while Organisms stays excluded')
ck(kept_z != default, 'the two selections genuinely differ')

print('\n=== 5. the sexes, which no tree can reach ===')
ck('Male' in default and 'Female' in default, 'excluded by default')
kept = vocabulary.resolve(mesh_csv_path=terms_csv, keep_sexes=True)[0]
ck('Male' not in kept and 'Female' not in kept, 'and kept when asked')
# The point of the switch: no choice of trees can do this on its own.
every_tree_kept = vocabulary.resolve(excluded_trees='', mesh_csv_path=terms_csv)[0]
ck('Male' in every_tree_kept,
   'keeping every tree does NOT keep them - only the switch does')

print('\n=== 6. hand-entered terms ===')
extra = vocabulary.resolve(mesh_csv_path=terms_csv,
                           extra_terms='Skin; Dermatitis, Allergic Contact')[0]
ck('Skin' in extra, 'a term listed by hand is excluded')
ck('Dermatitis, Allergic Contact' in extra,
   'a heading containing a comma survives as one term')
ck(vocabulary.parse_terms('Dermatitis, Allergic Contact')
   == ('Dermatitis, Allergic Contact',),
   'it is not split on the comma')

print('\n=== 7. with no tree column it falls back, and says so ===')
old_csv = os.path.join(tmp, 'old_terms.csv')
with open(old_csv, 'w', encoding='utf-8', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['DescriptorName', 'DescriptorUI', 'MeSH_stop_term'])
    w.writerow(['Skin', 'D000000', 'False'])
words, prov = vocabulary.resolve(excluded_trees='B;Z', mesh_csv_path=old_csv)
ck(len(words) > 1000, f'the shipped list is used, not an empty set: {len(words)}')
ck('shipped' in prov['source'], f'the source says so: {prov["source"]}')
ck('warning' in prov, 'and the mismatch is reported, not swallowed')
# The dangerous case is the quiet one: no stop words at all would look fine.
ck(len(vocabulary.resolve(mesh_csv_path=None)[0]) > 1000,
   'a missing file also falls back rather than excluding nothing')

print('\n=== 8. the trees the settings window lists ===')
ck(len(vocabulary.TREES) == 16, f'all sixteen: {len(vocabulary.TREES)}')
ck(vocabulary.DEFAULT_KEPT == ('A', 'C', 'D', 'G'),
   f'four kept by default: {vocabulary.DEFAULT_KEPT}')
ck(len(vocabulary.DEFAULT_EXCLUDED) == 12,
   f'twelve excluded: {len(vocabulary.DEFAULT_EXCLUDED)}')
ck(set(vocabulary.DEFAULT_KEPT) & set(vocabulary.DEFAULT_EXCLUDED) == set(),
   'and no tree is in both')
ck(vocabulary.parse_letters('BEF') == ('B', 'E', 'F'), 'letters with no separator')
ck(vocabulary.parse_letters('Q;9') == (), 'nonsense is dropped')
ck(vocabulary.format_letters('ZBE') == 'B;E;Z', 'stored in tree order')

print('\n=== 9. configure() reaches the modules that filter ===')
from mesh_aop import network, gt_network_validation                 # noqa: E402
vocabulary.configure({'Only This'})
ck(network.vocabulary.active() == frozenset({'Only This'}),
   'network sees the configured set')
ck(gt_network_validation.vocabulary.active() == frozenset({'Only This'}),
   'so does the benchmark')
vocabulary.configure(None)
ck(len(vocabulary.active()) > 1000, 'and None restores the shipped list')

print('\n=== 10. the guide describes the scheme in force ===')
md = guides.annotation_guide_markdown('/tmp/run_annotations.csv',
                                      strata_order='Alpha; Beta; Gamma')
ck('Alpha' in md and 'Beta' in md, 'it lists the configured strata')
ck('Stressor' not in md,
   'and NOT the AOP default when a different scheme is set')
ck(strata.COLUMN in md, f'it names the {strata.COLUMN!r} column')
ck('master annotations library' in md.lower(),
   'and warns about the shared library')
ck('permanent' in md.lower() or 'shared by every project' in md.lower(),
   'saying the merge is permanent')

default_md = guides.annotation_guide_markdown('/tmp/x.csv')
ck('Stressor' in default_md, 'the AOP seven are still the default')
ck('example, not a requirement' in default_md,
   'presented as an example rather than a rule')

print('\n=== 11. it is written as a page, not a text file ===')
out = os.path.join(tmp, 'run_annotations.csv')
open(out, 'w').close()
written = guides.write_annotation_guide(out, strata_order='Alpha; Beta')
ck(written and written.endswith('.html'), f'written as HTML: {written}')
body = open(written, encoding='utf-8').read()
ck(body.startswith('<!doctype html>'), 'a complete document')
ck('<style>' in body and '<script' not in body,
   'self-contained and inert')
ck('Alpha' in body, 'carrying this project\'s strata')
ck(guides.ANNOTATION_GUIDE_NAME in guides.ANNOTATION_GUIDE_NAMES
   and guides.ANNOTATION_GUIDE_NAME_FALLBACK in guides.ANNOTATION_GUIDE_NAMES,
   'both names are searched, so an older guide is not duplicated')

print('\n=== 12. the Stop words tab, and its sixteen checkboxes ===')
import tkinter as tk                                                # noqa: E402
from mesh_workbench import settings_schema as sschema               # noqa: E402
from mesh_workbench.app import Workbench                            # noqa: E402

names = [n for n, _f in sschema.TABS]
ck('Stop words' in names, f'the tab exists: {names}')
ck(names.index('Stop words') == names.index('Folders') + 1,
   f'directly after Folders: {names}')

app = Workbench(os.getcwd(), sys.executable)
app.update_idletasks()
app.show('settings')
app._tab('Stop words')
app.update()
app.update_idletasks()

var = app.vars.get('stop_words.excluded_trees')
ck(var is not None, 'the trees field has a variable like any other setting')

holder = app.widgets.get('stop_words.excluded_trees')
boxes = [w for w in (holder.winfo_children() if holder else [])
         if isinstance(w, tk.Checkbutton)]
ck(len(boxes) == 16, f'sixteen checkboxes are drawn: {len(boxes)}')
ck(all(b.winfo_ismapped() for b in boxes),
   'and every one of them is actually on screen')

labels = [str(b.cget('text')) for b in boxes]
ck(labels[0].startswith('A - '), f'labelled by letter and name: {labels[0]!r}')
ck(any('Geographicals' in t for t in labels), 'including Geographicals')

# The one that matters: ticking a box has to reach the stored string, and
# setting the string has to move the boxes. A widget that only did one of the
# two would look right on screen and save the wrong thing - or the reverse.
start = set(vocabulary.parse_letters(var.get()))
ck(start == set(vocabulary.DEFAULT_EXCLUDED),
   f'it opens on the twelve excluded by default: {sorted(start)}')
ck('A' not in start and 'Z' in start,
   'Anatomy unticked, Geographicals ticked - the way round the label reads')

anatomy = next(b for b in boxes if str(b.cget('text')).startswith('A - '))
anatomy.invoke()
app.update_idletasks()
ck('A' in vocabulary.parse_letters(var.get()),
   f'ticking Anatomy reaches the stored value: {var.get()!r}')
anatomy.invoke()
app.update_idletasks()
ck('A' not in vocabulary.parse_letters(var.get()),
   f'and un-ticking it takes it back out: {var.get()!r}')

var.set('B;Z')
app.update_idletasks()


def ticked():
    """Which boxes are on, read the way Tk itself reads them.

    NOT bool() of the raw Tcl value: after invoke() that value is the STRING
    '0', and bool('0') is True - which reports every box as ticked and makes
    this check pass on a widget that does not work.
    """
    out = {}
    for b, t in zip(boxes, labels):
        out[t.split(' - ')[0]] = bool(
            b.tk.getboolean(b.getvar(b.cget('variable'))))
    return out


state = ticked()
on = sorted(k for k, v in state.items() if v)
ck(on == ['B', 'Z'], f'loading a saved value ticks exactly those boxes: {on}')
ck(not state['A'], 'a box ticked earlier is cleared, not left behind')
ck(not state['E'], 'and the rest are unticked')

# What is actually written to mesh_config.json, which is the only thing the
# pipeline ever sees. A widget that looked right and saved something else
# would be the worst outcome here.
ck(var.get() == 'B;Z', f'the saved value is the selection: {var.get()!r}')
words_from_ui = vocabulary.resolve(excluded_trees=var.get(),
                                   mesh_csv_path=terms_csv)[0]
ck('Sweden' in words_from_ui,
   'and feeding it back through the resolver excludes Geographicals')
ck('Nurses' not in words_from_ui,
   'while Named Groups, now unticked, is kept')

# The tree block spans sixteen rows of a grid whose other fields are one row
# each. Get the rowspan wrong and it lands on top of the three settings below
# it - which looks like a rendering glitch and is actually two controls sharing
# one cell, where only the last one drawn can be clicked.
print('\n=== 12b. and it does not sit on top of the settings below it ===')
stack = []
for key in ('stop_words.excluded_trees', 'stop_words.keep_sexes',
            'stop_words.extra_terms', 'stop_words.rebuild'):
    wdg = app.widgets.get(key)
    ck(wdg is not None and wdg.winfo_ismapped(), f'{key} is on screen')
    if wdg is not None:
        wdg.update_idletasks()
        top = wdg.winfo_rooty() - app.winfo_rooty()
        stack.append((key, top, top + wdg.winfo_height()))

clashes = [(a[0], b[0]) for i, a in enumerate(stack) for b in stack[i + 1:]
           if a[1] < b[2] and b[1] < a[2]]
ck(not clashes, f'no two of them overlap', f'{clashes}')
ck(all(stack[i][2] <= stack[i + 1][1] for i in range(len(stack) - 1)),
   'and they stack in the order the schema lists them',
   f'{[(k, t, b) for k, t, b in stack]}')

print('\n=== 13. nothing user-facing still says "AOP level" ===')
blob = []
for _name, fields in sschema.TABS:
    for f in fields:
        blob += [str(f.label), str(f.what), str(f.deflt), str(f.note or '')]
text = ' '.join(blob)
ck('AOP level' not in text, 'the settings form does not')
ck('aop_level' not in text, 'nor does it name the old column')
ck('Strata order' in text, 'and the strata order is offered')

menu = app.nametowidget(app.cget('menu'))
found = []
for i in range(menu.index('end') + 1):
    if menu.type(i) == 'cascade':
        sub = app.nametowidget(menu.entrycget(i, 'menu'))
        for j in range(sub.index('end') + 1):
            if sub.type(j) != 'separator':
                found.append(str(sub.entrycget(j, 'label')))
ck(not any('AOP' in lbl for lbl in found), f'no menu entry does: {found}')
ck(any('strata' in lbl.lower() for lbl in found),
   f'and the strata entries are there: {[l for l in found if "strata" in l.lower()]}')

app.destroy()

import shutil                                                       # noqa: E402
shutil.rmtree(tmp, ignore_errors=True)

print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('all strata and vocabulary checks passed')
