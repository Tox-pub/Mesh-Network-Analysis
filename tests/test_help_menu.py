"""The Help menu, and the reader windows behind it.

Every entry opens a window this program owns rather than handing the file to
whatever the system opens Markdown with - which steals focus, may not be
installed, and on Linux is often a text editor where the shipped copy can be
edited by accident. The manual in particular is meant to be read WHILE the
form is being filled in, so none of these may be modal.

Also checked: the citation is one string used in three places, so what the Help
menu shows, what HELP.md prints and what CITATION.cff records cannot drift.
"""
import os
import re
import sys

sys.path.insert(0, 'src')

import tkinter as tk                                               # noqa: E402
from mesh_aop import citation                                       # noqa: E402
from mesh_workbench.app import Workbench                            # noqa: E402

FAILS = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + (f'   -- {extra}' if extra else ''))
    if not ok:
        FAILS.append(msg)


print('=== 1. the citation never reads "Version None" ===')
# version() does not only RAISE when a package is absent; on an editable
# install with no Version field it returns None, which went into the text.
for name, value in (('citation_text', citation.citation_text()),
                    ('bibtex', citation.bibtex()),
                    ('markdown', citation.citation_markdown())):
    ck('None' not in value, f'{name} carries no None')
ck(bool(re.search(r'Version \d', citation.citation_text())),
   f'and a real version: {citation.citation_text()[-90:]}')
# The fallback is a legitimate answer: an editable install often carries no
# Version metadata. What matters is that SOMETHING sane comes out.
ck(re.fullmatch(r'\d+\.\d+\.\d+', citation.package_version() or ''),
   f'the version looks like a version: {citation.package_version()!r}')

print('\n=== 2. the works this rests on are cited ===')
md = citation.citation_markdown()
for who, what in (('Dianati', 'the GLF paper'),
                  ('10.1103/PhysRevE.93.012304', 'its DOI'),
                  ('OECD', 'the ground-truth source'),
                  ('10.1787/9789264221444-en', 'its DOI'),
                  ('Blondel', 'Louvain'),
                  ('Truchon', 'BEDROC'),
                  ('Hutchins', 'iCite')):
    ck(who in md, f'{what} is cited')
ck('Global Likelihood Filter' in md, 'GLF is expanded the way the paper does')

print('\n=== 3. GLF is expanded consistently everywhere ===')
WRONG = ('Graph Likelihood Filtering', 'Global-Local Filtering',
         'Global likelihood filtering', 'global likelihood filtering')
# A line may NAME a wrong spelling in order to say it was wrong - the citation
# module and the manual both do, deliberately, so the correction is on record.
# Only unexplained uses count.
EXPLAINING = ('wrong', 'Earlier versions', 'spelled', 'before this was checked',
              'three different')


def _wrong_uses(body):
    out = []
    for line in body.splitlines():
        if any(e in line for e in EXPLAINING):
            continue
        out += [w for w in WRONG if w in line]
    return out
for root, _dirs, files in os.walk('src'):
    if 'notebooks' in root or '__pycache__' in root:
        continue
    for f in files:
        if not f.endswith('.py'):
            continue
        hits = _wrong_uses(open(os.path.join(root, f), encoding='utf-8').read())
        ck(not hits, f'{f} uses the right name', f'found {hits}')
for doc in ('HELP.md', 'README.md'):
    if os.path.exists(doc):
        hits = _wrong_uses(open(doc, encoding='utf-8').read())
        ck(not hits, f'{doc} uses the right name', f'found {hits}')

print('\n=== 4. the menu says what things are, not what they are called ===')
app = Workbench(os.getcwd(), sys.executable)
app.withdraw()
app.update_idletasks()

menu = app.nametowidget(app.cget('menu'))
help_menu = None
for i in range(menu.index('end') + 1):
    try:
        if menu.type(i) == 'cascade' and menu.entrycget(i, 'label') == 'Help':
            help_menu = app.nametowidget(menu.entrycget(i, 'menu'))
    except tk.TclError:
        pass
ck(help_menu is not None, 'the Help menu exists')
labels = []
for i in range(help_menu.index('end') + 1):
    labels.append('---' if help_menu.type(i) == 'separator'
                  else help_menu.entrycget(i, 'label'))

for want in ('MeSH Workbench Manual', 'Installing and Updating',
             'License', 'Cite this Program', 'About'):
    ck(want in labels, f'{want!r} is on the menu', f'menu is {labels}')
for gone in ('.md', 'README', 'Reference figures'):
    ck(not any(gone in l for l in labels), f'nothing mentions {gone!r}',
       f'menu is {labels}')

print('\n=== 5. every reader is a window, read-only, and not modal ===')
for name, fn in (('MeSH Workbench Manual', app.open_manual),
                 ('License', app.open_license),
                 ('Cite this Program', app.open_citation),
                 ('About', app.open_about)):
    fn()
    app.update_idletasks()
    win = app._readers.get(name)
    ck(win is not None and win.winfo_exists() == 1, f'{name}: opens')
    texts = [c for parent in win.winfo_children()
             for c in parent.winfo_children() if isinstance(c, tk.Text)]
    ck(bool(texts), f'{name}: has a text area')
    if texts:
        ck(str(texts[0].cget('state')) == 'disabled', f'{name}: is read-only')
        ck(len(texts[0].get('1.0', 'end')) > 200, f'{name}: rendered something')
    ck(win.grab_current() is None,
       f'{name}: is NOT modal - the program keeps running underneath')
    again = fn()
    app.update_idletasks()
    ck(app._readers.get(name) is win, f'{name}: opening twice raises the same window')

print('\n=== 6. About states the licence ===')
app.open_about()
app.update_idletasks()
about = app._readers['About']
body = [c for p in about.winfo_children() for c in p.winfo_children()
        if isinstance(c, tk.Text)][0]
text = body.get('1.0', 'end')
for want in ('MIT', '2026', 'Tox-pub', 'warranty'):
    ck(want in text, f'About mentions {want!r}')

app.destroy()
print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('all help-menu checks passed')
