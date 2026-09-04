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
# The Markdown documents now go to a browser, so the reader is their fallback
# rather than their normal path. It still has to work, and this is where that
# is checked - with the browser refused, so the suite never launches one. A
# test that opened the real browser would also be testing the machine.
import webbrowser                                                   # noqa: E402

_real_open = webbrowser.open
webbrowser.open = lambda *a, **k: False

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

print('\n=== 5b. when a browser IS available, the documents go to it ===')
asked = []
webbrowser.open = lambda url, *a, **k: (asked.append(url), True)[1]
for name, fn in (('MeSH Workbench Manual', app.open_manual),
                 ('Installing and Updating',
                  lambda: app.open_doc('INSTALL.md', 'Installing and Updating'))):
    app._readers.pop(name, None)
    fn()
    app.update_idletasks()
    ck(bool(asked), f'{name}: a browser was asked', f'asked {asked}')
    if asked:
        page = asked[-1]
        ck(page.startswith('file:') and page.endswith('.html'),
           f'{name}: for a local page, not a Markdown file', page)
        ck(name not in app._readers,
           f'{name}: and no reader window was opened as well')
    asked.clear()

# Opening one document has to write the others out too. HELP.md links to
# INSTALL.md and COMMAND-LINE.md; without the siblings on disk those links land
# on nothing, which is worse than the reader this replaced.
from pathlib import Path                                            # noqa: E402
from mesh_aop import paths as _paths                                # noqa: E402

folder = Path(_paths.user_root()) / 'manual'
app._readers.pop('MeSH Workbench Manual', None)
app.open_manual()
app.update_idletasks()
for stem in ('HELP', 'INSTALL', 'COMMAND-LINE', 'README'):
    page = folder / f'{stem}.html'
    ck(page.exists(), f'{stem}.html was written beside the manual', str(page))

manual = (folder / 'HELP.html').read_text(encoding='utf-8')
for target in ('INSTALL.html', 'COMMAND-LINE.html'):
    ck(f'href="{target}"' in manual, f'the manual links to {target}')
    ck((folder / target).exists(), f'and {target} is there to be opened')
ck('href="INSTALL.md"' not in manual, 'no link still points at the source file')

# License is not Markdown, so it has nothing to convert and stays in the reader.
asked.clear()
app._readers.pop('License', None)
app.open_license()
app.update_idletasks()
ck(not asked, 'License is not converted; it stays in the reader', f'asked {asked}')
ck('License' in app._readers, 'and its window still opens')

webbrowser.open = _real_open

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
