"""The manual opens as a page whose table of contents actually works.

A 900-line reference with a contents list is only useful if the links jump.
The reader window could not do that, so the documents are rendered to HTML and
handed to a browser, with the reader kept as the fallback.

The link check is the one that matters, and it is checked against the real
HELP.md rather than a fixture: every `](#anchor)` in the document must match an
`id=` the renderer emits. A slug rule that disagrees with the one the document
was written against fails silently, which is exactly the fault this replaces.
"""
import os
import re
import sys

sys.path.insert(0, 'src')

from mesh_aop import mdhtml                                           # noqa: E402

FAILS = []


def ck(ok, msg, extra=''):
    """`extra` says what went wrong, so it is only printed when something did."""
    print(('  [OK] ' if ok else '  [XX] ') + msg + ('' if ok or not extra
                                                    else f'   -- {extra}'))
    if not ok:
        FAILS.append(msg)


DOCS = ['HELP.md', 'INSTALL.md', 'README.md', 'COMMAND-LINE.md']

print('=== 1. every internal link resolves to a heading it emits ===')
for doc in DOCS:
    if not os.path.exists(doc):
        continue
    text = open(doc, encoding='utf-8').read()
    page = mdhtml.render_page(text, doc, doc)
    ids = set(re.findall(r'<h[1-6] id="([^"]+)"', page))
    wanted = set(re.findall(r'\]\(#([^)]+)\)', text))
    missing = sorted(wanted - ids)
    ck(not missing, f'{doc:<18} {len(wanted)} link(s) -> {len(ids)} heading(s)',
       f'unresolved: {missing[:6]}')

print('\n=== 2. the page stands alone ===')
page = mdhtml.render_page(open('HELP.md', encoding='utf-8').read(), 'Manual', 'HELP.md')
ck('<style>' in page, 'the stylesheet is inlined')
for bad in ('http://cdn', 'https://cdn', '<link rel="stylesheet"', '<script'):
    ck(bad not in page, f'nothing fetched or executed: no {bad!r}')
ck(page.startswith('<!doctype html>'), 'it is a complete document')
ck('<meta charset="utf-8">' in page, 'with an encoding declared')
ck('prefers-color-scheme' in page, 'and it follows the reader\'s theme')

print('\n=== 3. the structures these documents actually use ===')
SAMPLE = '''# Title

Intro **bold**, *italic*, `code`, and a [link](https://example.com).

## A Section

| Column | Meaning |
| :--- | ---: |
| `a` | first |
| b | second |

- one
- two
  - nested

1. first
2. second

> A quotation.

```
literal <not> markup
```

---
'''
out = mdhtml.to_html(SAMPLE)
for probe, what in (
    ('<h1 id="title">', 'h1 with an anchor'),
    ('<h2 id="a-section">', 'h2 with an anchor'),
    ('<strong>bold</strong>', 'bold'),
    ('<em>italic</em>', 'italic'),
    ('<code>code</code>', 'inline code'),
    ('<a href="https://example.com">link</a>', 'links'),
    ('<table>', 'tables'),
    ('<th style="text-align:left">Column</th>', 'table headers with alignment'),
    ('<td style="text-align:right">first</td>', 'right alignment from the rule row'),
    ('<ul>', 'bullet lists'),
    ('<ol>', 'numbered lists'),
    ('<blockquote>', 'block quotes'),
    ('<pre><code>', 'code blocks'),
    ('<hr>', 'rules'),
):
    ck(probe in out, what, f'missing {probe!r}')

ck('&lt;not&gt;' in out, 'markup inside a code block is escaped, not executed')

print('\n=== 3b. a code span is literal, whatever is inside it ===')
# The manual names a file pattern and a column pattern in one sentence. With
# code spans left in the emphasis pass, the asterisk closing the first paired
# with the one opening the second and italicised the words between them.
pair = mdhtml.to_html('Written to `*_mean_relevancy.db` as `score_*` columns.')
ck('<code>*_mean_relevancy.db</code>' in pair, 'the first span survives whole',
   pair)
ck('<code>score_*</code>' in pair, 'and so does the second', pair)
ck('<em>' not in pair, 'nothing between them was italicised', pair)
ck('*bold*' not in mdhtml.to_html('`**not bold**`'),
   'asterisks inside a span are not emphasis')
ck('<code>[a](b)</code>' in mdhtml.to_html('`[a](b)`'),
   'nor is a link written inside one')
ck('<a href="x"><code>y</code></a>' in mdhtml.to_html('[`y`](x)'),
   'but a span used AS a link label still links')

for doc in DOCS:
    if not os.path.exists(doc):
        continue
    page = mdhtml.render_page(open(doc, encoding='utf-8').read(), doc, doc)
    bad = re.findall(r'<code>[^<]*<(?!/code>)[^>]*>', page)
    ck(not bad, f'{doc:<18} no span in the real document is broken open',
       f'{len(bad)}: {bad[:3]}')

print('\n=== 3c. links between the documents point at the pages, not the source ===')
for url, want in (
    ('INSTALL.md', 'INSTALL.html'),
    ('COMMAND-LINE.md', 'COMMAND-LINE.html'),
    ('HELP.md#citation', 'HELP.html#citation'),
    ('INSTALL.md#macos', 'INSTALL.html#macos'),
    ('#ground-truth', '#ground-truth'),                 # within the page
    ('https://example.com/a.md', 'https://example.com/a.md'),   # not ours
    ('https://aopwiki.org/aops/40', 'https://aopwiki.org/aops/40'),
):
    got = mdhtml._local_href(url)
    ck(got == want, f'{url}  ->  {got}', f'expected {want}')

ck('<a href="INSTALL.html">' in mdhtml.to_html('see [INSTALL.md](INSTALL.md)'),
   'and the rewrite reaches the rendered link')

# The anchor a cross-document link asks for has to exist in THAT document.
pages = {d: mdhtml.render_page(open(d, encoding='utf-8').read(), d, d)
         for d in DOCS if os.path.exists(d)}
for doc, text in ((d, open(d, encoding='utf-8').read()) for d in pages):
    for target, frag in re.findall(r'\]\(([A-Za-z-]+\.md)#([^)]+)\)', text):
        if target not in pages:
            ck(False, f'{doc} links into {target}, which is not shipped')
            continue
        ck(f'id="{frag}"' in pages[target],
           f'{doc} -> {target}#{frag} lands on a real heading')

print('\n=== 4. slugs match the rule the documents were written to ===')
for heading, expect in (
    ('Ground Truth', 'ground-truth'),
    ('Validation & Benchmarking', 'validation--benchmarking'),
    ('The Run Ledger and PRISMA-Style Flow Overview',
     'the-run-ledger-and-prisma-style-flow-overview'),
    ('Data Acquisition & Prerequisites', 'data-acquisition--prerequisites'),
    ('The AOP Annotation Workflow (Biological Strata)',
     'the-aop-annotation-workflow-biological-strata'),
    ('`code` in a heading', 'code-in-a-heading'),
    ('**bold** heading', 'bold-heading'),
):
    got = mdhtml.slug(heading)
    ck(got == expect, f'{heading!r} -> {got!r}', f'expected {expect!r}')

print('\n=== 5. the application falls back rather than failing ===')
import inspect                                                      # noqa: E402
from mesh_workbench.app import Workbench                            # noqa: E402

src = inspect.getsource(Workbench.open_doc)
ck('_open_as_html' in src, 'the browser is tried first')
ck('self.show_reader(' in src, 'and the reader window is still reached')
ck(src.index('_open_as_html') < src.index('self.show_reader(heading'),
   'in that order')

fb = inspect.getsource(Workbench._open_as_html)
ck(fb.count('return False') >= 3,
   'every failure path returns False rather than raising')
ck('webbrowser' in fb, 'it asks for a browser, not the file-type handler')
ck('.md' in src, 'only Markdown is converted; LICENSE stays in the reader')

print('\n=== 6. a real document renders without loss ===')
text = open('HELP.md', encoding='utf-8').read()
page = mdhtml.render_page(text, 'Manual', 'HELP.md')
heads_md = len(re.findall(r'^#{1,6} ', text, re.M))
heads_html = len(re.findall(r'<h[1-6] id=', page))
ck(heads_html >= heads_md - 2,
   f'{heads_html} of {heads_md} headings survive (fenced examples excepted)')
tables_md = len(re.findall(r'^\|[^\n]*\n\s*\|?[\s:|-]+\|', text, re.M))
ck(page.count('<table>') >= tables_md - 1,
   f'{page.count("<table>")} of about {tables_md} tables rendered')
ck(len(page) > len(text), f'the page is not truncated: {len(page):,} chars')

print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('all manual-rendering checks passed')
