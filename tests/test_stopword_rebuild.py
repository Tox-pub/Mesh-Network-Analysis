"""The rebuild, end to end: MeSH XML -> terms CSV -> stop-word set.

This is the chain the Stop words tab depends on. If the tree column is not
written, or is written wrong, every checkbox silently falls back to the shipped
list - so it is checked against a real XML parse rather than a hand-made CSV.
"""
import csv
import os
import sys
import tempfile

sys.path.insert(0, 'src')
from mesh_aop import vocabulary                                     # noqa: E402
from mesh_aop.mesh_data_processor import process_raw_mesh_data      # noqa: E402

FAILS = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + ('' if ok or not extra
                                                    else f'   -- {extra}'))
    if not ok:
        FAILS.append(msg)


# A small but honest descriptor file: one term per tree, one in two trees, one
# with no tree numbers at all, and the two check tags.
DESCRIPTORS = [
    ('Skin', 'D012867', ['A17.360']),
    ('Escherichia coli', 'D004926', ['B03.440']),
    ('Dermatitis', 'D003872', ['C17.800.174']),
    ('Nickel', 'D009532', ['D01.268.556']),
    ('Biopsy', 'D001706', ['E01.370.225']),
    ('Anxiety', 'D001007', ['F01.470.132']),
    ('Apoptosis', 'D017209', ['G04.299.134']),
    ('Nurses', 'D009726', ['H02.403', 'M01.526.485']),
    ('Culture', 'D003468', ['I01.076']),
    ('Agriculture', 'D000383', ['J01.086']),
    ('Ethics', 'D004989', ['K01.316']),
    ('Databases, Factual', 'D016208', ['L01.313.500']),
    ('Sweden', 'D013548', ['Z01.542.248.819']),
    ('Review', 'D016454', ['V02.700']),
    ('Health Services', 'D006296', ['N02.421']),
    ('Inflammation', 'D007249', ['C15.378', 'G01.750']),
    ('Orphan Concept', 'D999999', []),
    ('Male', 'D008297', []),
    ('Female', 'D005260', []),
]

tmp = tempfile.mkdtemp()
xml_path = os.path.join(tmp, 'desc2025.xml')
with open(xml_path, 'w', encoding='utf-8') as fh:
    fh.write('<?xml version="1.0"?>\n<DescriptorRecordSet>\n')
    for name, ui, trees in DESCRIPTORS:
        fh.write('  <DescriptorRecord>\n')
        fh.write(f'    <DescriptorUI>{ui}</DescriptorUI>\n')
        fh.write(f'    <DescriptorName><String>{name}</String></DescriptorName>\n')
        if trees:
            fh.write('    <TreeNumberList>\n')
            for t in trees:
                fh.write(f'      <TreeNumber>{t}</TreeNumber>\n')
            fh.write('    </TreeNumberList>\n')
        fh.write('  </DescriptorRecord>\n')
    fh.write('</DescriptorRecordSet>\n')

out_csv = os.path.join(tmp, 'mesh_terms.csv')
out_py = os.path.join(tmp, 'mesh_stop_words.py')

print('=== 1. the rebuild runs and writes both files ===')
process_raw_mesh_data(xml_file=xml_path, output_csv=out_csv, output_py=out_py,
                      force_update=True)
ck(os.path.exists(out_csv), 'the terms CSV was written')
ck(os.path.exists(out_py), 'the stop-word module was written')

print('\n=== 2. the CSV carries the tree column, correctly filled ===')
with open(out_csv, encoding='utf-8', newline='') as fh:
    rows = {r['DescriptorName']: r for r in csv.DictReader(fh)}
ck(vocabulary.TREE_COLUMN in next(iter(rows.values())),
   f'{vocabulary.TREE_COLUMN} is a column')
ck(len(rows) == len(DESCRIPTORS), f'every descriptor is present: {len(rows)}')
for name, _ui, trees in DESCRIPTORS:
    want = ''.join(sorted({t[0] for t in trees}))
    got = rows.get(name, {}).get(vocabulary.TREE_COLUMN, None)
    ck(got == want, f'{name:<20} trees {got!r}', f'expected {want!r}')

print('\n=== 3. needs_tree_column no longer fires on the rebuilt file ===')
ck(not vocabulary.needs_tree_column(out_csv),
   'the rebuilt CSV is recognised as current')

print('\n=== 4. the default selection reproduces the old behaviour ===')
words, prov = vocabulary.resolve(mesh_csv_path=out_csv)
ck('shipped' not in prov['source'], f'it uses the CSV: {prov["source"]}')
for keep in ('Skin', 'Dermatitis', 'Nickel', 'Apoptosis', 'Inflammation'):
    ck(keep not in words, f'{keep} survives (A/C/D/G)')
for drop in ('Escherichia coli', 'Sweden', 'Review', 'Culture', 'Agriculture',
             'Ethics', 'Databases, Factual', 'Biopsy', 'Anxiety',
             'Health Services', 'Orphan Concept'):
    ck(drop in words, f'{drop} is excluded')
ck('Nurses' in words, 'Nurses (H + M, both excluded) is excluded')

print('\n=== 5. changing the selection changes the result, no XML needed ===')
w2, p2 = vocabulary.resolve(excluded_trees='B', mesh_csv_path=out_csv)
ck('Sweden' not in w2, 'un-ticking Geographicals keeps Sweden')
ck('Review' not in w2, 'and Publication Characteristics keeps Review')
ck('Escherichia coli' in w2, 'while Organisms stays excluded')
ck('warning' not in p2, 'and no fallback warning is raised')
ck('Orphan Concept' in w2, 'a descriptor in no tree is still excluded')

print('\n=== 6. the sexes, which the XML gives no tree numbers ===')
ck(rows['Male'][vocabulary.TREE_COLUMN] == '',
   'Male has no tree numbers in the XML, as expected')
ck('Male' in words and 'Female' in words, 'both excluded by default')
w3, _ = vocabulary.resolve(excluded_trees='', mesh_csv_path=out_csv,
                           keep_sexes=True)
ck('Male' not in w3 and 'Female' not in w3, 'and kept only by the switch')

print('\n=== 7. the generated module still imports ===')
ns = {}
exec(compile(open(out_py, encoding='utf-8').read(), out_py, 'exec'), ns)
ck('MESH_STOP_WORDS' in ns, 'it defines MESH_STOP_WORDS')
ck(isinstance(ns['MESH_STOP_WORDS'], list),
   f'as a list of {len(ns.get("MESH_STOP_WORDS", []))} terms')
ck('Sweden' in ns['MESH_STOP_WORDS'], 'containing the excluded terms')

import shutil                                                       # noqa: E402
shutil.rmtree(tmp, ignore_errors=True)

print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('the rebuild produces a terms file the tree checkboxes can actually use')
