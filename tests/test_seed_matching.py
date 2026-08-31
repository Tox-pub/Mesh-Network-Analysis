"""Asking secondary analysis for a node must not return a different node.

contributing_seeds is a ';'-separated list of exact MeSH headings. The queries
used to match it with LIKE '%term%', so a request for 'Skin' also returned
articles seeded by 'Skin Diseases', 'Skin Absorption' and 'Skin Tests'. Nothing
failed and nothing warned - the CSV simply contained the wrong articles, which
is the worst way for this to be wrong.

MeSH is full of headings that are prefixes of other headings, so this is the
normal case rather than a corner: Skin, Haptens, Allergens, Antigens and
Dermatitis all begin longer headings in the same vocabulary.
"""
import os
import sqlite3
import sys
import tempfile

sys.path.insert(0, 'src')

from mesh_aop.secondary_analysis import _fetch_relevance_rows, _seed_match  # noqa: E402

FAILS = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + (f'   -- {extra}' if extra and not ok else ''))
    if not ok:
        FAILS.append(msg)


ROWS = [
    (1,  'Skin'),
    (2,  'Haptens;Skin'),
    (3,  'Skin Diseases'),
    (4,  'Allergens;Skin Absorption'),
    (5,  'Dermatitis, Allergic Contact;Skin Tests'),
    (6,  'Haptens'),
    (7,  'Dermatitis, Allergic Contact;Skin'),
    (8,  'Dermatitis, Allergic Contact;Skin Diseases'),
    (9,  'Dermatitis, Allergic Contact, Systemic'),
    (10, 'Haptens, Incomplete'),
    (11, 'Skin_Test'),          # a literal underscore, LIKE's single-char wildcard
    (12, 'Skin 100%'),          # a literal percent, LIKE's any-run wildcard
]

box = tempfile.mkdtemp()
db = os.path.join(box, 'rel.db')
con = sqlite3.connect(db)
con.execute("CREATE TABLE article_relevance_scores (pmid INTEGER PRIMARY KEY, "
            "score_pagerank_centrality REAL, contributing_seeds TEXT, pub_date TEXT)")
for pmid, seeds in ROWS:
    con.execute("INSERT INTO article_relevance_scores VALUES (?,?,?,?)",
                (pmid, 1.0 / pmid, seeds, '2020'))
con.commit()

SEEDS = dict(ROWS)


def ask(*terms):
    return sorted(r[0] for r in _fetch_relevance_rows(con, *_seed_match(*terms)))


def expected(*terms):
    return sorted(p for p, s in ROWS if all(t in s.split(';') for t in terms))


print('=== 1. a node returns exactly the articles seeded by that heading ===')
for term in ('Skin', 'Haptens', 'Skin Diseases', 'Dermatitis, Allergic Contact',
             'Skin Tests', 'Allergens'):
    got, want = ask(term), expected(term)
    ck(got == want, f'{term!r:34} -> {got}', f'expected {want}')

print('\n=== 2. the headings that used to leak in are gone ===')
skin = ask('Skin')
for pmid in (3, 4, 5, 8):
    ck(pmid not in skin,
       f'pmid {pmid} ({SEEDS[pmid]!r}) is not returned for "Skin"')
hap = ask('Haptens')
ck(10 not in hap, f'pmid 10 ({SEEDS[10]!r}) is not returned for "Haptens"')
dac = ask('Dermatitis, Allergic Contact')
ck(9 not in dac, f'pmid 9 ({SEEDS[9]!r}) is not returned for that heading')

print('\n=== 3. nothing genuine was dropped ===')
ck(sorted(skin) == [1, 2, 7], f'"Skin" returns every article truly seeded by it: {skin}')
ck(sorted(hap) == [2, 6], f'"Haptens" returns {hap}')
ck(sorted(dac) == [5, 7, 8, 9] or sorted(dac) == [5, 7, 8],
   f'"Dermatitis, Allergic Contact" returns {dac}')

print('\n=== 4. edges need both headings on the same article ===')
for a, b in (('Skin', 'Dermatitis, Allergic Contact'),
             ('Haptens', 'Skin'),
             ('Allergens', 'Skin')):
    got, want = ask(a, b), expected(a, b)
    ck(got == want, f'{a!r} + {b!r} -> {got}', f'expected {want}')
ck(ask('Skin', 'Dermatitis, Allergic Contact') == [7],
   'the edge the user asked for returns only the article that has both')

print('\n=== 5. LIKE wildcards inside a heading are literal ===')
ck(ask('Skin_Test') == [11], f'underscore matches only itself: {ask("Skin_Test")}')
ck(ask('Skin Test') == [], f'"Skin Test" does not match "Skin_Test": {ask("Skin Test")}')
ck(ask('Skin 100%') == [12], f'percent matches only itself: {ask("Skin 100%")}')
ck(ask('Skin 100') == [], f'"Skin 100" does not match "Skin 100%": {ask("Skin 100")}')

print('\n=== 6. a heading that is in no article returns nothing, not everything ===')
for term in ('Kinase', 'kin', '', 'Ski'):
    ck(ask(term) == [], f'{term!r} returns nothing', f'got {ask(term)}')

print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('all seed-matching checks passed')
