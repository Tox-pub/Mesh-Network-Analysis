"""Secondary analysis has to work when no retrieval ever happened.

A run against the bundled reference corpus ships a finished network and no
article lists, so there is no cleaned citation database. Refusing outright was
wrong - the pipeline can generate what it needs, which is the whole point of a
demonstration - and scoring against zero citations is worse than refusing: it
drives Normalized_Citation to 1.0 for every article, clips every relevance score
to the same floor, and both impact scores come out identical, so the CSV looks
ranked and is noise.

So citations for the shortlist are fetched and cached. Checked here without
touching the network: the on-disk path is preferred, only genuinely unknown
articles are looked up, the fetch is bounded, the cache is written in the shape
the retrieval path writes, and a second query reads it instead of refetching.
"""
import os
import sqlite3
import sys
import tempfile

sys.path.insert(0, 'src')

import mesh_aop.secondary_analysis as sa                          # noqa: E402
from mesh_aop.secondary_analysis import (                          # noqa: E402
    _cache_citations, _ensure_citations, _require_databases, _calculate_impact_scores)

FAILS = []
CALLS = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + (f'   -- {extra}' if extra and not ok else ''))
    if not ok:
        FAILS.append(msg)


def fake_icite(pmids):
    """Stand in for the iCite API - the point is what gets ASKED for."""
    CALLS.append(list(pmids))
    return {p: {'cited_by': ';'.join(str(p * 10 + i) for i in range(3)), 'cites': ''}
            for p in pmids}


import mesh_aop.data_ops as data_ops                                # noqa: E402
data_ops.fetch_links_in_batches = fake_icite

box = tempfile.mkdtemp()

print('=== 1. the relevance DB is required; the cleaned one is not ===')
rel = os.path.join(box, 'P_mean_relevancy.db')
open(rel, 'w').close()
try:
    _require_databases(rel, os.path.join(box, 'nope.db'))
    ck(True, 'a missing cleaned database is no longer fatal')
except FileNotFoundError as e:
    ck(False, 'a missing cleaned database is no longer fatal', str(e)[:60])
try:
    _require_databases(os.path.join(box, 'gone.db'), None)
    ck(False, 'a missing relevance database still raises')
except FileNotFoundError as e:
    ck('relevance database' in str(e) and 'prefix' in str(e),
       'a missing relevance database still raises, and explains the prefix')

print('\n=== 2. nothing on disk: every article is fetched, then cached ===')
CALLS.clear()
clean = os.path.join(box, 'P_cleaned_pmids.db')
got = _ensure_citations([11, 22, 33], clean)
ck(sorted(got) == ['11', '22', '33'], f'all three returned: {sorted(got)}')
ck(got['11'] == '110;111;112', f"cited_by carried through: {got['11']}")
ck(os.path.exists(clean), 'a cleaned database was created')
ck(len(CALLS) == 1 and sorted(CALLS[0]) == [11, 22, 33],
   f'one batched call for the three unknowns: {CALLS}')

print('\n=== 3. the cache is the shape retrieval writes ===')
conn = sqlite3.connect(clean)
cols = [r[1] for r in conn.execute("PRAGMA table_info(pmids_table)")]
for c in ('pmid', 'generation', 'cited_by'):
    ck(c in cols, f'column {c!r} present', f'columns are {cols}')
rows = conn.execute("SELECT pmid, generation, cited_by FROM pmids_table ORDER BY pmid").fetchall()
ck(len(rows) == 3, f'three rows written: {len(rows)}')
ck(all(r[1] == 'P0' for r in rows), 'generation recorded')
conn.close()

print('\n=== 4. a second query reads the cache and fetches nothing ===')
CALLS.clear()
again = _ensure_citations([11, 22, 33], clean)
ck(CALLS == [], f'no fetch on the second call: {CALLS}')
ck(again['22'] == got['22'], 'same answer from disk')

print('\n=== 5. only the genuinely unknown ones are fetched ===')
CALLS.clear()
mixed = _ensure_citations([11, 22, 33, 44, 55], clean)
ck(len(CALLS) == 1 and sorted(CALLS[0]) == [44, 55],
   f'only the two new ones asked for: {CALLS}')
ck(len(mixed) == 5, f'all five returned: {len(mixed)}')

print('\n=== 6. the fetch is bounded ===')
CALLS.clear()
big = list(range(1000, 1000 + 500))
_ensure_citations(big, os.path.join(box, 'bounded.db'), fetch_budget=120)
asked = sum(len(c) for c in CALLS)
ck(asked == 120, f'budget of 120 respected, {asked} asked for')
ck(all(len(c) <= sa._ICITE_BATCH for c in CALLS),
   f'batched at {sa._ICITE_BATCH} per request, largest was {max(len(c) for c in CALLS)}')

print('\n=== 7. a retrieved corpus fetches nothing at all ===')
CALLS.clear()
own = os.path.join(box, 'own_cleaned.db')
c = sqlite3.connect(own)
c.execute("CREATE TABLE pmids_table (pmid INTEGER PRIMARY KEY, generation TEXT, "
          "cited_by TEXT, cites TEXT, mesh_terms TEXT)")
c.executemany("INSERT INTO pmids_table VALUES (?,'P0',?,'','')",
              [(p, f'{p}1;{p}2') for p in (7, 8, 9)])
c.commit(); c.close()
res = _ensure_citations([7, 8, 9], own)
ck(CALLS == [], 'an already-retrieved corpus makes no network call')
ck(res['8'] == '81;82', 'and reads its own citations')

print('\n=== 8. real citations keep the ranking from collapsing ===')
# what used to happen with no citation data at all
flat = _calculate_impact_scores(
    [(p, 0.5, 'Skin', None, 2015) for p in (1, 2, 3)], 'F1', 0.5, 10)
ck(len({r['F1_AIS'] for r in flat}) == 1,
   'with no citations every score is identical - the failure this avoids')
real = _calculate_impact_scores(
    [(1, 0.5, 'Skin', '10;11;12;13', 2015),
     (2, 0.5, 'Skin', '20', 2015),
     (3, 0.5, 'Skin', '', 2015)], 'F1', 0.5, 10)
ck(len({r['F1_AIS'] for r in real}) == 3,
   f"with citations the three separate: {sorted(round(r['F1_AIS'], 4) for r in real)}")
ck(real[0]['pmid'] == 1, f"the most-cited ranks first: pmid {real[0]['pmid']}")

print('\n=== 9. an uncited article never outranks a cited one ===')
# The floor standing in for "no citations" was log10(2)/max_log/2, which is
# ABOVE 1 as soon as the best article in the set manages under 0.41 citations a
# year. It then beat every real score, and the clip below dragged all the
# relevance scores up to it, so the ranking inverted and ARS stopped mattering.
for label, rows in (
    ('low-citation set (the regime that broke)',
     [(1, 0.5, 'S', '10;11;12;13', 2015), (2, 0.5, 'S', '20', 2015),
      (3, 0.5, 'S', '', 2015)]),
    ('well-cited set',
     [(1, 0.5, 'S', ';'.join(str(i) for i in range(400)), 2015),
      (2, 0.5, 'S', ';'.join(str(i) for i in range(40)), 2015),
      (3, 0.5, 'S', '', 2015)]),
    ('nothing cited at all',
     [(1, 0.9, 'S', '', 2015), (2, 0.5, 'S', '', 2015), (3, 0.1, 'S', '', 2015)]),
):
    out = _calculate_impact_scores(rows, 'F1', 0.5, 10)
    order = [r['pmid'] for r in out]
    uncited = [r for r in out if r['Incoming_Citations_Smoothed'] == 0]
    cited = [r for r in out if r['Incoming_Citations_Smoothed'] > 0]
    print(f'    {label}: order {order}')
    if cited and uncited:
        ck(max(u['Normalized_Citation'] for u in uncited)
           < min(c['Normalized_Citation'] for c in cited),
           '      uncited scores below every cited one')
        ck(order[0] != 3, '      the uncited article does not rank first')
    ck(all(r['Normalized_Citation'] <= 1.0 + 1e-9 for r in out),
       '      no normalised score exceeds 1.0')

# relevance has to survive the floor: equal citations, higher ARS wins
same = _calculate_impact_scores(
    [(1, 0.2, 'S', '10;11', 2015), (2, 0.8, 'S', '20;21', 2015)], 'F1', 0.5, 10)
ck(same[0]['pmid'] == 2,
   f"equal citations: the more relevant article wins (got {same[0]['pmid']})")

print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('all citation-fallback checks passed')
