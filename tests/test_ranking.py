"""The benchmark's ranking shortcut has to be the sort it replaced.

_positions used to rank the whole pool with a random tie-break key. It now
counts instead of sorting, which is several times cheaper but is only allowed
if it produces the same answer. "The same" means two different things and both
are checked here:

  no ties   the ranks must be IDENTICAL, whatever the seed. Nothing is random
            when no two scores are equal, so any disagreement is a bug.
  ties      the ranks are random by design, so the requirement is that they
            come from the same DISTRIBUTION as the sort would have drawn from.
            Checked by running both many times and comparing per-positive mean
            ranks, and the reported metrics, against Monte Carlo error.

Ties are not an edge case here: one of the thirteen scorers counts network
terms, so its scores are small integers and its tie groups hold millions.

The ground-truth set is also the user's to supply, so every check runs at 96
positives and again at 20,000. A method that is only correct for a small,
mostly-untied ground truth would be a trap for anyone bringing their own.
"""
import sys

import numpy as np

sys.path.insert(0, 'src')

from mesh_aop.validation_report import _positions, _metric_suite   # noqa: E402

FAILS = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + (f'   -- {extra}' if extra and not ok else ''))
    if not ok:
        FAILS.append(msg)


def by_sorting(scores, mask, is_pos, rng):
    """What _positions did before: rank the pool, random tie-break."""
    s = scores[mask]
    order = np.lexsort((rng.random(s.size), -s))
    rank = np.empty(s.size, dtype=np.int64)
    rank[order] = np.arange(1, s.size + 1)
    return rank[is_pos[mask]].astype(float), s.size


def positives(n, k, rng):
    m = np.zeros(n, bool)
    m[rng.choice(n, k, replace=False)] = True
    return m


N = 120_000
ALL = np.ones(N, bool)

print('=== 1. no ties: identical to the sort, for any seed ===')
rng = np.random.default_rng(0)
distinct = rng.permutation(N).astype(float)
for k in (1, 96, 20_000):
    ip = positives(N, k, rng)
    for seed_a, seed_b in ((1, 1), (1, 2), (7, 99)):
        a, na = by_sorting(distinct, ALL, ip, np.random.default_rng(seed_a))
        b, nb = _positions(distinct, ALL, ip, np.random.default_rng(seed_b))
        ck(np.array_equal(a, b) and na == nb,
           f'{k:>6,} positives, seeds {seed_a}/{seed_b}: ranks identical',
           f'max diff {np.abs(a - b).max():.0f}')

print('\n=== 2. ties: ranks stay inside the tie group, and never repeat ===')
rng = np.random.default_rng(1)
lumpy = rng.poisson(1.2, N).astype(float)          # the integer count scorer
vals, cnt = np.unique(lumpy, return_counts=True)
print(f'  ({vals.size} distinct scores over {N:,} articles; '
      f'largest tie group {cnt.max():,})')
for k in (96, 20_000):
    ip = positives(N, k, rng)
    r, n = _positions(lumpy, ALL, ip, np.random.default_rng(4))
    pv = lumpy[ip]
    lo = np.array([(lumpy > v).sum() + 1 for v in pv])
    hi = np.array([(lumpy >= v).sum() for v in pv])
    ck(bool(((r >= lo) & (r <= hi)).all()),
       f'{k:>6,} positives: every rank within its own tie group')
    ck(len(set(r.tolist())) == r.size,
       f'{k:>6,} positives: no two positives share a rank',
       f'{r.size - len(set(r.tolist()))} collisions')
    ck(n == N, f'{k:>6,} positives: pool size reported as {N:,}')

print('\n=== 3. ties: same distribution as the sort ===')
REPS = 250
for k in (96, 20_000):
    ip = positives(N, k, np.random.default_rng(5))
    old = np.empty((REPS, k))
    new = np.empty((REPS, k))
    for r in range(REPS):
        old[r], _ = by_sorting(lumpy, ALL, ip, np.random.default_rng(2_000 + r))
        new[r], _ = _positions(lumpy, ALL, ip, np.random.default_rng(9_000 + r))

    se = np.sqrt(old.var(0, ddof=1) / REPS + new.var(0, ddof=1) / REPS)
    z = np.where(se > 0, (old.mean(0) - new.mean(0)) / np.where(se > 0, se, 1), 0.0)
    # |z| > 4 is a 1-in-16,000 event per positive; at 20,000 positives roughly
    # one is expected, so the check is on the count, not on a single worst case.
    ck(int((np.abs(z) > 4).sum()) <= max(3, k // 5_000),
       f'{k:>6,} positives: mean ranks agree (max |z| {np.abs(z).max():.2f}, '
       f'{int((np.abs(z) > 4).sum())} over 4)')

    for metric in ('BEDROC', 'EF@0.01', 'MAP'):
        vo = np.array([_metric_suite(o, N)[metric] for o in old])
        vn = np.array([_metric_suite(x, N)[metric] for x in new])
        sd = np.sqrt(vo.var(ddof=1) / REPS + vn.var(ddof=1) / REPS)
        zz = abs(vo.mean() - vn.mean()) / sd if sd > 0 else 0.0
        ck(zz < 4, f'{k:>6,} positives: {metric} agrees (z {zz:.2f})',
           f'{vo.mean():.6f} vs {vn.mean():.6f}')

print('\n=== 4. the awkward shapes ===')
rng = np.random.default_rng(6)
none_pos = np.zeros(N, bool)
r, n = _positions(lumpy, ALL, none_pos, rng)
ck(r.size == 0 and n == N, 'no positives at all: empty ranks, pool size intact')

flat = np.zeros(N)
ip = positives(N, 96, rng)
r, _ = _positions(flat, ALL, ip, rng)
ck(bool(((r >= 1) & (r <= N)).all()) and len(set(r.tolist())) == 96,
   'every score identical: 96 distinct ranks spanning the pool')

everything = np.ones(N, bool)
r, n = _positions(lumpy, ALL, everything, np.random.default_rng(8))
ck(sorted(r.tolist()) == list(range(1, N + 1)),
   'every article a positive: ranks are exactly 1..N')

sub = np.zeros(N, bool)
sub[:5_000] = True
ip = np.zeros(N, bool)
ip[np.random.default_rng(9).choice(5_000, 40, replace=False)] = True
r, n = _positions(lumpy, sub, ip, np.random.default_rng(10))
ck(n == 5_000 and bool(((r >= 1) & (r <= 5_000)).all()),
   'masked frame: ranks are relative to the 5,000 articles in it')

small = np.array([3.0, 1.0, 3.0, 2.0])
ipx = np.array([True, False, False, False])
r, n = _positions(small, np.ones(4, bool), ipx, np.random.default_rng(11))
ck(n == 4 and r[0] in (1.0, 2.0), 'four articles: a tied top score ranks 1 or 2')

print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('all ranking checks passed')
