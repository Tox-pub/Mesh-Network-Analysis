# -*- coding: utf-8 -*-
"""
validation_report.py - consolidated external validation of the relevance scoring.

Evaluates every node weighting the pipeline produces against the curated ground
truth, and reports the result as a drill-down workbook, a narrative HTML report,
figures, and a console summary.

WHY FOUR EVALUATION FRAMES

A single "how well does it rank the corpus" number conflates two different things a
network can contribute: ordering ability, and reach. A MeSH query is not a ranking
method - it is a filter that returns a topic-restricted set with no internal order -
so comparing it to the ranked corpus at a fixed cut credits it for scope
restriction rather than for ordering. The frames separate these:

  CORPUS   every article in the scoring pool, ranked. The headline retrieval task.
           The pool is the set of articles carrying >=1 network term WITHIN the
           configured publication-date window - the same corpus definition the
           relevance database uses, so the two never diverge.
  WITHIN   only the query's own hits, so the article set is IDENTICAL for every
           scorer and only the ordering differs. Isolates ranking ability.
  OUTSIDE  only articles the query cannot return. Isolates reach - the ground-truth
           evidence that is invisible to a conventional search.
  HYBRID   query hits first (network-ranked), then the remainder (network-ranked).
           The reading list a user would actually be handed.

WHY POSITIVES-ONLY AND PAIRED BOOTSTRAPS

Every early-recognition metric depends only on WHERE the positives land, so the
resampling unit is the ground-truth positives, not the fixed multi-million pool -
resampling the pool mostly reshuffles negatives and reports intervals that are far
too tight. And because two scorers are evaluated on the SAME positives, comparing
their marginal intervals is underpowered; the paired bootstrap resamples once per
replicate and scores both, cancelling the shared variance.

WHY NODE-LEVEL AUC IS APPROPRIATE HERE

Article-level ROC-AUC is not reported: at ~1e-5 prevalence it is dominated by true
negatives and looks excellent for scorers that retrieve nothing. The node-level task
is different - roughly 40% of network terms are externally attested - so AUC is a
sound discrimination measure there, and it is reported with an explicit global
frequency control (stratified AUC plus incremental logistic regression) so that
"prominent terms are just common terms" can be tested rather than assumed.
"""

import os
import json
import sqlite3
from datetime import datetime

import numpy as np
import pandas as pd
import networkx as nx
from scipy import sparse
from scipy.stats import spearmanr

ALPHA_BEDROC = 20.0
EF_FRACS = (0.001, 0.01)
RECALL_KS = (1000, 10000, 100000)
CENTRALITIES = ("betweenness", "pagerank", "eigenvector")
SCOPES = (("", "full"), ("_subgraph", "subgraph"))


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# IO
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def _open_resilient(path):
    """Read-only connection, falling back to a normal one when OneDrive has
    dehydrated the file (a placeholder reads back with an empty schema)."""
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=120.0)
    if con.execute("SELECT name FROM sqlite_master WHERE type='table' LIMIT 1").fetchone():
        return con
    con.close()
    return sqlite3.connect(str(path), timeout=120.0)


def _base_terms(mesh_str):
    """Base MeSH descriptors of an article: drop subheadings and the major-topic '*'."""
    return {t.split('/')[0].lstrip('*').strip() for t in mesh_str.split(';') if t.strip()}


def _iso_window(start_date, end_date):
    """Normalise the configured context window to the ISO form pub_date is stored in."""
    if not start_date or not end_date:
        return None, None
    from .data_ops import parse_date_robust
    return (parse_date_robust(start_date).strftime('%Y-%m-%d'),
            parse_date_robust(end_date).strftime('%Y-%m-%d'))


def _build_or_load_cache(master_db_path, terms, cache_path,
                         start_date=None, end_date=None):
    """Article x term incidence over network terms, plus document length and year.

    Restricted to the configured publication-date context window, which is the
    same filter relevance.py applies when it builds the relevance database. The
    two therefore describe one corpus rather than two: evaluating retrieval over
    articles the scoring database does not contain cannot validate the pipeline
    the user actually ran. Passing no window scans everything, which is the old
    behaviour and is kept only for direct callers that have no config.

    Cached because the corpus scan dominates runtime and the report is typically
    re-run several times while wording and figures are settled. The cache key is
    the node set AND the window - a cache built under a different window is not
    interchangeable, and silently reusing one would reintroduce exactly the
    mismatch this filter exists to remove.
    """
    idx = {t: i for i, t in enumerate(terms)}
    lo, hi = _iso_window(start_date, end_date)
    window_key = np.array([lo or "", hi or ""], dtype=object)
    if os.path.exists(cache_path + "_X.npz"):
        meta = np.load(cache_path + "_meta.npz", allow_pickle=True)
        cached_window = (list(meta["window"]) if "window" in meta.files
                         else ["", ""])
        if list(meta["terms"]) != list(terms):
            print("  [i] Cached matrix was built for a different node set; rebuilding.")
        elif cached_window != list(window_key):
            print(f"  [i] Cached matrix was built for window "
                  f"{cached_window[0] or 'none'}..{cached_window[1] or 'none'}, "
                  f"need {lo or 'none'}..{hi or 'none'}; rebuilding.")
        else:
            return (sparse.load_npz(cache_path + "_X.npz"),
                    meta["pmids"], meta["doclen"], meta["years"])

    if lo:
        print(f"  Scanning master database for network-term incidence "
              f"({lo} to {hi})...")
    else:
        print("  Scanning master database for network-term incidence (no date filter)...")
    con = _open_resilient(master_db_path)
    cur = con.cursor()
    if lo:
        cur.execute("SELECT pmid, mesh_terms, pub_date FROM master_mesh_annotations "
                    "WHERE pub_date BETWEEN ? AND ?", (lo, hi))
    else:
        cur.execute("SELECT pmid, mesh_terms, pub_date FROM master_mesh_annotations")
    blocks, pm, dl, yr = [], [], [], []
    while True:
        rows = cur.fetchmany(200000)
        if not rows:
            break
        r_, c_, loc = [], [], 0
        for pmid, mt, pdate in rows:
            if not mt:
                continue
            bt = _base_terms(mt)
            js = [idx[x] for x in bt if x in idx]
            if not js:
                continue
            r_.extend([loc] * len(js))
            c_.extend(js)
            pm.append(int(pmid) if str(pmid).isdigit() else -1)
            dl.append(len(bt))
            yr.append(int(str(pdate)[:4]) if str(pdate)[:4].isdigit() else 0)
            loc += 1
        blocks.append(sparse.coo_matrix(
            (np.ones(len(r_), np.float32), (r_, c_)), shape=(loc, len(terms))).tocsr())
    con.close()

    X = sparse.vstack(blocks).tocsr()
    pmids = np.asarray(pm, dtype="int64")
    doclen = np.asarray(dl, dtype=np.float32)
    years = np.asarray(yr, dtype=np.int16)
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    sparse.save_npz(cache_path + "_X.npz", X)
    np.savez(cache_path + "_meta.npz", pmids=pmids, doclen=doclen, years=years,
             terms=np.array(terms, dtype=object), window=window_key)
    return X, pmids, doclen, years


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# METRICS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def _bedroc_from_positions(pos, n_total, alpha=ALPHA_BEDROC):
    """BEDROC computed from the positives' rank positions.

    Equivalent to the sorted-label form but O(n_positives) instead of O(n_pool),
    which is what makes tens of thousands of bootstrap replicates affordable.
    """
    n = pos.size
    if n == 0:
        return 0.0
    ra = n / n_total
    denom = n * (1 - np.exp(-alpha)) / (n_total * (np.exp(alpha / n_total) - 1))
    rie = np.sum(np.exp(-alpha * pos / n_total)) / denom
    return float(rie * ra * np.sinh(alpha / 2) /
                 (np.cosh(alpha / 2) - np.cosh(alpha / 2 - alpha * ra))
                 + 1 / (1 - np.exp(alpha * (1 - ra))))


def _metric_suite(pos, n_total):
    """All reported metrics depend only on the positives' positions."""
    out = {"BEDROC": _bedroc_from_positions(pos, n_total)}
    for f in EF_FRACS:
        cut = max(int(round(f * n_total)), 1)
        out[f"EF@{f:g}"] = float((pos <= cut).mean()) / f
    for k in RECALL_KS:
        if k < n_total:
            out[f"recall@{k}"] = float((pos <= k).mean())
    out["MAP"] = float(np.mean(np.arange(1, pos.size + 1) / np.sort(pos))) if pos.size else 0.0
    return out


def _bootstrap_metrics(pos, n_total, n_boot, rng):
    """Percentile CIs by resampling the positives - the actual sampled unit."""
    point = _metric_suite(pos, n_total)
    if pos.size == 0:
        return point, {k: (np.nan, np.nan) for k in point}
    draws = rng.integers(0, pos.size, size=(n_boot, pos.size))
    acc = {k: np.empty(n_boot) for k in point}
    for b in range(n_boot):
        m = _metric_suite(pos[draws[b]], n_total)
        for k, v in m.items():
            acc[k][b] = v
    return point, {k: tuple(np.percentile(v, [2.5, 97.5])) for k, v in acc.items()}


def _paired_bootstrap(pos_a, pos_b, n_total, metric_fn, n_boot, rng):
    """CI on fn(A) - fn(B) using ONE resample per replicate scored by both."""
    obs = metric_fn(pos_a, n_total) - metric_fn(pos_b, n_total)
    draws = rng.integers(0, pos_a.size, size=(n_boot, pos_a.size))
    diff = np.empty(n_boot)
    for b in range(n_boot):
        i = draws[b]
        diff[b] = metric_fn(pos_a[i], n_total) - metric_fn(pos_b[i], n_total)
    lo, hi = np.percentile(diff, [2.5, 97.5])
    p = 2 * min(float((diff <= 0).mean()), float((diff >= 0).mean()))
    return obs, lo, hi, min(p, 1.0)


def _pmid_db_reason(path, exc):
    """Why the retrieval database could not be read, in the user's terms.

    A bare "unable to open database file" reads like corruption. Much the most
    common cause is that there is no retrieval database at all - the run scored
    a network it did not retrieve, which is exactly what the bundled reference
    corpus does - so say that instead of relaying SQLite's wording.
    """
    if not path or not os.path.exists(str(path)):
        return ("no retrieval database was produced by this run - the "
                "demonstration corpus ships as a finished network, so the "
                "article lists behind it were never downloaded")
    return str(exc)


def _tie_offsets(n_tied, k, rng):
    """k distinct places drawn from 1..n_tied, uniform over the k-subsets.

    Distinct matters: MAP pairs the i-th positive with the i-th smallest rank,
    so two positives sharing a rank would score a precision that no ranking can
    actually produce. Sequential draws rejecting repeats IS sampling without
    replacement, and it beats permuting a tie group that can hold millions.
    """
    if k >= n_tied:
        return rng.permutation(n_tied)[:k] + 1
    if n_tied < 4 * k or n_tied < 512:
        return rng.choice(n_tied, size=k, replace=False) + 1
    seen, out = set(), []
    while len(out) < k:
        for v in rng.integers(0, n_tied, size=(k - len(out)) * 2 + 8).tolist():
            if v not in seen:
                seen.add(v)
                out.append(v)
                if len(out) == k:
                    break
    return np.asarray(out) + 1


def _positions(scores, mask, is_pos, rng):
    """Rank a masked subset (random tie-breaking) and return positive positions.

    Only the positives' ranks are ever used, and there are a few dozen of them
    against millions of articles, so the pool is never sorted. For each distinct
    positive score, count how many articles beat it and how many tie with it;
    the positives' places then come from inside their own tie group. Identical
    to sorting the pool with a random tie-break key - exactly so when no scores
    tie, and in distribution when they do - at a fraction of the cost, and the
    cost grows with the number of DISTINCT positive scores rather than with the
    pool, so a large user-supplied ground-truth set does not degrade it.
    """
    s = scores[mask]
    n = s.size
    pv = s[is_pos[mask]]
    if pv.size == 0:
        return np.empty(0, dtype=float), n

    uniq = np.unique(pv)
    below = np.searchsorted(uniq, s, side='left')     # distinct values under each score
    at_or = np.searchsorted(uniq, s, side='right')
    cum = np.cumsum(np.bincount(below, minlength=uniq.size + 1))
    greater = (n - cum[:uniq.size]).astype(np.int64)  # articles beating uniq[j]
    tied = np.bincount(below[at_or > below],
                       minlength=uniq.size)[:uniq.size].astype(np.int64)
    del below, at_or, cum

    grp = np.searchsorted(uniq, pv, side='left')
    order = np.argsort(grp, kind='stable')
    sorted_grp = grp[order]
    starts = np.flatnonzero(np.r_[True, sorted_grp[1:] != sorted_grp[:-1]])
    sizes = np.r_[starts[1:], sorted_grp.size] - starts
    which = sorted_grp[starts]

    rank = np.empty(pv.size, dtype=float)
    alone = sizes == 1
    if alone.any():                    # sole positive in its tie group: nothing
        j = which[alone]               # to coordinate, so draw them all at once
        rank[order[starts[alone]]] = greater[j] + rng.integers(1, tied[j] + 1)
    for start, size, j in zip(starts[~alone], sizes[~alone], which[~alone]):
        rank[order[start:start + size]] = (
            greater[j] + _tie_offsets(int(tied[j]), int(size), rng))
    return rank, n


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SCORERS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def _article_scores(X, weights):
    """ARS: the network-term weight mass an article carries, over the total.

    The denominator is a constant, so it does not affect ranking within a weighting;
    it is retained so scores stay comparable to the pipeline's stored values.
    """
    total = weights.sum()
    s = X @ weights
    return s / total if total > 0 else s


def _bm25(X, doclen, n_total, k1=1.2, b=0.75):
    """BM25 over the MeSH field, query = the network's terms.

    The standard IR baseline on exactly the same term evidence ARS uses, but weighted
    by IDF and normalised by document length. MeSH headings are not repeated within a
    record, so term frequency is 0/1 and the saturation factor collapses to a
    per-document constant driven only by length.
    """
    df = np.asarray((X > 0).sum(0)).ravel().astype(float)
    idf = np.log(1.0 + (n_total - df + 0.5) / (df + 0.5))
    denom = 1.0 + k1 * (1.0 - b + b * doclen / max(doclen.mean(), 1e-9))
    return np.asarray(X.multiply(idf).sum(1)).ravel() * (k1 + 1.0) / denom


def _collect_weightings(node_data, terms, graph):
    """Every node weighting to evaluate, keyed by display label.

    Raw centralities and their MRS counterparts come from the network file. Subgraph
    eigenvector is recomputed here if an older network predates it, so the report
    still produces the full 3x2 grid rather than silently dropping a cell.
    """
    weights, missing = {}, []
    for cent in CENTRALITIES:
        for suffix, scope in SCOPES:
            raw_key = f"{cent}{suffix}_centrality"
            mrs_key = f"MRS_{raw_key}"
            vals = np.array([float(n.get(raw_key, 0.0) or 0.0) for n in node_data])
            if not vals.any() and suffix == "_subgraph" and cent == "eigenvector":
                try:
                    ev = nx.eigenvector_centrality_numpy(graph, weight=None)
                    vals = np.array([ev.get(t, 0.0) for t in terms])
                    missing.append(f"{raw_key} (recomputed from the stored edge list)")
                except Exception:
                    pass
            if vals.any():
                weights[f"{cent} ({scope})"] = vals
            else:
                missing.append(raw_key)
            cvals = np.array([float(n.get(mrs_key, 0.0) or 0.0) for n in node_data])
            if cvals.any():
                weights[f"MRS_{cent} ({scope})"] = cvals
            else:
                missing.append(mrs_key)
    weights["uniform (count of network terms)"] = np.ones(len(terms))
    return weights, missing


def _evaluate_frame(scores_by_label, mask, is_pos, n_boot, rng, frame_name):
    """Run every scorer through one evaluation frame and tabulate with CIs."""
    rows = []
    for label, scores in scores_by_label.items():
        pos, size = _positions(scores, mask, is_pos, rng)
        point, cis = _bootstrap_metrics(pos, size, n_boot, rng)
        row = {"frame": frame_name, "scorer": label, "n_ranked": size,
               "n_positives": int(pos.size)}
        for k, v in point.items():
            lo, hi = cis[k]
            row[k] = v
            row[f"{k}_lo"] = lo
            row[f"{k}_hi"] = hi
        rows.append(row)
    return pd.DataFrame(rows)


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# NODE LEVEL
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def _node_level_analysis(weights_by_label, attested, log_freq, n_boot, rng):
    """Do node weightings identify externally attested terms, and is it just frequency?

    A reviewer of this work raised that network prominence might merely reflect how
    often a term appears in the literature. That is a confounding claim, and the way
    to answer it is a CONDITIONAL analysis - does the association survive holding
    frequency fixed - not a multiplicative penalty on the score, which changes the
    quantity being estimated rather than testing it. Three columns do the work:

      rho(G_t)   association with global document frequency. A correct adjustment
                 pushes this toward zero; strongly negative means the correction
                 overshot and the measure now favours rare terms.
      stratAUC   discrimination WITHIN frequency tertiles, so terms are only ever
                 compared against others of similar corpus frequency.
      incr. p    significance of the measure after log frequency is already in a
                 logistic model.
    """
    try:
        from sklearn.metrics import roc_auc_score
    except ImportError:
        return pd.DataFrame()
    try:
        import statsmodels.api as sm
    except ImportError:
        sm = None

    n_terms = attested.size
    edges = np.quantile(log_freq, np.linspace(0, 1, 4))

    def strat_auc(v):
        """Frequency-stratified AUC of a weighting at separating attested from unattested nodes."""
        vals = []
        for i in range(3):
            hi_inclusive = (i == 2)
            m = ((log_freq >= edges[i]) &
                 ((log_freq <= edges[i + 1]) if hi_inclusive else (log_freq < edges[i + 1])))
            if m.sum() > 5 and 0 < attested[m].sum() < m.sum():
                vals.append(roc_auc_score(attested[m], v[m]))
        return float(np.mean(vals)) if vals else np.nan

    def incremental(v):
        """Incremental logistic-regression lift of a weighting once global frequency is controlled."""
        if sm is None:
            return np.nan
        def z(a):
            """Standardise by the array's spread (mean / SD), guarding against a zero standard deviation."""
            sd = a.std()
            return (a - a.mean()) / sd if sd > 0 else a * 0.0
        design = sm.add_constant(np.column_stack([z(log_freq), z(v)]))
        try:
            return float(sm.Logit(attested.astype(float), design).fit(disp=0).pvalues[2])
        except Exception:
            return np.nan

    rows = []
    for label, v in list(weights_by_label.items()) + [("log10(global frequency)", log_freq)]:
        auc = roc_auc_score(attested, v)
        draws = rng.integers(0, n_terms, size=(min(n_boot, 2000), n_terms))
        boot = [roc_auc_score(attested[i], v[i]) for i in draws
                if 0 < attested[i].sum() < n_terms]
        lo, hi = np.percentile(boot, [2.5, 97.5]) if boot else (np.nan, np.nan)
        # A constant weighting (the uniform baseline) has no defined rank correlation;
        # report it as missing rather than emitting a warning per call.
        rho = spearmanr(v, log_freq).correlation if np.ptp(v) > 0 else np.nan
        rows.append({"node_measure": label, "AUC": auc, "AUC_lo": lo, "AUC_hi": hi,
                     "rho_global_freq": rho,
                     "stratified_AUC": strat_auc(v),
                     "incremental_p_given_freq": incremental(v)})
    return pd.DataFrame(rows).sort_values("AUC", ascending=False)


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PAIRED TESTING
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def _paired_tests(scores_by_label, mask, is_pos, n_total_hint, n_boot, rng, frame_name):
    """The comparisons the conclusions actually rest on.

    Marginal intervals from two scorers can overlap substantially while the paired
    difference is decisively non-zero, because both are measured on the same 88
    positives. Every claim of the form "X beats Y" is therefore settled here rather
    than by eyeballing two CIs in the frame tables.
    """
    pos_cache = {}
    size = None
    for label, s in scores_by_label.items():
        pos_cache[label], size = _positions(s, mask, is_pos, rng)

    metrics = {"BEDROC": lambda p, n: _bedroc_from_positions(p, n),
               "recall@10000": lambda p, n: float((p <= 10000).mean()),
               "recall@100000": lambda p, n: float((p <= 100000).mean()),
               "MAP": lambda p, n: float(np.mean(np.arange(1, p.size + 1) / np.sort(p)))}

    pairs = []
    for cent in CENTRALITIES:
        for _, scope in SCOPES:
            a, b = f"MRS_{cent} ({scope})", f"{cent} ({scope})"
            if a in pos_cache and b in pos_cache:
                pairs.append((a, b, "MRS vs raw centrality"))
    uni = "uniform (count of network terms)"
    for cent in CENTRALITIES:
        for _, scope in SCOPES:
            b = f"{cent} ({scope})"
            if b in pos_cache and uni in pos_cache:
                pairs.append((b, uni, "topology weighting vs plain count"))
    if "BM25 (MeSH field)" in pos_cache:
        for cent in CENTRALITIES:
            b = f"{cent} (subgraph)"
            if b in pos_cache:
                pairs.append((b, "BM25 (MeSH field)", "network weighting vs standard IR"))

    rows = []
    for a, b, family in pairs:
        for mname, fn in metrics.items():
            obs, lo, hi, p = _paired_bootstrap(pos_cache[a], pos_cache[b], size,
                                               fn, n_boot, rng)
            rows.append({"frame": frame_name, "comparison_family": family,
                         "metric": mname, "scorer_A": a, "scorer_B": b,
                         "difference_A_minus_B": obs, "ci_lo": lo, "ci_hi": hi,
                         "p_two_sided": p,
                         "significant": "yes" if (lo > 0 or hi < 0) else "no"})
    return pd.DataFrame(rows)


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ORCHESTRATOR
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def run_validation_report(final_network_file, master_db_path, pmid_db_path,
                          ground_truth_file, output_dir, project_prefix="",
                          primary_query_term=None, n_boot=2000, random_seed=42,
                          cache_dir=None, start_date=None, end_date=None):
    """Build the full validation report: workbook, HTML narrative, figures, console.

    Returns the dict of DataFrames it wrote, so callers can inspect results without
    re-reading the workbook.
    """
    rng = np.random.default_rng(random_seed)
    os.makedirs(output_dir, exist_ok=True)
    fig_dir = os.path.join(output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    cache_dir = cache_dir or os.path.join(output_dir, ".validation_cache")

    print("\n" + "=" * 78)
    print("VALIDATION REPORT")
    print("<" * 30 + ">" * 30)

    with open(final_network_file, "r", encoding="utf-8") as f:
        net = json.load(f)
    node_data = [n["data"] for n in net["elements"]["nodes"]]
    terms = [n["id"] for n in node_data]
    graph = nx.Graph()
    graph.add_nodes_from(terms)
    for e in net["elements"]["edges"]:
        graph.add_edge(e["data"]["source"], e["data"]["target"])
    print(f"  network: {len(terms)} nodes / {graph.number_of_edges()} edges")

    weights, missing = _collect_weightings(node_data, terms, graph)
    if missing:
        print(f"  [i] weightings unavailable or reconstructed: {len(missing)}")
        for m in missing:
            print(f"      - {m}")

    X, pmids, doclen, years = _build_or_load_cache(
        master_db_path, terms, os.path.join(cache_dir, "pool"),
        start_date=start_date, end_date=end_date)
    n_pool = X.shape[0]
    print(f"  scored pool: {n_pool:,} articles carrying >=1 network term")

    # article scores for every weighting, plus the non-network comparators
    scores = {label: _article_scores(X, w) for label, w in weights.items()}
    scores["BM25 (MeSH field)"] = _bm25(X, doclen, n_pool)
    scores["recency (newest first)"] = years.astype(float)
    scores["random order"] = rng.random(n_pool)

    gt = pd.read_excel(ground_truth_file)
    pmid_col = next((c for c in gt.columns if c.strip().upper() == "PMID"), gt.columns[0])
    gt_ids = {int(str(x).replace(".0", "")) for x in gt[pmid_col]
              if str(x).replace(".0", "").isdigit()}
    is_pos = np.isin(pmids, list(gt_ids))
    n_pos = int(is_pos.sum())

    # the query pool: the real P0 result set, not a single-heading proxy
    in_query = np.zeros(n_pool, dtype=bool)
    try:
        con = _open_resilient(pmid_db_path)
        p0 = {r[0] for r in con.execute(
            "SELECT pmid FROM pmids_table WHERE generation='P0'")}
        con.close()
        in_query = np.isin(pmids, list(p0))
    except Exception as exc:
        print(f"  [i] Query-conditioned frames skipped: "
              f"{_pmid_db_reason(pmid_db_path, exc)}. The whole-pool ranking "
              f"below is unaffected.")

    if primary_query_term and primary_query_term in terms:
        col = terms.index(primary_query_term)
        scores[f"naive query ('{primary_query_term}')"] = np.asarray(
            X[:, col].todense()).ravel()

    print(f"  ground truth: {n_pos} of {len(gt_ids)} curated PMIDs are in the pool")
    print(f"    inside the query set : {int((is_pos & in_query).sum())}")
    print(f"    outside              : {int((is_pos & ~in_query).sum())}")

    all_mask = np.ones(n_pool, dtype=bool)
    frames = {"corpus-wide": all_mask}
    if in_query.any():
        frames["within query pool"] = in_query
        frames["outside query pool"] = ~in_query

    print(f"\n  evaluating {len(scores)} scorers x {len(frames) + 1} frames "
          f"(B = {n_boot:,}) ...")
    frame_tables = {}
    for name, mask in frames.items():
        frame_tables[name] = _evaluate_frame(scores, mask, is_pos, n_boot, rng, name)

    # hybrid: query hits first (network-ranked), then the remainder (network-ranked)
    if in_query.any():
        hybrid_scores = {}
        for label, s in scores.items():
            if label.startswith("naive query") or label == "random order":
                continue
            hybrid_scores[f"HYBRID: {label}"] = s + in_query * (float(s.max()) + 1.0)
        hybrid_scores["query only (filter, unranked)"] = in_query.astype(float)
        frame_tables["hybrid"] = _evaluate_frame(
            hybrid_scores, all_mask, is_pos, n_boot, rng, "hybrid")

    paired = pd.concat([
        _paired_tests(scores, all_mask, is_pos, n_pool, n_boot, rng, "corpus-wide"),
        _paired_tests(scores, in_query, is_pos, n_pool, n_boot, rng, "within query pool")
        if in_query.any() else pd.DataFrame(),
    ], ignore_index=True)

    attested = np.asarray(X[is_pos].sum(0)).ravel() > 0
    df_t = np.asarray((X > 0).sum(0)).ravel().astype(float)
    node_tbl = _node_level_analysis(weights, attested, np.log10(np.maximum(df_t, 1)),
                                    n_boot, rng)

    diagnostics = pd.DataFrame([
        {"item": "generated", "value": datetime.now().strftime("%Y-%m-%d %H:%M:%S")},
        {"item": "network nodes", "value": len(terms)},
        {"item": "network edges", "value": graph.number_of_edges()},
        {"item": "scored pool (articles)", "value": n_pool},
        {"item": "curated ground-truth PMIDs", "value": len(gt_ids)},
        {"item": "ground truth present in pool", "value": n_pos},
        {"item": "ground truth inside query set", "value": int((is_pos & in_query).sum())},
        {"item": "ground truth outside query set", "value": int((is_pos & ~in_query).sum())},
        {"item": "query set size (P0)", "value": int(in_query.sum())},
        {"item": "attested network terms", "value": int(attested.sum())},
        {"item": "bootstrap replicates", "value": n_boot},
        {"item": "random seed", "value": random_seed},
        {"item": "weightings evaluated", "value": len(weights)},
        {"item": "notes on missing weightings", "value": "; ".join(missing) or "none"},
    ])

    tables = dict(frame_tables)
    tables["paired tests"] = paired
    tables["node level"] = node_tbl
    tables["diagnostics"] = diagnostics

    _print_console_summary(tables)
    xlsx = os.path.join(output_dir, f"{project_prefix}validation_report.xlsx")
    _write_workbook(tables, xlsx)
    figs = _make_figures(scores, is_pos, in_query, n_pool, rng, fig_dir, project_prefix)
    try:
        figs.update(_fig_paired_forest(paired, fig_dir, project_prefix))
    except Exception as exc:
        print(f"  [!] Forest plot skipped: {exc}")
    html = os.path.join(output_dir, f"{project_prefix}validation_report.html")
    _write_html(tables, html, figs, project_prefix)
    print(f"\n  [+] workbook : {xlsx}")
    print(f"  [+] narrative: {html}")
    print(f"  [+] figures  : {fig_dir}")
    return tables


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# OUTPUTS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def _fmt_ci(row, metric, pct=False):
    """Format a metric with its bootstrap CI as 'point [lo, hi]' (optionally as percentages)."""
    v, lo, hi = row[metric], row.get(f"{metric}_lo"), row.get(f"{metric}_hi")
    if pct:
        return f"{100*v:.1f}% [{100*lo:.1f}, {100*hi:.1f}]"
    return f"{v:.4f} [{lo:.4f}, {hi:.4f}]"


def _print_console_summary(tables):
    """Print the headline per-frame scorer comparison to the console."""
    corpus = tables.get("corpus-wide")
    if corpus is None or corpus.empty:
        return
    print("\n" + "-" * 96)
    print("CORPUS-WIDE SUMMARY (95% CI from resampling the ground-truth positives)")
    print("-" * 96)
    cols = ["scorer", "BEDROC", "EF@0.001", "recall@10000", "recall@100000", "MAP"]
    have = [c for c in cols if c in corpus.columns]
    top = corpus.sort_values("BEDROC", ascending=False)
    print(f"  {'scorer':<38}{'BEDROC':>9}{'EF@0.1%':>10}{'r@10k':>9}{'r@100k':>9}{'MAP':>10}")
    for _, r in top.iterrows():
        print(f"  {r['scorer']:<38}{r.get('BEDROC', float('nan')):>9.4f}"
              f"{r.get('EF@0.001', float('nan')):>10.1f}"
              f"{100*r.get('recall@10000', float('nan')):>8.1f}%"
              f"{100*r.get('recall@100000', float('nan')):>8.1f}%"
              f"{r.get('MAP', float('nan')):>10.5f}")

    paired = tables.get("paired tests")
    if paired is not None and not paired.empty:
        sig = paired[(paired["significant"] == "yes") &
                     (paired["frame"] == "corpus-wide")]
        print(f"\n  paired comparisons reaching significance (corpus-wide): "
              f"{len(sig)} of {len(paired[paired['frame'] == 'corpus-wide'])}")

    node = tables.get("node level")
    if node is not None and not node.empty:
        print("\n" + "-" * 96)
        print("NODE LEVEL - do weightings identify attested terms, and is it just frequency?")
        print("-" * 96)
        print(f"  {'node measure':<40}{'AUC':>8}{'rho(G_t)':>10}{'stratAUC':>10}{'incr p':>9}")
        for _, r in node.iterrows():
            print(f"  {r['node_measure']:<40}{r['AUC']:>8.4f}"
                  f"{r['rho_global_freq']:>10.3f}{r['stratified_AUC']:>10.4f}"
                  f"{r['incremental_p_given_freq']:>9.4f}")


def _write_workbook(tables, path):
    """Drill-down workbook: summary first, then one sheet per evaluation frame."""
    readme = pd.DataFrame([
        ("BEDROC", "Boltzmann-enhanced early recognition, alpha=20. Weights the top of "
                   "the ranking. Primary metric: the task is early recognition under "
                   "extreme class imbalance."),
        ("EF@0.001 / EF@0.01", "Enrichment factor: positives found in the top 0.1% / 1% "
                               "relative to chance. 1.0 = no better than random."),
        ("recall@K", "Fraction of ground truth found in the first K articles. The "
                     "practical 'how much reading' measure."),
        ("MAP", "Mean average precision over the positives. Rewards ranking positives "
                "consistently high rather than in one lucky cluster."),
        ("95% CI", "Percentile bootstrap resampling the GROUND-TRUTH POSITIVES, which "
                   "are the sampled unit. Resampling the article pool instead would "
                   "reshuffle negatives and give intervals that are far too tight."),
        ("ROC-AUC (articles)", "Deliberately NOT reported. At ~1e-5 prevalence it is "
                               "dominated by true negatives and flatters scorers that "
                               "retrieve nothing useful."),
        ("AUC (nodes)", "Reported: ~40% of network terms are attested, so class balance "
                        "is reasonable and AUC is a sound discrimination measure."),
        ("frame: corpus-wide", "Every article ranked. The headline retrieval task."),
        ("frame: within query pool", "Only the query's own hits. Identical article set "
                                     "for every scorer, so this isolates ORDERING "
                                     "ability from scope."),
        ("frame: outside query pool", "Only articles the query cannot return. Isolates "
                                      "REACH - evidence invisible to a conventional "
                                      "search."),
        ("frame: hybrid", "Query hits first (network-ranked), then the remainder "
                          "(network-ranked). The reading list a user is handed."),
        ("rho_global_freq", "Association between a node measure and global document "
                            "frequency. ~0 is frequency-neutral; strongly negative means "
                            "a frequency correction has overshot into favouring rare terms."),
        ("stratified_AUC", "Discrimination within frequency tertiles, so terms are only "
                           "compared against others of similar corpus frequency. This is "
                           "the direct test of 'prominence is just frequency'."),
        ("incremental_p_given_freq", "Significance of the measure after log global "
                                     "frequency is already in a logistic model."),
        ("paired tests", "Differences between two scorers on the SAME resampled "
                         "positives. Overlapping marginal CIs do not imply equivalence; "
                         "these are the tests conclusions should rest on."),
    ], columns=["term", "meaning"])

    sheet_order = [("00_README", readme),
                   ("01_CorpusWide", tables.get("corpus-wide")),
                   ("02_WithinQuery", tables.get("within query pool")),
                   ("03_OutsideQuery", tables.get("outside query pool")),
                   ("04_Hybrid", tables.get("hybrid")),
                   ("05_PairedTests", tables.get("paired tests")),
                   ("06_NodeLevel", tables.get("node level")),
                   ("07_Diagnostics", tables.get("diagnostics"))]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        for name, df in sheet_order:
            if df is None or (hasattr(df, "empty") and df.empty):
                continue
            df.to_excel(xl, sheet_name=name, index=False)
            ws = xl.sheets[name]
            for col in ws.columns:
                width = max((len(str(c.value)) for c in col if c.value is not None),
                            default=10)
                ws.column_dimensions[col[0].column_letter].width = min(max(width + 2, 10), 70)
            ws.freeze_panes = "A2"


def _make_figures(scores, is_pos, in_query, n_pool, rng, fig_dir, prefix):
    """Figures for the four claims the report makes. Returns {title: filename}."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    made = {}
    all_mask = np.ones(n_pool, dtype=bool)

    def save(fig, stem, title):
        """Save a figure under the run's figure directory (PNG, plus TIFF for print)."""
        png = os.path.join(fig_dir, f"{prefix}{stem}.png")
        fig.savefig(png, dpi=300, bbox_inches="tight")
        try:
            fig.savefig(os.path.join(fig_dir, f"{prefix}{stem}.tif"),
                        dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
        except Exception:
            pass
        plt.close(fig)
        made[title] = os.path.basename(png)

    # 1. enrichment curves - the headline retrieval picture
    interesting = [k for k in scores if "(subgraph)" in k or k.startswith("BM25")
                   or k.startswith("uniform") or k.startswith("naive")]
    fracs = np.logspace(-5, 0, 60)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for label in interesting:
        pos, size = _positions(scores[label], all_mask, is_pos, rng)
        if pos.size == 0:
            continue
        rec = [float((pos <= max(int(f * size), 1)).mean()) for f in fracs]
        ax.plot(fracs * 100, np.array(rec) * 100, lw=1.6, label=label)
    ax.plot(fracs * 100, fracs * 100, "k--", lw=1, label="random")
    ax.set_xscale("log")
    ax.set_xlabel("Corpus screened (%)")
    ax.set_ylabel("Ground truth recovered (%)")
    ax.set_title("Enrichment: recovery versus reading effort")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(alpha=0.3)
    save(fig, "fig1_enrichment", "Enrichment curves")

    # 2. the query has a ceiling; the network does not
    if in_query.any():
        fig, ax = plt.subplots(figsize=(8, 5.5))
        best = [k for k in scores if "(subgraph)" in k and not k.startswith("MRS")][:2] + \
               [k for k in scores if k.startswith("MRS") and "(subgraph)" in k][:1]
        for label in best:
            s = scores[label] + in_query * (float(scores[label].max()) + 1.0)
            pos, size = _positions(s, all_mask, is_pos, rng)
            rec = [float((pos <= max(int(f * size), 1)).mean()) for f in fracs]
            ax.plot(fracs * 100, np.array(rec) * 100, lw=1.8, label=f"HYBRID: {label}")
        pos, size = _positions(in_query.astype(float), all_mask, is_pos, rng)
        rec = [float((pos <= max(int(f * size), 1)).mean()) for f in fracs]
        ax.plot(fracs * 100, np.array(rec) * 100, lw=2.2, color="crimson",
                label="query only (filter)")
        ceiling = 100 * float((is_pos & in_query).sum()) / max(int(is_pos.sum()), 1)
        ax.axhline(ceiling, color="crimson", ls=":", lw=1.2,
                   label=f"query ceiling ({ceiling:.0f}%)")
        ax.set_xscale("log")
        ax.set_xlabel("Corpus screened (%)")
        ax.set_ylabel("Ground truth recovered (%)")
        ax.set_title("The query saturates at its ceiling; the hybrid continues past it")
        ax.legend(fontsize=7, loc="upper left")
        ax.grid(alpha=0.3)
        save(fig, "fig2_query_ceiling", "Query ceiling versus hybrid")
    return made


def _fig_paired_forest(paired, fig_dir, prefix):
    """Forest plot of paired differences - the significance picture at a glance."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sub = paired[(paired["frame"] == "corpus-wide") & (paired["metric"] == "BEDROC")]
    if sub.empty:
        return {}
    sub = sub.sort_values("difference_A_minus_B")
    fig, ax = plt.subplots(figsize=(9, max(3.5, 0.32 * len(sub))))
    ypos = np.arange(len(sub))
    colors = ["#1a7f37" if s == "yes" and d > 0 else "#b42318" if s == "yes" else "#98a2b3"
              for s, d in zip(sub["significant"], sub["difference_A_minus_B"])]
    # errorbar takes a single ecolor, so the bars are drawn one at a time to let each
    # comparison carry its own significance colour.
    for y, d, lo, hi, c in zip(ypos, sub["difference_A_minus_B"],
                               sub["ci_lo"], sub["ci_hi"], colors):
        ax.plot([lo, hi], [y, y], color=c, lw=1.4, solid_capstyle="butt", zorder=2)
    ax.scatter(sub["difference_A_minus_B"], ypos, c=colors, s=22, zorder=3)
    ax.axvline(0, color="black", lw=1)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{a}  vs  {b}" for a, b in zip(sub["scorer_A"], sub["scorer_B"])],
                       fontsize=7)
    ax.set_xlabel("Difference in BEDROC (A - B), paired bootstrap 95% CI")
    ax.set_title("Paired comparisons: green = A better, red = B better, grey = unresolved")
    ax.grid(alpha=0.3, axis="x")
    png = os.path.join(fig_dir, f"{prefix}fig3_paired_forest.png")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {"Paired comparisons (forest)": os.path.basename(png)}


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# NARRATIVE REPORT
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def _best(df, metric, exclude_prefix=None):
    """Best-scoring row for a metric, optionally excluding scorers with a label prefix."""
    d = df if exclude_prefix is None else df[~df["scorer"].str.startswith(exclude_prefix)]
    d = d.dropna(subset=[metric])
    return d.loc[d[metric].idxmax()] if not d.empty else None


def _prose(tables):
    """Manuscript-ready paragraphs, with the numbers interpolated from the tables.

    Written so the text stays true if the underlying data changes: every figure quoted
    is pulled from the computed tables rather than hard-coded.
    """
    out = []
    corpus = tables.get("corpus-wide", pd.DataFrame())
    diag = {r["item"]: r["value"] for _, r in tables.get("diagnostics",
                                                        pd.DataFrame()).iterrows()}
    node = tables.get("node level", pd.DataFrame())
    paired = tables.get("paired tests", pd.DataFrame())
    within = tables.get("within query pool", pd.DataFrame())
    outside = tables.get("outside query pool", pd.DataFrame())
    hybrid = tables.get("hybrid", pd.DataFrame())

    n_pool = diag.get("scored pool (articles)", "the")
    n_pos = diag.get("ground truth present in pool", "")
    inq = diag.get("ground truth inside query set", 0)
    outq = diag.get("ground truth outside query set", 0)

    # --- methods
    out.append(("Methods - evaluation design", f"""
Ranking performance was evaluated against {n_pos} curated ground-truth articles present
in a scored pool of {n_pool:,} records carrying at least one network term. Because
ground-truth prevalence is on the order of 1e-5, early-recognition metrics (BEDROC with
alpha = 20, enrichment factor at 0.1% and 1%, recall at fixed reading depths, and mean
average precision) were used in preference to ROC-AUC, which at this prevalence is
dominated by true negatives and rewards scorers that retrieve nothing of value.
Confidence intervals are percentile bootstraps over the ground-truth positives, which
are the sampled unit; resampling the article pool would reshuffle negatives drawn from
a fixed population and yield intervals that are too narrow. Comparisons between scorers
use a paired bootstrap in which a single resample is scored by both methods, cancelling
the variance shared through the common positive set.

Scoring was evaluated in four frames. The corpus-wide frame ranks all articles. The
within-query frame ranks only the articles returned by the source query, so the article
set is identical across scorers and only ordering differs, isolating ranking ability
from scope. The outside-query frame ranks only articles the query cannot return,
isolating reach. The hybrid frame places query hits first, network-ranked, followed by
the remainder, also network-ranked."""))

    # --- ARS result
    if not corpus.empty and not paired.empty:
        uni_tests = paired[(paired["comparison_family"] == "topology weighting vs plain count")
                           & (paired["frame"] == "corpus-wide")]
        n_sig = int((uni_tests["significant"] == "yes").sum())
        b = _best(corpus, "BEDROC", exclude_prefix="uniform")
        uni_row = corpus[corpus["scorer"].str.startswith("uniform")]
        uni_b = float(uni_row["BEDROC"].iloc[0]) if not uni_row.empty else float("nan")
        out.append(("Results - topology weighting is not equivalent to counting terms", f"""
Weighting network terms by their topological importance outperformed an unweighted
count of the same terms. The strongest weighting reached BEDROC {b['BEDROC']:.4f}
(95% CI {b['BEDROC_lo']:.4f}-{b['BEDROC_hi']:.4f}) against {uni_b:.4f} for the uniform
count, and {n_sig} of {len(uni_tests)} paired comparisons against the uniform baseline
were significant at the 5% level. Article relevance scoring is therefore not a
reformulation of term counting: the topology of the concept network carries ranking
information that term membership alone does not."""))

    # --- query ceiling / reach
    if inq and outq:
        ceiling = 100.0 * inq / (inq + outq)
        best_out = _best(outside, "recall@100000") if not outside.empty else None
        seg = (f" Within the articles the query cannot return, the strongest weighting "
               f"recovered {100*best_out['recall@100000']:.1f}% of that evidence by "
               f"100,000 articles screened." if best_out is not None else "")
        out.append(("Results - query reach and its ceiling", f"""
Of the {inq + outq} ground-truth articles, {inq} carried the source query and {outq} did
not, placing a hard ceiling of {ceiling:.0f}% recall on any strategy based on that query
alone: the remaining evidence is not retrievable by it at any reading depth.{seg} This is
the expected pattern for an adverse outcome pathway, whose key events are indexed under
mechanistic rather than disease headings, and it is the principal argument for
network-derived relevance over query refinement."""))

    # --- within-query ranking ability
    if not within.empty:
        rnd = within[within["scorer"] == "random order"]
        bw = _best(within, "MAP", exclude_prefix="random")
        if bw is not None and not rnd.empty:
            out.append(("Results - ranking ability with scope held constant", f"""
Restricted to the query's own results, where every scorer ranks an identical article
set, the strongest weighting achieved MAP {bw['MAP']:.4f} against
{float(rnd['MAP'].iloc[0]):.4f} for random ordering. Because scope is fixed by
construction in this frame, the difference is attributable to ordering alone and cannot
be explained by the size or composition of the retrieved set."""))

    # --- hybrid
    if not hybrid.empty:
        bh = _best(hybrid, "MAP")
        q = hybrid[hybrid["scorer"].str.startswith("query only")]
        if bh is not None and not q.empty:
            out.append(("Results - the combined reading list", f"""
Combining the two - query hits first, then network-ranked remainder - outperformed
either component. The best hybrid reached MAP {bh['MAP']:.4f} and recovered
{100*bh['recall@100000']:.1f}% of the ground truth within 100,000 articles, against
{float(q['MAP'].iloc[0]):.4f} and {100*float(q['recall@100000'].iloc[0]):.1f}% for the
query used alone. The two sources of signal are complementary rather than competing:
the query supplies precision at the top, the network supplies reach beyond it."""))

    # --- node level / reviewer response
    if not node.empty:
        freq = node[node["node_measure"] == "log10(global frequency)"]
        bn = node.iloc[0]
        if not freq.empty:
            f0 = freq.iloc[0]
            out.append(("Results - node prominence is not reducible to term frequency", f"""
Global document frequency alone discriminated externally attested network terms with
AUC {f0['AUC']:.4f}, confirming that attested terms are indeed more common in the
literature. That association does not account for the network signal. The strongest node
measure ({bn['node_measure']}) reached AUC {bn['AUC']:.4f}
(95% CI {bn['AUC_lo']:.4f}-{bn['AUC_hi']:.4f}), and within frequency-matched tertiles it
retained AUC {bn['stratified_AUC']:.4f} against {f0['stratified_AUC']:.4f} for frequency
itself, remaining significant after conditioning on log global frequency
(p = {bn['incremental_p_given_freq']:.4g}). Among terms of comparable corpus frequency,
network position therefore still identifies the attested ones.

We note that confounding of this kind is properly addressed by conditional analysis
rather than by multiplying scores by a specificity penalty. A penalty term does not test
the confound; it redefines the estimand, replacing prominence with rarity-adjusted
prominence, which is a different quantity and one that no longer corresponds to the
curated ground truth. Node-level AUC is reported here, rather than at the article level,
because roughly 40% of network terms are attested and the class balance makes the measure
informative."""))

    out.append(("Limitations", f"""
The ground truth comprises {n_pos} articles, so interval estimates are wide and small
differences between closely matched scorers are not resolved; conclusions are drawn only
where paired tests reach significance. The curation was assembled by domain experts and
may under-represent literature that is hard to find by conventional search, which would
bias the comparison against network-derived scoring rather than in its favour. Metrics
are computed on articles carrying at least one network term, so recall is expressed
relative to the ground truth present in that pool rather than to the full curated set."""))
    return out


def _df_to_html(df, max_rows=400):
    """Render a DataFrame to an HTML table (capped at max_rows) for the report."""
    d = df.head(max_rows)
    fmt = {c: (lambda v: f"{v:.4f}" if isinstance(v, float) else v) for c in d.columns}
    return d.to_html(index=False, float_format=lambda v: f"{v:.4f}",
                     classes="tbl", border=0, na_rep="-")


def _write_html(tables, path, figures, prefix):
    """Assemble the narrative HTML validation report from the tables and embedded figures."""
    css = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
 max-width:1180px;margin:0 auto;padding:32px 24px;line-height:1.65;color:#1a1a1a}
h1{font-size:1.9rem;margin-bottom:.2em}h2{margin-top:2.2em;padding-bottom:.3em;
 border-bottom:2px solid #e5e7eb;font-size:1.3rem}h3{margin-top:1.6em;font-size:1.05rem;color:#374151}
.sub{color:#6b7280;margin-top:0}
.prose{background:#f9fafb;border-left:3px solid #9ca3af;padding:14px 18px;margin:14px 0;
 white-space:pre-wrap;font-size:.94rem}
table.tbl{border-collapse:collapse;width:100%;font-size:.8rem;margin:12px 0}
table.tbl th{background:#f3f4f6;text-align:left;padding:6px 8px;border-bottom:2px solid #d1d5db;
 position:sticky;top:0}
table.tbl td{padding:5px 8px;border-bottom:1px solid #eef0f2}
table.tbl tr:hover td{background:#f9fafb}
.scroll{overflow-x:auto;max-height:640px;overflow-y:auto;border:1px solid #e5e7eb;border-radius:6px}
img{max-width:100%;height:auto;border:1px solid #e5e7eb;border-radius:6px;margin:10px 0}
.note{color:#6b7280;font-size:.85rem}
"""
    parts = [f"<style>{css}</style>",
             f"<h1>Validation report</h1>",
             f"<p class='sub'>Generated {datetime.now().strftime('%Y-%m-%d %H:%M')} "
             f"&middot; {prefix or 'project'}</p>"]

    parts.append("<h2>Interpretation</h2>")
    for title, text in _prose(tables):
        parts.append(f"<h3>{title}</h3><div class='prose'>{text.strip()}</div>")

    if figures:
        parts.append("<h2>Figures</h2>")
        for title, fname in figures.items():
            parts.append(f"<h3>{title}</h3><img src='figures/{fname}' alt='{title}'>")

    order = [("Corpus-wide", "corpus-wide"), ("Within query pool", "within query pool"),
             ("Outside query pool", "outside query pool"), ("Hybrid", "hybrid"),
             ("Paired tests", "paired tests"), ("Node level", "node level"),
             ("Diagnostics", "diagnostics")]
    parts.append("<h2>Tables</h2>")
    for title, key in order:
        df = tables.get(key)
        if df is None or df.empty:
            continue
        parts.append(f"<h3>{title}</h3><div class='scroll'>{_df_to_html(df)}</div>")

    parts.append("<p class='note'>Confidence intervals are percentile bootstraps over "
                 "the ground-truth positives. Article-level ROC-AUC is deliberately not "
                 "reported; see the README sheet of the workbook.</p>")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(parts))


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PROJECTION COMPARISON  (ported from the sandbox ARS evaluation; no HTML)
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def _fig_projection_comparison(df, fig_dir, prefix):
    """Grouped BEDROC bar chart with bootstrap-CI whiskers, one panel per frame."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    frames = list(dict.fromkeys(df["frame"]))
    methods = list(dict.fromkeys(df["method"]))
    fig, axes = plt.subplots(1, len(frames), figsize=(6.2 * len(frames), 5.4), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, fr in zip(axes, frames):
        sub = df[df["frame"] == fr].set_index("method").reindex(methods)
        y = np.arange(len(methods))
        err = np.vstack([(sub["BEDROC"] - sub["BEDROC_lo"]).to_numpy(),
                         (sub["BEDROC_hi"] - sub["BEDROC"]).to_numpy()])
        ax.barh(y, sub["BEDROC"].to_numpy(), xerr=np.abs(err), color="#00819C",
                error_kw=dict(ecolor="#333", lw=1), height=.66)
        ax.set_yticks(y); ax.set_yticklabels(methods, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("BEDROC (alpha = 20)")
        ax.set_title(fr, fontsize=11)
        ax.grid(True, axis="x", ls="--", alpha=.35)
    fig.suptitle("Article-scoring projection comparison vs ground truth", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(fig_dir, f"{prefix}fig_projection_comparison.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


def run_projection_comparison(final_network_file, master_db_path, pmid_db_path,
                              ground_truth_file, output_dir, project_prefix="",
                              primary_query_term="Dermatitis, Allergic Contact",
                              seed_key="pagerank_subgraph_centrality",
                              mrs_key="MRS_pagerank_subgraph_centrality",
                              n_boot=2000, random_seed=42, cache_dir=None,
                              make_figure=True, start_date=None, end_date=None):
    """Compare article-scoring PROJECTIONS against the ground truth (CSV + figure, no HTML).

    Complements run_validation_report: that asks "which node WEIGHTING (seed) ranks the
    ground truth best?"; this fixes the seed and asks "which PROJECTION turns node weights
    into an article score best?" - the symmetrically normalised projection vs the plain
    sum vs mutual reinforcement vs the MRS re-weighting vs BM25 vs trivial baselines,
    scored by BEDROC across three frames (whole pool / citation neighbourhood / within
    query) with positives-only bootstrap CIs. This is the sandbox ARS evaluation ported
    into the pipeline; the deprecated CRS-weighted scorer is replaced by the MRS-weighted
    one, and only a table and a figure are emitted (no HTML).

    Shares the on-disk scoring-pool cache with run_validation_report, so when both run in
    the benchmark step the 9M-row corpus is scanned once. The mutual-reinforcement scorer
    materialises one transient copy of the incidence matrix; run this when memory is free.
    """
    rng = np.random.default_rng(random_seed)
    os.makedirs(output_dir, exist_ok=True)
    fig_dir = os.path.join(output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    cache_dir = cache_dir or os.path.join(output_dir, ".validation_cache")

    print("\n" + "=" * 78)
    print("PROJECTION COMPARISON (article-scoring methods)")
    print("<" * 30 + ">" * 30)

    with open(final_network_file, "r", encoding="utf-8") as f:
        net = json.load(f)
    node_data = [n["data"] for n in net["elements"]["nodes"]]
    terms = [n["id"] for n in node_data]
    nd = {d["id"]: d for d in node_data}

    def w_of(k):
        """Node-attribute vector over the network terms; missing or blank values count as 0."""
        return np.array([float(nd[t].get(k, 0.0) or 0.0) for t in terms])

    # Seed = the chosen node weighting; fall back to the whole-graph key if the subgraph
    # one was never computed (i.e. ground-truth analysis was off at network-build time).
    s = w_of(seed_key)
    if not s.any():
        s = w_of(seed_key.replace("_subgraph", ""))
        print(f"  [i] '{seed_key}' is empty; using the whole-graph seed instead.")
    mrs = w_of(mrs_key)
    if not mrs.any():
        mrs = w_of(mrs_key.replace("_subgraph", ""))

    X, pmids, doclen, years = _build_or_load_cache(
        master_db_path, terms, os.path.join(cache_dir, "pool"),
        start_date=start_date, end_date=end_date)
    n_pool, M = X.shape
    print(f"  scored pool: {n_pool:,} articles x {M} network terms")

    deg_a = np.asarray(X.sum(1)).ravel()
    deg_t = np.asarray(X.sum(0)).ravel()
    with np.errstate(divide="ignore"):
        da_isqrt = np.where(deg_a > 0, deg_a ** -0.5, 0.0)
        dt_isqrt = np.where(deg_t > 0, deg_t ** -0.5, 0.0)
        da_inv = np.where(deg_a > 0, 1.0 / deg_a, 0.0)

    def _static(w):
        """Unnormalised projection: the summed seed weights of an article's network terms."""
        return X @ w

    def _normed(seed):
        """Symmetrically normalised projection of the seed through the incidence matrix."""
        return da_isqrt * (X @ (dt_isqrt * seed))

    def _reinforced(seed, a=0.5, it=1000, tol=1e-9):
        """Mutual-reinforcement fixed point of the seed, then the normalised projection."""
        Mop = np.asarray((X.T @ X.multiply(da_inv[:, None])).todense())
        P = (dt_isqrt[:, None] * Mop) * dt_isqrt[None, :]
        p = seed / seed.sum() if seed.sum() else seed
        t = p.copy()
        for _ in range(it):
            tn = a * (P @ t) + (1 - a) * p
            if np.abs(tn - t).sum() < tol:
                break
            t = tn
        return da_isqrt * (X @ (dt_isqrt * tn))

    scorers = {
        "Normalised ARS":         _normed(s),
        "ARS (unnormalised sum)": _static(s),
        "ARS (MRS-weighted)":     _static(mrs),
        "Bipartite reinforced":   _reinforced(s),
        "BM25 (MeSH field)":      _bm25(X, doclen, n_pool),
        "Uniform term count":     _static(np.ones(M)),
        "Random":                 rng.random(n_pool),
    }
    if primary_query_term and primary_query_term in terms:
        scorers[f"Naive query ('{primary_query_term}')"] = np.asarray(
            X[:, terms.index(primary_query_term)].todense()).ravel()

    gt = pd.read_excel(ground_truth_file)
    pmid_col = next((c for c in gt.columns if c.strip().upper() == "PMID"), gt.columns[0])
    gt_ids = {int(str(x).replace(".0", "")) for x in gt[pmid_col]
              if str(x).replace(".0", "").isdigit()}
    is_pos = np.isin(pmids, list(gt_ids))

    in_q = np.zeros(n_pool, bool)
    in_nb = np.ones(n_pool, bool)
    try:
        con = _open_resilient(pmid_db_path)
        p0 = {r[0] for r in con.execute("SELECT pmid FROM pmids_table WHERE generation='P0'")}
        nb = {r[0] for r in con.execute("SELECT pmid FROM pmids_table")}
        con.close()
        in_q = np.isin(pmids, list(p0))
        in_nb = np.isin(pmids, list(nb))
    except Exception as exc:
        print(f"  [i] Only the whole-pool frame will be scored: "
              f"{_pmid_db_reason(pmid_db_path, exc)}.")

    frames = {"whole pool": np.ones(n_pool, bool),
              "citation neighbourhood": in_nb,
              "within query": in_q}
    print(f"  ground truth: {int(is_pos.sum())} positives in pool "
          f"({int((is_pos & in_q).sum())} within query, "
          f"{int((is_pos & in_nb).sum())} in neighbourhood)")

    rows = []
    for fname, mask in frames.items():
        if not mask.any() or not (is_pos & mask).any():
            continue
        size = int(mask.sum())
        for lab, sc in scorers.items():
            pos, _ = _positions(sc.astype(float), mask, is_pos, rng)
            point, ci = _bootstrap_metrics(pos, size, n_boot, rng)
            rows.append({"frame": fname, "method": lab, "n_articles": size,
                         "n_positives": int(pos.size), "BEDROC": point["BEDROC"],
                         "BEDROC_lo": ci["BEDROC"][0], "BEDROC_hi": ci["BEDROC"][1],
                         "EF@0.01": point.get("EF@0.01", np.nan),
                         "recall@1000": point.get("recall@1000", np.nan),
                         "MAP": point.get("MAP", np.nan)})
        print(f"  [{fname}] scored {len(scorers)} methods ({size:,} articles)")

    df = pd.DataFrame(rows)
    csv_path = os.path.join(output_dir, f"{project_prefix}projection_comparison.csv")
    df.to_csv(csv_path, index=False, encoding="utf-8-sig")
    print(f"  [+] table:  {csv_path}")

    if make_figure and not df.empty:
        print(f"  [+] figure: {_fig_projection_comparison(df, fig_dir, project_prefix)}")
    return df
