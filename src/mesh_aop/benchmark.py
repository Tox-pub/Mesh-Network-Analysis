# -*- coding: utf-8 -*-
"""
benchmark.py

Validation and benchmarking of the MeSH AOP network against an external,
curated ground-truth citation set (e.g. the OECD AOP reference bibliography).

This module treats the evaluation as an *extreme-imbalance early-retrieval*
problem (a handful of known-relevant PMIDs hidden in a multi-million article
candidate pool) rather than a balanced binary-classification problem. It is
built around three principles:

  1. VALIDATE THE GROUND TRUTH FIRST. External citation->PMID resolution is
     noisy; mis-resolved PMIDs silently corrupt every downstream metric. We
     flag temporally-impossible PMIDs (a year-vs-PMID plausibility check) and
     quarantine them before scoring.

  2. USE METRICS THAT ONLY NEED THE POSITIONS OF KNOWN POSITIVES. The
     unjudged (non-ground-truth) articles are NOT confirmed negatives, so
     precision / PR-AUC over the full pool are pessimistically biased. Recall@K,
     MAP, NDCG, Enrichment Factor and BEDROC only ask "where do the known
     positives rank?" and are therefore robust to incomplete judgments. ROC-AUC
     and PR-AUC are still reported, but clearly demoted to secondary diagnostics.

  3. QUANTIFY UNCERTAINTY AND COMPARE TO BASELINES. With ~70 positives every
     point estimate has a wide confidence interval, so all headline metrics are
     reported with bootstrap CIs and benchmarked against a permutation null
     (random ranking) and simple structural baselines.

Author note: this module is intentionally self-contained and defensive so it can
be developed and reviewed without a live database connection. Run it via
`mesh-pipeline --step benchmark` or import `run_benchmark` directly.
"""

import os
import re
import json
import math
import sqlite3
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd

# scikit-learn provides well-tested AUC / NDCG implementations.
from sklearn.metrics import roc_auc_score, average_precision_score, ndcg_score

# Matplotlib is a hard dependency of the package but rendering is optional here;
# fall back gracefully so a headless / broken display never aborts a benchmark.
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _HAS_MPL = True
except Exception:  # pragma: no cover - environment dependent
    _HAS_MPL = False


# =============================================================================
# 1. GROUND-TRUTH VALIDATION
# =============================================================================

# Approximate (year -> largest PMID assigned by ~end of that year) anchors.
# PMIDs are assigned roughly chronologically, so a PMID far larger than the
# anchor for a citation's publication year is almost certainly a resolution
# error. These are deliberately generous; they only need to catch gross errors
# (e.g. a 2008 paper resolved to a 34-million PMID). Update as PubMed grows.
_PMID_YEAR_ANCHORS = {
    1975: 200_000,
    1985: 6_000_000,
    1995: 8_500_000,
    2000: 11_000_000,
    2005: 16_000_000,
    2010: 21_000_000,
    2015: 26_500_000,
    2020: 33_000_000,
    2023: 38_000_000,
    2025: 40_500_000,
    2027: 43_000_000,
}


def normalize_pmid(value) -> str:
    """Coerce any PMID representation to a canonical digit-only string.

    Critically handles the pandas float-coercion trap: a PMID column containing
    any non-numeric token is read as ``object``/``float``, so ``str(19033392.0)``
    becomes ``'19033392.0'`` and silently fails to match the database's
    ``'19033392'``. We strip to digits only, which fixes both the ``.0`` suffix
    and stray whitespace. Returns ``''`` for unusable values.
    """
    if value is None:
        return ""
    s = str(value).strip()
    if not s or s.upper() == "NOT_FOUND" or s.lower() in {"nan", "none"}:
        return ""
    # Handle a trailing ".0" from float coercion, then keep digits only.
    if re.fullmatch(r"\d+\.0+", s):
        s = s.split(".")[0]
    digits = re.sub(r"\D", "", s)
    return digits


def _max_plausible_pmid(year: int) -> float:
    """Linear-interpolated upper bound on a plausible PMID for a given year."""
    years = sorted(_PMID_YEAR_ANCHORS)
    if year <= years[0]:
        return _PMID_YEAR_ANCHORS[years[0]]
    if year >= years[-1]:
        return _PMID_YEAR_ANCHORS[years[-1]]
    for lo, hi in zip(years, years[1:]):
        if lo <= year <= hi:
            frac = (year - lo) / (hi - lo)
            return _PMID_YEAR_ANCHORS[lo] + frac * (_PMID_YEAR_ANCHORS[hi] - _PMID_YEAR_ANCHORS[lo])
    return _PMID_YEAR_ANCHORS[years[-1]]


def _extract_citation_year(text: str) -> int:
    """Best-effort extraction of a 4-digit publication year from a citation string."""
    if not isinstance(text, str):
        return None
    candidates = [int(y) for y in re.findall(r"\b(19\d{2}|20\d{2})\b", text)]
    # Citations sometimes contain volume/page numbers that look like years; the
    # publication year is almost always the smallest plausible 19xx/20xx token.
    candidates = [y for y in candidates if 1950 <= y <= datetime.now().year + 1]
    return min(candidates) if candidates else None


# Filenames the benchmark recognizes as a ground-truth drop-in, searched in
# order within the active raw directory (data/raw when use_reference_data is
# False, data/reference_raw when True). Users supplying their own ground truth
# simply place one of these in data/raw/ -- no config edit required.
KNOWN_GT_FILENAMES = (
    "ground_truth_pmids.csv",
    "ground_truth_pmids.txt",
    "ground_truth.csv",
    "ground_truth.txt",
    "oecd_resolved_citations.csv",  # bundled reference set
)

# Header aliases for auto-detecting columns in arbitrary user files.
_PMID_COL_ALIASES = ("pmid", "pmids", "pubmed_id", "pubmed id", "pubmedid", "pm_id", "id")
_REF_COL_ALIASES = ("raw_reference", "reference", "references", "citation", "cite", "title", "text")


def _read_text_any_encoding(path: Path) -> str:
    for enc in ("utf-8", "cp1252", "latin-1"):
        try:
            return path.read_text(encoding=enc)
        except (UnicodeDecodeError, UnicodeError):
            continue
    return path.read_text(encoding="utf-8", errors="replace")


def _read_csv_any_encoding(path, **kwargs) -> pd.DataFrame:
    for enc in ("utf-8", "cp1252", "latin-1"):
        try:
            return pd.read_csv(path, encoding=enc, **kwargs)
        except (UnicodeDecodeError, UnicodeError):
            continue
    return pd.read_csv(path, encoding="utf-8", encoding_errors="replace", **kwargs)


def resolve_ground_truth_path(raw_dir, configured_name: str = "") -> str:
    """Locate a ground-truth file in `raw_dir`.

    If `configured_name` is set, only that file is used (returns None if absent
    so the caller can emit a precise error). Otherwise the known default
    filenames are tried in order. Returns an absolute path string or None.
    """
    raw_dir = Path(raw_dir)
    if configured_name:
        p = raw_dir / configured_name
        return str(p) if p.exists() else None
    for name in KNOWN_GT_FILENAMES:
        p = raw_dir / name
        if p.exists():
            return str(p)
    return None


def _load_gt_dataframe(path: str) -> pd.DataFrame:
    """Load a ground-truth file in whatever shape the user provided.

    Accepts, in decreasing structure:
      * the bundled `Raw_Reference;PMID` semicolon CSV,
      * any CSV/TSV with a recognizable PMID column (and optional reference text),
      * a single column of PMIDs (with or without a header),
      * a plain `.txt` list with one PMID per line.

    Always returns a DataFrame guaranteed to contain a usable PMID column;
    a reference-text column is included when present (enables the year check).
    """
    p = Path(path)
    if p.suffix.lower() in (".txt", ""):
        lines = [ln.strip() for ln in _read_text_any_encoding(p).splitlines() if ln.strip()]
        # Drop a non-numeric header line (e.g. "PMID") if present.
        if lines and not any(ch.isdigit() for ch in lines[0]):
            lines = lines[1:]
        return pd.DataFrame({"Raw_Reference": lines, "PMID": lines})

    # Delimited text: sniff the separator from the header line.
    first_line = _read_text_any_encoding(p).splitlines()[:1]
    first_line = first_line[0] if first_line else ""
    sep = ";" if ";" in first_line else ("\t" if "\t" in first_line else ",")
    df = _read_csv_any_encoding(p, sep=sep, dtype=str)

    # A single column whose header is itself a PMID means there was no header row.
    if df.shape[1] == 1 and re.fullmatch(r"\d+(\.0+)?", str(df.columns[0]).strip() or ""):
        df = _read_csv_any_encoding(p, sep=sep, dtype=str, header=None, names=["PMID"])
    return df


def _identify_gt_columns(df: pd.DataFrame):
    """Return (pmid_col, ref_col) by header alias, then by digit-density fallback."""
    lut = {str(c).strip().lower(): c for c in df.columns}
    pmid_col = next((lut[a] for a in _PMID_COL_ALIASES if a in lut), None)
    ref_col = next((lut[a] for a in _REF_COL_ALIASES if a in lut), None)
    if pmid_col is None:
        if df.shape[1] == 1:
            pmid_col = df.columns[0]
        else:
            # Pick the column whose cells are most PMID-like (a 4-9 digit run).
            best, best_frac = None, -1.0
            for c in df.columns:
                vals = df[c].dropna().astype(str)
                if len(vals) == 0:
                    continue
                frac = vals.map(lambda x: bool(re.search(r"\b\d{4,9}\b", x))).mean()
                if frac > best_frac:
                    best, best_frac = c, frac
            pmid_col = best
    return pmid_col, ref_col


def validate_ground_truth(path: str, grace_years: int = 3,
                          min_year_for_check: int = 1990,
                          pmid_col: str = None, ref_col: str = None) -> dict:
    """Load, normalize and sanity-check a ground-truth file (format-flexible).

    The file may be the bundled `Raw_Reference;PMID` CSV, any CSV/TSV with a
    PMID column, a single column of PMIDs, or a `.txt` list (see
    `_load_gt_dataframe`). Column names are auto-detected unless overridden.

    Two complementary plausibility tests flag likely citation->PMID resolution
    errors for human review (we flag, never auto-delete):

      1. ABSOLUTE CEILING: a PMID larger than the largest plausible PMID
         assigned to date cannot exist yet, regardless of the cited year.
      2. CHRONOLOGICAL: for papers published in `min_year_for_check` or later
         (where new PMIDs track publication time), a PMID far exceeding the
         anchor for that year is almost certainly a mis-resolution.

    The chronological test needs reference text (to read the year) and is
    skipped for bare PMID lists and for pre-1990 papers (MEDLINE backfilled
    older literature, so a 1965 paper can legitimately carry a mid-range PMID;
    catching those needs online re-verification, not a magnitude heuristic).

    Returns a dict with:
        resolved    -> set[str] of clean, plausible PMIDs (the scoring target)
        valid_df    -> DataFrame of kept rows
        not_found   -> DataFrame of rows whose citation never resolved to a PMID
        suspect_df  -> DataFrame of rows flagged as implausible
        n_raw, n_resolved, n_not_found, n_suspect, n_duplicate, n_clean
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Ground-truth file not found: {path}")

    df = _load_gt_dataframe(path)
    auto_pmid, auto_ref = _identify_gt_columns(df)
    pmid_col = pmid_col or auto_pmid
    ref_col = ref_col or auto_ref
    if pmid_col is None or pmid_col not in df.columns:
        raise KeyError(f"Could not find a PMID column in {os.path.basename(path)}; "
                       f"columns present: {list(df.columns)}. Expected one of "
                       f"{_PMID_COL_ALIASES}, a single PMID column, or a .txt list.")
    n_raw = len(df)

    df["_pmid_norm"] = df[pmid_col].map(normalize_pmid)
    not_found = df[df["_pmid_norm"] == ""].copy()
    resolved_rows = df[df["_pmid_norm"] != ""].copy()

    # Year-vs-PMID plausibility flag.
    ref_series = resolved_rows[ref_col] if ref_col in resolved_rows.columns else pd.Series("", index=resolved_rows.index)
    resolved_rows["_cite_year"] = ref_series.map(_extract_citation_year)

    # Ceiling = the largest PMID plausibly assigned by now. Based on the current
    # year (not the far-future anchor) so impossibly-large PMIDs are caught while
    # genuine just-published papers are not. Keep _PMID_YEAR_ANCHORS current.
    absolute_ceiling = _max_plausible_pmid(datetime.now().year)

    def _is_suspect(row):
        try:
            pmid_int = int(row["_pmid_norm"])
        except ValueError:
            return True
        # Test 1: cannot exist yet, independent of citation year.
        if pmid_int > absolute_ceiling:
            return True
        # Test 2: chronological, only where PMIDs reliably track time.
        year = row["_cite_year"]
        if year is None or year < min_year_for_check:
            return False
        return pmid_int > _max_plausible_pmid(year + grace_years)

    resolved_rows["_suspect"] = resolved_rows.apply(_is_suspect, axis=1)
    suspect_df = resolved_rows[resolved_rows["_suspect"]].copy()
    valid_df = resolved_rows[~resolved_rows["_suspect"]].copy()

    # Deduplicate on normalized PMID (a reference list can cite the same paper twice).
    before = len(valid_df)
    valid_df = valid_df.drop_duplicates(subset="_pmid_norm")
    n_duplicate = before - len(valid_df)

    resolved = set(valid_df["_pmid_norm"])
    return {
        "resolved": resolved,
        "valid_df": valid_df,
        "not_found": not_found,
        "suspect_df": suspect_df,
        "n_raw": n_raw,
        "n_resolved": int((df["_pmid_norm"] != "").sum()),
        "n_not_found": len(not_found),
        "n_suspect": len(suspect_df),
        "n_duplicate": n_duplicate,
        "n_clean": len(resolved),
    }


# =============================================================================
# 2. RANKING METRICS (robust to incomplete relevance judgments)
# =============================================================================
# Convention: every metric operates on a label array already SORTED by score
# (descending). `labels[i] == 1` means the article at rank i+1 is a known
# positive. These metrics depend only on the positions of positives and never
# assume the unlabeled articles are true negatives.

def _positive_ranks(sorted_labels: np.ndarray) -> np.ndarray:
    """1-indexed ranks of the positive items in a score-sorted label array."""
    return np.flatnonzero(sorted_labels == 1) + 1


def recall_at_k(sorted_labels: np.ndarray, k: int, total_pos: int) -> float:
    if total_pos <= 0:
        return float("nan")
    k = min(k, len(sorted_labels))
    return float(sorted_labels[:k].sum()) / total_pos


def precision_at_k(sorted_labels: np.ndarray, k: int) -> float:
    """Precision over the top-K. NOTE: biased low under incomplete judgments
    because unjudged articles in the top-K are counted as non-relevant."""
    k = min(k, len(sorted_labels))
    if k == 0:
        return float("nan")
    return float(sorted_labels[:k].sum()) / k


def average_precision(sorted_labels: np.ndarray, total_pos: int) -> float:
    """Mean Average Precision for a single query (a.k.a. AP)."""
    if total_pos <= 0:
        return float("nan")
    hits = 0
    cum = 0.0
    for i, lab in enumerate(sorted_labels, start=1):
        if lab == 1:
            hits += 1
            cum += hits / i
    return cum / total_pos


def r_precision(sorted_labels: np.ndarray, total_pos: int) -> float:
    """Precision at rank R, where R = number of positives."""
    if total_pos <= 0:
        return float("nan")
    r = min(total_pos, len(sorted_labels))
    return float(sorted_labels[:r].sum()) / total_pos


def enrichment_factor(sorted_labels: np.ndarray, frac: float, total_pos: int) -> float:
    """Enrichment factor at the top `frac` of the ranking.

    EF = (positives recovered in top fraction) / (positives expected by chance).
    EF = 1.0 is random; EF = 10 means the top slice is 10x richer in positives
    than the pool. Standard in virtual screening, which has the same shape as
    this problem (few 'actives' in a huge library).
    """
    n = len(sorted_labels)
    if total_pos <= 0 or n == 0:
        return float("nan")
    top_n = max(1, int(math.ceil(frac * n)))
    found = float(sorted_labels[:top_n].sum())
    expected = total_pos * (top_n / n)
    return found / expected if expected > 0 else float("nan")


def bedroc(sorted_labels: np.ndarray, alpha: float = 20.0) -> float:
    """Boltzmann-Enhanced Discrimination of ROC (Truchon & Bayly, 2007).

    An early-recognition metric in [0, 1] that exponentially weights top ranks;
    alpha=20 means ~80% of the score comes from the top ~8% of the ranking.
    BEDROC ~ 0.5 corresponds to random for a balanced set; the random baseline
    for a sparse set is reported separately via the permutation null, so treat
    this as a relative score and compare against that null.
    """
    n_total = len(sorted_labels)
    ranks = _positive_ranks(sorted_labels)
    n = len(ranks)
    if n == 0 or n_total == 0 or n == n_total:
        return float("nan")
    ra = n / n_total
    s = np.sum(np.exp(-alpha * ranks / n_total))
    rie = s / ((n / n_total) * (1 - math.exp(-alpha)) / (math.exp(alpha / n_total) - 1))
    factor = (ra * math.sinh(alpha / 2)) / (math.cosh(alpha / 2) - math.cosh(alpha / 2 - alpha * ra))
    return float(rie * factor + 1.0 / (1 - math.exp(alpha * (1 - ra))))


# =============================================================================
# 3. SCORING UTILITIES (tie-aware sorting, bootstrap, permutation null)
# =============================================================================

def _sorted_labels(scores: np.ndarray, labels: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Sort labels by score descending, breaking ties RANDOMLY.

    The relevance scores contain many exact ties (articles sharing the same set
    of contributing seeds get identical scores). A naive stable sort would then
    order tied articles by their original DB order, biasing every rank metric.
    We inject a random secondary key so ties are resolved fairly and the metric
    reflects expected performance over arbitrary tie orderings.
    """
    jitter = rng.random(len(scores))
    order = np.lexsort((jitter, -scores))  # primary: -scores (desc); secondary: random
    return labels[order]


def compute_metric_suite(scores: np.ndarray, labels: np.ndarray, k_values, ef_fracs,
                         ndcg_ks, alpha: float, rng: np.random.Generator) -> dict:
    """Compute the full ranking-metric suite for one scorer on one pool."""
    total_pos = int(labels.sum())
    sl = _sorted_labels(scores, labels, rng)
    out = {
        "n_pool": int(len(labels)),
        "n_positives": total_pos,
        "prevalence": total_pos / len(labels) if len(labels) else float("nan"),
        "MAP": average_precision(sl, total_pos),
        "R_precision": r_precision(sl, total_pos),
        "BEDROC": bedroc(sl, alpha=alpha),
    }
    # AUCs are secondary diagnostics; guard against degenerate single-class pools.
    if 0 < total_pos < len(labels):
        out["ROC_AUC"] = float(roc_auc_score(labels, scores))
        out["PR_AUC"] = float(average_precision_score(labels, scores))
        out["PR_AUC_lift_vs_random"] = out["PR_AUC"] / out["prevalence"] if out["prevalence"] else float("nan")
        for k in ndcg_ks:
            out[f"NDCG@{k}"] = float(ndcg_score(labels[None, :], scores[None, :], k=min(k, len(labels))))
    out["recall_at_k"] = {int(k): recall_at_k(sl, k, total_pos) for k in k_values}
    out["precision_at_k"] = {int(k): precision_at_k(sl, k) for k in k_values}
    out["enrichment_factor"] = {f"top_{f:g}": enrichment_factor(sl, f, total_pos) for f in ef_fracs}
    return out


def bootstrap_ci(scores: np.ndarray, labels: np.ndarray, metric_fn, n_boot: int,
                 rng: np.random.Generator, ci: float = 0.95) -> dict:
    """Percentile bootstrap CI for a scalar metric by resampling the pool rows.

    Resampling rows with replacement is a nonparametric bootstrap of the
    candidate pool; it propagates the substantial uncertainty that comes from
    having only a few dozen positives.
    """
    n = len(labels)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        try:
            v = metric_fn(scores[idx], labels[idx])
        except Exception:
            v = float("nan")
        if v is not None and not (isinstance(v, float) and math.isnan(v)):
            vals.append(v)
    if not vals:
        return {"point": float("nan"), "lo": float("nan"), "hi": float("nan")}
    lo = (1 - ci) / 2 * 100
    hi = (1 + ci) / 2 * 100
    return {
        "point": float(np.mean(vals)),
        "lo": float(np.percentile(vals, lo)),
        "hi": float(np.percentile(vals, hi)),
        "n_boot": len(vals),
    }


def permutation_null(labels: np.ndarray, metric_on_ranked, n_perm: int,
                     rng: np.random.Generator, observed: float) -> dict:
    """Random-ranking null distribution for a metric, plus the observed lift.

    A randomly shuffled label array IS a random ranking, so evaluating the
    metric on shuffled labels samples its null distribution. `metric_on_ranked`
    maps a label array (already in ranked order) to the metric value. `observed`
    is the REAL scorer's value (computed elsewhere); we report its lift over the
    null mean and a one-sided empirical p-value.
    """
    null_vals = np.empty(n_perm)
    perm_labels = labels.copy()
    for i in range(n_perm):
        rng.shuffle(perm_labels)
        null_vals[i] = metric_on_ranked(perm_labels)
    null_mean = float(np.nanmean(null_vals))
    p_emp = float((np.sum(null_vals >= observed) + 1) / (n_perm + 1))
    lift = float(observed / null_mean) if null_mean not in (0.0,) and not math.isnan(null_mean) else float("nan")
    return {
        "observed": float(observed),
        "null_mean": null_mean,
        "lift": lift,
        "p_empirical": p_emp,
    }


# =============================================================================
# 4. DATA LOADING
# =============================================================================

def load_relevance_scores(db_path: str,
                          table: str = "article_relevance_scores") -> pd.DataFrame:
    """Load the per-article relevance scores written by relevance.py."""
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Relevance database not found: {db_path}")
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        df = pd.read_sql_query(
            f"SELECT pmid, score_betweenness_centrality, score_eigenvector_centrality, "
            f"contributing_seeds FROM {table}",
            conn,
        )
    finally:
        conn.close()
    df["pmid"] = df["pmid"].map(normalize_pmid)
    df = df[df["pmid"] != ""].copy()
    # Number of network nodes (seed terms) an article bridges -> a free
    # structural baseline scorer and the basis of the edge-bridging filter.
    df["n_seeds"] = df["contributing_seeds"].fillna("").str.count(";") + 1
    df.loc[df["contributing_seeds"].fillna("") == "", "n_seeds"] = 0
    return df


# =============================================================================
# 5. ORCHESTRATION
# =============================================================================

def run_benchmark(resolved_csv_path: str, relevance_db_path: str,
                  output_dir: str, file_prefix: str = "benchmark",
                  primary_node: str = "Dermatitis, Allergic Contact",
                  negative_control_csv: str = None,
                  k_values=(100, 500, 1000, 5000, 10000, 50000),
                  ef_fracs=(0.001, 0.01, 0.05, 0.10),
                  ndcg_ks=(1000, 10000),
                  bedroc_alpha: float = 20.0,
                  n_boot: int = 2000,
                  n_perm: int = 2000,
                  random_seed: int = 42,
                  make_figures: bool = True) -> dict:
    """End-to-end validation + benchmark. Persists a JSON report and figures.

    Scorers evaluated on the common candidate pool:
        * betweenness  -> score_betweenness_centrality
        * eigenvector  -> score_eigenvector_centrality
        * n_seeds      -> structural baseline (# of bridged network nodes)
    plus a permutation (random-ranking) null for every metric.
    """
    os.makedirs(output_dir, exist_ok=True)
    rng = np.random.default_rng(random_seed)
    report = {"meta": {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "relevance_db": os.path.basename(relevance_db_path),
        "ground_truth": os.path.basename(resolved_csv_path),
        "primary_node": primary_node,
        "random_seed": random_seed,
        "bedroc_alpha": bedroc_alpha,
        "n_boot": n_boot,
        "n_perm": n_perm,
    }}

    print("\n" + "=" * 80)
    print("<<< MeSH AOP NETWORK :: GROUND-TRUTH VALIDATION & BENCHMARK >>>")
    print("=" * 80)

    # ---- 5.1 Validate the ground truth -------------------------------------
    gt = validate_ground_truth(resolved_csv_path)
    print("\n--- 1. GROUND-TRUTH VALIDATION ---")
    print(f"  Raw reference rows ............ {gt['n_raw']:,}")
    print(f"  Resolved to a PMID ............ {gt['n_resolved']:,}")
    print(f"  Unresolved (NOT_FOUND) ........ {gt['n_not_found']:,}   <- retrieval ceiling loss")
    print(f"  Flagged temporally implausible  {gt['n_suspect']:,}   <- quarantined")
    print(f"  Duplicate PMIDs removed ....... {gt['n_duplicate']:,}")
    print(f"  CLEAN positives used .......... {gt['n_clean']:,}")

    # Persist the quarantine so a human can adjudicate the flagged rows.
    if gt["n_suspect"] > 0:
        qpath = os.path.join(output_dir, f"{file_prefix}_quarantined_pmids.csv")
        cols = [c for c in ["Raw_Reference", "PMID", "_pmid_norm", "_cite_year"] if c in gt["suspect_df"].columns]
        gt["suspect_df"][cols].to_csv(qpath, sep=";", index=False)
        print(f"  [!] Wrote {gt['n_suspect']} suspect rows to {os.path.basename(qpath)} for manual review.")
    report["validation"] = {k: gt[k] for k in
                            ("n_raw", "n_resolved", "n_not_found", "n_suspect", "n_duplicate", "n_clean")}

    target = gt["resolved"]
    if not target:
        raise RuntimeError("No clean ground-truth PMIDs remain after validation; aborting.")

    # ---- 5.2 Load scores and decompose coverage ----------------------------
    print(f"\n  Loading relevance scores from {os.path.basename(relevance_db_path)} ...")
    scores_df = load_relevance_scores(relevance_db_path)
    pool_pmids = set(scores_df["pmid"])
    scores_df["y_true"] = scores_df["pmid"].isin(target).astype(int)
    scores_df["has_primary"] = scores_df["contributing_seeds"].str.contains(
        primary_node, na=False, regex=False).astype(int)

    in_pool = target & pool_pmids
    naive_pmids = set(scores_df.loc[scores_df["has_primary"] == 1, "pmid"])
    topology_exclusive = in_pool - naive_pmids

    print("\n--- 2. COVERAGE / RETRIEVABILITY (ceiling, independent of ranking) ---")
    print(f"  Candidate pool size ........... {len(pool_pmids):,}")
    print(f"  Clean ground-truth positives .. {len(target):,}")
    print(f"  ...present in pool (ceiling) .. {len(in_pool):,}  ({100*len(in_pool)/len(target):.1f}%)")
    print(f"  ...also under primary node .... {len(in_pool & naive_pmids):,}")
    print(f"  ...topology-exclusive ......... {len(topology_exclusive):,}  "
          f"(found WITHOUT the primary disease node)")
    report["coverage"] = {
        "pool_size": len(pool_pmids),
        "clean_positives": len(target),
        "positives_in_pool": len(in_pool),
        "retrievability_ceiling": len(in_pool) / len(target),
        "under_primary_node": len(in_pool & naive_pmids),
        "topology_exclusive": len(topology_exclusive),
    }

    # ---- 5.3 Ranking metrics per scorer, with CIs and null -----------------
    scorer_cols = {
        "betweenness": "score_betweenness_centrality",
        "eigenvector": "score_eigenvector_centrality",
        "n_seeds (baseline)": "n_seeds",
    }
    labels = scores_df["y_true"].to_numpy()

    print("\n--- 3. RANKING PERFORMANCE ON FULL POOL (primary metrics) ---")
    print("  (Recall@K / MAP / EF / BEDROC need only positive positions; "
          "ROC/PR-AUC are secondary.)")
    report["ranking"] = {}
    for name, col in scorer_cols.items():
        scores = scores_df[col].to_numpy(dtype=float)
        suite = compute_metric_suite(scores, labels, k_values, ef_fracs,
                                     ndcg_ks, bedroc_alpha, rng)

        # Bootstrap CIs for two headline metrics (MAP and EF@1%).
        ci_map = bootstrap_ci(scores, labels,
                              lambda s, l: average_precision(_sorted_labels(s, l, rng), int(l.sum())),
                              n_boot, rng)
        ci_ef = bootstrap_ci(scores, labels,
                             lambda s, l: enrichment_factor(_sorted_labels(s, l, rng), 0.01, int(l.sum())),
                             n_boot, rng)
        # Permutation null / lift for MAP (observed = the real scorer's MAP).
        null_map = permutation_null(labels,
                                    lambda l: average_precision(l, int(l.sum())),
                                    n_perm, rng, observed=suite["MAP"])
        suite["MAP_ci95"] = [ci_map["lo"], ci_map["hi"]]
        suite["EF_top_0.01_ci95"] = [ci_ef["lo"], ci_ef["hi"]]
        suite["MAP_vs_random_lift"] = null_map["lift"]
        suite["MAP_p_empirical"] = null_map["p_empirical"]
        report["ranking"][name] = suite

        print(f"\n  [{name}]  scorer = {col}")
        print(f"    MAP ............ {suite['MAP']:.4f}  "
              f"(95% CI {ci_map['lo']:.4f}-{ci_map['hi']:.4f}; "
              f"{null_map['lift']:.1f}x random, p={null_map['p_empirical']:.3g})")
        print(f"    R-precision .... {suite['R_precision']:.4f}")
        print(f"    BEDROC(a={bedroc_alpha:g}) {suite['BEDROC']:.4f}")
        if "ROC_AUC" in suite:
            print(f"    ROC-AUC ....... {suite['ROC_AUC']:.4f}   PR-AUC {suite['PR_AUC']:.5f} "
                  f"({suite['PR_AUC_lift_vs_random']:.0f}x baseline)  [secondary]")
        ef = suite["enrichment_factor"]
        print("    Enrichment ..... " + "  ".join(f"{k}:{v:.1f}x" for k, v in ef.items()))
        rk = suite["recall_at_k"]
        print("    Recall@K ....... " + "  ".join(
            f"@{k}:{100*v:.1f}%" for k, v in rk.items() if not math.isnan(v)))

    # ---- 5.4 Optional negative control -------------------------------------
    if negative_control_csv and os.path.exists(negative_control_csv):
        print("\n--- 4. NEGATIVE CONTROL (unrelated ground truth should score ~random) ---")
        nc = validate_ground_truth(negative_control_csv)
        nc_target = nc["resolved"]
        nc_labels = scores_df["pmid"].isin(nc_target).to_numpy().astype(int)
        if 0 < nc_labels.sum() < len(nc_labels):
            s = scores_df["score_betweenness_centrality"].to_numpy(dtype=float)
            nc_map = average_precision(_sorted_labels(s, nc_labels, rng), int(nc_labels.sum()))
            nc_null = permutation_null(nc_labels,
                                       lambda l: average_precision(l, int(l.sum())),
                                       n_perm, rng, observed=nc_map)
            print(f"  Negative-control MAP {nc_map:.4f}  (lift {nc_null['lift']:.2f}x, "
                  f"p={nc_null['p_empirical']:.3g}) -> should be ~1x / non-significant")
            report["negative_control"] = {"MAP": nc_map, "lift": nc_null["lift"],
                                          "p_empirical": nc_null["p_empirical"],
                                          "n_positives": int(nc_labels.sum())}
        else:
            print("  [!] Negative-control set has no usable overlap with the pool; skipped.")

    # ---- 5.5 Figures -------------------------------------------------------
    if make_figures and _HAS_MPL:
        try:
            _plot_enrichment_curves(scores_df, scorer_cols, output_dir, file_prefix, rng)
            report["meta"]["figures"] = f"{file_prefix}_enrichment.png"
            print(f"\n  [+] Wrote enrichment figure to {file_prefix}_enrichment.png")
        except Exception as e:  # pragma: no cover
            print(f"  [!] Figure generation failed: {e}")

    # ---- 5.6 Persist report ------------------------------------------------
    json_path = os.path.join(output_dir, f"{file_prefix}_results.json")
    with open(json_path, "w") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\n  [+] Wrote machine-readable results to {os.path.basename(json_path)}")
    print("=" * 80 + "\n")
    return report


def _plot_enrichment_curves(scores_df, scorer_cols, output_dir, file_prefix, rng):
    """Cumulative recall (enrichment) curve: recall vs fraction of pool screened."""
    labels = scores_df["y_true"].to_numpy()
    total_pos = int(labels.sum())
    n = len(labels)
    fig, ax = plt.subplots(figsize=(7, 5))
    fracs = np.linspace(0.0001, 1.0, 300)
    for name, col in scorer_cols.items():
        sl = _sorted_labels(scores_df[col].to_numpy(dtype=float), labels, rng)
        cum = np.cumsum(sl)
        recall = [cum[min(int(math.ceil(f * n)), n) - 1] / total_pos for f in fracs]
        ax.plot(fracs * 100, np.array(recall) * 100, label=name, lw=1.8)
    ax.plot([0, 100], [0, 100], "k--", lw=1, label="random")
    ax.set_xlabel("Fraction of pool screened (%)")
    ax.set_ylabel("Ground-truth recall (%)")
    ax.set_title("Enrichment: ground-truth recovery vs ranking depth")
    ax.set_xscale("log")
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{file_prefix}_enrichment.png"), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    # Direct-run convenience for reference data. The CLI step resolves paths via
    # MeshConfig instead; see cli.py --step benchmark.
    run_benchmark(
        resolved_csv_path="data/reference_raw/oecd_resolved_citations.csv",
        relevance_db_path="data/reference_processed/DAC_Mesh_contextual_relevance.db",
        output_dir="results",
        file_prefix="DAC_Mesh_benchmark",
    )
