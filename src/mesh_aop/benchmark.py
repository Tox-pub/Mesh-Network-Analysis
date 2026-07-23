# -*- coding: utf-8 -*-
"""
benchmark.py

Validates and benchmarks the MeSH AOP network against an external, curated
ground-truth citation set (e.g. the OECD AOP reference bibliography).

Treats the evaluation as an extreme-imbalance early-retrieval problem (a handful
of known-relevant PMIDs in a multi-million article pool): validate the ground
truth first, score with metrics that only need positive positions (robust to
incomplete judgments), and quantify uncertainty against baselines.
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

from sklearn.metrics import roc_auc_score, average_precision_score, ndcg_score

from .baseline_manager import _keep_on_device

# Matplotlib rendering is optional; fall back so a headless display never aborts.
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _HAS_MPL = True
except Exception:
    _HAS_MPL = False


def _open_readonly_resilient(db_path):
    """Open a database read-only, resilient to OneDrive dehydration.

    Pins the file and opens read-only; if the read-only view exposes no tables -
    a dehydrated OneDrive placeholder reads back empty - it reopens with a normal
    connection, which forces OneDrive to hydrate the real file. No writes are
    issued, and an already-hydrated file keeps the read-only path unchanged.
    """
    _keep_on_device(db_path)
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True, timeout=120.0)
    if conn.execute("SELECT name FROM sqlite_master WHERE type='table' LIMIT 1").fetchone():
        return conn
    conn.close()
    return sqlite3.connect(str(db_path), timeout=120.0)


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 1. GROUND-TRUTH VALIDATION
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Filenames recognized as a ground-truth drop-in, searched in order within the
# active raw directory (data/raw, or data/reference_raw when using reference data).
# Optional extra scorer: ARS weighted by PageRank recomputed on the filtered
# consensus subgraph instead of the full co-occurrence graph. Only present when
# relevance ran with benchmark.run_ground_truth_analysis enabled; everything that
# touches it degrades gracefully when the column is absent.
OPTIONAL_SUBGRAPH_PR_COL = "score_pagerank_subgraph_centrality"

KNOWN_GT_FILENAMES = (
    "ground_truth_pmids.csv",
    "ground_truth_pmids.txt",
    "ground_truth.csv",
    "ground_truth.txt",
    "oecd_resolved_citations.csv",
)

# Header aliases for auto-detecting columns in arbitrary user files.
_PMID_COL_ALIASES = ("pmid", "pmids", "pubmed_id", "pubmed id", "pubmedid", "pm_id", "id")
_REF_COL_ALIASES = ("raw_reference", "reference", "references", "citation", "cite", "title", "text")


def normalize_pmid(value) -> str:
    """Coerce any PMID to a canonical digit-only string (handles the pandas float '123.0' trap)."""
    if value is None:
        return ""

    s = str(value).strip()
    if not s or s.upper() == "NOT_FOUND" or s.lower() in {"nan", "none"}:
        return ""

    # Strip a trailing ".0" from float coercion, then keep digits only.
    if re.fullmatch(r"\d+\.0+", s):
        s = s.split(".")[0]

    return re.sub(r"\D", "", s)


def _read_text_any_encoding(path: Path) -> str:
    """Read a text file, trying UTF-8 then common Windows code pages."""
    for enc in ("utf-8", "cp1252", "latin-1"):
        try:
            return path.read_text(encoding=enc)
        except (UnicodeDecodeError, UnicodeError):
            continue

    return path.read_text(encoding="utf-8", errors="replace")


def _read_csv_any_encoding(path, **kwargs) -> pd.DataFrame:
    """Read a CSV, trying UTF-8 then common Windows code pages."""
    for enc in ("utf-8", "cp1252", "latin-1"):
        try:
            return pd.read_csv(path, encoding=enc, **kwargs)
        except (UnicodeDecodeError, UnicodeError):
            continue

    return pd.read_csv(path, encoding="utf-8", encoding_errors="replace", **kwargs)


def resolve_ground_truth_path(raw_dir, configured_name: str = "", root=None) -> str:
    """Locate a ground-truth file; honors a configured name, else tries the known defaults.

    A configured name may be a bare filename (looked up in raw_dir), an absolute path,
    or a path relative to the project root (e.g. 'data/reference_processed/...') when
    ``root`` is given -- so the curated OECD ground truth can live outside data/raw.
    """
    raw_dir = Path(raw_dir)

    if configured_name:
        candidates = [Path(configured_name)]                     # absolute or cwd-relative
        if root:
            candidates.append(Path(root) / configured_name)      # project-root relative
        candidates.append(raw_dir / configured_name)             # bare filename in raw_dir
        for p in candidates:
            if p.exists():
                return str(p)
        return None

    for name in KNOWN_GT_FILENAMES:
        p = raw_dir / name
        if p.exists():
            return str(p)

    return None


def _load_gt_dataframe(path: str) -> pd.DataFrame:
    """Load a ground-truth file in any shape (Excel, semicolon CSV, PMID column, single column, or .txt list)."""
    p = Path(path)

    # Excel workbook (e.g. the curated OECD ground truth): read the first sheet.
    if p.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(p, dtype=str)

    # Plain text list: one token per line, with an optional non-numeric header.
    if p.suffix.lower() in (".txt", ""):
        lines = [ln.strip() for ln in _read_text_any_encoding(p).splitlines() if ln.strip()]
        if lines and not any(ch.isdigit() for ch in lines[0]):
            lines = lines[1:]
        return pd.DataFrame({"Raw_Reference": lines, "PMID": lines})

    # Delimited text: sniff the separator from the header line.
    first_line = _read_text_any_encoding(p).splitlines()[:1]
    first_line = first_line[0] if first_line else ""
    sep = ";" if ";" in first_line else ("\t" if "\t" in first_line else ",")

    df = _read_csv_any_encoding(p, sep=sep, dtype=str)

    # A lone column whose header is itself a PMID means there was no header row.
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
            # Pick the column whose cells most look like PMIDs (a 4-9 digit run).
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


def validate_ground_truth(path: str, pmid_col: str = None, ref_col: str = None) -> dict:
    """Load and normalize a ground-truth file, raising clear errors on missing, empty or malformed input.

    Returns the clean PMID set, the unresolved (NOT_FOUND/blank) rows, and summary
    counts. PMID correctness is not guessed here - that is the job of the upstream
    citation re-query; this function only loads the file robustly and reports it.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Ground-truth file not found: {path}\n"
            f"    Place a file named one of {KNOWN_GT_FILENAMES} in the active raw "
            f"directory, or set benchmark.ground_truth_csv in mesh_config.json."
        )

    try:
        df = _load_gt_dataframe(path)
    except Exception as e:
        raise ValueError(
            f"Could not read ground-truth file '{os.path.basename(path)}': {e}. "
            f"Expected a CSV/TSV, a single column of PMIDs, or a .txt list."
        ) from e

    if df is None or len(df) == 0:
        raise ValueError(f"Ground-truth file '{os.path.basename(path)}' is empty (no rows).")

    auto_pmid, auto_ref = _identify_gt_columns(df)
    pmid_col = pmid_col or auto_pmid
    ref_col = ref_col or auto_ref

    if pmid_col is None or pmid_col not in df.columns:
        raise KeyError(
            f"No PMID column found in '{os.path.basename(path)}' (columns: {list(df.columns)}). "
            f"Expected one of {_PMID_COL_ALIASES}, a single column of PMIDs, or a .txt list. "
            f"If the file is delimited, check that the delimiter is ';', ',', or a tab."
        )

    n_raw = len(df)
    df["_pmid_norm"] = df[pmid_col].map(normalize_pmid)
    not_found = df[df["_pmid_norm"] == ""].copy()
    valid_df = df[df["_pmid_norm"] != ""].copy()

    # Deduplicate on the normalized PMID (a reference list can cite a paper twice).
    before = len(valid_df)
    valid_df = valid_df.drop_duplicates(subset="_pmid_norm")
    n_duplicate = before - len(valid_df)

    resolved = set(valid_df["_pmid_norm"])
    if not resolved:
        raise ValueError(
            f"No usable PMIDs parsed from '{os.path.basename(path)}' "
            f"(column '{pmid_col}', {n_raw:,} rows). Check the delimiter and that the "
            f"PMID column actually contains numeric identifiers."
        )

    return {
        "resolved": resolved,
        "valid_df": valid_df,
        "not_found": not_found,
        "n_raw": n_raw,
        "n_resolved": int((df["_pmid_norm"] != "").sum()),
        "n_not_found": len(not_found),
        "n_duplicate": n_duplicate,
        "n_clean": len(resolved),
    }


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 2. RANKING METRICS (robust to incomplete relevance judgments)
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Every metric takes a label array already SORTED by score (descending), where
# labels[i] == 1 marks a known positive at rank i+1. They depend only on positive
# positions and never assume the unlabeled articles are true negatives.

def _positive_ranks(sorted_labels: np.ndarray) -> np.ndarray:
    """1-indexed ranks of the positives in a score-sorted label array."""
    return np.flatnonzero(sorted_labels == 1) + 1


def recall_at_k(sorted_labels: np.ndarray, k: int, total_pos: int) -> float:
    """Fraction of positives recovered within the top K."""
    if total_pos <= 0:
        return float("nan")

    k = min(k, len(sorted_labels))
    return float(sorted_labels[:k].sum()) / total_pos


def precision_at_k(sorted_labels: np.ndarray, k: int) -> float:
    """Precision over the top K (biased low under incomplete judgments)."""
    k = min(k, len(sorted_labels))
    if k == 0:
        return float("nan")

    return float(sorted_labels[:k].sum()) / k


def average_precision(sorted_labels: np.ndarray, total_pos: int) -> float:
    """Average Precision (MAP for a single query)."""
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
    """Precision at rank R, where R is the number of positives."""
    if total_pos <= 0:
        return float("nan")

    r = min(total_pos, len(sorted_labels))
    return float(sorted_labels[:r].sum()) / total_pos


def enrichment_factor(sorted_labels: np.ndarray, frac: float, total_pos: int) -> float:
    """Fold-enrichment of positives in the top `frac` of the ranking (1.0 == random)."""
    n = len(sorted_labels)
    if total_pos <= 0 or n == 0:
        return float("nan")

    top_n = max(1, int(math.ceil(frac * n)))
    found = float(sorted_labels[:top_n].sum())
    expected = total_pos * (top_n / n)

    return found / expected if expected > 0 else float("nan")


def bedroc(sorted_labels: np.ndarray, alpha: float = 20.0) -> float:
    """Boltzmann-Enhanced Discrimination of ROC, an early-recognition score (Truchon & Bayly, 2007)."""
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


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 3. SCORING UTILITIES (tie-aware sorting, bootstrap, permutation null)
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def _sorted_labels(scores: np.ndarray, labels: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Sort labels by score (descending), breaking the many exact ties RANDOMLY to avoid rank bias."""
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
        if out["prevalence"]:
            out["PR_AUC_lift_vs_random"] = out["PR_AUC"] / out["prevalence"]
        else:
            out["PR_AUC_lift_vs_random"] = float("nan")
        for k in ndcg_ks:
            out[f"NDCG@{k}"] = float(
                ndcg_score(labels[None, :], scores[None, :], k=min(k, len(labels)))
            )

    out["recall_at_k"] = {
        int(k): recall_at_k(sl, k, total_pos)
        for k in k_values
    }
    out["precision_at_k"] = {
        int(k): precision_at_k(sl, k)
        for k in k_values
    }
    out["enrichment_factor"] = {
        f"top_{f:g}": enrichment_factor(sl, f, total_pos)
        for f in ef_fracs
    }

    return out


def bootstrap_ci(scores: np.ndarray, labels: np.ndarray, metric_fn, n_boot: int,
                 rng: np.random.Generator, ci: float = 0.95) -> dict:
    """Percentile bootstrap CI for a scalar metric by resampling the pool rows with replacement."""
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


def bootstrap_ci_multi(scores: np.ndarray, labels: np.ndarray, metric_fns: dict,
                       n_boot: int, rng: np.random.Generator, ci: float = 0.95) -> dict:
    """Percentile bootstrap CIs for several metrics off a SHARED set of resamples.

    Resampling the pool and re-sorting it dominates the cost (the sort is O(n log n)
    on the full candidate pool), so each replicate is sorted once and every metric is
    read off that same sorted copy. Computing the statistics in separate passes
    repeated the sort per statistic for no statistical gain: sharing replicates
    leaves each metric's marginal percentile interval unchanged.

    `metric_fns` maps a name to a callable taking (sorted_labels, total_positives).
    """
    n = len(labels)
    vals = {name: [] for name in metric_fns}

    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        s, l = scores[idx], labels[idx]
        total_pos = int(l.sum())
        sl = _sorted_labels(s, l, rng)          # the expensive step - once per replicate
        for name, fn in metric_fns.items():
            try:
                v = fn(sl, total_pos)
            except Exception:
                v = float("nan")
            if v is not None and not (isinstance(v, float) and math.isnan(v)):
                vals[name].append(v)

    lo_q = (1 - ci) / 2 * 100
    hi_q = (1 + ci) / 2 * 100
    out = {}
    for name, v in vals.items():
        if not v:
            out[name] = {"point": float("nan"), "lo": float("nan"), "hi": float("nan")}
        else:
            out[name] = {
                "point": float(np.mean(v)),
                "lo": float(np.percentile(v, lo_q)),
                "hi": float(np.percentile(v, hi_q)),
                "n_boot": len(v),
            }
    return out


def permutation_null(labels: np.ndarray, metric_on_ranked, n_perm: int,
                     rng: np.random.Generator, observed: float) -> dict:
    """Random-ranking null for a metric (shuffled labels), reporting the observed value's lift and p."""
    null_vals = np.empty(n_perm)
    perm_labels = labels.copy()

    for i in range(n_perm):
        rng.shuffle(perm_labels)
        null_vals[i] = metric_on_ranked(perm_labels)

    null_mean = float(np.nanmean(null_vals))
    p_emp = float((np.sum(null_vals >= observed) + 1) / (n_perm + 1))

    if null_mean not in (0.0,) and not math.isnan(null_mean):
        lift = float(observed / null_mean)
    else:
        lift = float("nan")

    return {
        "observed": float(observed),
        "null_mean": null_mean,
        "lift": lift,
        "p_empirical": p_emp,
    }


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 4. DATA LOADING
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def load_relevance_scores(db_path: str,
                          table: str = "article_relevance_scores") -> pd.DataFrame:
    """Load the per-article relevance scores written by relevance.py."""
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Relevance database not found: {db_path}")

    conn = _open_readonly_resilient(db_path)
    try:
        # The subgraph-PageRank score is optional: it only exists when relevance was
        # run with that weighting enabled, so select it only if the column is there.
        present = {r[1] for r in conn.execute(f"PRAGMA table_info({table})")}
        cols = ["pmid", "score_betweenness_centrality", "score_pagerank_centrality",
                "contributing_seeds"]
        if OPTIONAL_SUBGRAPH_PR_COL in present:
            cols.insert(3, OPTIONAL_SUBGRAPH_PR_COL)
        df = pd.read_sql_query(f"SELECT {', '.join(cols)} FROM {table}", conn)
    finally:
        conn.close()

    df["pmid"] = df["pmid"].map(normalize_pmid)
    df = df[df["pmid"] != ""].copy()

    # Number of network nodes (seed terms) an article bridges -> a free structural
    # baseline scorer and the basis of the edge-bridging filter.
    df["n_seeds"] = df["contributing_seeds"].fillna("").str.count(";") + 1
    df.loc[df["contributing_seeds"].fillna("") == "", "n_seeds"] = 0

    return df


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 5. ORCHESTRATION
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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
    """Run the end-to-end validation and benchmark, persisting a JSON report and an enrichment figure."""
    os.makedirs(output_dir, exist_ok=True)
    rng = np.random.default_rng(random_seed)

    report = {
        "meta": {
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "relevance_db": os.path.basename(relevance_db_path),
            "ground_truth": os.path.basename(resolved_csv_path),
            "primary_node": primary_node,
            "random_seed": random_seed,
            "bedroc_alpha": bedroc_alpha,
            "n_boot": n_boot,
            "n_perm": n_perm,
        }
    }

    print("\n" + "<"*30 + ">"*30)
    print("<<< MeSH AOP NETWORK :: GROUND-TRUTH VALIDATION & BENCHMARK >>>")
    print("<"*30 + ">"*30)

    # <<< 5.1 Validate the ground truth >>>
    gt = validate_ground_truth(resolved_csv_path)

    print("\n<<< 1. GROUND-TRUTH VALIDATION >>>")
    print(f"  Raw reference rows ............ {gt['n_raw']:,}")
    print(f"  Resolved to a PMID ............ {gt['n_resolved']:,}")
    print(f"  Unresolved (NOT_FOUND/blank) .. {gt['n_not_found']:,}   (retrieval ceiling loss)")
    print(f"  Duplicate PMIDs removed ....... {gt['n_duplicate']:,}")
    print(f"  CLEAN positives used .......... {gt['n_clean']:,}")

    report["validation"] = {
        k: gt[k]
        for k in ("n_raw", "n_resolved", "n_not_found", "n_duplicate", "n_clean")
    }

    target = gt["resolved"]
    if not target:
        raise RuntimeError("No clean ground-truth PMIDs remain after validation; aborting.")

    # <<< 5.2 Load scores and decompose coverage >>>
    print(f"\n  Loading relevance scores from {os.path.basename(relevance_db_path)} ...")
    scores_df = load_relevance_scores(relevance_db_path)

    pool_pmids = set(scores_df["pmid"])
    scores_df["y_true"] = scores_df["pmid"].isin(target).astype(int)
    scores_df["has_primary"] = scores_df["contributing_seeds"].str.contains(
        primary_node, na=False, regex=False
    ).astype(int)

    # The seed string is only needed to derive n_seeds and has_primary (both done).
    # On a ~9M-row pool it is a multi-GB object column, so release it before the
    # resampling phase allocates its own per-scorer arrays.
    scores_df.drop(columns=["contributing_seeds"], inplace=True)

    in_pool = target & pool_pmids
    naive_pmids = set(scores_df.loc[scores_df["has_primary"] == 1, "pmid"])
    topology_exclusive = in_pool - naive_pmids

    print("\n<<< 2. COVERAGE / RETRIEVABILITY (ceiling, independent of ranking) >>>")
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

    # <<< 5.3 Ranking metrics per scorer, with CIs and null >>>
    scorer_cols = {
        "betweenness": "score_betweenness_centrality",
        "pagerank": "score_pagerank_centrality",
    }
    # Scored head-to-head against whole-corpus PageRank under identical conditions
    # when available; simply absent otherwise.
    if OPTIONAL_SUBGRAPH_PR_COL in scores_df.columns:
        scorer_cols["pagerank (subgraph)"] = OPTIONAL_SUBGRAPH_PR_COL
    scorer_cols["n_seeds (baseline)"] = "n_seeds"
    labels = scores_df["y_true"].to_numpy()

    n_pool = int(len(labels))
    n_pos = int(labels.sum())
    prevalence = (n_pos / n_pool) if n_pool else float("nan")

    print("\n<<< 3. RANKING PERFORMANCE ON FULL POOL (primary metrics) >>>")
    print(f"  Pool {n_pool:,} | positives {n_pos:,} | prevalence {prevalence:.6%}")
    if prevalence < 1e-3:
        # Report the imbalance actually measured rather than assuming it: ROC-AUC stays
        # flattering when negatives outnumber positives this heavily, so the
        # early-recognition metrics are the ones to read.
        print("  Extreme class imbalance -> lead with the early-recognition metrics")
        print("  (Recall@K / MAP / EF / BEDROC): they depend only on WHERE the positives rank.")
        print("  ROC-AUC can look strong here even when precision is ~0, so ROC/PR-AUC are")
        print("  reported as SECONDARY. See README -> 'Benchmark metrics' for why.")
    else:
        print("  ROC/PR-AUC are informative at this prevalence; the early-recognition")
        print("  metrics (Recall@K / MAP / EF / BEDROC) remain the primary ranking view.")

    report["ranking"] = {}
    for name, col in scorer_cols.items():
        scores = scores_df[col].to_numpy(dtype=float)
        suite = compute_metric_suite(
            scores, labels, k_values, ef_fracs, ndcg_ks, bedroc_alpha, rng
        )

        # Named helpers (rather than inline lambdas) for the resampling metrics.
        # Both statistics are read off one shared set of bootstrap replicates: the
        # per-replicate sort of the full pool is the dominant cost, so drawing
        # separate replicate sets for MAP and EF doubled the work for no gain.
        def _map_stat(sorted_labels, total_pos):
            return average_precision(sorted_labels, total_pos)

        def _ef1_stat(sorted_labels, total_pos):
            return enrichment_factor(sorted_labels, 0.01, total_pos)

        def _map_ranked(l):
            return average_precision(l, int(l.sum()))

        cis = bootstrap_ci_multi(
            scores, labels, {"MAP": _map_stat, "EF": _ef1_stat}, n_boot, rng
        )
        ci_map, ci_ef = cis["MAP"], cis["EF"]
        null_map = permutation_null(labels, _map_ranked, n_perm, rng, observed=suite["MAP"])

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

        ef_parts = [f"{k}:{v:.1f}x" for k, v in suite["enrichment_factor"].items()]
        print("    Enrichment ..... " + "  ".join(ef_parts))

        rk_parts = [
            f"@{k}:{100*v:.1f}%"
            for k, v in suite["recall_at_k"].items()
            if not math.isnan(v)
        ]
        print("    Recall@K ....... " + "  ".join(rk_parts))

    # <<< 5.4 Optional negative control >>>
    if negative_control_csv and os.path.exists(negative_control_csv):
        print("\n<<< 4. NEGATIVE CONTROL (unrelated ground truth should score ~random) >>>")
        # The control is optional: a missing/malformed file warns and is skipped
        # rather than aborting the whole benchmark.
        try:
            nc = validate_ground_truth(negative_control_csv)
        except (FileNotFoundError, ValueError, KeyError) as e:
            print(f"  [!] Skipping negative control - could not read it: {e}")
            nc = None

        nc_target = nc["resolved"] if nc else set()
        nc_labels = scores_df["pmid"].isin(nc_target).to_numpy().astype(int)

        if 0 < nc_labels.sum() < len(nc_labels):
            s = scores_df["score_betweenness_centrality"].to_numpy(dtype=float)
            nc_map = average_precision(_sorted_labels(s, nc_labels, rng), int(nc_labels.sum()))

            def _nc_map_ranked(l):
                return average_precision(l, int(l.sum()))

            nc_null = permutation_null(nc_labels, _nc_map_ranked, n_perm, rng, observed=nc_map)

            print(f"  Negative-control MAP {nc_map:.4f}  (lift {nc_null['lift']:.2f}x, "
                  f"p={nc_null['p_empirical']:.3g}) -> should be ~1x / non-significant")
            report["negative_control"] = {
                "MAP": nc_map,
                "lift": nc_null["lift"],
                "p_empirical": nc_null["p_empirical"],
                "n_positives": int(nc_labels.sum()),
            }
        else:
            print("  [!] Negative-control set has no usable overlap with the pool; skipped.")

    # <<< 5.5 Figures >>>
    if make_figures and _HAS_MPL:
        try:
            _plot_enrichment_curves(scores_df, scorer_cols, output_dir, file_prefix, rng)
            report["meta"]["figures"] = f"{file_prefix}_enrichment.png"
            print(f"\n  [+] Wrote enrichment figure to {file_prefix}_enrichment.png")
        except Exception as e:
            print(f"  [!] Figure generation failed: {e}")

    # <<< 5.6 Persist report >>>
    json_path = os.path.join(output_dir, f"{file_prefix}_results.json")
    with open(json_path, "w") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\n  [+] Wrote machine-readable results to {os.path.basename(json_path)}")

    print("<"*30 + ">"*30 + "\n")
    return report


def _plot_enrichment_curves(scores_df, scorer_cols, output_dir, file_prefix, rng):
    """Plot the cumulative-recall (enrichment) curve: recall vs fraction of pool screened."""
    labels = scores_df["y_true"].to_numpy()
    total_pos = int(labels.sum())
    n = len(labels)

    fig, ax = plt.subplots(figsize=(7, 5))
    fracs = np.linspace(0.0001, 1.0, 300)

    for name, col in scorer_cols.items():
        sl = _sorted_labels(scores_df[col].to_numpy(dtype=float), labels, rng)
        cum = np.cumsum(sl)
        recall = [
            cum[min(int(math.ceil(f * n)), n) - 1] / total_pos
            for f in fracs
        ]
        ax.plot(fracs * 100, np.array(recall) * 100, label=name, lw=1.8)

    # Random ranking recovers positives in proportion to the fraction screened.
    # Evaluated on the same grid as the curves: a two-point line starting at x=0 is
    # invalid on the log axis below and renders as a flat line at 100% recall.
    ax.plot(fracs * 100, fracs * 100, "k--", lw=1, label="random")
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
    # Direct-run convenience for reference data; the CLI step resolves paths via
    # MeshConfig instead (see cli.py --step benchmark).
    run_benchmark(
        resolved_csv_path="data/reference_raw/oecd_resolved_citations.csv",
        relevance_db_path="data/reference_processed/DAC_Mesh_contextual_relevance.db",
        output_dir="results",
        file_prefix="DAC_Mesh_benchmark",
    )
