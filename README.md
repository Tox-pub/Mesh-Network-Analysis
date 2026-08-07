# Medical Subject Headings (MeSH) Adverse Outcome Pathway (AOP) Network Pipeline

## Overview

This repository contains a comprehensive computational pipeline designed to construct, filter, and analyze knowledge graphs representing Adverse Outcome Pathways (**AOPs**) and biological flows. By leveraging the NCBI Entrez API and the complete offline NLM PubMed Baseline, the pipeline extracts primary literature associated with specific Medical Subject Headings (**MeSH**), maps multi-generational citation topologies, and calculates semantic co-occurrence networks.

The system connects **Stressors** (e.g., chemicals) to **Adverse Outcomes** (e.g., diseases) through biological intermediates. It utilizes Global-Local Filtering (GLF) and Simulated Annealing (SA) to optimize subgraph density, Louvain heuristics for community detection, and Mean Relevancy Scoring (MRS) to rank nodes and edges based on their impact within the global corpus of literature.

---

## Repository Structure

The package assumes and enforces the following directory architecture. The required MeSH XML input file must be placed in the raw data directory prior to execution.

```text
Mesh-Network-Analysis-Main-Library/
│
├── data/                               # Data storage
│   ├── raw/                            # Inputs for a run
│   │   ├── aop_annotations_master.csv  # Ships w/ repo: AOP strata dictionary (pre-seeded; grows each run)
│   │   ├── desc2025.xml                # Auto-downloaded from NLM if missing (or place manually); not in repo
│   │   ├── ground_truth_pmids.csv      # Optional, you place this: YOUR benchmark set (see "Supplying Your Own Ground Truth")
│   │   ├── master_mesh_database.db     # Auto-generated: offline PubMed corpus (Step 0)
│   │   ├── pubmed_baseline/            # Auto-downloaded: NLM Baseline XMLs (~40GB, Step 0)
│   │   └── pubmed_updates/             # Auto-downloaded: NLM Daily Update XMLs (optional)
│   ├── processed/                      # Auto-generated: pipeline databases and JSONs (starts empty)
│   ├── reference_raw/                  # Ships w/ repo: bundled reference inputs
│   │   └── oecd_resolved_citations.csv # OECD AOP-40 citation->PMID table (the bundled ground-truth source)
│   └── reference_processed/            # Ships w/ repo: curated OECD ground-truth set + bundled reference network
│
├── results/                            # Output artifacts (auto-generated)
│   ├── figures/                        # High-resolution pipeline plots (.png, .tif, .html)
│   ├── benchmark/                      # All --step benchmark outputs (ranking + ground-truth)
│   │   └── validation/                 # Node-weighting + projection report
│   ├── logs/                           # System logs and failed fetch records
│   ├── *_run_annotations.csv           # Run-specific AOP annotation templates
│   ├── *_Top_Network_Articles.csv      # Secondary analysis exports
│   └── *_export.xlsx                   # Exported full network tables
│
├── src/
│   ├── mesh_aop/                       # Core Python package modules
│   │   ├── __init__.py
│   │   ├── baseline_manager.py         # Multi-core MapReduce ETL for the Master Database
│   │   ├── benchmark.py                # Ground-truth validation & performance benchmarking
│   │   ├── check_env.py                # System environment & dependency verification
│   │   ├── cli.py                      # Orchestrator and CLI entry point
│   │   ├── config_parser.py            # Two-tier configuration engine
│   │   ├── data_ops.py                 # SQLite and NCBI Entrez querying
│   │   ├── gt_network_validation.py    # Node/edge convergent ground-truth validation
│   │   ├── mesh_data_processor.py      # Unified XML extraction and stop-word generation
│   │   ├── mesh_stop_words.py          # Auto-generated MeSH stop-word set
│   │   ├── network.py                  # NetworkX assembly, filtering, and centrality
│   │   ├── node2vec_embedding.py       # Vendored Node2Vec embedding (removes the node2vec dep)
│   │   ├── relevance.py                # Mean Relevancy Scoring (Semantic Re-ranking)
│   │   ├── secondary_analysis.py       # Metadata hydration and targeted graph querying
│   │   ├── stats.py                    # GLF/SA mathematical models and graph statistics
│   │   ├── validation_report.py        # Consolidated node-weighting + projection evaluation
│   │   ├── viz.py                      # Matplotlib, Seaborn, and Plotly graphics
│   │   └── wizard.py                   # Interactive configuration module
│   └── mesh_aop_notebooks/             # Jupyter notebook equivalents of each module
│       └── *.ipynb                     # One notebook per module for interactive exploration
│
├── environment.yml                     # Mamba/Conda cross-platform dependency resolution
├── pyproject.toml                      # Modern Python package specification
├── mesh_config.json                    # Runtime user config (auto-generated; git-ignored, not in repo)
├── LICENSE                             # MIT License
└── README.md                           # This document


```

---

## Data Acquisition & Prerequisites

### 1. The MeSH XML File (automatic)

The NLM has officially discontinued the MeSH ASCII format as of 2026. This pipeline now utilizes the computational gold-standard **MeSH XML format**.

* **Automatic download (default):** You do **not** need to fetch this file by hand. When Step 1 builds the MeSH support files, it checks `data/raw/` for the current descriptor XML and, if it is missing or a newer annual release has appeared, downloads the latest `descYYYY.xml` straight from NLM into `data/raw/` (verifying it is the genuine ~300 MB file, not an error page). The year advances automatically — `desc2025.xml` today, `desc2026.xml` once NLM publishes it.
* **Manual placement (optional / offline):** If you prefer, or if the build machine has no internet, download `descYYYY.xml` yourself from the [NLM MeSH Data Distribution page](https://nlmpubs.nlm.nih.gov/projects/mesh/) and drop it into `data/raw/`. An existing local copy is always reused, and is the offline fallback when NLM cannot be reached.

### 2. Internet Connectivity & Disk Budget

Internet access is required only for the steps that talk to NCBI; everything else runs against your local databases:

| Pipeline step | Internet? | When |
|---|---|---|
| Step 0 — Master DB build (baseline FTP) | **Yes** | First run only; again only to refresh to a newer yearly baseline or apply daily updates |
| Step 2 — Article collection (Entrez + citation links) | **Yes** | **Every** time you build a new query's citation database |
| Step 3.5 — Secondary metadata hydration | **Yes** | Whenever you export top articles / run `--step secondary` |
| Step 1 — MeSH support-file build | Only if XML absent | Downloads `descYYYY.xml` from NLM when it is missing or superseded; fully offline once the file is present |
| Steps 3, 4, benchmark | No | Run fully offline against the local databases |

Once the master database is built, everyday analysis is offline — you only reconnect to build a **new query's** database (Step 2, which fetches the P0 cohort plus its incoming/outgoing citation links) or to update the master corpus. These baseline and daily-update archives are the official NLM/NCBI PubMed releases, distributed at the [NCBI PubMed Data Distribution page](https://pubmed.ncbi.nlm.nih.gov/download/).

**Disk budget.** Set aside **~60–80 GB free** to build the master database the first time: the pipeline downloads the full set of NLM baseline `.xml.gz` archives (~40 GB) and expands them into the ~7 GB SQLite master database, with working headroom during extraction. This one-time bulk download plus local SQLite lookups is dramatically faster and more reliable than issuing millions of live Entrez queries per run, and is precisely what lets Steps 1/3/4/benchmark run fully offline. Once the master database is built and verified, you may **delete the `pubmed_baseline/` archives** to reclaim the ~40 GB — provided you keep `master_mesh_database.db` and do not intend to run local daily-update ingestion (which re-reads those archives).

**Daily updates (optional).** Between annual baselines, NLM publishes daily update archives (`pubmed_updates/`). Applying them keeps the master corpus current with the newest PMIDs, but re-runs the multi-core ETL and requires the update archives to be present. A fresh yearly baseline supersedes accumulated daily updates, so for most analyses the annual baseline alone is sufficient — enable daily updates only if you specifically need very recent publications.

---

## Environment Setup & Installation

This pipeline is computationally intensive and relies heavily on numerical arrays and graph operations. Due to the complexity of the underlying C-extensions in libraries like SciPy and NumPy, we strongly recommend using Mamba/Micromamba to resolve dependencies and ensure cross-platform compatibility.

### 1. System Requirements

* **Python:** **3.11–3.13** (`requires-python = ">=3.11,<3.14"`; validated against the 3.13 stack). Python 3.13 **is** supported — the Node2Vec embedding is vendored in `node2vec_embedding.py`, so the `node2vec` package, which pinned `numpy<2.0` and had no 3.13 wheel, is no longer a dependency. Releases newer than 3.13 are untested.
* **Memory:** Minimum 16GB RAM; 32GB+ highly recommended for Step 0 (Database Compilation) and networks exceeding 1 citation generation.
* **Storage:** 100GB+ free space (The NLM Baseline XMLs and resulting SQLite Master Database expand rapidly).
* **Windows path length:** Create the virtual environment at a **short path** (e.g. `C:\Users\<you>\mesh_env`), *not* inside a deeply nested or OneDrive-synced project folder. Some dependencies (e.g. `statsmodels`) ship very long filenames that overflow the Windows 260-character `MAX_PATH` limit and abort the install. See **Troubleshooting**.

### 2. Installation (Standard Pip/Venv)

If Mamba is unavailable, use a standard **Python 3.11/3.12** virtual environment. This relies on `pyproject.toml` to resolve dependencies via pip. The steps differ by OS — macOS/Linux are straightforward; Windows needs two workarounds (path length and script/exe blocking).

#### macOS / Linux

```bash
cd path/Mesh-Network-Analysis-Main-Library
python3.12 -m venv mesh_env          # or python3.11
source mesh_env/bin/activate
pip install --upgrade pip
pip install -e .
mesh-check-env --auto                 # entry-point commands work directly here
```

On macOS, if you don't have Python 3.12, install it with `brew install python@3.12` (Homebrew) or from python.org. **None of the Windows caveats below apply on macOS/Linux** — activation works, the `mesh-pipeline` / `mesh-check-env` commands run directly, and there's no path-length limit.

#### Windows

Two machine-level guardrails (common on corporate/managed laptops) break the Linux-style flow, so the procedure is different:

* **Path length:** create the venv at a **short path outside** the project — the deep OneDrive path plus `statsmodels`' long filenames overflow the 260-character `MAX_PATH` limit and abort the install.
* **Script/exe blocking:** `Activate.ps1` is often blocked by execution policy, and the generated `mesh-pipeline.exe` / `mesh-check-env.exe` launchers are often blocked by security software. **You don't need either** — call the venv's `python.exe` by full path and use the module form (`-m mesh_aop.cli`). `python.exe` itself is not blocked.

```powershell
cd "C:\path\to\Mesh-Network-Analysis-Main-Library"

# Create the venv OUTSIDE the project, at a short path (Python 3.12)
py -3.12 -m venv "$env:USERPROFILE\mesh_env"

# Install using the venv's python by full path (no activation needed)
& "$env:USERPROFILE\mesh_env\Scripts\python.exe" -m pip install --upgrade pip
& "$env:USERPROFILE\mesh_env\Scripts\python.exe" -m pip install -e .
& "$env:USERPROFILE\mesh_env\Scripts\python.exe" -m mesh_aop.check_env --auto
```

Run every pipeline command the same way — full path to `python.exe`, module form, which bypasses the blocked launcher:

```powershell
& "$env:USERPROFILE\mesh_env\Scripts\python.exe" -m mesh_aop.cli --step all --interactive
```

*(Prefer typing plain `python` / `mesh-pipeline`? Run `Set-ExecutionPolicy -ExecutionPolicy Bypass -Scope Process -Force` then `& "$env:USERPROFILE\mesh_env\Scripts\Activate.ps1"` to activate. But the full-path form above needs no policy change and is unaffected by corporate lockdowns.)*

*Note: The environment checker resolves the `community` / `python-louvain` namespace collision automatically — see Troubleshooting for the other Windows-specific errors.*

### 3. Alternative Installation: via Mamba/Micromamba

This approach uses the provided `environment.yml` to fetch pre-compiled binaries via `conda-forge`, bypassing common OS-level compilation errors.

```bash
# Navigate to the repository root
cd path/Mesh-Network-Analysis-Main-Library

# Create the environment and resolve dependencies
mamba env create -f environment.yml

# Activate the environment
mamba activate mesh_aop_network

# Run the verification script to ensure successful installation and provision OS-level rendering libraries
mesh-check-env --auto

```

### 4. Verification

Regardless of the installation method chosen, verify the package is correctly linked to your PATH by calling the command line interface:

```bash
mesh-pipeline --help
mesh-pipeline --version

```

---

## Execution Guide

The pipeline is entirely modular and controlled via a terminal interface. Configuration is handled by an interactive command-line wizard, allowing users to modify runtime parameters safely without touching source code.

> **Invocation by platform.** The examples below use the `mesh-pipeline` command, which works on **macOS/Linux** (and on Windows after activating the venv). On **Windows**, if activation or the `.exe` launcher is blocked, use the equivalent module form with the venv's Python by full path — it behaves identically:
> ```powershell
> & "$env:USERPROFILE\mesh_env\Scripts\python.exe" -m mesh_aop.cli --step all --interactive
> ```
> i.e. replace `mesh-pipeline` with `& "$env:USERPROFILE\mesh_env\Scripts\python.exe" -m mesh_aop.cli` in any command. Always run from the project root so it finds `mesh_config.json` and the `data/` folders.

### CLI Flags

| Flag | Description |
|------|-------------|
| `--step <name>` | Which pipeline segment to run (`all`, `process`, `data_ops`, `network`, `secondary`, `viz`, `benchmark`). Defaults to `all`. |
| `--interactive` | Launches the interactive wizard before execution. |
| `--config <path>` | Path to a custom config JSON. Defaults to `mesh_config.json` in the current directory. |
| `--readme` | Opens this documentation file in your default OS viewer. |
| `-v` / `--version` | Prints the installed package version and exits. |

### Running the Complete Pipeline

To construct a network from the ground up, execute the `all` step. The `--interactive` flag invokes the wizard.

```bash
mesh-pipeline --step all --interactive

```

### Running Individual Modules

If upstream dependencies are already built, specific modules can be executed in isolation.

* **Step 0 & 1:** `mesh-pipeline --step process --interactive` (Database Compilation & MeSH processing)
* **Step 2:** `mesh-pipeline --step data_ops --interactive` (Entrez API Collection)
* **Step 3:** `mesh-pipeline --step network --interactive` (Topology & Filtering)
* **Step 3.5:** `mesh-pipeline --step secondary --interactive` (Targeted Export Analysis)
* **Step 4:** `mesh-pipeline --step viz --interactive` (Biological Figure Generation)
* **Step 5:** `mesh-pipeline --step benchmark` (Ground-Truth Validation & Performance Benchmarking)

---

## Configuration Wizard Parameter Glossary

The interactive wizard is categorized into discrete blocks. Below is the scientific and computational rationale for each tunable variable.

### 1. Control Flags & Directories

* **Use Reference Data:** If `True`, bypasses API downloading and utilizes a static, pre-curated reference dataset.
* **Pause for Annotation (AFK Mode):** If `False` (Default), the pipeline operates in AFK Mode. It will run uninterrupted from start to finish, automatically assigning 'Unassigned' to all biological levels. If `True`, the pipeline will safely pause after Step 3 to allow the user to manually annotate the network before rendering the final biological visualizations.
* **Custom Prefix:** The naming convention prepended to all output files (e.g., `DAC_Mesh`).

### 2. Master Database Status (Step 0 ETL)

The wizard actively probes your local Master SQLite Database for corruption, completion status, and age (checking if a new yearly baseline is available).

* **Compile PubMed Baseline / Daily Updates:** Initiates a multi-core MapReduce extraction of the NLM XMLs into the local cache.

### 3. NCBI Credentials

* **Entrez Email & API Key:** Registration with NCBI allows for 10 API requests per second. Without a key, requests are hard-limited to 3 per second, roughly tripling Step 2 processing time. Sign in or create a free account at the [NCBI/NLM account portal](https://account.ncbi.nlm.nih.gov/) to obtain a key (found under *Account Settings → API Key Management*).

### 4. Search Parameters

* **MeSH Search Term:** The PubMed query used to retrieve the base (P0) cohort. This string is sent **verbatim** to NCBI as the P0 `esearch` query, so any valid PubMed syntax is accepted — MeSH tags (`[Mesh]`, `[Majr]`), field tags (`[tiab]`, `[au]`), quoted phrases, and the boolean operators `AND` / `OR` / `NOT` with parentheses for grouping. Examples:
  * Single MeSH heading — `Dermatitis, Allergic Contact [Mesh]`
  * MeSH plus a free-text phrase — `"Dermatitis, Allergic Contact"[Mesh] OR "skin sensitization"[tiab]`
  * Boolean combination — `(Dermatitis, Allergic Contact[Mesh]) AND (Haptens[Mesh])`

  A few rules keep the query valid and the downstream steps consistent:
  * Spell MeSH headings exactly as they appear in the [MeSH Browser](https://meshb.nlm.nih.gov/) — internal commas are fine here (the term is a single value, not a CSV field). A misspelled heading silently degrades to a free-text search.
  * Wrap any multi-word heading or phrase in double quotes, e.g. `"Dermatitis, Allergic Contact"[Mesh]`.
  * Boolean operators must be **uppercase** (`AND`, `OR`, `NOT`); lowercase is treated as an ordinary search word.
  * The benchmark's `primary_node` defaults to this term's MeSH node. If you use a compound or free-text query rather than a single MeSH heading, set `benchmark.primary_node` explicitly so the "naive query" baseline and the topology-exclusive split resolve to a real node.
* **Start / End Date:** Constrains the temporal boundaries of the initial P0 PubMed search.
* **Citation Generations:** Controls the depth of the citation scrape. The value counts *levels including P0*, matching the wizard label (`P0=1, P0+G1=2, ...`).
  * `1`: Only the base parental generation (P0) articles.
  * `2`: P0 articles + all articles they cite + all articles that cite them (G1).
  * `3`: The above + a further citation hop (G2), and so on.
  * *Warning:* Values $\ge 2$ result in exponential data growth.
* **Update MeSH Support Files:** If `True`, forces re-extraction of the MeSH term list and stop-word set from the XML file even if cached outputs already exist.



### 5. Analysis Parameters

* **Calculate Full Centrality (Boolean):** Controls **how betweenness is computed, not whether centrality is computed at all**. If `True`, betweenness is calculated **exactly** over every node pair. If `False` (Default), it is **estimated** from a sample of `betweenness_k_samples` source nodes, which is dramatically faster and bounds memory on large graphs.
  * Eigenvector and PageRank centrality are computed **either way** and are unaffected by this flag, so Article Relevance Scores (ARS) and Mean Relevancy Scores (MRS) are produced normally in both modes. The only difference is that betweenness — and therefore the betweenness-weighted ARS/MRS — is a sampled estimate rather than an exact value. Report this in your methods if you leave it `False`.
* **Betweenness K-Samples:** Heuristic sampling limit for Centrality calculation. Lower values increase speed but reduce precision.
* **Context Start / End Date:** Temporal constraints applied exclusively to Step 3 Mean Relevancy Scoring, allowing the simulation of historical network states.
* **Random Seed:** Integer seed passed to NetworkX and scikit-learn for reproducible t-SNE projections and community detection. Default: `42`.
* **Eigenvector Max Iterations / Tolerance:** Power-iteration convergence controls for Eigenvector centrality. Increase `eigenvector_max_iter` (default `1000`) if convergence warnings appear on very large or sparse graphs.

### 6. Network & Simulation Parameters

* **Lambda Value:** The distance-decay factor for generational node weighting: a node's generation weight is $W = e^{-\lambda d}$, where $d$ is the citation-generation distance from the P0 seed set ($d = 0$ for P0, $1$ for G1, …). Larger $\lambda$ penalizes distant generations more steeply.
* **Node Weight Factors:** A four-component dictionary combined into a single per-node `adjusted_node_weight` attribute recorded on the final network. It is a **user-tunable node-importance metric, not a driver of the GLF/SA filtering** — consensus filtering is driven by edge co-occurrence strength, independent of these factors. The attribute is reported on the final graph (and offered as one of the candidate weightings in the validation step), so it is available as a base metric for your own downstream analysis. The four keys and their defaults are:

  | Key | Default | Role |
  |-----|---------|------|
  | `centrality` | `0.45` | Combined betweenness + eigenvector centrality contribution |
  | `article_rank` | `0.15` | PageRank-style article influence |
  | `rank_median_cit` | `0.20` | Median incoming citations across linked articles |
  | `rank_total_cit` | `0.20` | Total incoming citation volume |

  The four values must sum to `1.0`.

* **Target Edges (`target_num_edges`):** The target subgraph size handed to the GLF and SA optimizers. Each selects a subgraph of about this many edges; their intersection (the consensus, reduced to its LCC) is typically somewhat smaller. *(Distinct from the `target_edges` in Secondary Analysis.)*
* **GLF / SA Iterations:** Monte Carlo search and thermal cooling steps for the optimization heuristics.
* **SA Temperature Start:** The initial thermal energy of the Simulated Annealing system (default `5000.0`). Higher values allow larger disruptive jumps early in the search.
* **SA Cooling Rate:** The multiplicative decay applied to temperature each iteration (default `0.999995`). Values closer to `1.0` cool more slowly and explore more widely at the cost of run time.

### 7. Secondary Analysis Parameters

Executes highly targeted queries against the finalized network to extract specific publications for manual review.

* **Export Top Articles:** If `True` (default), always exports the highest-scoring network-wide articles at the end of Step 3 without requiring `--step secondary`.
* **Export Limit:** Maximum number of articles returned per query (default `500`).
* **Exclude Review Articles:** Filters out broad review articles to isolate primary literature.
* **Target Nodes:** Evaluates the literature density of specific nodes. **Must be semicolon-separated** (e.g., `Skin; T-Lymphocytes`).
* **Target Edges (`target_edges`):** Evaluates literature specifically linking two concepts. Separate the two node names with a **space-padded hyphen** (` - `), and separate multiple edge queries with semicolons: `NodeA - NodeB; NodeC - NodeD`. *(Distinct from the numeric `target_num_edges`)*
* **Sort Metric (`sort_metric`):** Selects the ranking engine — `Linear` (compensatory weighted average) or `F1` (penalizing harmonic mean). Default `F1`.

  **1. Linear (Weighted Additive Model)**

  * **Mechanism:** Calculates the arithmetic weighted average of the normalized Article Relevance Score (ARS) and the normalized citation rate.
  * **Behavior:** Compensatory. A high score in one metric can offset a low score in the other, governed by the user-defined `linear_weight_ars`.
  * **Formula:** $(ARS \times w) + (Cit \times (1 - w))$

  **2. F1 (Harmonic Mean)**

  * **Mechanism:** Calculates the strict harmonic mean of the normalized ARS and the normalized citation rate.
  * **Behavior:** Penalizing. The final output skews heavily toward the lower of the two input values. An article must possess both high topological relevance and high community impact to achieve a high score.
  * **Formula:** $2 \times \frac{ARS \times Cit}{ARS + Cit}$
* **ARS Weight (`linear_weight_ars`):** Weight given to the normalized ARS in the `Linear` metric (0–1.0); the remaining `1 − w` weights the normalized **citation rate** (citations per year, which corrects the age bias that otherwise favors older papers purely for having accumulated citations longer). Default `0.5`.
* **Compare Multiple Networks (`compare_networks`):** Off by default. When enabled, secondary analysis runs a **node-overlap comparison** across a set of saved networks you list in **Networks to Compare** (`comparison_networks`) — a comma-separated, quote-wrapped list such as `"Ex_Graph_1.json","Ex_Graph_2.graphml"`. Bare names are resolved against the processed data folder (or `data/reference_processed/` when `Use Reference Data` is on); explicit paths and non-JSON formats (`.graphml`, `.gml`, `.gexf`, …) are also accepted. Missing files trigger a warning naming where the pipeline looked. It writes a membership matrix, a pairwise intersection/Jaccard table, and an overlap figure to `results/`.

### 8. Benchmark Parameters

Controls the optional `--step benchmark` evaluation (see the **Validation & Benchmarking** section below). Configured under the `benchmark` block of `mesh_config.json`.

* **`ground_truth_csv`:** Filename of your ground-truth PMID set inside the active raw directory. Leave **empty** (default) to auto-detect a recognized filename; set it to pin a specific file.
* **`negative_control_csv`:** Optional filename of an *unrelated* ground-truth set used as a specificity check (it should score near random). Empty disables the control.
* **`primary_node`:** The base disease/seed node used to separate "naive" articles (those carrying the primary node) from topology-exclusive hits. Defaults to the search term's MeSH node.
* **`n_boot` / `n_perm`:** Bootstrap resamples and permutation iterations behind the confidence intervals and the random-ranking null. Default `25` each. Each bootstrap resample re-ranks the full article pool, so runtime scales linearly with both this value and the size of that pool. Measured on a ~9-million-article pool:

  | `n_boot` / `n_perm` | Approx. runtime |
  | --- | --- |
  | `25` *(default)* | ~15 min |
  | `50` | ~30 min |
  | `100` | ~1 hour |
  | `200` | ~2 hours |

  Increasing the count tightens the bootstrap confidence intervals and stabilizes the permutation p-value, but runtime grows proportionally; decreasing it runs faster at the cost of noisier, less reliable estimates. The default (`25`) is a practical balance. Interval precision is ultimately bounded by the number of ground-truth positives rather than by `n_boot`, so values well beyond `200` give diminishing returns.
* **`run_ground_truth_analysis` (Boolean):** Master switch for the whole step. **If unset it follows `Use Reference Data`** — on when you are running against the bundled reference corpus (which the bundled ground truth describes), off when you are running your own data (where that ground truth would not apply). Set it explicitly to override. The interactive wizard prompts for it first in this section.
  When enabled, this **also** causes betweenness **and** PageRank to be recomputed on the filtered consensus subgraph, in addition to the whole-corpus versions, and scored as extra benchmark scorers (`betweenness (subgraph)`, `pagerank (subgraph)`). Whole-corpus centrality measures importance across the entire literature, where generic high-degree MeSH terms dominate; subgraph centrality measures it within the curated concept space. Computing both for both algorithms makes centrality **scope** and centrality **type** a full 2×2 rather than confounding them. Note the subgraph betweenness is computed **exactly** (the subgraph is small), unlike the whole-graph estimate. The `n_seeds` baseline is a uniform weight of 1 per node and therefore scope-invariant — it is the single control for all four. Because this decision changes which centralities are computed, it is read at the **network** step, so the wizard offers the flag during `--step network` as well; the benchmark scores whichever `score_*` columns relevance actually produced.
* **`run_network_validation` (Boolean):** If `True` (default), also runs the node/edge convergent validation described under **Validation & Benchmarking**. Set `False` to run only the article ranking benchmark.
* **`run_projection_comparison` (Boolean):** If `True` (default), also runs the article-scoring **projection comparison** — with the node seed fixed, it scores the alternative ways of turning node weights into an article score (normalised ARS, unnormalised sum, MRS-weighted, bipartite-reinforced, BM25, uniform, random, naive query) by BEDROC across three frames with positives-only bootstrap CIs, writing `{prefix}_projection_comparison.csv` and a figure to `results/validation/`. Set `False` to skip it.
* **`network_validation_weight_key`:** Node attribute used as the "network weight" when correlating a node's ground-truth prominence against its importance. Default `MRS_pagerank_centrality`; set to `MRS_betweenness_centrality` to compare against the betweenness weighting instead. The wizard also lets you choose any raw or MRS centrality (betweenness / pagerank / eigenvector, whole-graph or subgraph).
* **`min_articles_per_node`:** Minimum number of ground-truth articles a term must appear in to become a node of the ground-truth co-occurrence network (default `2`, which suppresses singleton noise).
* **`background_pool_size`:** Number of randomly sampled articles used to estimate corpus base rates and to build the permutation nulls (default `50000`). Larger is more precise but slower.

---

## The AOP Annotation Workflow (Biological Strata)

**Critical Concept:** The pipeline algorithms can identify the statistical relationships between MeSH terms, but they cannot automatically determine if a term represents a "Molecular", "Cellular", or "Tissue" level event. To generate accurate biological Sankey flows, human context must be provided.

To streamline this, the pipeline utilizes a **Semicolon-Delimited Master Dictionary System**.
*Note: Because MeSH terms frequently contain commas (e.g., `Dermatitis, Allergic Contact`), standard CSV comma-delimiters will corrupt the data. All annotation files in this pipeline exclusively use semicolons (`;`).*

### How to Annotate Your Network:

1. **Enable Pausing:** In the Wizard, ensure `Pause for Annotation` is set to `True`.
2. **Run the Pipeline:** Let the pipeline run. It will execute Steps 0 through 3 (Network Construction), output AOP-independent figures (like distributions and convergence trajectories), and then pause — after the node set is final but before any biological Sankey/alluvial figures are rendered (Step 4).
3. **Open the Run Template:** Navigate to the `results/` directory and open your run-specific template: `[PREFIX]_run_annotations.csv`.
4. **Assign Strata:** This file contains every surviving node in your network. It automatically pulls any known assignments from your Master Dictionary. For any term listed as `Unassigned`, replace the text with one of the following 7 strata. *(Unsure how to categorize a term? Look up its official scope note and hierarchy in the [NLM MeSH Browser](https://meshb.nlm.nih.gov/) — the official MeSH description is the best guide for determining the correct AOP stratum.)*
* `Stressor` - External stimuli that initiate a biological reaction (e.g., `UV Rays`, `Chemicals`)
* `Molecular` - Gene, protein, or receptor level events (e.g., `Receptors, Antigen`)
* `Cellular` - Cellular level events (e.g., `Chemotaxis`, `Apoptosis`, `T-Cell`)
* `Tissue` - Events localized to a subsection of tissue (e.g., `Necrosis`, `Skin Absorption`)
* `Organ` - Organ level names or events (e.g., `Liver`, `Skin`)
* `Adverse Outcome` - High-level disease outcomes and disease (e.g., `Drug Hypersensitivity`)
* `Uncategorized` - Broad biological or non-biological terms that do not fall into a distinct AOP stratum.


5. **Save the File:** Save the file, ensuring it remains **semicolon-delimited**.
6. **Resume the Pipeline:** Run `mesh-pipeline --step viz`. The pipeline will ask if you want to sync your new annotations to the Master Dictionary for future runs, and then generate the final biological figures.

### What "syncing to the Master Dictionary" does

When you resume (`--step viz`), the pipeline offers to write your run's annotations back to the **Master Dictionary** (`data/raw/aop_annotations_master.csv`). Saying yes merges every `term → stratum` assignment from this run into that master file:

* The assignments become **persistent in your copy of the repository** — every future run pre-fills these terms automatically, so you never re-annotate the same MeSH heading twice, and you can re-export your curated strata at any time from the master file.
* The Master Dictionary ships **partially pre-seeded**: a subset of terms relevant to the bundled Dermatitis, Allergic Contact (DAC / AOP 40) query are already annotated from prior curation. **Most MeSH terms are not annotated** and will appear as `Unassigned` in your run template — you must assign those yourself for the biological Sankey/alluvial flows to be meaningful. Anything left `Unassigned` (or run in AFK mode) is treated as `Uncategorized`.

---

## Output Artifacts

Upon successful completion of the pipeline, the following critical files are generated:

### Data & Network Artifacts (`data/processed/`)

* `master_mesh_database.db`: A persistent offline cache of all PubMed IDs and MeSH annotations.
* `*_pmids.db` / `*_cleaned_pmids.db`: The raw and de-duplicated SQLite schemas holding the hierarchical citation linkage of your extracted articles.
* `mesh_terms.csv`: The MeSH vocabulary extracted from the XML in Step 1.
* `*_full_network_data.json`: The raw, unfiltered NetworkX graphical representation.
* `*_glf_optimal_subgraph.json` / `*_sa_optimal_subgraph.json`: The two independent optimizer solutions whose intersection forms the consensus subgraph.
* `*_optimization_history.json`: The GLF/SA score trajectory (feeds the optimization-trajectory figure).
* `*_consensus_lcc_network.json`: The optimized graph containing the intersection of the GLF and SA algorithms, reduced to its Largest Connected Component (LCC).
* `*_final_network_with_relevance.json`: The fully annotated terminal graph, populated with semantic Mean Relevancy Scores (MRS) and Louvain Community classifications.
* `*_mean_relevancy.db`: The per-article Mean Relevancy Scoring database — the `score_*` columns behind ARS/MRS and the secondary/benchmark steps.

### Analytical Exports (`results/`)

* `*_export.xlsx`: A tabular summary of the final nodes and edges.
* `*_Top_Network_Articles.csv`: The highest-scoring primary literature driving the network's structure.
* `*_run_annotations.csv`: The per-run biological-strata template you edit during the annotation pause.
* `*_network_overlap_membership.csv` / `*_network_overlap_matrix.csv` / `*_Network_Overlap.png`: node-overlap comparison across the networks you named — produced only when **Compare Multiple Networks** is enabled.
* **Benchmark & validation outputs** (all under `results/benchmark/`): everything the `--step benchmark` step produces is consolidated in one folder — the article-ranking benchmark, the node/edge convergent validation with its figures, and the node-weighting/projection report (nested in `results/benchmark/validation/`). See the **Validation & Benchmarking** section below for the complete list (`*_benchmark_results.json`, `*_benchmark_enrichment.png`, `*_benchmark_quarantined_pmids.csv`, `*_gt_network_validation.xlsx`, `*_gt_cooccurrence_network.json`, the `*_GT_*.png` figures, and `validation/*_validation_report.xlsx`/`.html` + `validation/*_projection_comparison.csv`).
* `logs/`: run logs and failed-fetch records.
* **Figures (`results/figures/`)**:
* **Figure 1:** Edge weight distribution (Power law analysis) to assess network topology.
* **Figure 2:** Optimization Trajectory (GLF vs SA convergence).
* **Figure 3:** Community Composition bar charts detailing biological strata makeup per cluster.
* **Figure 4:** MRS Centrality correlations (Betweenness vs PageRank).
* **Figure 5:** t-SNE projection of the network colored by Louvain community.
* **Figure 6:** AOP Alluvial/Sankey flows (The primary visualization connecting Stressors to Outcomes). *(Note: Provided as interactive `.html` files for deep pathway inspection).*
* **Figure 7:** Dumbbell plots assessing shift in topological vs semantic relevance.



---

## Validation & Benchmarking

The `--step benchmark` module evaluates how well the network's per-article relevance scores recover an external, curated set of "ground-truth" PMIDs (e.g. the bibliography of an authoritative OECD Adverse Outcome Pathway document). It is designed for the realistic situation here — a few dozen known-relevant papers hidden in a multi-million article candidate pool — and follows three principles: **validate the ground truth first, use metrics robust to incomplete relevance judgments, and quantify uncertainty against baselines.**

```bash
mesh-pipeline --step benchmark
```

### The Bundled OECD Ground Truth

The repository ships a curated positive set derived from the reference bibliography of the **OECD Adverse Outcome Pathway for skin sensitisation (AOP 40)**. It is the default target of `benchmark.ground_truth_csv`.

| File | Location | Contents |
| --- | --- | --- |
| `oecd_ground_truth_curated.xlsx` | `data/reference_processed/` | **96 curated positives** — primary research articles with resolved PMIDs |
| `oecd_ground_truth_exclusions.xlsx` | `data/reference_processed/` | **22 excluded citations**, each with an explicit reason (book chapters, OECD guidance documents, and records that are not indexed primary articles) |
| `oecd_ground_truth.bib` | `data/reference_processed/` | Typed BibTeX of the parsed citations (`@article`, `@incollection`, `@techreport`, `@misc`) |
| `oecd_resolved_citations.csv` | `data/reference_raw/` | The raw citation → PMID resolution table the curation was built from |

**How it was built.** Citations were parsed from the OECD reference list, resolved to PMIDs via Entrez `esearch` on author + year + volume + first page (rather than free-text matching, which mis-resolves), and then manually adjudicated. Non-article references and records lacking MeSH indexing were moved to the exclusions file rather than silently dropped, so the **96 + 22 = 118** entries fully account for the source bibliography — the exclusions are auditable, not hidden losses.

**Why the exclusions matter.** Anything without MeSH indexing is unreachable by a MeSH co-occurrence network *by construction*. Keeping those citations in the positive set would depress every metric for a reason unrelated to ranking quality, so they are excluded explicitly and counted separately.

### Supplying Your Own Ground Truth

The benchmark reads its ground-truth file from the **active raw directory**, which follows the `Use Reference Data` flag:

* `Use Reference Data = True`  → reads from `data/reference_raw/` (the bundled OECD set).
* `Use Reference Data = False` → reads from **`data/raw/`** — drop your own file here.

To benchmark against your own citations, save a file into `data/raw/` using one of these recognized names (searched in order; no config edit required):

```text
ground_truth_pmids.csv
ground_truth_pmids.txt
ground_truth.csv
ground_truth.txt
```

Alternatively, pin an explicit filename via the `benchmark.ground_truth_csv` config key.

**Accepted formats** (auto-detected, with delimiter and encoding sniffing):

* An **`.xlsx` workbook** with a `PMID` column (the format of the bundled curated set).
* A CSV/TSV with a **`PMID`** column (aliases such as `pmids`, `pubmed_id`, `id` are recognized) and an optional reference/citation column.
* The bundled `Raw_Reference;PMID` semicolon-delimited format.
* A **single column of PMIDs**, with or without a header.
* A plain **`.txt` list** with one PMID per line.

If no ground-truth file is found, the step prints the exact directory it searched, the accepted filenames, and the accepted formats.

### What It Reports

1. **Ground-truth validation.** PMIDs are digit-normalized (fixing the common `19033392.0` float-coercion mismatch), de-duplicated, and split into resolved vs. unresolved (`NOT_FOUND`). Each resolved PMID is sanity-checked for plausibility: an *absolute ceiling* test rejects PMIDs too large to exist yet, and a *chronological* test flags post-1990 papers whose PMID is implausibly large for the cited year (a hallmark of citation→PMID mis-resolution). **Flagged rows are quarantined to a review CSV — never deleted automatically.** When a reference/citation column is present the year check is active; for bare PMID lists only the absolute-ceiling test applies.
2. **Coverage vs. ranking (kept separate).** *Retrievability ceiling* — how many ground-truth PMIDs are even present in the candidate pool — is reported independently of *ranking quality*, so a low score is correctly attributed to either missing coverage or poor ranking rather than conflated.
3. **Ranking metrics suited to extreme imbalance & incomplete judgments.** Primary metrics depend only on the positions of known positives and never assume unjudged articles are negatives: **Recall@K, MAP, R-precision, NDCG, Enrichment Factor (top 1 %/5 %)** and **BEDROC** (early recognition). ROC-AUC and PR-AUC are reported but explicitly demoted to secondary diagnostics, because under ~10⁻⁵ prevalence ROC-AUC is over-optimistic and PR-AUC is biased low by treating unjudged articles as negatives.
4. **Uncertainty and baselines.** Every headline metric carries a **bootstrap confidence interval**, and a **permutation null** reports lift-over-random with an empirical p-value. A structural baseline (number of bridged network nodes) and an optional **negative-control** ground-truth set (which should score near random) guard against the network merely surfacing generically popular papers.

### Outputs

Written to `results/benchmark/`:

* `*_benchmark_results.json` — all metrics, confidence intervals, lift, and coverage figures, with the run's seed and configuration for reproducibility.
* `*_benchmark_quarantined_pmids.csv` — ground-truth rows flagged as implausible, for manual adjudication.
* `*_benchmark_enrichment.png` — cumulative recall (enrichment) curve: ground-truth recovery vs. fraction of the pool screened.
* `validation/*_validation_report.xlsx` / `.html` + figures — the consolidated **node-weighting** comparison across four framings (corpus / within-query / outside / hybrid) with BEDROC and paired bootstraps (from `run_network_validation`).
* `validation/*_projection_comparison.csv` + figure — the article-scoring **projection** comparison (normalised ARS vs unnormalised sum vs MRS-weighted vs bipartite vs BM25/baselines) by BEDROC across three frames (from `run_projection_comparison`).

> **Note on incomplete judgments:** the non-ground-truth articles are *unjudged*, not confirmed negatives. This is why recall- and rank-based metrics (which only need positive positions) are the headline numbers, while precision/PR-AUC are interpreted with caution. For the most rigorous citation→PMID verification, cross-check publication years online via NCBI Entrez (not enabled by default, as it requires network access).

### Node/Edge Convergent Validation

`--step benchmark` also runs a **structural** check (`gt_network_validation.py`) that asks a different question from the article ranking above: *is the network's own vocabulary and wiring independently reproduced by the ground-truth literature?* A network could rank articles well simply by naming popular terms — recovering its **edges** from an independent corpus is far harder to achieve by chance.

The ground-truth articles' MeSH terms are pulled from the master database and filtered through the **same node-eligibility stop-word list used to build the network**, so stop-listed headings (check tags, organisms, geographics) can never appear as spurious "near misses" — they were never eligible to be nodes in the first place. An audit sheet records exactly which stop terms were removed and asserts that none leaked into the node or miss lists.

Both levels are calibrated against a **permutation null of random article draws of the same size** — essential, because two corpora about the same topic will overlap somewhat by construction, so raw overlap alone is not evidence:

* **Node level** — how many network nodes are attested in the ground truth, their enrichment over corpus base rates, and the Spearman correlation between a node's ground-truth prominence and its network weight.
* **Edge level** — how many of the network's edges reappear as ground-truth co-occurrences.
* **Misses** — ground-truth-frequent eligible terms that are *not* network nodes, which localize where the network under-covers the domain.
* **Node-weighting comparison** — scores every available node weighting (each `MRS_*`, each raw centrality, `adjusted_node_weight`) against external ground-truth prominence, to test whether the MRS transformation adds node-importance information **over** the centrality it is built from. Because an MRS is derived from its centrality the two are collinear, so a plain comparison of their correlations is not sufficient; the table therefore reports a **partial correlation controlling for the raw centrality** and a **paired bootstrap CI on the difference**. A difference interval spanning zero means the transformation is not distinguishable from the centrality underneath it. Three criteria are reported side by side: correlation with GT article frequency, correlation with **enrichment** (base-rate corrected, which removes any residual frequency component and is therefore the stricter test), and an AUC for ranking attested nodes above unattested ones across *all* nodes.

It runs *before* the ranking benchmark (it takes minutes rather than tens of minutes), so a long benchmark never blocks the structural result. Disable it with `benchmark.run_network_validation = false`.

**Outputs** (to `results/benchmark/`):

* `*_gt_network_validation.xlsx` — sheets: `summary`, `stopword_audit`, `network_nodes`, `GT_terms`, `GT_misses`, `network_edge_validation`, `node_weighting_comparison`
* `*_gt_cooccurrence_network.json` — the ground-truth co-occurrence graph in **Cytoscape.js** format, every node and edge annotated with its overlap status (`shared` vs `GT_only_miss`, `recovered` vs `GT_only`), enrichment, and counts
* `*_GT_Node_Validation.png`, `*_GT_Cooccurrence_Network.png`, `*_GT_Permutation_Nulls.png`

---

## Jupyter Notebook Interface

Every module in `src/mesh_aop/` has a corresponding Jupyter notebook in `src/mesh_aop_notebooks/`. These notebooks mirror the source modules cell-by-cell and are intended for:

* **Interactive exploration** — step through the pipeline one cell at a time and inspect intermediate data structures.
* **Prototyping** — experiment with individual functions (e.g., tweak GLF parameters and re-run only the filtering step) without triggering the full CLI orchestration.
* **Debugging** — isolate a specific module and inspect its inputs and outputs in a notebook environment.

Each notebook contains its own `mesh_config.json` and `environment.yml` references so it can be run independently from the `src/mesh_aop_notebooks/` directory if needed.

---

## Programmatic API Usage

The package exposes a clean Python API through its `__init__.py`. All pipeline functions can be imported and called directly without using the CLI, which is useful for embedding the analysis within a larger workflow or Jupyter-based research pipeline.

```python
from mesh_aop import (
    MeshConfig,
    process_raw_mesh_data,
    run_initial_data_collection,
    run_network_construction,
    run_consensus_filtering_and_lcc,
    run_community_detection,
    run_mean_relevancy_scoring,
    get_top_network_articles,
    plot_sankey_alluvial,
)

# Load config (merges factory defaults with your local mesh_config.json)
config = MeshConfig(config_path="mesh_config.json")

# Step 1 – Extract MeSH terms from XML
process_raw_mesh_data(
    xml_file=config.files['mesh_xml'],
    output_csv=config.files['mesh_terms_csv'],
    output_py=config.files['mesh_stopwords_py'],
)

# Step 2 – Collect articles from NCBI Entrez
run_initial_data_collection(
    search_term_param=config.get('search_parameters', 'search_term'),
    start_date_str=config.get('search_parameters', 'start_date'),
    end_date_str=config.get('search_parameters', 'end_date'),
    generations_n_param=config.get('search_parameters', 'generations_n'),
    db_path=config.files['pmids_db'],
    entrez_email=config.get('credentials', 'entrez_email'),
    entrez_api_key=config.get('credentials', 'entrez_api_key'),
)

# Step 3 – Build and filter the network
run_network_construction(db_path_param=config.files['cleaned_db'],
                         output_json_path=config.files['full_network'], ...)
run_consensus_filtering_and_lcc(...)
run_community_detection(network_file_path=config.files['consensus_lcc'], ...)
run_mean_relevancy_scoring(...)
```

The full list of exported symbols is defined in `src/mesh_aop/__init__.py`.

---

## Troubleshooting

### Python 3.13: NumPy tries to compile and fails

Symptom: during `pip install -e .` you see meson/`Unknown compiler` errors building NumPy, e.g. `Could not find ... vswhere.exe`, or `ResolutionImpossible` mentioning `numpy`.

Cause: an earlier revision depended on the `node2vec` package, which pinned `numpy<2.0`; that NumPy ships **no prebuilt wheel for Python 3.13**, so pip fell back to compiling it from C source — which needs an MSVC compiler you likely don't have.

Fix: **this no longer applies.** The Node2Vec embedding is vendored in `node2vec_embedding.py` and the package dependency was removed, so Python 3.13 installs cleanly. If you still hit this, you are installing from an out-of-date checkout or a stale copy of `pyproject.toml`: pull the current revision, delete the virtual environment, and reinstall.

### Windows `MAX_PATH` / "No such file or directory" during install

Symptom: `pip install -e .` aborts with `OSError: [Errno 2] No such file or directory: '...\statsmodels\tsa\vector_ar\tests\JMulTi_results\...txt'` and a hint about *"Windows Long Path support"*.

Cause: the Windows 260-character path limit. A deeply nested project location (especially under `OneDrive - <Org>\Documents\...`) plus a venv plus `statsmodels`' long test filenames exceeds 260 characters, so the file write fails and pip rolls back the whole install.

Fix (no admin needed): create the venv at a **short path outside** the project, then install from the repo root — only the deep dependency files need the short location; the editable package just links back to your source:

```powershell
py -3.12 -m venv "$env:USERPROFILE\mesh_env"     # e.g. C:\Users\you\mesh_env
& "$env:USERPROFILE\mesh_env\Scripts\python.exe" -m pip install -e .
```

(Alternative, if you have admin rights: enable long paths via `HKLM\SYSTEM\CurrentControlSet\Control\FileSystem\LongPathsEnabled = 1` and reboot.)

### "Activate.ps1 cannot be loaded ... running scripts is disabled"

Symptom: activating the virtual environment fails with `running scripts is disabled on this system` (a `PSSecurityException` / `UnauthorizedAccess`).

Cause: Windows PowerShell's **execution policy** (often `Restricted` by default) blocks the `Activate.ps1` script. This only blocks the activation *script* — `python.exe` itself is unaffected. `mesh-check-env` detects this and prints the same guidance.

Fix — enable local scripts once (no admin needed; works in all future windows):

```powershell
Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned
& "$env:USERPROFILE\mesh_env\Scripts\Activate.ps1"
```

If a corporate policy blocks even `CurrentUser`, set it per-window instead: `Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass -Force`. Or skip activation entirely and call the venv Python by full path (needs no policy change):

```powershell
& "$env:USERPROFILE\mesh_env\Scripts\python.exe" -m mesh_aop.cli --step all --interactive
```

### "Access is denied" running `mesh-pipeline` / `mesh-check-env`

Symptom: the console-script launchers fail with *Access is denied* / `ApplicationFailedException`, even though `python.exe` itself works.

Cause: corporate security/EDR software on managed machines often blocks freshly-created, unsigned `.exe` launcher stubs that pip generates in the venv's `Scripts` folder.

Fix: invoke the tools as modules, which bypasses the launcher exes entirely:

```powershell
python -m mesh_aop.cli --step all --interactive
python -m mesh_aop.check_env --auto
```

### `community` / `python-louvain` Namespace Collision

The `python-louvain` package installs itself under the name `community`, which conflicts with an unrelated package also named `community` on PyPI. If you see `ImportError: cannot import name 'community_louvain' from 'community'`, run:

```bash
pip uninstall community python-louvain -y
pip install python-louvain
```

The `mesh-check-env --auto` script detects and resolves this automatically.

### RAM Exhaustion on Large Networks

Networks with `generations_n >= 2` can load several million records into RAM during centrality calculation. If the process is killed by the OS:

1. Set `calculate_full_centrality: false` in `mesh_config.json` so betweenness is estimated from a node sample rather than computed exactly.
2. Reduce `betweenness_k_samples` (default `1000`) to lower the heuristic sampling budget.
3. Reduce `target_num_edges` so GLF/SA operate on a smaller subgraph.

Setting `calculate_full_centrality: false` does **not** disable centrality or prevent MRS from being calculated — eigenvector and PageRank are computed regardless, and only betweenness becomes a sampled estimate (see the Analysis Parameters section above).

### SQLite Lock Errors on Network-Attached Storage

The `baseline_manager` uses a verified safe-transfer architecture specifically to handle write failures on NAS and cloud-synced directories (OneDrive, Dropbox, etc.). If you still encounter lock errors, ensure no other process (e.g., a cloud-sync agent) has the `.db` file open, then re-run `mesh-pipeline --step process`.

### Convergence Warnings for Eigenvector Centrality

On sparse or disconnected graphs, the power-iteration solver may not converge within the default 1000 iterations. Increase `eigenvector_max_iter` in `mesh_config.json` or lower `eigenvector_tol` to relax the stopping criterion.

---

## Citation

If you use this code or methodology, please cite:

> *[To be updated after publication]*
> [](https://doi.org/10.5281/zenodo.18662959)

## License

[MIT License]
