# Medical Subject Headings (MeSH) Adverse Outcome Pathway (AOP) Network Pipeline

## Overview

This repository contains a comprehensive computational pipeline designed to construct, filter, and analyze knowledge graphs representing Adverse Outcome Pathways (**AOPs**) and biological flows. By leveraging the NCBI Entrez API and the complete offline NLM PubMed Baseline, the pipeline extracts primary literature associated with specific Medical Subject Headings (**MeSH**), maps multi-generational citation topologies, and calculates semantic co-occurrence networks.

The system connects **Stressors** (e.g., chemicals) to **Adverse Outcomes** (e.g., diseases) through biological intermediates. It utilizes Global-Local Filtering (GLF) and Simulated Annealing (SA) to optimize subgraph density, Louvain heuristics for community detection, and Contextual Relevance Scoring (CRS) to rank nodes and edges based on their impact within the global corpus of literature.

---

## Repository Structure

The package assumes and enforces the following directory architecture. The required MeSH XML input file must be placed in the raw data directory prior to execution.

```text
Mesh-Network-Analysis-Main-Library/
│
├── data/                               # Data storage (auto-generated)
│   ├── raw/                            # User-defined target inputs
│   │   ├── desc2025.xml                # Required: MeSH XML source
│   │   ├── aop_annotations_master.csv  # Auto-generated Master Dictionary
│   │   ├── pubmed_baseline/            # Auto-downloaded NLM Baseline XMLs (~40GB)
│   │   └── pubmed_updates/             # Auto-downloaded NLM Daily Update XMLs
│   ├── processed/                      # Target pipeline databases and JSONs
│   ├── reference_raw/                  # Curated reference inputs
│   └── reference_processed/            # Curated reference outputs
│
├── results/                            # Output artifacts (auto-generated)
│   ├── figures/                        # High-resolution plots (.png, .tif, .html)
│   ├── logs/                           # System logs and failed fetch records
│   ├── *_run_annotations.csv           # Run-specific AOP annotation templates
│   ├── *_relevance_*.csv               # Secondary analysis exports
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
│   │   ├── mesh_data_processor.py      # Unified XML extraction and stop-word generation
│   │   ├── mesh_stop_words.py          # Auto-generated MeSH stop-word set
│   │   ├── network.py                  # NetworkX assembly, filtering, and centrality
│   │   ├── relevance.py                # Contextual Relevance Scoring (Semantic Re-ranking)
│   │   ├── secondary_analysis.py       # Metadata hydration and targeted graph querying
│   │   ├── stats.py                    # GLF/SA mathematical models and graph statistics
│   │   ├── viz.py                      # Matplotlib, Seaborn, and Plotly graphics
│   │   └── wizard.py                   # Interactive configuration module
│   └── mesh_aop_notebooks/             # Jupyter notebook equivalents of each module
│       └── *.ipynb                     # One notebook per module for interactive exploration
│
├── environment.yml                     # Mamba/Conda cross-platform dependency resolution
├── mesh_config.json                    # Auto-generated user configuration file
├── pyproject.toml                      # Modern Python package specification
└── README.md                           # This document


```

---

## Data Acquisition & Prerequisites

### 1. Acquiring the MeSH XML File

The NLM has officially discontinued the MeSH ASCII format as of 2026. This pipeline now utilizes the computational gold-standard **MeSH XML format**.

* **Download:** Navigate to the [NLM MeSH Data Distribution page](https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/) and download `desc2025.xml` (or the most current yearly release).
* **Placement:** Place this file directly into your `data/raw/` directory.

### 2. Internet Connectivity Requirement

**Active internet connectivity is strictly required for the first execution.** The pipeline must connect to the NCBI FTP servers to download the full PubMed Baseline (thousands of `.xml.gz` files) and compile the ~30-million record local master database. Subsequent analytical runs can be performed completely offline. These baseline and daily-update archives are the official NLM/NCBI annual PubMed releases, distributed at the [NCBI PubMed Data Distribution page](https://pubmed.ncbi.nlm.nih.gov/download/).

---

## Environment Setup & Installation

This pipeline is computationally intensive and relies heavily on numerical arrays and graph operations. Due to the complexity of the underlying C-extensions in libraries like SciPy and NumPy, we strongly recommend using Mamba/Micromamba to resolve dependencies and ensure cross-platform compatibility.

### 1. System Requirements

* **Python:** **3.11 or 3.12 only** (enforced by `requires-python = ">=3.11,<3.13"`). Python 3.13+ is *not* supported: `node2vec` requires `numpy<2.0`, and `numpy<2.0` has no prebuilt wheel for 3.13, so pip tries to compile NumPy from source and fails without a C/C++ compiler. See **Troubleshooting** if you only have 3.13.
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

### 2. Alternative Installation: via Mamba/Micromamba

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

* **Entrez Email & API Key:** Registration with NCBI allows for 10 API requests per second. Without a key, requests are hard-limited to 3 per second, increasing Step 2 processing time exponentially. Sign in or create a free account at the [NCBI/NLM account portal](https://account.ncbi.nlm.nih.gov/) to obtain a key (found under *Account Settings → API Key Management*).

### 4. Search Parameters

* **MeSH Search Term:** The primary indexing term used to retrieve the base (P0) cohort of articles from PubMed (e.g., `Dermatitis, Allergic Contact [Mesh]`).
* **Start / End Date:** Constrains the temporal boundaries of the initial P0 PubMed search.
* **Citation Generations:** Controls the depth of the citation scrape.
  * `0`: Only the base parental generation (P0) articles.
  * `1`: P0 articles + all articles they cite + all articles that cite them (G1).
  * *Warning:* Values $\ge 1$ result in exponential data growth.
* **Update MeSH Support Files:** If `True`, forces re-extraction of the MeSH term list and stop-word set from the XML file even if cached outputs already exist.



### 5. Analysis Parameters

* **Calculate Full Centrality (Boolean):** If `True` (Default), calculates Eigenvector and Betweenness centralities (using the K-Samples parameter to estimate Betweenness for speed).
  * **[!] WARNING:** If set to `False`, the pipeline skips this heavy graph math to prevent RAM exhaustion on massive networks. You will receive a "Bare Bones" network based purely on co-occurrence. Advanced downstream metrics, including Article Relevance Scores (ARS) and Contextual Relevance Scores (CRS), **cannot and will not be calculated**.
* **Betweenness K-Samples:** Heuristic sampling limit for Centrality calculation. Lower values increase speed but reduce precision.
* **Context Start / End Date:** Temporal constraints applied exclusively to Step 3 Contextual Relevance Scoring, allowing the simulation of historical network states.
* **Random Seed:** Integer seed passed to NetworkX and scikit-learn for reproducible t-SNE projections and community detection. Default: `42`.
* **Eigenvector Max Iterations / Tolerance:** Power-iteration convergence controls for Eigenvector centrality. Increase `eigenvector_max_iter` (default `1000`) if convergence warnings appear on very large or sparse graphs.

### 6. Network & Simulation Parameters

* **Lambda Value:** The distance penalty decay factor applied to generational node weighting ($W = e^{-\lambda d}$).
* **Node Weight Factors:** A four-component weighting dictionary that controls how a node's composite importance score is assembled prior to GLF/SA filtering. The four keys and their defaults are:

  | Key | Default | Role |
  |-----|---------|------|
  | `centrality` | `0.45` | Combined betweenness + eigenvector centrality contribution |
  | `article_rank` | `0.15` | PageRank-style article influence |
  | `rank_median_cit` | `0.20` | Median incoming citations across linked articles |
  | `rank_total_cit` | `0.20` | Total incoming citation volume |

  The four values must sum to `1.0`.

* **Target Edges:** The hard threshold for the final consensus subgraph size. The optimization algorithms will prune the graph until exactly this number of edges remains.
* **GLF / SA Iterations:** Monte Carlo search and thermal cooling steps for the optimization heuristics.
* **SA Temperature Start:** The initial thermal energy of the Simulated Annealing system (default `5000.0`). Higher values allow larger disruptive jumps early in the search.
* **SA Cooling Rate:** The multiplicative decay applied to temperature each iteration (default `0.999995`). Values closer to `1.0` cool more slowly and explore more widely at the cost of run time.

### 7. Secondary Analysis Parameters

Executes highly targeted queries against the finalized network to extract specific publications for manual review.

* **Export Top Articles:** If `True` (default), always exports the highest-scoring network-wide articles at the end of Step 3 without requiring `--step secondary`.
* **Export Limit:** Maximum number of articles returned per query (default `500`).
* **Exclude Review Articles:** Filters out broad review articles to isolate primary literature.
* **Target Nodes:** Evaluates the literature density of specific nodes. **Must be semicolon-separated** (e.g., `Skin; T-Lymphocytes`).
* **Target Edges:** Evaluates literature specifically linking two concepts. Format with a dash-space-dash: `NodeA - NodeB; NodeC - NodeD`.
* **Sort Metric:** The alternative options for the `sort_metric` parameter is **`F1`** or **`Linear`**.

  **1. Linear (Weighted Additive Model)**

  * **Mechanism:** Calculates the arithmetic weighted average of the normalized Article Relevance Score (ARS) and the normalized Citation Score.
  * **Behavior:** Compensatory. A high score in one metric can offset a low score in the other, governed by the user-defined `linear_weight_ars`.
  * **Formula:** $(ARS \times w) + (Cit \times (1 - w))$

  **2. F1 (Harmonic Mean)**

  * **Mechanism:** Calculates the strict harmonic mean of the normalized ARS and the normalized Citation Score.
  * **Behavior:** Penalizing. The final output skews heavily toward the lower of the two input values. An article must possess both high topological relevance and high community impact to achieve a high score.
  * **Formula:** $2 \times \frac{ARS \times Cit}{ARS + Cit}$
* **ARS Weight (`linear_weight_ars`):** The weight for **`Linear`** sorting metric of the ARS scores per article (0–1.0), with the remaining weight applied to incoming citations (*article popularity*). Default `0.5`.

### 8. Benchmark Parameters

Controls the optional `--step benchmark` evaluation (see the **Validation & Benchmarking** section below). Configured under the `benchmark` block of `mesh_config.json`.

* **`ground_truth_csv`:** Filename of your ground-truth PMID set inside the active raw directory. Leave **empty** (default) to auto-detect a recognized filename; set it to pin a specific file.
* **`negative_control_csv`:** Optional filename of an *unrelated* ground-truth set used as a specificity check (it should score near random). Empty disables the control.
* **`primary_node`:** The base disease/seed node used to separate "naive" articles (those carrying the primary node) from topology-exclusive hits. Defaults to the search term's MeSH node.
* **`n_boot` / `n_perm`:** Bootstrap resamples and permutation iterations for confidence intervals and the random-ranking null. Default `100` each, chosen to keep the step tractable: **every bootstrap resample re-sorts the entire article pool**, so cost scales with pool size (~4 s per iteration on a ~9 M-row pool, run twice per scorer across three scorers). Approximate total runtime for `--step benchmark`:

  | `n_boot` / `n_perm` | Runtime | Use |
  | --- | --- | --- |
  | `50` | ~20 min | fast smoke test |
  | `100` *(default)* | ~40 min | routine runs |
  | `200` | ~1.3 h | solid intervals |
  | `1000`–`2000` | ~6–13 h | publication-grade |

  Raise for tighter intervals, but note the confidence-interval width is ultimately limited by **how many ground-truth positives you have**, not by `n_boot`.

---

## The AOP Annotation Workflow (Biological Strata)

**Critical Concept:** The pipeline algorithms can identify the statistical relationships between MeSH terms, but they cannot automatically determine if a term represents a "Molecular", "Cellular", or "Tissue" level event. To generate accurate biological Sankey flows, human context must be provided.

To streamline this, the pipeline utilizes a **Semicolon-Delimited Master Dictionary System**.
*Note: Because MeSH terms frequently contain commas (e.g., `Dermatitis, Allergic Contact`), standard CSV comma-delimiters will corrupt the data. All annotation files in this pipeline exclusively use semicolons (`;`).*

### How to Annotate Your Network:

1. **Enable Pausing:** In the Wizard, ensure `Pause for Annotation` is set to `True`.
2. **Run the Pipeline:** Let the pipeline run. It will execute Steps 0 through 3 (Network Construction), output AOP-independent figures (like distributions and convergence trajectories), and then pause.
3. **Open the Run Template:** Navigate to the `results/` directory and open your run-specific template: `[PREFIX]_run_annotations.csv`.
4. **Assign Strata:** This file contains every surviving node in your network. It automatically pulls any known assignments from your Master Dictionary. For any term listed as `Unassigned`, replace the text with one of the following 7 strata. *(Unsure how to categorize a term? Look up its official scope note and hierarchy in the [NLM MeSH Browser](https://meshb.nlm.nih.gov/) — the official MeSH description is the best guide for determining the correct AOP stratum.)*
* `Stressor` - External stimuli that initiate a biological reaction (e.g., `UV Rays`, `Chemicals`)
* `Molecular` - Gene, protein, or receptor level events (e.g., `Receptors, Antigen, T-Cell`)
* `Cellular` - Cellular level events (e.g., `Chemotaxis`, `Apoptosis`)
* `Tissue` - Events localized to a subsection of tissue (e.g., `Necrosis`)
* `Organ` - Organ level names or events (e.g., `Liver`, `Skin`)
* `Adverse Outcome` - High-level disease outcomes (e.g., `Drug Hypersensitivity`)
* `Uncategorized` - Broad biological terms that do not fall into a distinct AOP stratum.


5. **Save the File:** Save the file, ensuring it remains **semicolon-delimited**.
6. **Resume the Pipeline:** Run `mesh-pipeline --step viz`. The pipeline will ask if you want to sync your new annotations to the Master Dictionary for future runs, and then generate the final biological figures.

---

## Output Artifacts

Upon successful completion of the pipeline, the following critical files are generated:

### Data & Network Artifacts (`data/processed/`)

* `master_mesh_database.db`: A persistent offline cache of all PubMed IDs and MeSH annotations.
* `*_cleaned_pmids.db`: The SQLite schema containing the hierarchical linkage of your extracted articles.
* `*_full_network_data.json`: The raw, unfiltered NetworkX graphical representation.
* `*_consensus_lcc_network.json`: The optimized graph containing the intersection of the GLF and SA algorithms, reduced to its Largest Connected Component (LCC).
* `*_final_network_with_relevance.json`: The fully annotated terminal graph, populated with semantic Contextual Relevance Scores (CRS) and Louvain Community classifications.

### Analytical Exports (`results/`)

* `*_export.xlsx`: A tabular summary of the final nodes and edges.
* `*_Top_Network_Articles.csv`: The highest-scoring primary literature driving the network's structure.
* **Figures (`results/figures/`)**:
* **Figure 1:** Edge weight distribution (Power law analysis) to assess network topology.
* **Figure 2:** Optimization Trajectory (GLF vs SA convergence).
* **Figure 3:** Community Composition bar charts detailing biological strata makeup per cluster.
* **Figure 4:** CRS Centrality correlations (Betweenness vs PageRank).
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

Written to `results/`:

* `*_benchmark_results.json` — all metrics, confidence intervals, lift, and coverage figures, with the run's seed and configuration for reproducibility.
* `*_benchmark_quarantined_pmids.csv` — ground-truth rows flagged as implausible, for manual adjudication.
* `*_benchmark_enrichment.png` — cumulative recall (enrichment) curve: ground-truth recovery vs. fraction of the pool screened.

> **Note on incomplete judgments:** the non-ground-truth articles are *unjudged*, not confirmed negatives. This is why recall- and rank-based metrics (which only need positive positions) are the headline numbers, while precision/PR-AUC are interpreted with caution. For the most rigorous citation→PMID verification, cross-check publication years online via NCBI Entrez (not enabled by default, as it requires network access).

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
    run_contextual_relevance_scoring,
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
run_contextual_relevance_scoring(...)
```

The full list of exported symbols is defined in `src/mesh_aop/__init__.py`.

---

## Troubleshooting

### Python 3.13: NumPy tries to compile and fails

Symptom: during `pip install -e .` you see meson/`Unknown compiler` errors building NumPy, e.g. `Could not find ... vswhere.exe`, or `ResolutionImpossible` mentioning `numpy`.

Cause: `node2vec` requires `numpy<2.0`, and `numpy<2.0` ships **no prebuilt wheel for Python 3.13**, so pip falls back to compiling NumPy from C source — which needs an MSVC compiler you likely don't have.

Fix: use **Python 3.11 or 3.12**. On Windows, install 3.12 and recreate the venv with `py -3.12 -m venv ...`. (Any 3.12.x patch release works.)

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

1. Set `calculate_full_centrality: false` in `mesh_config.json` to skip Eigenvector and Betweenness calculations.
2. Reduce `betweenness_k_samples` (default `1000`) to lower the heuristic sampling budget.
3. Reduce `target_num_edges` so GLF/SA operate on a smaller subgraph.

Note that disabling full centrality prevents CRS scores from being calculated (see the Analysis Parameters section above).

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