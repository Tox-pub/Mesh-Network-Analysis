# MeSH Workbench — Help and reference

How the pipeline works, what every setting does, and what it produces. For
getting the program onto a machine see [INSTALL.md](INSTALL.md); for a short
overview see the [README](README.md).


## Overview

This repository contains a comprehensive computational pipeline designed to construct, filter, and analyze knowledge graphs representing Adverse Outcome Pathways (**AOPs**) and biological flows. By leveraging the NCBI Entrez API and the complete offline NLM PubMed Baseline, the pipeline extracts primary literature associated with specific Medical Subject Headings (**MeSH**), maps multi-generational citation topologies, and calculates semantic co-occurrence networks.

The system connects **Stressors** (e.g., chemicals) to **Adverse Outcomes** (e.g., diseases) through biological intermediates. It utilizes Global-Local Filtering (GLF) and Simulated Annealing (SA) to optimize subgraph density, Louvain heuristics for community detection, and Mean Relevancy Scoring (MRS) to rank nodes and edges based on their impact within the global corpus of literature.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Data Acquisition & Prerequisites](#data-acquisition--prerequisites)
  - [The MeSH XML File (automatic)](#1-the-mesh-xml-file-automatic)
  - [Internet Connectivity & Disk Budget](#2-internet-connectivity--disk-budget)
- [User-Provided Files — Quick Reference](#user-provided-files--quick-reference)
- [Setting up a development environment](#setting-up-a-development-environment)
  - [Requirements](#requirements)
  - [Install from source](#install-from-source)
- [Execution Guide](#execution-guide)
  - [CLI Flags](#cli-flags)
  - [Running the Complete Pipeline](#running-the-complete-pipeline)
  - [Running Individual Modules](#running-individual-modules)
- [Configuration Wizard Parameter Glossary](#configuration-wizard-parameter-glossary)
  - [1. Control Flags & Directories](#1-control-flags--directories)
  - [2. Master Database Status](#2-master-database-status-step-0-etl)
  - [3. NCBI Credentials](#3-ncbi-credentials)
  - [4. Search Parameters](#4-search-parameters)
  - [5. Analysis Parameters](#5-analysis-parameters)
  - [6. Network & Simulation Parameters](#6-network--simulation-parameters)
  - [7. Secondary Analysis Parameters](#7-secondary-analysis-parameters)
  - [8. Benchmark Parameters](#8-benchmark-parameters)
- [The AOP Annotation Workflow](#the-aop-annotation-workflow-biological-strata)
  - [How to Annotate Your Network](#how-to-annotate-your-network)
  - [Syncing to the Master Dictionary](#what-syncing-to-the-master-dictionary-does)
- [Output Artifacts](#output-artifacts)
- [The Workbench Window](#the-workbench-window)
- [Controlling a Run](#controlling-a-run)
- [When Files Go Wrong](#when-files-go-wrong)
- [The Run Ledger & PRISMA Report](#the-run-ledger--prisma-report)
- [Validation & Benchmarking](#validation--benchmarking)
  - [The Bundled OECD Ground Truth](#the-bundled-oecd-ground-truth)
  - [Supplying Your Own Ground Truth](#supplying-your-own-ground-truth)
  - [What It Reports](#what-it-reports)
  - [Node/Edge Convergent Validation](#nodeedge-convergent-validation)
- [Jupyter Notebook Interface](#jupyter-notebook-interface)
- [Programmatic API Usage](#programmatic-api-usage)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

---

## Repository Structure

The package assumes and enforces the following directory architecture. The required MeSH XML input file must be placed in the raw data directory prior to execution.

```text
Mesh-Network-Analysis/
│
├── data/                               # Data storage
│   ├── raw/                            # Inputs for a run
│   │   ├── aop_annotations_master.csv  # Ships w/ repo: AOP strata dictionary (pre-seeded; grows each run)
│   │   ├── desc2025.xml                # Auto-downloaded from NLM if missing (or place manually); not in repo
│   │   ├── ground_truth_pmids.template.csv # Ships w/ repo: copy+fill for your own benchmark set
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

## User-Provided Files — Quick Reference

Everything a user supplies, where it goes, and how the pipeline picks it up. Most inputs are automatic or optional: only the MeSH XML (auto-downloaded) is always needed, plus a ground-truth list **if** you benchmark your own data.

| File | Where it goes | How it's picked up | Format / structure |
|---|---|---|---|
| **MeSH descriptor XML** | `data/raw/descYYYY.xml` | Auto-downloaded from NLM if missing (§1); or place manually | Unmodified NLM MeSH descriptor XML |
| **Ground-truth PMIDs** (benchmark) | `data/raw/` (own data), or a path | Auto-detected by name, or `benchmark.ground_truth_csv`; requires `run_ground_truth_analysis = true` | `PMID` column required — see the template and **Supplying Your Own Ground Truth** |
| **Negative-control PMIDs** (optional) | `data/raw/`, or a path | `benchmark.negative_control_csv` (filename or path) | Same structure as ground truth |
| **Comparison networks** (optional) | `data/processed/`, or a path | `comparison_networks` list when `compare_networks` is on | Pipeline network JSON (Cytoscape) or a networkx-readable graph — **reuse `*_consensus_lcc_network.json` outputs; do not hand-author** |
| **AOP strata annotations** | `results/*_run_annotations.csv` (generated) | You edit the generated template during the Step-3 pause | Semicolon-delimited; assign one of the 7 strata |
| **Entrez credentials** | wizard, or environment variables | `MESH_ENTREZ_EMAIL` / `MESH_ENTREZ_API_KEY`, or the wizard | Email + NCBI API key |

**Templates** for the files you create yourself ship next to where they belong (e.g. [`data/raw/ground_truth_pmids.template.csv`](data/raw/ground_truth_pmids.template.csv)) — copy, rename, and fill. Everything else in `data/processed/`, the network JSONs, the relevance databases, and all figures is **generated** — never hand-authored.

---

## Setting up a development environment

Installing the **application** is covered in [INSTALL.md](INSTALL.md) — that is
the right document for anyone who simply wants to run it, and none of the below
is needed. What follows is for working on the source.

### Requirements

* **Python 3.11–3.13** (`requires-python = ">=3.11,<3.14"`). 3.13 is supported:
  the Node2Vec embedding is vendored in `node2vec_embedding.py`, so the
  `node2vec` package — which pinned `numpy<2.0` and had no 3.13 wheel — is no
  longer a dependency. Newer releases are untested.
* **Memory** — 16 GB minimum, 32 GB or more for the database build and for
  networks past one citation generation.
* **Storage** — 100 GB+ free. The NLM archives and the SQLite master database
  expand quickly.
* **tkinter**, for the desktop application. It ships with the python.org
  installer, but is a separate package on most Linux distributions
  (`python3-tk`) and is absent from the Microsoft Store build of Python.

### Install from source

```bash
git clone https://github.com/Tox-pub/Mesh-Network-Analysis.git
cd Mesh-Network-Analysis
python -m venv ~/mesh_env
~/mesh_env/bin/python -m pip install -e .
```

On **Windows**, create the environment at a short path such as
`C:\Users\<you>\mesh_env` — not inside a deeply nested or cloud-synced folder.
Several dependencies ship very long filenames that overflow the 260-character
`MAX_PATH` limit and abort the install part-way through. See **Troubleshooting**.

Verify the entry points resolve:

```bash
mesh-pipeline --version
mesh-check-env
```

`mesh-check-env` also reports missing OS-level rendering libraries.

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

* **Use Bundled Reference Data (demonstration only):** Analyses the finished network that ships with the program instead of building one from your own search. See [Using the bundled reference data](#using-the-bundled-reference-data) below for what it does and does not change.
* **Pause for Annotation (AFK Mode):** If `False` (Default), the pipeline operates in AFK Mode. It will run uninterrupted from start to finish, automatically assigning 'Unassigned' to all biological levels. If `True`, the pipeline will safely pause after Step 3 to allow the user to manually annotate the network before rendering the final biological visualizations.
* **Custom Prefix:** The naming convention prepended to all output files (e.g., `DAC_Mesh`).

#### Using the bundled reference data

**This is for demonstration, not for research.** A first real run means a long PubMed download and a full rebuild before you see a single figure. Ticking this box instead analyses a network that already exists — the allergic contact dermatitis network published with this software, already built and scored — so the figures, the workflow report and the benchmark come out in minutes. Use it to judge whether the outputs are what you want, and to learn where everything lands on disk, before committing to a retrieval of your own.

What you **cannot** do with it is change the corpus. The articles behind that network are fixed, the retrieval that produced them is not repeated and cannot be varied, and any finding in it is already published — it is not yours. Untick the box and run your own search for that.

| | |
| :--- | :--- |
| **Set for you, and greyed out** | Search term, both date windows, citation depth, random seed, benchmark primary node — anything else would label the figures with a query that did not make them. |
| **Still yours** | Figure resolution and formats, which figures are drawn, your folders, your credentials, your project prefix. |
| **Does not run** | Retrieval (there is nothing to fetch), network construction (it ships built), and secondary analysis (it needs the cleaned citation database, which only retrieval produces). |
| **Your settings** | Untouched. Untick the box and your own search term and dates come back as you left them. |
| **File naming** | Outputs are prefixed `Reference_` so they can never be confused with your own. |

The reference networks are **copied into your own networks folder** the first time you tick the box, and everything downstream reads them from there. They are yours to open, edit or delete like any other result. Copies already in that folder are left alone, so your edits survive later runs; delete one and the pristine original returns from the program folder.

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

> **Target nodes and edges do not run as part of a full pipeline run.** Each one hydrates hundreds of records from NCBI, which would add a long tail to every run, so they wait until you ask for them: choose **secondary — targeted queries & exports** from the step list and press Run. Everything they need is on disk once a full run has finished, so this can be done days later as long as the project prefix is the same. The Secondary tab in the application carries the same note. (Export Top Articles is not affected — that one *does* run in a full pipeline.)

**Where to get the names to type.** A target must be a MeSH heading that is in your network, spelled exactly as the network spells it — `Dermatitis, Allergic Contact`, not `allergic contact dermatitis`. Two places give you the real list:

* The **network JSON** in your networks folder — the easiest way to see what you are choosing between. Open it in **Cytoscape** (*File → Import → Network from File*) and the whole graph is drawn for you; click any node or edge and its exact id appears in the table below, ready to copy. Gephi and the `networkx` library read these files too. The **Open network folder** button on the Results screen takes you there.
* The **network Excel export** in your results folder (written by *Export network to Excel* on the same tab) — one row per node and per edge, so you can sort by centrality and copy the names out.

**How an article is matched.** `relevance.py` records, for every article, which of your network's nodes appear among the MeSH headings it was indexed with (base headings — `Skin/drug effects` is stored as `Skin`). A target node returns the articles carrying that heading; a target edge returns the articles carrying **both**. Matching is on the **whole heading**, so `Skin` returns articles indexed under Skin and never those under `Skin Diseases`, `Skin Absorption` or `Skin Tests`. A name that is not in the network finds nothing, so what you type is checked before the query runs and the nearest matches are offered.

* **Export Top Articles:** If `True` (default), always exports the highest-scoring network-wide articles at the end of Step 3 without requiring `--step secondary`.
* **Export Limit:** Maximum number of articles returned per query (default `500`).
* **Exclude Review Articles:** Filters out broad review articles to isolate primary literature.
* **Target Nodes:** The articles behind specific nodes. **Must be semicolon-separated** (e.g., `Skin; T-Lymphocytes`).
* **Target Edges (`target_edges`):** The articles behind a specific relationship — those indexed under both headings. Separate the two node names with a **space-padded hyphen** (` - `), and separate multiple edge queries with semicolons: `NodeA - NodeB; NodeC - NodeD`. *(Distinct from the numeric `target_num_edges`)*
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

* **`ground_truth_csv`:** Ground-truth PMID set. Leave **empty** (default) to auto-detect a recognized filename in `data/raw/`; set it to pin a specific filename or path. When `Use Reference Data` is on, the bundled OECD set is used automatically. See **Supplying Your Own Ground Truth** for the required structure and the template.
* **`negative_control_csv`:** Optional *unrelated* ground-truth set used as a specificity check (it should score near random). Same location, resolution, and structure as `ground_truth_csv` (drop it in `data/raw/` or give a path). Empty (default) disables the control.
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
* `*_pmids.db`: What retrieval downloaded — every PMID reached, which citation generation it came from, and its `cited_by` list.
* `*_cleaned_pmids.db`: **Belongs to one search, not to the installation.** It is a copy of `*_pmids.db` with the citation generations at or beyond your generation depth deleted, then annotated with MeSH headings from the master database. Because it depends on the search term, the date window and the generation depth, it is rebuilt per project prefix rather than shared — it is *not* a derivation of the master MeSH database and does not belong in the Databases screen.

  It is created automatically in **Step 2 (Retrieval & Database Operations)** of a full run, and it is a hard prerequisite of network construction, secondary analysis and the benchmark. If you delete it, the next full run rebuilds it from `*_pmids.db` without re-downloading anything.

  The **bundled reference corpus has no cleaned database and cannot have one**: it ships as a finished network, and the article lists behind it were never downloaded. That is why secondary analysis is refused under *Use bundled reference data* rather than producing empty output — and it is not a path problem that copying files could fix.
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
* `*_run_ledger.csv`: Every count the run produced, at every stage — records per citation generation, MeSH annotations screened, what each optimiser kept, what the consensus and the LCC discarded. Semicolon-delimited (MeSH headings contain commas), one ledger per file prefix, with the timestamp each figure was recorded. See **The Run Ledger & PRISMA Report** below.
* `*_prisma_flow_report.txt` and `figures/*_prisma_flow.*`: A PRISMA-style flow of records through the pipeline, in text and as a figure — the overview of the whole search, suitable for a methods section or a supplementary figure.
* `*_network_overlap_membership.csv` / `*_network_overlap_matrix.csv` / `*_Network_Overlap.png`: node-overlap comparison across the networks you named — produced only when **Compare Multiple Networks** is enabled.
* **Benchmark & validation outputs** (all under `results/benchmark/`): the `--step benchmark` step asks several different questions of the same network, so its folder is grouped by question rather than left as one pile of files:

  | Folder | What is in it |
  | --- | --- |
  | `inputs/` | The ground truth the run actually used, kept under the project prefix (`*_ground_truth.*`, and `*_negative_control.*` when one is configured). A benchmark number cannot be read without knowing what it was scored against, and the original lives wherever you put it — this copy is the one the results belong to. |
  | `ranking/` | The article-ranking benchmark: `*_benchmark_results.json`, `*_benchmark_quarantined_pmids.csv`, and `figures/*_benchmark_enrichment.png`. |
  | `ranking_validation/` | Every node weighting scored across the evaluation frames: `*_validation_report.xlsx` / `.html`, `*_projection_comparison.csv`, and `figures/`. |
  | `network_validation/` | The node/edge convergent validation — is the network's vocabulary and wiring reproduced? `*_gt_network_validation.xlsx`, `*_gt_cooccurrence_network.json`, and `figures/*_GT_*.png`. |

  See the **Validation & Benchmarking** section below for what each file contains.
* `logs/`: run logs and failed-fetch records.
* **Figures (`results/figures/`)**:
* **Figure 1:** Edge weight distribution (Power law analysis) to assess network topology.
* **Figure 2:** Optimization Trajectory (GLF vs SA convergence).
* **Figure 3:** Community Composition bar charts detailing biological strata makeup per cluster.
* **Figure 4:** MRS Centrality correlations (Betweenness vs PageRank).
* **Figure 5:** t-SNE projection of the network colored by Louvain community.
* **Figure 6:** AOP Alluvial/Sankey flows (The primary visualization connecting Stressors to Outcomes). *(Note: Provided as interactive `.html` files for deep pathway inspection).*



---

## The Workbench Window

### Data setup

Lists what is on disk **at the paths the pipeline actually uses** — the configured data folder, and the current project prefix. Each row shows Present/Missing, its size, and the actions available:

| Row | What it is |
| :--- | :--- |
| Master annotation database | Built from the archive. Required for every run. Hours to rebuild. |
| PubMed baseline archive | The yearly snapshot, ~44 GB. Only needed to build the database. |
| PubMed daily updates | Records published since the baseline snapshot. |
| MeSH descriptor file | Defines the stop-word vocabulary. Fetched automatically. |
| Retrieved PMIDs / Citation database / Relevance database | Per-project, named with your prefix. |

**Also fetch the daily update files when building** is a checkbox on this screen, so you can see whether it is on without committing to a build first.

Every Build, Rebuild, Download and Delete here asks you to **type `REBUILD` or `DELETE`** before it proceeds. These actions either destroy something that took hours to produce or start something that will take hours to finish, and a yes/no box is one mis-aimed click. Deleting a database also removes its `-wal`, `-shm` and health sidecars — a stale write-ahead log left beside a rebuilt database is worse than useless, because SQLite will try to replay it.

### Settings

Nine tabs. **Search**, **Folders** and **Credentials** come first, because the query, where things are written, and which NCBI account retrieves them are what a new user has to settle before anything else works.

The **Search** tab carries a standing note on the project prefix: how to resume an interrupted run by keeping it, and how to preserve an earlier analysis by changing it. The description pane at the foot of the window shows what the focused control does, its default, and the caveat that matters.

### Before a run starts

Pressing **Run** checks whether this prefix has already produced output for the steps about to execute. If it has, you are shown exactly what would be replaced, with sizes, and offered the chance to change the prefix instead. **If the prefix has not been used for those steps, nothing is asked** — a fresh project never meets a dialog it has no reason to read.

### Results

Opens on the **run overview**: the PRISMA-style account of the whole search, in a scrollable, selectable box — records retrieved per citation generation, what the MeSH screen excluded and why, what each optimiser kept, what the consensus and the largest connected component discarded, and what was finally included. Buttons open the full report and the flow diagram; the file listing sits underneath.

### Menus

- **Run** — any single step.
- **Database** — the Data setup screen.
- **Results** — the run overview.
- **Tools** — Annotate AOP levels, Show the annotation guide (both greyed until a run has produced something to annotate), Check and repair files, Uninstall.
- **Help** — this manual, the installation notes, the read-me, the reference-figure provenance, and the annotation instructions.

---

## Controlling a Run

### Pause

**Pause** asks the run to stop at its next safe point — between pipeline stages, and between batches inside the long retrieval loops. It is not instant, and the window says so: on a long step it can take a few minutes, and the status line reads *Pausing* until the run actually reports back, at which point it changes to *Paused at a safe point*. **Resume** lets it continue from exactly where it stopped.

This is deliberately cooperative rather than a process suspend, for two reasons:

- **A suspend is not reliable.** On a managed Windows machine the suspend calls return success, the per-thread suspend counts come back as expected, and the process carries on computing — the security layer virtualises them. A Pause button that silently does nothing is worse than none.
- **A suspend is not safe.** It freezes the process wherever it happens to be, which here is often the middle of a multi-gigabyte SQLite transaction. Left suspended long enough for a laptop to sleep, it produces precisely the corruption the repair tool exists to clean up.

Time spent paused is excluded from the elapsed clock, so reported durations remain honest.

### Abort

**Abort run** stops the step for good; its work is lost. The confirmation says how long the run has been going, because that is the fact that changes the decision. Completed steps are kept — a later run resumes from the last one that finished rather than starting over.

If the run reaches a checkpoint before the process is killed, it exits tidily and reports *Stopped at a safe point*, meaning nothing was left half-written.

From a terminal, the same controls are files: create `pause.flag` in the directory named by `MESH_CONTROL_DIR` to pause, delete it to resume, create `abort.flag` to stop.

---

## When Files Go Wrong

A run that is interrupted — a machine that slept, a disk that filled, a sync client caught mid-write — can leave a file that looks complete and is not. It has the right name and a plausible size, and it fails much later, inside a step that had no way to know its input was rubbish.

### Tools → Check and repair files

Opens every artefact the pipeline depends on and reports its condition:

| Status | Meaning |
| :--- | :--- |
| **ok** | Opens and reads cleanly |
| *(blank)* | Not built yet — normal, and not a problem |
| **EMPTY** | Present but zero bytes: a write that never started |
| **CORRUPT** | Opening or parsing it fails outright |
| **SUSPECT** | Opens, but the contents do not add up |
| **ORPHAN** | A temporary or part-file no finished run would leave |

Definitely-broken files are pre-ticked; merely suspect ones are listed unticked, because deleting something odd is not a call to make on your behalf. The dialog warns before removing anything that costs hours to rebuild, and afterwards names the step to resume from. Sidecar files go with their parent — a stale SQLite write-ahead log left beside a rebuilt database is worse than useless, because SQLite will try to replay it.

Nothing listed is your own work. Every one of these files is machinery the pipeline rebuilds; your results are never touched.

From a terminal:

```bash
python -m mesh_aop.cli --check-files
```

Add `--repair-files` to delete what is broken, and `--deep-check` for SQLite's full integrity check (thorough, and slow on a large database).

### The master annotation database

This one is treated differently because it is the only artefact whose replacement costs a night rather than a coffee break — roughly 44 GB downloaded and several hours compiled.

- After a successful build, a **health record** (`master_mesh_database.db.health.json`) stores its size, modification time and row count. Later damage is then detectable in seconds instead of by re-scanning eight gigabytes.
- Before any step that reads it, the pipeline runs a **pre-flight check**: open it, structural `quick_check`, confirm the table has rows, compare against the health record. A corrupt database stops the run *before* it spends hours on work that depends on it; a merely suspect one warns and continues.
- **Rebuild instructions are written beside the database itself**, as `HOW TO REBUILD THIS DATABASE - read me.txt`. A build log is not where anyone looks a month later when the file will not open.

---

## The Run Ledger & PRISMA Report

Every run writes two provenance artefacts into `results/`, without being asked and regardless of which `--step` was run.

### `*_run_ledger.csv` — the counts

The pipeline computes a great many numbers, prints them once, and throws them away; recovering one afterwards used to mean re-running the stage that produced it. The ledger keeps them. It is a semicolon-delimited table — semicolons because MeSH headings contain commas as a matter of course (*Dermatitis, Allergic Contact*) — with five columns:

| Column | Meaning |
| :--- | :--- |
| `stage` | `run`, `identification`, `screening`, `pruning`, `included`, `unipartite` — in pipeline order |
| `metric` | What was counted |
| `value` | The count, or the setting |
| `recorded` | When that row was written |
| `note` | A plain-language gloss, so the file reads without this manual |

There is **one ledger per file prefix**, and recording **merges**: re-running the network step rewrites the `pruning` and `included` rows and leaves `identification` alone. That is why every row carries its own timestamp — a ledger whose stages were written weeks apart is legible rather than misleading, and the text report calls the mixed dates out explicitly when they differ.

Counts come from whichever source is most trustworthy, in this order: the stage that just ran hands back its own numbers; failing that, the ledger reuses what it recorded last time; failing that, it reads the artefacts on disk. The middle step matters — a resumed run would otherwise re-parse a several-hundred-megabyte network JSON to recover a number it had already written down.

### `*_prisma_flow_report.txt` and `figures/*_prisma_flow.*` — the overview

The same counts, laid out as a PRISMA 2020-style flow: **Identification → Screening → Pruning → Included**, with what was excluded at each step, and why, set beside it. The pipeline has genuinely the same shape as a systematic review — a seeded PubMed search expanded over citation generations, an annotation screen, two optimisers pruning the graph, and a largest connected component that everything downstream reports on — so it is drawn the same way.

The figure is written at the dpi and in the formats set under **Figures** in the settings (default 300 dpi, JPEG + TIFF). A stage the run never reached is left out rather than drawn as zero: *not run* and *found nothing* are different claims, and a reader cannot tell them apart from a box of zeroes.

---

## Validation & Benchmarking

The `--step benchmark` module evaluates how well the network's per-article relevance scores recover an external, curated set of "ground-truth" PMIDs (e.g. the bibliography of an authoritative OECD Adverse Outcome Pathway document). It is designed for the realistic situation here — a few dozen known-relevant papers hidden in a multi-million article candidate pool — and follows three principles: **validate the ground truth first, use metrics robust to incomplete relevance judgments, and quantify uncertainty against baselines.**

```bash
mesh-pipeline --step benchmark
```

### The Bundled OECD Ground Truth

The repository ships a curated positive set derived from the reference bibliography of the **OECD Adverse Outcome Pathway for skin sensitisation (AOP 40)**. It is the default target of `benchmark.ground_truth_csv`.

**Source document.** The bibliography is taken from:

> OECD (2014), *The Adverse Outcome Pathway for Skin Sensitisation Initiated by Covalent Binding to Proteins*, OECD Series on Testing and Assessment, No. 168, OECD Publishing, Paris. <https://doi.org/10.1787/9789264221444-en>

[Direct PDF](https://www.oecd.org/content/dam/oecd/en/publications/reports/2014/09/the-adverse-outcome-pathway-for-skin-sensitisation-initiated-by-covalent-binding-to-proteins_g1g48567/9789264221444-en.pdf) · [DOI landing page](https://doi.org/10.1787/9789264221444-en)

Cite this document, not this repository, as the origin of the ground-truth set: the curated files below are a PMID resolution of its reference list, not an independent literature selection. If you replace the ground truth with your own, cite whatever source your positives came from in its place.

| File | Location | Contents |
| --- | --- | --- |
| `oecd_ground_truth_curated.xlsx` | `data/reference_processed/` | **96 curated positives** — primary research articles with resolved PMIDs |
| `oecd_ground_truth_exclusions.xlsx` | `data/reference_processed/` | **22 excluded citations**, each with an explicit reason (book chapters, OECD guidance documents, and records that are not indexed primary articles) |
| `oecd_ground_truth.bib` | `data/reference_processed/` | Typed BibTeX of the parsed citations (`@article`, `@incollection`, `@techreport`, `@misc`) |
| `oecd_resolved_citations.csv` | `data/reference_raw/` | The raw citation → PMID resolution table the curation was built from |

**How it was built.** Citations were parsed from the OECD reference list, resolved to PMIDs via Entrez `esearch` on author + year + volume + first page (rather than free-text matching, which mis-resolves), and then manually adjudicated. Non-article references and records lacking MeSH indexing were moved to the exclusions file rather than silently dropped, so the **96 + 22 = 118** entries fully account for the source bibliography — the exclusions are auditable, not hidden losses.

**Why the exclusions matter.** Anything without MeSH indexing is unreachable by a MeSH co-occurrence network *by construction*. Keeping those citations in the positive set would depress every metric for a reason unrelated to ranking quality, so they are excluded explicitly and counted separately.

### Supplying Your Own Ground Truth

Two steps to benchmark against your own citations:

1. **Enable it.** Set `benchmark.run_ground_truth_analysis = true`. It defaults to *off* for your own data and *on* only for the bundled reference corpus, so this switch is required.
2. **Provide the file.** Drop a file into **`data/raw/`** using one of these recognized names — it is auto-detected in order, **no config edit needed**:

   ```text
   ground_truth_pmids.csv
   ground_truth_pmids.txt
   ground_truth.csv
   ground_truth.txt
   ```

   Or point `benchmark.ground_truth_csv` at any filename or path. When `Use Reference Data` is on, the bundled OECD set is substituted automatically — your own file and the bundled one never shadow each other.

**Required structure.** A ready-to-fill template ships at [`data/raw/ground_truth_pmids.template.csv`](data/raw/ground_truth_pmids.template.csv) — copy it to one of the names above and replace the example rows. The only hard requirement is a **`PMID` column**:

| Column | Required? | Notes |
|---|---|---|
| `PMID` | **Yes** | One PubMed ID per row. Aliases `pmids`, `pubmed_id`, `id` are also accepted. Use integer digits only — `19033392`, **not** `19033392.0`. |
| `Reference` (or any 2nd column) | Optional | Free-text citation string. When present it enables the publication-year sanity check that flags citation→PMID mis-resolution. Quote any value containing a comma. |

**Accepted file shapes** (auto-detected, with delimiter/encoding sniffing):

* An **`.xlsx` workbook** with a `PMID` column (the format of the bundled curated set).
* A **CSV/TSV** with a `PMID` column and an optional reference/citation column.
* The bundled `Raw_Reference;PMID` **semicolon-delimited** format.
* A **single column of PMIDs**, with or without a header.
* A plain **`.txt` list**, one PMID per line (simplest — no delimiter concerns).

If no ground-truth file is found, the step prints the exact directory it searched, the accepted filenames, and these formats.

**Negative control (optional).** `benchmark.negative_control_csv` takes an *unrelated* ground-truth set (which should score near random) as a specificity check. It uses the **same location, name-or-path resolution, and structure** as the ground truth above — drop it in `data/raw/` (or give a path) and set the config key to its filename; the template applies unchanged. Empty (default) disables it.

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

Archived release: [10.5281/zenodo.18662959](https://doi.org/10.5281/zenodo.18662959)

**If you also use the bundled ground truth**, cite its source document as well — the curated PMID set is a resolution of this bibliography, not an independent selection:

> OECD (2014), *The Adverse Outcome Pathway for Skin Sensitisation Initiated by Covalent Binding to Proteins*, OECD Series on Testing and Assessment, No. 168, OECD Publishing, Paris. <https://doi.org/10.1787/9789264221444-en>

See [The Bundled OECD Ground Truth](#the-bundled-oecd-ground-truth) for how the set was derived from it.

## License

[MIT License]
