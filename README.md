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
│   └── mesh_aop/                       # Core Python package modules
│       ├── __init__.py
│       ├── baseline_manager.py         # Multi-core MapReduce ETL for the Master Database
│       ├── cli.py                      # Orchestrator and CLI entry point
│       ├── config_parser.py            # Two-tier configuration engine
│       ├── data_ops.py                 # SQLite and NCBI Entrez querying
│       ├── mesh_data_processor.py      # Unified XML extraction and stop-word generation
│       ├── network.py                  # NetworkX assembly, filtering, and centrality
│       ├── relevance.py                # Contextual Relevance Scoring (Semantic Re-ranking)
│       ├── secondary_analysis.py       # Metadata hydration and targeted graph querying
│       ├── viz.py                      # Matplotlib, Seaborn, and Plotly graphics
│       └── wizard.py                   # Interactive configuration module
│
├── check_env.py                        # System environment & dependency verification script
├── environment.yml                     # Mamba/Conda cross-platform dependency resolution
├── mesh_config.json                    # Auto-generated user configuration file
├── pyproject.toml                      # Modern Python package specification
└── README.md                           # This document


```

---

## Data Acquisition & Prerequisites

### 1. Acquiring the MeSH XML File

The NLM has officially discontinued the MeSH ASCII format as of 2026. This pipeline now utilizes the computational gold-standard **MeSH XML format**.

* **Download:** Navigate to the NLM MeSH Data Distribution Page and download `desc2025.xml` (or the most current yearly release).
* **Placement:** Place this file directly into your `data/raw/` directory.

### 2. Internet Connectivity Requirement

**Active internet connectivity is strictly required for the first execution.** The pipeline must connect to the NCBI FTP servers to download the full PubMed Baseline (thousands of `.xml.gz` files) and compile the ~30-million record local master database. Subsequent analytical runs can be performed completely offline.

---

## Environment Setup & Installation

This pipeline is computationally intensive and relies heavily on numerical arrays and graph operations. Due to the complexity of the underlying C-extensions in libraries like SciPy and NumPy, we strongly recommend using Mamba/Micromamba to resolve dependencies and ensure cross-platform compatibility.

### 1. System Requirements

* **Python:** Version 3.11 or greater is required. (Python 3.13+ may cause compilation conflicts).
* **Memory:** Minimum 16GB RAM; 32GB+ highly recommended for Step 0 (Database Compilation) and networks exceeding 1 citation generation.
* **Storage:** 100GB+ free space (The NLM Baseline XMLs and resulting SQLite Master Database expand rapidly).

### 2. Installation (Standard Pip/Venv)

If Mamba is unavailable, you may use a standard Python virtual environment. This approach relies on the `pyproject.toml` file to resolve dependencies via pip.

```bash
# Navigate to the repository root
cd path/Mesh-Network-Analysis-Main-Library

# Create and activate the virtual environment
python3 -m venv mesh_env
source mesh_env/bin/activate  # On Windows use: mesh_env\Scripts\activate

# Run the Environment Checker (Automatically initiates 'pip install -e .')
python check_env.py --auto

```

*Note: The environment checker specifically resolves namespace collisions between the generic `community` package and the required `python-louvain` package during pip installations.*

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
python check_env.py --auto

```

### 4. Verification

Regardless of the installation method chosen, verify the package is correctly linked to your PATH by calling the command line interface:

```bash
mesh-pipeline --help

```

---

## Execution Guide

The pipeline is entirely modular and controlled via a terminal interface. Configuration is handled by an interactive command-line wizard, allowing users to modify runtime parameters safely without touching source code.

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

* **Entrez Email & API Key:** Registration with NCBI allows for 10 API requests per second. Without a key, requests are hard-limited to 3 per second, increasing Step 2 processing time exponentially.

### 4. Search Parameters

* **MeSH Search Term:** The primary indexing term used to retrieve the base (P0) cohort of articles from PubMed (e.g., `Dermatitis, Allergic Contact [Mesh]`).
* **Start / End Date:** Constrains the temporal boundaries of the initial P0 PubMed search.
* **Citation Generations:** Controls the depth of the citation scrape.
* `0`: Only the base parental generation (P0) articles.
* `1`: P0 articles + all articles they cite + all articles that cite them (G1).
* *Warning:* Values $\ge 1$ result in exponential data growth.



### 5. Analysis Parameters

* **Calculate Full Centrality (Boolean):** If `True` (Default), calculates Eigenvector and Betweenness centralities (using the K-Samples parameter to estimate Betweenness for speed).
  * **[!] WARNING:** If set to `False`, the pipeline skips this heavy graph math to prevent RAM exhaustion on massive networks. You will receive a "Bare Bones" network based purely on co-occurrence. Advanced downstream metrics, including Article Relevance Scores (ARS) and Contextual Relevance Scores (CRS), **cannot and will not be calculated**.
* **Betweenness K-Samples:** Heuristic sampling limit for Centrality calculation. Lower values increase speed but reduce precision.
* **Context Start / End Date:** Temporal constraints applied exclusively to Step 3 Contextual Relevance Scoring, allowing the simulation of historical network states.

### 6. Network & Simulation Parameters

* **Lambda Value:** The distance penalty decay factor applied to generational node weighting ($W = e^{-\lambda d}$).
* **Target Edges:** The hard threshold for the final consensus subgraph size. The optimization algorithms will prune the graph until exactly this number of edges remains.
* **GLF / SA Iterations:** Monte Carlo search and thermal cooling steps for the optimization heuristics.

### 7. Secondary Analysis Parameters

Executes highly targeted queries against the finalized network to extract specific publications for manual review.

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
* **ARS Weight:** The weight for **`Linear`** sorting metric of the ARS scores per article (0-1.0) with the other percentage made up by incoming citations (*article popularity*).

---

## The AOP Annotation Workflow (Biological Strata)

**Critical Concept:** The pipeline algorithms can identify the statistical relationships between MeSH terms, but they cannot automatically determine if a term represents a "Molecular", "Cellular", or "Tissue" level event. To generate accurate biological Sankey flows, human context must be provided.

To streamline this, the pipeline utilizes a **Semicolon-Delimited Master Dictionary System**.
*Note: Because MeSH terms frequently contain commas (e.g., `Dermatitis, Allergic Contact`), standard CSV comma-delimiters will corrupt the data. All annotation files in this pipeline exclusively use semicolons (`;`).*

### How to Annotate Your Network:

1. **Enable Pausing:** In the Wizard, ensure `Pause for Annotation` is set to `True`.
2. **Run the Pipeline:** Let the pipeline run. It will execute Steps 0 through 3 (Network Construction), output AOP-independent figures (like distributions and convergence trajectories), and then pause.
3. **Open the Run Template:** Navigate to the `results/` directory and open your run-specific template: `[PREFIX]_run_annotations.csv`.
4. **Assign Strata:** This file contains every surviving node in your network. It automatically pulls any known assignments from your Master Dictionary. For any term listed as `Unassigned`, replace the text with one of the following 7 strata:
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
* **Figure 4:** CRS Centrality correlations (Betweenness vs Eigenvector).
* **Figure 5:** t-SNE projection of the network colored by Louvain community.
* **Figure 6:** AOP Alluvial/Sankey flows (The primary visualization connecting Stressors to Outcomes). *(Note: Provided as interactive `.html` files for deep pathway inspection).*
* **Figure 7:** Dumbbell plots assessing shift in topological vs semantic relevance.



---

## Citation

If you use this code or methodology, please cite:

> *[To be updated after publication]*
> [](https://doi.org/10.5281/zenodo.18662959)

## License

[MIT License]
