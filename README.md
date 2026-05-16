# MeSH Network Analysis to Identify Biological Link Plausibility

## Overview
This repository contains a computational pipeline designed to mimic **Adverse Outcome Pathways (AOPs)** biological flows and networks using **Medical Subject Headings (MeSH)** co-occurrence networks.

By analyzing the co-occurrence of MeSH terms across millions of PubMed articles, this tool builds a weighted network graph that connects **Stressors** (e.g., chemicals) to **Adverse Outcomes** (e.g., diseases) through biological intermediates. The pipeline utilizes graph theory, optimization algorithms, Global Likelihood Filter (GLF) & Simulated Annealing (SA), and a novel **Article Relevance Score (ARS)** and **Contextual Relevance Score (CRS)** to filter noise and identify the most biologically relevant pathways contained within the network of interest.

**Key Capabilities:**
* **Automated Data Mining:** Fetches and processes PubMed PMIDs, citations, and associated MeSH terms via the Entrez API.
* **Network Construction:** Builds rank-normalized co-occurrence networks from MeSH terms.
* **Optimization:** Uses GLF and SA to extract optimal subgraphs and form a Consensus Model (CM).
* **Contextual Scoring:** Calculates the Article Relevancy Scores (ARS) and Contextual Relevancy Scores (CRS) to rank nodes and edges to determine their impact on the final network within the context of the global corpus.
* **Visualization:** Generates publication-ready figures (Sankey diagrams, t-SNE, Alluvial flows).

---

## Project Structure

```text
Mesh-Network_Analysis/
│
├── scripts/
│   └── python/
│       ├── config.py                 # Central configuration (API keys, search terms)
│       ├── run_pipeline.py           # Orchestrator script to run the full workflow
│       ├── mesh_stop_words.py        # List of excluded terms (Geography, Publication types)
│       ├── mesh_data_processor.py    # Utilities to process raw MeSH binary files
│       ├── master_mesh_network.py    # Core script: Data fetching, network building, simulation
│       ├── secondary_analysis.py     # Exports data to Excel & runs node/edge queries
│       └── figures.py                # Generates all visualization figures
│
├── results/
│   ├── figures/                      # Final high-res TIFF/SVG/JPEG output
│   └── logs/                         # Execution logs and error reports
│
├── data/
│   ├── raw/                           # Input location for NEW analyses
│   │   ├── aop_annotations_master.csv # (Manual Input) Users fill this after Step 3
│   │   ├── master_mesh_database.db    # Local SQL database of PMIDs (Large file)
│   │   ├── d2025.bin                  # Raw MeSH ASCII descriptors (from NLM) (.7z format, must be unzipped)
│   │   └── 20250301_marc_full2025.bin # Raw MeSH MARC binary (from NLM) (.7z format, must be unzipped)
│   │
│   ├── processed/                     # Output location for NEW analyses (Auto-generated)
│   │
│   ├── reference_raw/                 # Read-only inputs for the "Allergic Contact Dermatitis" reference case
│   │
│   └── reference_processed/           # Outputs for the reference case
│
├── requirements.txt                   # Python dependencies
└── README.md                          # This file
```

---

## Initialization:

### 1. Clone Repository
```bash
git clone https://github.com/YourRepo/Mesh-Network-Analysis.git
cd Mesh-Network-Analysis
```

### 2. Install Dependencies
Ensure you have a modern version of Python (version 3.8 or higher is required).
```bash
pip install -r requirements.txt
```

### 3. Unzip and Setup Stop Words
The raw MeSH files must be unzipped before running the pipeline.
**Note:** If you are running the reference analysis, you can skip this step.
1. Navigate to `data/raw/`
2. Unzip `d2025.bin.7z` to `d2025.bin`
3. Unzip `20250301_marc_full2025.bin.7z` to `20250301_marc_full2025.bin`
If newer versions of these files are required, they can be downloaded from the National Library of Medicine (NLM):
- MARC: https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/meshmarc/
- ASCII: https://nlmpubs.nlm.nih.gov/projects/mesh/.asciimesh/

---

## Quick Start: Running the Reference Analysis
*This will reproduce the "Dermatitis, Allergic Contact" network described in the publication.*

### 1. Configuration
Open `scripts/python/config.py` and ensure the reference dataset flag is set to `True`:
```python
USE_REFERENCE_DATA = True
```

### 2. Execution
Run the pipeline:
```bash
python scripts/python/run_pipeline.py
```
*The pipeline will detect the existing processed files in `data/reference_processed/` and automatically skip heavy computation steps, jumping directly to figure generation.*

---

## Running a New Analysis (Fresh Network Analysis)

To analyze a new search term (e.g., "Liver Cirrhosis"), follow this workflow:

### 1. Configuration
Open `scripts/python/config.py` and update the following settings:
* `USE_REFERENCE_DATA = False`
* `SEARCH_TERM = "Your Search Query"`
* `ENTREZ_EMAIL` and `ENTREZ_API_KEY`: Required for PubMed access. You can obtain these by creating an account at NCBI (https://www.ncbi.nlm.nih.gov/myncbi/).
* Adjust all other configuration settings as necessary.

### 2. Execution
Execute the pipeline. It will process raw data, build the network, and export the initial results.
```bash
python scripts/python/run_pipeline.py
```
* **Step 1:** Checks MeSH raw data organization and creates the stop word list if it is not present.
* **Step 2:** Scrapes PubMed, builds the network, runs subgraph optimization, and calculates ARS and CRS. (This is a computationally intensive step).
* **Step 3:** Exports final network files, databases, and Excel files to the `results/` and `data/processed/` directories.
* **Step 4:** Attempts to generate figures. Note that full figure utility requires manual biological strata assignment (see below).

---

## Manual Biological Strata Assignment
**Critical:** The pipeline cannot automatically determine if a MeSH term represents a "Molecular", "Cellular", or "Tissue" level event. You must provide this biological context manually for new analyses.

**Instructions:**
1. After running the pipeline (specifically Step 3), navigate to the `results/` directory and open the newly generated Excel file (e.g., `[FILE_PREFIX]_final_network_with_relevance_export.xlsx`).
2. Open the **Nodes** sheet and review the extracted terms.
3. Open `data/raw/aop_annotations_master.csv`.
4. Copy the nodes from the generated Excel file to the `mesh_terms` column in the CSV.
5. Assign each term one of the following 7 Biological Strata in the `aop_level` column:
    * `Stressor`: External stimuli that initiate a biological reaction (e.g., "UV Rays").
    * `Molecular`: Gene or protein level names or events.
    * `Cellular`: Cellular level events (e.g., "Chemotaxis").
    * `Tissue`: Events localized to a subsection of tissue (e.g., "Necrosis").
    * `Organ`: Organ level names or events (e.g., "Liver").
    * `Adverse Outcome`: High-level disease outcomes (e.g., "Drug Hypersensitivity").
    * `Uncategorized`: Terms that do not fall into one of the above categories.
6. Save the updated `aop_annotations_master.csv` file.
7. Re-run the pipeline (or execute `scripts/python/figures.py` directly) to generate the final visualizations with your new annotations applied.

---

## Script Descriptions

### `run_pipeline.py`
A master pipeline orchestrator set up for convenience. It checks your environment and configuration, then executes the individual scripts in the correct sequential order, handling error checking and timing.

### `config.py`
This is the master configuration file. It contains parameters that dictate the behavior of the analysis. Notable options include:

* `USE_REFERENCE_DATA`: Determines if the scripts utilize the pre-computed reference dataset for "Dermatitis, Allergic Contact".
* `CUSTOM_FILE_PREFIX`: Sets the file header for all generated files (e.g., `{FILE_PREFIX}_full_network_data.json`).
* `ENTREZ_EMAIL` & `ENTREZ_API_KEY`: Required credentials for PubMed API access. See the NCBI documentation for more details.
* `SEARCH_TERM`, `START_DATE`, `END_DATE`: Defines the MEDLINE search query. The script concatenates these, so do not add dates directly into the `SEARCH_TERM`.
* `GENERATIONS_N`: The number of citation generations to retrieve. 0 or 1 retrieves only direct articles. Values >= 2 pull extended citation networks, finding more distant relationships at a significantly increased computational cost.
* `UPDATE_MESH_SUPPORT_FILES`: Set to `True` to force the regeneration of `mesh_stop_words.py` using the `.bin` files in `data/raw/`. Useful when updating to a newer MeSH release.

### Standard Execution Pipeline

#### 1. `mesh_data_processor.py`
* **Usage:** Executes automatically if `mesh_stop_words.py` is missing or if `UPDATE_MESH_SUPPORT_FILES = True` in `config.py`.
* **Action:** Parses the binary MARC and ASCII files to extract MeSH descriptors, identifies tree numbers, and generates the `mesh_stop_words.py` exclusion list based on biological categories of interest.

#### 2. `master_mesh_network.py`
* **Usage:** The computational core of the pipeline.
* **Action:** Queries PubMed, builds a rank-normalized co-occurrence graph, runs GLF and SA optimization algorithms, extracts the Largest Connected Component (LCC), performs Louvain community detection, and calculates Article Relevance Scores (ARS) and Contextual Relevance Scores (CRS).

#### 3. `secondary_analysis.py`
* **Usage:** Prepares data for manual review and visualization.
* **Action:** Converts the complex JSON network output into an Excel file for manual biological strata annotation. Provides utilities for running specific queries on nodes or edges to retrieve contributing PMIDs.

#### 4. `figures.py`
* **Usage:** Generates high-resolution visualizations in the `results/figures/` directory.
* **Action:** Produces edge weight distributions, optimization trajectories, community composition charts, t-SNE projections, AOP Alluvial flows, and centrality comparison plots. Requires `node2vec` for dendrogram generation.

---

## Notes & Troubleshooting

* **`master_mesh_database.db`:** This SQLite database stores fetched PMIDs to prevent redundant downloads. It can grow extremely large (>6GB). If deleted, the script will rebuild it automatically, but the initial run will be significantly slower.
* **Node2Vec:** The dendrogram generation in `figures.py` requires the `node2vec` package. If it is not installed or fails to import, the dendrogram figure will be skipped gracefully with a warning.
* **Memory Usage:** Network construction (Step 2) can be highly memory-intensive for broad search queries (>10^6 articles) or when using multiple citation generations. A minimum of 16GB RAM is recommended for large datasets.
* **Environment Configuration:** The pipeline expects to be run from a directory where `scripts/python/config.py` is accessible (typically the repository root). Ensure you are executing scripts using paths like `python scripts/python/run_pipeline.py`.

---

## Citation
If you use this code or methodology, please cite:
> *[To be updated after publication]*
> [![DOI](https://zenodo.org/badge/1145152319.svg)](https://doi.org/10.5281/zenodo.18662959)

## License
[MIT License]
