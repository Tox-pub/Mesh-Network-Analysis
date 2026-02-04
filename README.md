# MeSH Network Analysis for AOP Discovery

## Overview
This repository contains a computational pipeline designed to reconstruct **Adverse Outcome Pathways (AOPs)** using **Medical Subject Headings (MeSH)** co-occurrence networks.

By analyzing the co-occurrence of MeSH terms across thousands of PubMed articles, this tool builds a weighted network graph that connects **Stressors** (e.g., chemicals) to **Adverse Outcomes** (diseases) through biological intermediates. The pipeline utilizes graph theory, optimization algorithms (GLF & Simulated Annealing), and a novel **Contextual Relevance Score (CRS)** to filter noise and identify the most biologically relevant pathways.

**Key Capabilities:**
* **Automated Data Mining:** Fetches and processes PubMed citations via Entrez API.
* **Network Construction:** Builds rank-normalized co-occurrence networks from MeSH terms.
* **Optimization:** Uses **Graph Likelihood Filter (GLF)** and **Simulated Annealing (SA)** to extract optimal subgraphs.
* **Contextual Scoring:** Calculates CRS to rank nodes based on their bridging centrality between research domains.
* **Visualization:** Generates publication-ready figures (Sankey diagrams, t-SNE, Alluvial flows).

---

## 📂 Project Structure

```text
Project-Root/
│
├── scripts/
│   └── python/
│       ├── config.py                 # Central configuration (API keys, search terms)
│       ├── run_pipeline.py           # Orchestrator script to run the full workflow
│       ├── mesh_stop_words.py        # List of excluded terms (Geography, Publication types)
│       ├── mesh_data_processor.py    # Utilities to process raw MeSH binary files
│       ├── master_mesh_network.py    # CORE SCRIPT: Data fetching, network building, simulation
│       ├── secondary_analysis.py     # Exports data to Excel & runs node/edge queries
│       └── figures.py                # Generates all visualization figures
│
├── results/
│   ├── figures/                      # Final high-res TIFF/SVG/PNG output
│   └── logs/                         # Execution logs and error reports
│
├── data/
│   ├── raw/                          # Input location for NEW analyses
│   │   ├── aop_annotations_master.csv # (Manual Input) Users fill this after Step 3
│   │   ├── master_mesh_database.db    # Local SQL database of PMIDs (Large file)
│   │   ├── d2025.bin                  # Raw MeSH ASCII descriptors (from NLM)
│   │   └── 20250301_marc_full2025.bin # Raw MeSH MARC binary (from NLM)
│   │
│   ├── processed/                    # Output location for NEW analyses (Auto-generated)
│   │
│   ├── reference_raw/                # Read-only inputs for the "Dermatitis" reference case
│   │   └── [Reference DBs and Annotation files]
│   │
│   └── reference_processed/          # "Golden Copy" outputs for the reference case
│       └── [Pre-computed Networks, JSONs, and CSVs]
│
├── requirements.txt                  # Python dependencies
└── README.md                         # This file

---

## 🚀 Quick Start: Running the Reference Analysis
*Goal: Reproduce the "Dermatitis, Allergic Contact" network described in the publication.*

1.  Open `scripts/python/config.py`.
2.  Ensure the reference flag is set to **True**:
    ```python
    USE_REFERENCE_DATA = True
    ```
3.  Run the pipeline:
    ```bash
    python run_pipeline.py
    ```
    *The pipeline will detect the existing processed files in `data/reference_processed/` and automatically skip heavy computation steps, jumping directly to figure generation.*

---

## 🧪 Running a New Analysis (Fresh Data)

To analyze a **new search term** (e.g., "Liver Cirrhosis"), follow this 4-step workflow:

### 1. Configuration
Open `scripts/python/config.py` and update the following:
* `USE_REFERENCE_DATA = False`
* `SEARCH_TERM = "Your Search Query [Mesh]"`
* `ENTREZ_EMAIL` and `ENTREZ_API_KEY` (Required for PubMed access).

### 2. Execution (Steps 1-3)
Run the pipeline. It will process raw data, build the network, and export the initial results.
```bash
python run_pipeline.py


### Section 2: Script Descriptions

```markdown
* **Step 1:** Checks MeSH raw data.
* **Step 2:** Scrapes PubMed, builds the network, and runs GLF/SA optimization (Computationally intensive).
* **Step 3:** Exports the final network to an Excel file in `results/`.

### 3. ⚠️ Manual Biological Strata Assignment
**Crucial Step:** The pipeline cannot automatically determine if a MeSH term represents a "Molecular", "Cellular", or "Tissue" level event. You must provide this context.

1.  Go to the `results/` folder and open the newly generated Excel file:
    * `[FILE_PREFIX]_final_network_with_relevance_export.xlsx`
2.  Open the **Nodes** sheet.
3.  Review the high-scoring nodes.
4.  Open `data/raw/aop_annotations_master.csv`.
5.  Add your new terms and assign them one of the following **7 Biological Strata**:
    * `Stressor`
    * `Molecular`
    * `Cellular`
    * `Tissue`
    * `Organ`
    * `Adverse Outcome`
    * `Uncategorized`
6.  Save the CSV.

### 4. Generate Figures
Run the pipeline again (or just the figures script) to generate the final visualizations with your new annotations applied.
---

## 📜 Script Descriptions

### `run_pipeline.py`
The master orchestrator. It checks your configuration and executes the following scripts in the correct order, handling error checking and timing.

### 1. `mesh_data_processor.py`
* **Input:** Raw NLM files (`.bin` MARC and ASCII).
* **Action:** Extracts unique MeSH descriptors and generates the `mesh_stop_words.py` file based on tree numbers (excluding Geography, Publication Types, etc.).
* **Usage:** Only runs if stop words need updating or `UPDATE_MESH_SUPPORT_FILES = True` in config.

### 2. `master_mesh_network.py`
The computational core of the project.
* **Data Collection:** Queries PubMed/Entrez for the search term.
* **Network Building:** Constructs a co-occurrence graph where nodes are MeSH terms and edges are weighted by frequency (Rank-Normalized).
* **Optimization:** Runs **GLF (Graph Likelihood Filter)** and **SA (Simulated Annealing)** to identify the most statistically significant subgraph.
* **LCC Extraction:** Extracts the Largest Connected Component.
* **Community Detection:** Applies the Louvain algorithm.
* **CRS Calculation:** Computes Contextual Relevance Scores (Betweenness & Eigenvector centrality within the subgraph).

### 3. `secondary_analysis.py`
* **Action:** Converts the complex JSON network output into a user-friendly Excel file (`_export.xlsx`).
* **Analysis:** Can be configured to run specific queries on single Nodes or Edges to retrieve the exact PMIDs contributing to those connections.

### 4. `figures.py`
Generates high-resolution visualizations found in the `results/figures/` folder:
* **Figure 1:** Edge weight distribution (Power law analysis).
* **Figure 2:** Optimization Trajectory (GLF vs SA convergence).
* **Figure 3:** Community Composition bar charts.
* **Figure 5:** t-SNE projection of the network.
* **Figure 6:** AOP Alluvial/Sankey flow (The primary AOP visualization).
* **Figure 7/8:** Centrality comparisons (Dumbbell plots and Scatter panels).

---

## 📝 Notes & Troubleshooting

* **`master_mesh_database.db`**: This SQLite database stores fetched PMIDs to prevent re-downloading millions of citations. It can grow large (>6GB). If you delete it, the script will rebuild it, but the first run will be significantly slower.
* **Node2Vec**: The dendrogram generation in `figures.py` requires the `node2vec` library. If not installed, that specific figure is skipped gracefully.
* **Memory Usage**: Step 2 (Network Construction) can be memory intensive for broad search terms (>50,000 articles). Ensure you have at least 16GB RAM for large datasets.

---

## Citation
If you use this code or methodology, please cite:
> *[Insert Citation Here]*

## License
[MIT License / Your License Choice]
