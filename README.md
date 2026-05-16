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
│   ├── figures/                      # Final high-res TIFF/SVG/PNG output
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

*Note: The script currently assumes `run_pipeline.py` is being run from the main directory so `config.py` knows where to look.*
* Execute `run_pipeline.py` or `mesh_data_processor.py`.

---

## Quick Start: Running the Reference Analysis
*This will reproduce the "Dermatitis, Allergic Contact" network described in the publication.*

### 1.  Open `scripts/python/config.py` and configure user defined settings.
### 2.  Ensure the reference dataset flag is set to **True**:
```python
USE_REFERENCE_DATA = True
```
### 3.  Run the pipeline:
```bash
python run_pipeline.py
```
*The pipeline will detect the existing processed files in `data/reference_processed/` and automatically skip heavy computation steps, jumping directly to figure generation.*

---

## Running a New Analysis (Fresh Network Analysis)

To analyze a **new search term** (e.g., "Liver Cirrhosis"), follow this workflow:

### 1. Configuration
Open `scripts/python/config.py` and update the following settings:
* `USE_REFERENCE_DATA = False`
* `SEARCH_TERM = "Your Search Query"`
* `ENTREZ_EMAIL` and `ENTREZ_API_KEY`: Required for PubMed access. You can obtain these by creating an account at NCBI (https://www.ncbi.nlm.nih.gov/myncbi/).
* Adjust all other configuration settings to your liking.

### 2. Execution
You can choose to run the pipeline at this step or run the scripts in order of execution independently. It will process raw data, build the network, and export the initial results.
```bash
python run_pipeline.py
```
* **Step 1:** Checks MeSH raw data organization and creates the stop word list if it is not present.
* **Step 2:** Scrapes PubMed, builds the network, runs subgraph optimization, and calculates ARS and CRS (Computationally intensive step).
* **Step 3:** Exports final network files, databases, and Excel files to the `results/`, `processed/`, and `raw/` directories.

---

## Script Descriptions

### `run_pipeline.py`
A master pipeline orchestrator set up for convenience. It checks your configuration and environment, then executes the following scripts in the correct order, handling error checking and timing.

### `config.py`
This is the master operator control file. There are several options here the user should check and edit depending on their intended analysis.
Some notable options are as follows:

```bash
USE_REFERENCE_DATA
```
This determines if the `run_pipeline.py` and specifically `master_mesh_network.py` uses the reference dataset for search string `Dermatitis, Allergic Contact [Mesh] AND 1950/01/01[EDAT] : 2025/01/01[EDAT]`.

```bash
CUSTOM_FILE_PREFIX
```
This will be the file header for all files generated across all scripts e.g. `{FILE_PREFIX}_full_network_data.json`. For the reference dataset `DAC` is used as the `FILE_PREFIX`.

```bash
ENTREZ_EMAIL
ENTREZ_API_KEY
```
Required for PubMed API access. Login is located at: `https://www.ncbi.nlm.nih.gov/myncbi/`. More information can be found at `https://www.ncbi.nlm.nih.gov/books/NBK25501/`.

```bash
SEARCH_TERM
START_DATE
END_DATE
```
This is the desired MEDLINE search string; it can be any string that can be put into PubMed `https://pubmed.ncbi.nlm.nih.gov/`. This is read in the script as:
```bash
{SEARCH_TERM} AND {START_DATE}[EDAT] : {END_DATE}[EDAT]
```
So there is no need to add dates to the `SEARCH_TERM` argument.

```bash
GENERATIONS_N
```
This is the number of generations the `master_mesh_network.py` pulls for incoming and outgoing citations under the premise that scientific literature will reference related works.
0 or 1 leads to just the articles in the initial search string returning. >= 2 generations pulls more possible relevant relationships between terms but at a significant computational cost.

```bash
UPDATE_MESH_SUPPORT_FILES
```
This script relies on `20250301_marc_full2025.bin` and `d2025.bin` to create the stop words in `mesh_stop_words.py`. This means post 01/01/2025 they may not contain new tree assignments or new MeSH terms created after this date. This can be easily updated by downloading the files from `https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/meshmarc/` and `https://nlmpubs.nlm.nih.gov/projects/mesh/.asciimesh/` respectively, replacing these in `data/raw`, and running `mesh_data_processor.py` or `run_pipeline.py` with:
    ```python
    UPDATE_MESH_SUPPORT_FILES = True
    ```

## Standard Execution Pipeline

### 1. `mesh_data_processor.py`
* **Usage:** Only runs if `mesh_stop_words.py` is missing or `UPDATE_MESH_SUPPORT_FILES = True` in `config.py`, allowing the user to update quickly. The created `mesh_terms.csv` is a manual diagnostic tool to check valid terms listed.
    * **Input:** Raw NLM files (`.bin` MARC and ASCII).
    * **Action:** Parses the binary MARC file to extract all unique MeSH descriptors and UIs into a CSV.
    * Parses the ASCII descriptor file to identify tree numbers.
    * Generates the `mesh_stop_words.py` file by filtering out terms that fall outside the biological categories of interest (e.g., Geography, Publication Types, etc.).

### 2. `master_mesh_network.py`
* **Usage:** The computational core. Gathers and parses MeSH terms, creates network JSONs, and calculates centrality metrics, contextual relevance scores, etc.
    * **Data Collection:** Queries PubMed/Entrez for the search term defined in `config.py` and retrieves all relevant PMIDs.
    * **Network Construction:** Builds a rank-normalized co-occurrence graph where nodes are MeSH terms and edges are weighted by the frequency of co-occurrence across the article set.
    * **Optimization:** Runs two competing optimization algorithms, Global Likelihood Filter (GLF) and Simulated Annealing (SA), to identify the most statistically significant subgraph by removing noise edges, and further finding the consensus of these two models.
    * **LCC Extraction:** Extracts the Largest Connected Component (LCC) from the optimized consensus subgraph.
    * **Community Detection:** Applies the Louvain algorithm to identify clusters of related terms.
    * **ARS Calculation:** Computes Article Relevance Scores (ARS) by utilizing Betweenness and Eigenvector centrality to determine the article's influence within the subgraph.
    * **CRS Calculation:** Computes Contextual Relevance Scores (CRS) by using the cumulative ARS per node and weighting the quantity of evidence in the full literature corpus.

### 3. `secondary_analysis.py`
* **Usage:** Allows for several post-computational analyses on the final filtered network created in `master_mesh_network.py`. This includes an Excel file.
    * **Action:** Converts the complex JSON final network output from the master script into an Excel file (`_export.xlsx`) with separate sheets for Nodes and Edges that can be used to update the biological annotations in `aop_annotations_master.csv`.
    * **Analysis:** Can be configured to run specific queries on single Nodes or Edges to retrieve the exact PMIDs contributing to those connections, aiding in the manual verification of pathways and literature identification.

### 4. `figures.py`
* **Usage:** Generates high-resolution visualizations found in the `results/figures/` folder.

**[!] WARNING:** If running a new analysis, you will need to **manually assign biological strata** in `aop_annotations_master.csv` using nodes from `{FILE_PREFIX}_final_network_with_relevance_export.xlsx` to get the full utility out of all figure metrics.
  * **Figure 1:** Edge weight distribution (Power law analysis) to assess network topology.
  * **Figure 2:** Optimization Trajectory (GLF vs SA convergence) to validate the filtering process.
  * **Figure 3:** Community Composition bar charts showing the biological makeup of each cluster.
  * **Figure 5:** t-SNE projection of the network colored by Louvain community.
  * **Figure 6:** AOP Alluvial/Sankey flow (The primary visualization connecting Stressors to Outcomes).
  * **Figure 7/8:** Centrality comparisons (Dumbbell plots and Scatter panels) comparing CRS to raw centrality metrics.

---

### Manual Biological Strata Assignment
**Critical:** The pipeline cannot automatically determine if a MeSH term represents a "Molecular", "Cellular", or "Tissue" level event. You must provide this biological context.
*Manually annotating thousands of MeSH terms by hand was not within the scope of the project, but if you require this, please let us know.*

**Instructions:**
1.  Navigate to the `results/` folder and open the newly generated Excel file:
    * `[FILE_PREFIX]_final_network_with_relevance_export.xlsx`
2.  Open the **Nodes** sheet.
3.  Review the extracted nodes.
4.  Open `data/raw/aop_annotations_master.csv` and copy nodes from the Excel file to the `mesh_terms` column.
5.  Add your new terms and assign them one of the following **7 Biological Strata**:
    * `Stressor`: Usually external stimuli that initiates a biological reaction (e.g., `UV Rays`).
    * `Molecular`: Usually gene or protein level names or events.
    * `Cellular`: Cellular level events (e.g., `Chemotaxis`).
    * `Tissue`: Events localized to a subsection of tissue (e.g., `Necrosis`).
    * `Organ`: Organ level names or events (e.g., `Liver`).
    * `Adverse Outcome`: High level disease outcomes (e.g., `Drug Hypersensitivity`).
    * `Uncategorized`: Terms that do not fall into one of the above categories.
6.  Save the CSV as comma delimited.
7.  Run the pipeline again (or just `figures.py`) to generate the final visualizations with your new annotations applied.

---

## Notes & Troubleshooting

* **`master_mesh_database.db`**: This SQLite database stores fetched PMIDs to prevent re-downloading millions of citations. It can grow extremely large (>6GB). If you delete it, the script will rebuild it, but the first run will be significantly slower.
* **Node2Vec**: The dendrogram generation in `figures.py` requires the `node2vec` library. If it is not installed or fails to import, the dendrogram figure will be skipped gracefully with a warning in the console.
* **Memory Usage**: Step 2 (Network Construction) can be highly memory intensive for broad search terms (>10^6 articles) or when using multiple citation generations. Ensure you have at least 16GB RAM for large datasets.
* **Environment Configuration**: Ensure you are running `run_pipeline.py` from the main project directory. The pipeline expects to find `scripts/python/config.py` correctly located.

---

## Citation
If you use this code or methodology, please cite:
> *[To be updated after publication]*
> [![DOI](https://zenodo.org/badge/1145152319.svg)](https://doi.org/10.5281/zenodo.18662959)

## License
[MIT License]
