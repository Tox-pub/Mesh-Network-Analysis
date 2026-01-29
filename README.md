# Mesh-Network-Analysis

This repository contains the source code and datasets for the manuscript titled:
**"[Insert Your Paper Title Here]"**

## Repository Structure

The code is organized sequentially to reproduce the analysis pipeline:

1.  **`scripts/01_network_analysis_pipeline.py`**:
    * Processes the raw data found in `data/raw/`.
    * Performs GLF and Simulated Annealing optimization.
    * Outputs graph structures (JSON) to `data/processed/`.

2.  **`scripts/02_generate_figures.py`**:
    * Loads the processed data from `data/processed/`.
    * Generates the Sankey diagrams, Dendrograms, and Distribution plots.
    * Saves outputs to `results/figures/`.

## Data Availability

* **Raw Data**: Located in `data/raw/`.
* **Intermediate Networks**: The optimized subgraphs (GLF/SA) are available in `data/processed/` to allow for immediate reproduction of figures.

## Instructions for Reproduction

### 1. Environment Setup
To replicate the environment, install the required dependencies:

```bash
pip install -r requirements.txt
