# -*- coding: utf-8 -*-
"""
Config.py
"""

import os
import sys
from datetime import date

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CONTROL FLAGS (USER SETTINGS)
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Set to TRUE to generate figures using the pre-curated (published) network data.
# Set to FALSE to use the new data and generate networks from scratch.
USE_REFERENCE_DATA = True
# Set desired file prefix for naming scheme of data
CUSTOM_FILE_PREFIX = "DAC_Mesh" # "DAC_Mesh" is used for reference data

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 1. USER PARAMETERS (EDIT THESE)
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# --- A. Entrez Credentials (REQUIRED) ---
# PLEASE REPLACE THESE BEFORE RUNNING (Required if USE_REFERENCE_DATA = False)
ENTREZ_EMAIL   = #"YOUR_EMAIL_HERE"
ENTREZ_API_KEY = #"YOUR_API_KEY_HERE"

# --- B. Search Parameters ---
# This concatinates the SEARCH_TERM to the DATES
SEARCH_TERM   = "Dermatitis, Allergic Contact [Mesh]"
START_DATE    = "1950/01/01"
END_DATE      = "2025/01/01"
# Default = 1 (i.e. no children) Generations of incoming and outgoing citations to expand pool.
GENERATIONS_N = 1
# Set to TRUE if updating stop words with mesh_data_processor.py (FALSE = Default).
UPDATE_MESH_SUPPORT_FILES = False

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 2. DYNAMIC PARAMETERS (EDIT THESE IF YOU FEEL LIKE IT)
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# --- C. Network Analysis Parameters ---
# [!] WARNING: Setting this to false turns off CRS, betweeeness, Eigenvector scoring
# Set to FALSE to significantly speeds up Step 4
CALCULATE_FULL_CENTRALITY = True # (Default == TRUE)
# A seed for random sampling.
RANDOM_SEED = 42
# The number of nodes to sample for betweenness centrality approximation.
BETWEENNESS_K_SAMPLES = 1000  # Lower this for faster betweenness calc. (1^6=edges & k=1000 is approx 80 min)

# --- D. Sub-network simulation Parameters (optimized for edges n=500)---
TARGET_NUM_EDGES = 500        # The target number of edges for the optimal subgraphs.
GLF_ITERATIONS   = 5_000_000  # Number of iterations for the GLF simulation.
SA_ITERATIONS    = 5_000_000  # Number of iterations for the Simulated Annealing simulation.
SA_TEMP_START    = 5000.0     # Initial temperature for SA
SA_COOLING_RATE  = 0.999995   # Cooling rate for SA.

# --- E. Contextual Relevance ---
#This allows for a refined look at relevant windows of relevance.
CONTEXT_START_DATE = START_DATE
# Note: You can change this to a fixed date if needed, e.g., "2025/01/01"
CONTEXT_END_DATE   = date.today().strftime("%Y/%m/%d")

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 3. DYNAMIC PROJECT ROOT DETECTION
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if "__file__" in globals():
    CONFIG_DIR = os.path.dirname(os.path.abspath(__file__))
else:
    CONFIG_DIR = os.path.abspath(os.getcwd())

# We are in 'Repo/scripts/python', so Root is TWO levels up
PROJECT_ROOT = os.path.abspath(os.path.join(CONFIG_DIR, '..', '..'))
# Fallback check
if not os.path.exists(os.path.join(PROJECT_ROOT, 'data')):
    PROJECT_ROOT = os.path.abspath(os.path.join(CONFIG_DIR, '..'))

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 4. DIRECTORY STRUCTURE
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')
RAW_DIR  = os.path.join(DATA_DIR, 'raw')
PROCESSED_DIR = os.path.join(DATA_DIR, 'processed')
REF_PROCESSED_DIR = os.path.join(DATA_DIR, 'reference_processed')
REF_RAW_DIR = os.path.join(DATA_DIR, 'reference_raw')

RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')
FIGURES_DIR = os.path.join(RESULTS_DIR, 'figures')
LOG_DIR     = os.path.join(RESULTS_DIR, 'logs')

os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 5. FILE NAMES & PATHS
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Adjust to your preferred file name output and determine source directory
if USE_REFERENCE_DATA:
    FILE_PREFIX = "DAC_Mesh"
    SOURCE_DIR = REF_PROCESSED_DIR
    SOURCE_RAW_DIR = REF_RAW_DIR
else:
    FILE_PREFIX = CUSTOM_FILE_PREFIX
    SOURCE_DIR = PROCESSED_DIR
    SOURCE_RAW_DIR = RAW_DIR

FILES = {
    # Raw Inputs  (Keep master_db location the same between runs for fastest results)
    "annotations":  os.path.join(SOURCE_RAW_DIR, "aop_annotations_master.csv"),
    "master_db":    os.path.join(SOURCE_RAW_DIR, "master_mesh_database.db"),
    "pmids_db":     os.path.join(SOURCE_RAW_DIR, f"{FILE_PREFIX}_pmids.db"),

    # Processed Outputs (Switchable Source)
    "full_network": os.path.join(SOURCE_DIR, f"{FILE_PREFIX}_full_network_data.json"),
    "glf_subgraph":  os.path.join(SOURCE_DIR, f"{FILE_PREFIX}_glf_optimal_subgraph.json"),
    "sa_subgraph":   os.path.join(SOURCE_DIR, f"{FILE_PREFIX}_sa_optimal_subgraph.json"),
    "consensus_lcc": os.path.join(SOURCE_DIR, f"{FILE_PREFIX}_consensus_lcc_network.json"),
    "final_network": os.path.join(SOURCE_DIR, f"{FILE_PREFIX}_final_network_with_relevance.json"),

    # Processed DBs
    "cleaned_db":    os.path.join(SOURCE_DIR, f"{FILE_PREFIX}_cleaned_pmids.db"),
    "relevance_db":  os.path.join(SOURCE_DIR, f"{FILE_PREFIX}_contextual_relevance.db"),

    # Logs
    "failed_mesh":   os.path.join(LOG_DIR, f"{FILE_PREFIX}_failed_mesh_fetches.tsv"),
    "empty_mesh":    os.path.join(LOG_DIR, f"{FILE_PREFIX}_empty_mesh_pmids.tsv"),
    "log_file":      os.path.join(LOG_DIR, f"{FILE_PREFIX}_processing_errors.log"),

    # MeSH Inputs - Update these to get the latest MeSH terms. (Download these from https://nlmpubs.nlm.nih.gov/projects/mesh/)
    "mesh_marc":  os.path.join(SOURCE_RAW_DIR, "20250301_marc_full2025.bin"), # https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/meshmarc/
    "mesh_ascii": os.path.join(SOURCE_RAW_DIR, "d2025.bin"), # https://nlmpubs.nlm.nih.gov/projects/mesh/.asciimesh/

    # New MeSH Outputs
    "mesh_terms_csv":    os.path.join(SOURCE_DIR, "mesh_terms.csv"),
    "mesh_stopwords_py": os.path.join(PROJECT_ROOT, "scripts", "python", "mesh_stop_words.py"), # Saves directly to scripts folder
}
