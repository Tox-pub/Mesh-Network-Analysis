# MeSH Workbench — Command-Line Guide

This document is how to run the pipeline from a shell, and how to work on the
source. It ships with the source package.

**Most people do not need it.** MeSH Workbench is a desktop application, and
everything the pipeline does can be driven from its window. If that is what you
want, see [INSTALL.md](INSTALL.md) to install it and
[HELP.md](HELP.md) for the manual — HELP.md is also what the application's
**Help → MeSH Workbench Manual** entry opens.

Use this document if you are scripting runs, working on the code, or running on
a machine with no desktop.

**Every setting named here is the same setting the application shows.** The
config keys below (`search_parameters.search_term`, `benchmark.primary_node`,
and so on) are the dotted paths into `mesh_config.json`; the application writes
the same file. The parameter glossary in [HELP.md](HELP.md) documents what each
one means, once, for both.

---

## Contents

- [Running the pipeline from a downloaded build](#running-the-pipeline-from-a-downloaded-build)
- [Setting up a development environment](#setting-up-a-development-environment)
- [Execution Guide](#execution-guide)
- [Repository Structure](#repository-structure)
- [Jupyter Notebook Interface](#jupyter-notebook-interface)
- [Programmatic API Usage](#programmatic-api-usage)
- [Where to find everything else](#where-to-find-everything-else)

## Running the pipeline from a downloaded build

Every downloadable build carries its own Python. The command line requires no
system Python, no virtual environment and no `pip install`. Nothing in the
section after this one applies to a downloaded build; it covers installing from
source.

Each build provides a launcher named `mesh-pipeline`, placed beside the
application in the program folder. It accepts every flag listed under
[CLI Flags](#cli-flags) and passes them through unchanged.

### Windows

Extract the portable zip, or install the `.msi`. Open a shell in the program
folder and run:

```powershell
mesh-pipeline.bat --step viz
mesh-pipeline.bat --step all --interactive
```

From another directory, give the full path in quotes — the folder name contains
a space:

```powershell
& "C:\path\to\MeshWorkbench\mesh-pipeline.bat" --step viz
```

**If the launcher does not run**, call the module directly. The command is
identical in effect; the launcher adds nothing but the path lookup:

```powershell
& "C:\path\to\MeshWorkbench\python\python.exe" -m mesh_aop.cli --step viz
```

That second form also applies where policy blocks `.bat` execution, since it
invokes the bundled interpreter rather than a script.

### Linux and macOS

Extract the tarball, change into the extracted folder, and run:

```bash
./mesh-pipeline --step viz
./mesh-pipeline --step all --interactive
```

The first invocation unpacks the bundled wheels, which takes about a minute and
requires no network. Subsequent invocations start immediately.

**If the launcher reports "Permission denied"**, the executable bit was lost in
transit — extracting on Windows or copying through a FAT32 or exFAT volume does
this. Restore it:

```bash
chmod +x mesh-pipeline "MeSH Workbench" mesh-uninstall python/bin/python3.12
```

**If the launcher still does not run**, call the module directly. Run the
launcher once first, or the bundled libraries will not have been unpacked:

```bash
PYTHONPATH="$PWD/app/src" ./python/bin/python3.12 -m mesh_aop.cli --step viz
```

### Startup time

The pipeline loads NumPy, SciPy, pandas, scikit-learn, statsmodels and
matplotlib before it does anything, so every invocation carries a fixed cost of
roughly four seconds on an internal disk.

Run the program from an internal disk. Measured on the same build, a command
taking four seconds from an internal SSD took over two minutes from an external
drive — the libraries total about 460 MB, and that volume is read on every
start. Extracting to a USB stick is the usual cause of an apparently frozen
command line.

---

## Setting up a development environment

Installing the **application** is covered in [INSTALL.md](INSTALL.md) — that is
the right document for anyone who simply wants to run it, and none of the below
is needed. What follows is for working on the source.

### Requirements

* **Python 3.11–3.13** (`requires-python = ">=3.11,<3.14"`).
* **Memory** — 16 GB is the ideal minimum; 32 GB or more for the database
  build and for networks past one citation generation. It will run on less.
* **Storage** — 100 GB+ free. The NLM baseline archives are about 50 GB and
  the SQLite master database it builds is about 10 GB. Once that database is
  built and verified you can delete `data/raw/pubmed_baseline/` and reclaim
  the 50 GB, provided you do not intend to apply daily updates later — those
  re-read the archives.

### Install from source

This installs the **command-line pipeline** from source. To install the
application instead, see [INSTALL.md](INSTALL.md); the application is described
under [The Workbench Window](HELP.md#the-workbench-window) below.

The commands assume **PowerShell** on Windows or **bash** on Linux and macOS.
Adjust the paths and the activation command for another shell.

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

## Repository Structure

The package assumes and enforces the following directory architecture.

```text
Mesh-Network-Analysis/
│
├── data/                               # Data storage
│   ├── raw/                            # Inputs for a run
│   │   ├── aop_annotations_master.csv  # Ships w/ repo: AOP strata dictionary (pre-seeded; grows each run)
│   │   ├── desc2025.xml                # Auto-downloaded from NLM if missing (or place manually); not in repo
│   │   ├── ground_truth_pmids.template.csv # Ships w/ repo: copy+fill for your own benchmark set
│   │   ├── ground_truth_pmids.csv      # Optional, you place this: YOUR benchmark set (see "Ground Truth")
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
│   ├── *_run_annotations.csv           # Run-specific strata annotation templates
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
│   │   ├── vocabulary.py               # Which MeSH trees an analysis may see
│   │   ├── strata.py                   # The annotation scheme and its order
│   │   ├── mdhtml.py                   # Renders the shipped documents as HTML
│   │   ├── network.py                  # NetworkX assembly, filtering, and centrality
│   │   ├── node2vec_embedding.py       # Node2Vec embedding used by the dendrogram figure
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

---

## Where to find everything else

This document covers only how to drive the pipeline from a shell. Everything
about what the pipeline does, and what each setting means, is in
[HELP.md](HELP.md):

| Looking for | See |
| :--- | :--- |
| What every setting does | **Configuration Wizard Parameter Glossary** |
| How ARS and MRS are calculated | **How articles and terms are scored** |
| Prerequisites and disk budget | **Data Acquisition & Prerequisites** |
| Assigning strata | **The Annotation Workflow (Strata)** |
| What the run produces | **Output Artifacts** |
| Ground truth and benchmarking | **Ground Truth**, **Validation & Benchmarking** |
| Damaged files | **When Files Go Wrong** |
| Citing this program | **Citation** |
