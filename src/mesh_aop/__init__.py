# -*- coding: utf-8 -*-
"""
MeSH AOP Network Analysis - package initializer and public API.

This package builds and analyzes knowledge graphs that connect chemical
stressors to adverse outcomes through biological intermediates, using the
co-occurrence of Medical Subject Headings (MeSH) across the PubMed literature.

The pipeline runs in stages: process the MeSH vocabulary, collect article and
citation data from NCBI Entrez, construct and statistically filter the
co-occurrence network (GLF/SA optimization, Louvain communities, mean relevancy
scoring), validate against ground truth, and visualize the result as
Adverse Outcome Pathway flows.

This module re-exports the primary function of each stage as a flat, importable
API so the pipeline can be driven programmatically as well as from the CLI.

The re-exports are resolved on first use rather than at import time. Importing
them eagerly pulled numpy, scipy, pandas, gensim, scikit-learn, statsmodels,
matplotlib, plotly, Biopython and seaborn into the process - 2,812 modules and
about 16 seconds - before anything had been asked for. That cost fell on every
consumer equally, including ones that touch none of it: the desktop application
imports this package only to look up where its configuration file lives, and
paid 16 seconds at startup for the privilege.

Nothing about the public API changes. `from mesh_aop import run_benchmark` and
`mesh_aop.run_benchmark` both work exactly as before; the submodule is imported
the first time the name is actually read, and cached in the module globals so
later reads cost nothing.
"""

from importlib import import_module

__version__ = "3.2.10"

# Public name -> the submodule that defines it. Keeping this as data is what
# makes the deferral possible; it is also the whole public surface in one place.
_EXPORTS = {
    # --- 1. Configuration ---
    "MeshConfig": "config_parser",
    # --- 2. Data Processing (MARC/ASCII) ---
    "process_raw_mesh_data": "mesh_data_processor",
    "ensure_mesh_descriptor_xml": "mesh_data_processor",
    # --- 3. Database Operations & Entrez Fetching ---
    "run_initial_data_collection": "data_ops",
    "clean_database": "data_ops",
    "populate_master_mesh_database": "data_ops",
    "update_db_with_mesh_batch": "data_ops",
    # --- 4. Network Construction ---
    "run_network_construction": "network",
    "run_consensus_filtering_and_lcc": "network",
    "run_community_detection": "network",

    # --- Run provenance and the PRISMA flow report ---
    "RunLedger": "run_ledger",
    "open_ledger": "run_ledger",
    "ledger_path": "run_ledger",
    "write_prisma_report": "prisma",
    "build_flow": "prisma",
    # --- 5. Filtering & Relevance Scoring ---
    "run_mean_relevancy_scoring": "relevance",
    # --- 6. Secondary Analysis & Export ---
    "analyze_node_relevancy": "secondary_analysis",
    "analyze_edge_relevancy": "secondary_analysis",
    "get_top_network_articles": "secondary_analysis",
    "convert_network_json_to_excel": "secondary_analysis",
    "run_network_overlay_comparison": "secondary_analysis",
    # --- 6b. Validation & Benchmarking ---
    "run_benchmark": "benchmark",
    "validate_ground_truth": "benchmark",
    "resolve_ground_truth_path": "benchmark",
    "normalize_pmid": "benchmark",
    "run_gt_network_validation": "gt_network_validation",
    # --- 7. Visualization ---
    "analyze_dispersion": "viz",
    "plot_cooccurrance_distribution": "viz",
    "run_optimization_comparison": "viz",
    "plot_louvain_community_bars": "viz",
    "plot_tsne_louvain_overlap": "viz",
    "plot_sankey_alluvial": "viz",
    "plot_dendrogram": "viz",
    "generate_filtering_summary_sankeys": "viz",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name):
    """Import the submodule that provides `name`, the first time it is read.

    PEP 562. This is what `from mesh_aop import run_benchmark` falls back to
    when the name is not already a module global, so the flat API is unchanged.
    """
    module = _EXPORTS.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(f".{module}", __name__), name)
    globals()[name] = value          # subsequent reads skip this path entirely
    return value


def __dir__():
    return sorted(list(globals()) + __all__)
