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
"""

__version__ = "2.2.0"

# --- 1. Configuration ---
from .config_parser import MeshConfig

# --- 2. Data Processing (MARC/ASCII) ---
from .mesh_data_processor import process_raw_mesh_data

# --- 3. Database Operations & Entrez Fetching ---
from .data_ops import (
    run_initial_data_collection,
    clean_database,
    populate_master_mesh_database,
    update_db_with_mesh_batch
)

# --- 4. Network Construction ---
from .network import (
    run_network_construction,
    run_consensus_filtering_and_lcc,
    run_community_detection
)

# --- 5. Filtering & Relevance Scoring ---
from .relevance import (
    run_mean_relevancy_scoring
)

# --- 6. Secondary Analysis & Export ---
from .secondary_analysis import (
    analyze_node_relevancy,
    analyze_edge_relevancy,
    get_top_network_articles,
    convert_network_json_to_excel,
    run_network_overlay_comparison
)

# --- 6b. Validation & Benchmarking ---
from .benchmark import (
    run_benchmark,
    validate_ground_truth,
    resolve_ground_truth_path,
    normalize_pmid
)
from .gt_network_validation import run_gt_network_validation

# --- 7. Visualization ---
from .viz import (
    analyze_dispersion,
    plot_cooccurrance_distribution,
    run_optimization_comparison,
    plot_louvain_community_bars,
    plot_joint_plot,
    plot_tsne_louvain_overlap,
    plot_sankey_alluvial,
    plot_dumbell_plot,
    plot_scatter_panels,
    plot_dendrogram,
    generate_filtering_summary_sankeys
)

# Define the public API of the package
__all__ = [
    "MeshConfig",
    "process_raw_mesh_data",
    "run_initial_data_collection",
    "clean_database",
    "populate_master_mesh_database",
    "update_db_with_mesh_batch",
    "run_network_construction",
    "run_community_detection",
    "run_consensus_filtering_and_lcc",
    "run_mean_relevancy_scoring",
    "analyze_node_relevancy",
    "analyze_edge_relevancy",
    "get_top_network_articles",
    "convert_network_json_to_excel",
    "run_network_overlay_comparison",
    "run_benchmark",
    "validate_ground_truth",
    "resolve_ground_truth_path",
    "normalize_pmid",
    "run_gt_network_validation",
    "analyze_dispersion",
    "plot_cooccurrance_distribution",
    "run_optimization_comparison",
    "plot_louvain_community_bars",
    "plot_joint_plot",
    "plot_tsne_louvain_overlap",
    "plot_sankey_alluvial",
    "plot_dumbell_plot",
    "plot_scatter_panels",
    "plot_dendrogram",
    "generate_filtering_summary_sankeys"
]