# -*- coding: utf-8 -*-
"""
cli.py - command-line orchestrator for the MeSH AOP Network pipeline (`mesh-pipeline`).

Parses command-line arguments and runs the pipeline either end to end or one
stage at a time, wiring together configuration, MeSH processing, NCBI data
collection, network construction, relevance scoring, secondary analysis,
benchmarking, and visualization.

Stages are selected with `--step` (process, data_ops, network, secondary, viz,
benchmark, or all); `--interactive` launches the configuration wizard first.
The orchestrator also enforces per-step prerequisites, manages the optional AOP
annotation pause, and generates the run-specific annotation templates.

Usage:
    mesh-pipeline --step all
    mesh-pipeline --step process --interactive
    mesh-pipeline --step viz
    mesh-pipeline --step benchmark
"""

import os
import sys
import time
import argparse
import json
import subprocess
import pandas as pd
from importlib.metadata import version, PackageNotFoundError

from .config_parser import MeshConfig
from .wizard import run_interactive_wizard
from .mesh_data_processor import process_raw_mesh_data
from .baseline_manager import PubMedBaselineManager
from .data_ops import (
    run_initial_data_collection,
    clean_database,
    populate_master_mesh_database,
    update_db_with_mesh_batch
)
from .network import run_network_construction, run_consensus_filtering_and_lcc, run_community_detection
from .relevance import run_contextual_relevance_scoring
from .secondary_analysis import convert_network_json_to_excel, analyze_node_relevancy, analyze_edge_relevancy, get_top_network_articles
from .benchmark import run_benchmark, resolve_ground_truth_path, KNOWN_GT_FILENAMES
from .viz import (
    load_and_prepare_data, load_full_raw_data, analyze_dispersion,
    plot_cooccurrance_distribution, run_optimization_comparison,
    plot_louvain_community_bars, plot_joint_plot, plot_tsne_louvain_overlap,
    plot_sankey_alluvial, plot_dumbell_plot, plot_scatter_panels, plot_dendrogram
)

try:
    from .mesh_stop_words import MESH_STOP_WORDS
except ImportError:
    MESH_STOP_WORDS = set()

def get_version() -> str:
    """Dynamically retrieves the package version established in pyproject.toml."""
    try:
        return version("mesh_aop_network")
    except PackageNotFoundError:
        return "Unknown (Package not installed via pip/mamba)"

def open_readme():
    """Locates and opens the README.md file using the default OS viewer."""
    # Attempt to locate the README one directory up from the src/ folder
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    readme_path = os.path.join(base_dir, 'README.md')

    # Fallback to current working directory if not found
    if not os.path.exists(readme_path):
        readme_path = 'README.md'

    if not os.path.exists(readme_path):
        print("[!] Error: README.md could not be located in the project root.")
        return

    print(f"Opening documentation: {readme_path}")

    # Cross-platform execution
    try:
        if sys.platform == "win32":
            os.startfile(readme_path)
        elif sys.platform == "darwin":
            subprocess.call(["open", readme_path])
        else:
            subprocess.call(["xdg-open", readme_path])
    except Exception as e:
        print(f"[!] Error opening file: {e}")

def _is_valid_db(db_path, table_name="pmids_table"):
    """Return True if the SQLite file exists, is non-empty, and contains the expected table."""
    if not os.path.exists(db_path) or os.path.getsize(db_path) == 0:
        return False
    try:
        import sqlite3
        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        cur.execute(f"SELECT 1 FROM sqlite_master WHERE type='table' AND name='{table_name}'")
        valid = cur.fetchone() is not None
        conn.close()
        return valid
    except Exception:
        return False

def _ensure_prerequisites(required_files: dict, step_name: str):
    """Abort with a clear message listing any required input files that are missing for a step."""
    missing = []
    for label, filepath in required_files.items():
        if not os.path.exists(filepath):
            missing.append((label, filepath))

    if missing:
        print("\n" + "<"*30 + ">"*30)
        print(f"[CRITICAL ERROR] Missing Prerequisites for {step_name}")
        print("<"*30 + ">"*30)
        print("The following required files were not found:")
        for label, filepath in missing:
            print(f"  - {label}: {filepath}")
        sys.exit(1)

def _initialize_master_annotations(mesh_csv_path: str, master_anno_path: str):
    """Pre-populates the Master Annotations Dictionary with all valid MeSH terms from the XML extraction."""
    if os.path.exists(master_anno_path):
        return
    print("\n<<< Initializing Master AOP Annotations Library >>>")
    try:
        df = pd.read_csv(mesh_csv_path)
        term_col = 'mesh_term' if 'mesh_term' in df.columns else df.columns[0]
        all_terms = df[term_col].dropna().astype(str).unique()

        # Filter out stop words so we only have valid biological/chemical terms
        valid_terms = [t for t in all_terms if t.lower() not in {s.lower() for s in MESH_STOP_WORDS}]

        anno_df = pd.DataFrame({'mesh_term': valid_terms, 'aop_level': ['Unassigned'] * len(valid_terms)})
        os.makedirs(os.path.dirname(master_anno_path), exist_ok=True)
        anno_df.to_csv(master_anno_path, sep=';', index=False)
        print(f"  [+] Created Master Annotations Dictionary with {len(valid_terms):,} valid MeSH terms.")
    except Exception as e:
        print(f"  [!] Failed to initialize master annotations: {e}")

def _generate_run_annotations(json_path: str, master_anno_path: str, run_anno_path: str):
    """Creates a run-specific annotation template populated from the Master Dictionary."""
    print("\n<<< Generating Run-Specific Annotations Template >>>")
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
        nodes = [n['data']['id'] for n in data.get('elements', {}).get('nodes', [])]

        if not nodes:
            print("  [!] No nodes found in network to annotate.")
            return

        master_dict = {}
        if os.path.exists(master_anno_path):
            try:
                master_df = pd.read_csv(master_anno_path, sep=';')
                master_df.columns = [c.lower() for c in master_df.columns]
                master_dict = dict(zip(master_df['mesh_term'], master_df['aop_level']))
            except Exception:
                pass

        run_data = []
        for n in nodes:
            # Pull from master if it exists and isn't unassigned, otherwise default to Unassigned
            lvl = master_dict.get(n, 'Unassigned')
            run_data.append({'mesh_term': n, 'aop_level': lvl})

        run_df = pd.DataFrame(run_data)
        run_df.to_csv(run_anno_path, sep=';', index=False)
        print(f"  [+] Created run-specific annotation file: {os.path.basename(run_anno_path)}")

    except Exception as e:
        print(f"  [!] Failed to generate run annotations: {e}")

def _sync_run_to_master(run_anno_path: str, master_anno_path: str, is_afk: bool):
    """Asks the user if they want to update the Master Dictionary with their new annotations."""
    if not os.path.exists(run_anno_path) or not os.path.exists(master_anno_path):
        return

    if is_afk:
        print("  [+] AFK Mode: Skipping Master Annotation Sync prompt to prevent pipeline pause.")
        return

    try:
        run_df = pd.read_csv(run_anno_path, sep=';')
        unassigned_count = (run_df['aop_level'].str.strip().str.lower() == 'unassigned').sum()
    except Exception:
        unassigned_count = 0

    print("\n<<< Master Annotation Sync >>>")
    print("  [!] WARNING: This action will permanently update your Master Annotations Library.")
    if unassigned_count > 0:
        print(f"  [!] CAUTION: Your run-specific file contains {unassigned_count} 'Unassigned' term(s).")
        print("      (Note: 'Unassigned' terms are safely ignored and will NOT overwrite existing Master assignments).")

    ans = input("\n  [?] Do you want to sync the AOP levels from your run-specific file back to the Master Library? [y/n/Enter to skip]: ").strip().lower()
    if ans in ['y', 'yes']:
        try:
            master_df = pd.read_csv(master_anno_path, sep=';')
            run_dict = dict(zip(run_df['mesh_term'], run_df['aop_level']))

            def update_lvl(row):
                term = row['mesh_term']
                if term in run_dict and run_dict[term].strip().lower() != 'unassigned':
                    return run_dict[term]
                return row['aop_level']

            master_df['aop_level'] = master_df.apply(update_lvl, axis=1)
            existing_terms = set(master_df['mesh_term'])
            new_rows = [r for _, r in run_df.iterrows() if r['mesh_term'] not in existing_terms]
            if new_rows:
                master_df = pd.concat([master_df, pd.DataFrame(new_rows)], ignore_index=True)

            master_df.to_csv(master_anno_path, sep=';', index=False)
            print("  [+] Master Annotations Library updated successfully.")
        except Exception as e:
            print(f"  [!] Failed to update Master Annotations Library: {e}")

def main():
    """Parse command-line arguments and run the requested pipeline stage(s)."""
    parser = argparse.ArgumentParser(
        prog="mesh-pipeline",
        description=(
            "Medical Subject Headings (MeSH) Adverse Outcome Pathway (AOP) Network Pipeline\n"
            "Orchestrates knowledge graph generation, semantic filtering, and visualization."
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-v', '--version', action='version',
                        version=f"%(prog)s v{get_version()}",
                        help="Print the current package version and exit.")

    parser.add_argument('--readme', action='store_true',
                        help="Open the comprehensive README.md documentation in the default system viewer.")

    parser.add_argument('--step', choices=['all', 'process', 'data_ops', 'network', 'secondary', 'viz', 'benchmark'], default='all',
                        help=(
                            "Specify which segment of the pipeline to execute:\n"
                            "  all        : Run the entire pipeline sequentially (Default).\n"
                            "  process    : Step 1 - Process MeSH XML and initialize Master Annotations.\n"
                            "  data_ops   : Step 2 - Collect Entrez data and build SQLite PMIDs database.\n"
                            "  network    : Step 3 - Construct network, filter via GLF/SA, and calculate CRS.\n"
                            "  secondary  : Step 3.5 - Execute targeted node/edge queries and metadata hydration.\n"
                            "  viz        : Step 4 - Generate biological plots, joint distributions, and HTML networks.\n"
                            "  benchmark  : Step 5 - Validate ground-truth citations and benchmark ranking performance."
                        ))

    parser.add_argument('--config', default='mesh_config.json',
                        help="Path to the configuration JSON file. Defaults to 'mesh_config.json'.")

    parser.add_argument('--interactive', action='store_true',
                        help="Launch the interactive command-line wizard to modify parameters before execution.")

    args = parser.parse_args()

    # Intercept Documentation Request
    if args.readme:
        open_readme()
        sys.exit(0)

    print("<"*30 + ">"*30)
    print("STARTING MESH NETWORK ANALYSIS PIPELINE")
    print("<"*30 + ">"*30)

    try:
        config = MeshConfig(config_path=args.config)
    except Exception as e:
        print(f"\n[!] Configuration Error: {e}")
        sys.exit(1)

    if args.interactive:
        reset_triggered = run_interactive_wizard(config, args.step)
        if reset_triggered:
            config.reset()
            print("\n<<< Pipeline halted. Factory defaults have been restored. Run again to execute. >>>\n")
            sys.exit(0)
        else:
            config.save()

    total_start = time.time()

    opt_hist_path = os.path.join(os.path.dirname(config.files['full_network']), f"{config.prefix}_optimization_history.json")
    run_anno_path = os.path.join(config.results_dir, f"{config.prefix}_run_annotations.csv")

    try:
        # <<< STEP 0: Offline Master DB Generation >>>
        if config.params.get('_run_baseline_etl', False):
            print("\n>>> STARTING: Step 0 - PubMed Baseline ETL (Master Database Generation)")

            if config.params.get('_delete_corrupt_db', False):
                db_path = str(config.files['master_db'])
                if os.path.exists(db_path):
                    # A corrupt master DB can't be reused and a multi-GB copy would
                    # only waste (OneDrive-synced) space, so delete it and rebuild.
                    os.remove(db_path)
                    print(f"  [!] Corrupt master DB deleted; rebuilding from scratch.")

            baseline_mgr = PubMedBaselineManager(
                raw_data_dir=config.active_raw_dir,
                master_db_path=config.files['master_db']
            )
            baseline_mgr.run_downloads(
                include_baseline=not config.params.get('_skip_baseline_download', False),
                include_updates=config.params.get('_run_baseline_updates', False)
            )
            baseline_mgr.compile_database()

        # <<< STEP 1: MeSH Support Files Processing >>>
        if args.step in ['all', 'process']:
            print("\n>>> STARTING: Step 1 - MeSH Raw Data Processing")

            xml_file_path = config.files.get('mesh_xml', os.path.join(config.active_raw_dir, 'desc2025.xml'))

            if config.get('search_parameters', 'update_mesh_support_files') or not os.path.exists(config.files['mesh_terms_csv']):
                _ensure_prerequisites({"MeSH XML": xml_file_path}, "Step 1 (Raw Data Processing)")

            process_raw_mesh_data(
                xml_file=xml_file_path,
                output_csv=config.files['mesh_terms_csv'],
                output_py=config.files['mesh_stopwords_py'],
                force_update=config.get('search_parameters', 'update_mesh_support_files')
            )

            # Pre-populate Master Annotations Library
            _initialize_master_annotations(config.files['mesh_terms_csv'], config.files['annotations'])

        # <<< STEP 2: Entrez Data Collection >>>
        if args.step in ['all', 'data_ops']:
            print("\n>>> STARTING: Step 2 - Entrez Data Collection & Database Operations")

            if not _is_valid_db(config.files['pmids_db'], "pmids_table"):
                run_initial_data_collection(
                    search_term_param=config.get('search_parameters', 'search_term'),
                    start_date_str=config.get('search_parameters', 'start_date'),
                    end_date_str=config.get('search_parameters', 'end_date'),
                    generations_n_param=config.get('search_parameters', 'generations_n'),
                    db_path=config.files['pmids_db'],
                    entrez_email=config.get('credentials', 'entrez_email'),
                    entrez_api_key=config.get('credentials', 'entrez_api_key')
                )
            if not _is_valid_db(config.files['cleaned_db'], "pmids_table"):
                clean_database(
                    original_db=config.files['pmids_db'],
                    cleaned_db=config.files['cleaned_db'],
                    start_gen_to_remove_label=f"G{config.get('search_parameters', 'generations_n')}"
                )

            populate_master_mesh_database(
                source_pmids_input=config.files['cleaned_db'],
                master_db_path=config.files['master_db'],
                entrez_email=config.get('credentials', 'entrez_email'),
                entrez_api_key=config.get('credentials', 'entrez_api_key'),
                failed_log=config.files['failed_mesh'],
                empty_log=config.files['empty_mesh']
            )

            update_db_with_mesh_batch(
                local_db_path=config.files['cleaned_db'],
                master_db_path=config.files['master_db']
            )

        # <<< STEP 3: Network Construction & Analysis >>>
        if args.step in ['all', 'network']:
            print("\n>>> STARTING: Step 3 - Network Construction & Statistical Filtering")

            _ensure_prerequisites({"Cleaned Database": config.files['cleaned_db']}, "Step 3 (Network Construction)")

            if not os.path.exists(config.files['full_network']):
                run_network_construction(
                db_path_param=config.files['cleaned_db'],
                output_json_path=config.files['full_network'],
                lambda_val=config.get('network_parameters', 'lambda_val'),
                node_weight_factors=config.params.get('network_parameters', {}).get('node_weight_factors', {'centrality': 0.45, 'article_rank': 0.15, 'rank_median_cit': 0.2, 'rank_total_cit': 0.2}),
                run_full_centrality=config.get('analysis_parameters', 'calculate_full_centrality'),
                betweenness_k_samples=config.get('analysis_parameters', 'betweenness_k_samples'),
                random_seed=config.get('analysis_parameters', 'random_seed'),
                eigenvector_max_iter=config.get('analysis_parameters', 'eigenvector_max_iter'),
                eigenvector_tol=config.get('analysis_parameters', 'eigenvector_tol')
            )

            if not os.path.exists(config.files['consensus_lcc']):
                run_consensus_filtering_and_lcc(
                    input_json_path=config.files['full_network'],
                    glf_output_path=config.files['glf_subgraph'],
                    sa_output_path=config.files['sa_subgraph'],
                    final_lcc_output_path=config.files['consensus_lcc'],
                    history_output_path=opt_hist_path,
                    target_num_edges=config.get('simulation_parameters', 'target_num_edges'),
                    glf_iterations=config.get('simulation_parameters', 'glf_iterations'),
                    sa_iterations=config.get('simulation_parameters', 'sa_iterations'),
                    sa_initial_temp=config.get('simulation_parameters', 'sa_temp_start'),
                    sa_cooling_rate=config.get('simulation_parameters', 'sa_cooling_rate'),
                    random_seed=config.get('analysis_parameters', 'random_seed') or 42
                )

            run_step_6 = True
            if os.path.exists(config.files['consensus_lcc']):
                try:
                    with open(config.files['consensus_lcc'], 'r') as f:
                        nodes = json.load(f).get('elements', {}).get('nodes', [])
                        if nodes and nodes[0].get('data', {}).get('filtered_louvain_community_id') is not None:
                            run_step_6 = False
                except Exception:
                    pass

            if run_step_6:
                run_community_detection(
                    network_file_path=config.files['consensus_lcc'],
                    random_seed=config.get('analysis_parameters', 'random_seed')
                )

            if not os.path.exists(config.files['final_network']):
                run_contextual_relevance_scoring(
                    input_nodes_file=config.files['consensus_lcc'],
                    output_nodes_file=config.files['final_network'],
                    master_db_path=config.files['master_db'],
                    relevance_db_path=config.files['relevance_db'],
                    id_key='id',
                    weight_key_1="betweenness_centrality", final_key_1="CRS_betweenness_centrality",
                    weight_key_2="eigenvector_centrality", final_key_2="CRS_eigenvector_centrality",
                    start_date_param=config.get('analysis_parameters', 'context_start_date'),
                    end_date_param=config.get('analysis_parameters', 'context_end_date'),
                    entrez_email=config.get('credentials', 'entrez_email'),
                    entrez_api_key=config.get('credentials', 'entrez_api_key'),
                    calculate_full_centrality=config.get('analysis_parameters', 'calculate_full_centrality')
                )

            excel_path = os.path.join(config.results_dir, f"{config.prefix}_export.xlsx")
            convert_network_json_to_excel(config.files['final_network'], excel_path)

            # Generate Run-Specific Annotations Template
            _generate_run_annotations(config.files['final_network'], config.files['annotations'], run_anno_path)

            # <<< AOP-INDEPENDENT VISUALIZATIONS >>>
            print("\n<<< Generating Topological & Distribution Figures (AOP-Independent) >>>")
            _, edge_df, _ = load_and_prepare_data(config.files['final_network'], run_anno_path)
            raw_edge_df, raw_all_edges = load_full_raw_data(config.files['full_network'])

            if raw_edge_df is not None:
                analyze_dispersion(raw_edge_df, str(config.figures_dir), config.prefix)
                plot_cooccurrance_distribution(raw_edge_df, edge_df, str(config.figures_dir), config.prefix)

            if raw_all_edges and os.path.exists(opt_hist_path):
                run_optimization_comparison(
                    history_json_path=opt_hist_path,
                    output_dir=str(config.figures_dir),
                    file_prefix=config.prefix
                )

        # <<< STEP 3.5: Secondary Analysis >>>
        sec_params = config.params.get('secondary_analysis', {})

        # Safely determine if we should run secondary analysis based on the execution mode
        should_run_sec = True if args.step == 'secondary' else sec_params.get('export_top_articles', True)

        if args.step in ['all', 'secondary'] and should_run_sec:
            print("\n>>> STARTING: Secondary Analysis Operations")
            _ensure_prerequisites({
                "Relevance Database": config.files['relevance_db'],
                "Cleaned Database": config.files['cleaned_db']
            }, "Secondary Analysis")

            api_email = config.get('credentials', 'entrez_email')
            api_key = config.get('credentials', 'entrez_api_key')
            exclude_rev = sec_params.get('exclude_reviews', True)

            # --- DUAL-ENGINE METRIC PARAMETERS ---
            export_limit = sec_params.get('export_limit', 500)
            sort_metric = sec_params.get('sort_metric', 'Linear')
            linear_weight_ars = sec_params.get('linear_weight_ars', 0.5)

            if sec_params.get('export_top_articles', True):
                get_top_network_articles(
                    db_path=config.files['relevance_db'],
                    cleaned_db_path=config.files['cleaned_db'],
                    results_dir=str(config.results_dir),
                    file_prefix=config.prefix, limit=export_limit, exclude_reviews=exclude_rev,
                    entrez_email=api_email, entrez_api_key=api_key,
                    sort_metric=sort_metric, linear_weight_ars=linear_weight_ars
                )
            if args.step == 'secondary':
                t_nodes = sec_params.get('target_nodes', '')
                if t_nodes:
                    for node_name in [n.strip() for n in t_nodes.split(';') if n.strip()]:
                        analyze_node_relevancy(
                            node_name=node_name, db_path=config.files['relevance_db'],
                            cleaned_db_path=config.files['cleaned_db'],
                            results_dir=str(config.results_dir),
                            file_prefix=config.prefix, limit=export_limit, exclude_reviews=exclude_rev,
                            entrez_email=api_email, entrez_api_key=api_key,
                            sort_metric=sort_metric, linear_weight_ars=linear_weight_ars
                        )

                t_edges = sec_params.get('target_edges', '')
                if t_edges:
                    for edge_str in [e.strip() for e in t_edges.split(';') if e.strip()]:
                        parts = [p.strip() for p in edge_str.split(' - ')]
                        if len(parts) == 2:
                            analyze_edge_relevancy(
                                node1=parts[0], node2=parts[1], db_path=config.files['relevance_db'],
                                cleaned_db_path=config.files['cleaned_db'],
                                results_dir=str(config.results_dir),
                                file_prefix=config.prefix, limit=export_limit, exclude_reviews=exclude_rev,
                                entrez_email=api_email, entrez_api_key=api_key,
                                sort_metric=sort_metric, linear_weight_ars=linear_weight_ars
                            )

        if args.step == 'all':
            if config.params.get('control_flags', {}).get('pause_for_annotation', False):
                print(f"\n[!] Pipeline paused. AOP level pipeline functions not completed.")
                print(f"    Annotations for MeSH terms in LCC need to checked and edited manually for new runs")
                print(f"    Annotate {os.path.basename(run_anno_path)} and then run mesh-pipeline --step viz")
                sys.exit(0)
            else:
                print(f"AFK Override Active. Proceeding to visualization with 'Unassigned' AOP levels...")

        # <<< STEP 4: Visualization (AOP-Dependent) >>>
        if args.step in ['all', 'viz']:
            print("\n>>> STARTING: Step 4 - Biological Figure Generation")

            _ensure_prerequisites({
                "Final Network JSON": config.files['final_network'],
                "Run-Specific Annotations CSV": run_anno_path
            }, "Step 4 (Biological Figure Generation)")

            is_afk = not config.params.get('control_flags', {}).get('pause_for_annotation', False)
            _sync_run_to_master(run_anno_path, config.files['annotations'], is_afk=is_afk)

            node_df, edge_df, G = load_and_prepare_data(config.files['final_network'], run_anno_path)

            plot_louvain_community_bars(node_df, str(config.figures_dir), config.prefix)
            plot_joint_plot(node_df, str(config.figures_dir), config.prefix)
            plot_tsne_louvain_overlap(node_df, G, str(config.figures_dir), config.prefix)
            plot_sankey_alluvial(G, node_df, str(config.figures_dir), config.prefix)
            plot_dumbell_plot(node_df, str(config.figures_dir), config.prefix)
            plot_scatter_panels(node_df, str(config.figures_dir), config.prefix)
            plot_dendrogram(G, node_df, str(config.figures_dir), config.prefix)

        # <<< STEP 5: Ground-Truth Validation & Benchmarking >>>
        if args.step == 'benchmark':
            print("\n>>> STARTING: Step 5 - Ground-Truth Validation & Benchmarking")

            bench_params = config.params.get('benchmark', {})
            configured_gt = bench_params.get('ground_truth_csv', '')

            # Resolve the ground-truth file from the ACTIVE raw directory. This is
            # data/reference_raw/ when 'Use Reference Data' is True, and data/raw/
            # when False -- so users supply their own PMID list by dropping a file
            # named e.g. 'ground_truth_pmids.csv' into data/raw/.
            gt_path = resolve_ground_truth_path(config.active_raw_dir, configured_gt)
            if not gt_path:
                print("\n" + "<"*30 + ">"*30)
                print("[CRITICAL ERROR] No ground-truth file found for Benchmarking")
                print("<"*30 + ">"*30)
                print(f"  Searched in: {config.active_raw_dir}")
                if configured_gt:
                    print(f"  Configured filename '{configured_gt}' was not present there.")
                else:
                    print("  Place a ground-truth file there using one of these names:")
                    for name in KNOWN_GT_FILENAMES:
                        print(f"    - {name}")
                print("\n  Accepted formats: a CSV/TSV with a 'PMID' column (optionally a")
                print("  reference/citation column for date validation), a single column of")
                print("  PMIDs, or a plain .txt list with one PMID per line.")
                sys.exit(1)

            nc_filename = bench_params.get('negative_control_csv', '')
            nc_path = resolve_ground_truth_path(config.active_raw_dir, nc_filename) if nc_filename else None

            _ensure_prerequisites({
                "Relevance Database": config.files['relevance_db']
            }, "Step 5 (Benchmarking)")

            print(f"  [+] Using ground-truth file: {os.path.basename(gt_path)}")

            try:
                run_benchmark(
                    resolved_csv_path=gt_path,
                    relevance_db_path=str(config.files['relevance_db']),
                    output_dir=str(config.results_dir),
                    file_prefix=f"{config.prefix}_benchmark",
                    primary_node=bench_params.get('primary_node', 'Dermatitis, Allergic Contact'),
                    negative_control_csv=nc_path,
                    random_seed=config.get('analysis_parameters', 'random_seed') or 42,
                    n_boot=bench_params.get('n_boot', 2000),
                    n_perm=bench_params.get('n_perm', 2000)
                )
            except (ValueError, KeyError) as e:
                print("\n" + "<"*30 + ">"*30)
                print("[BENCHMARK ERROR] The ground-truth file could not be parsed")
                print("<"*30 + ">"*30)
                print(f"  File: {gt_path}")
                print(f"  Reason: {e}")
                print("\n  Accepted formats: a CSV/TSV with a 'PMID' column, a single")
                print("  column of PMIDs, or a plain .txt list with one PMID per line.")
                sys.exit(1)

    except Exception as e:
        print("\n" + "<"*30 + ">"*30)
        print(f"[CRITICAL ERROR] Pipeline Execution Failed: {e}")
        print("<"*30 + ">"*30)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\n\n[!] PIPELINE CANCELLED BY USER (Ctrl+C).")
        sys.exit(130)

    total_time = time.time() - total_start
    print("\n" + "<"*30 + ">"*30)
    print(f"PIPELINE COMPLETED SUCCESSFULLY in {total_time/60:.2f} minutes.")
    print("<"*30 + ">"*30)

if __name__ == "__main__":
    main()