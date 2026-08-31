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
from pathlib import Path
from importlib.metadata import version, PackageNotFoundError

from .config_parser import MeshConfig
from .wizard import run_interactive_wizard
from .mesh_data_processor import process_raw_mesh_data, ensure_mesh_descriptor_xml
from .baseline_manager import PubMedBaselineManager
from .data_ops import (
    run_initial_data_collection,
    clean_database,
    populate_master_mesh_database,
    update_db_with_mesh_batch
)
from .network import run_network_construction, run_consensus_filtering_and_lcc, run_community_detection
from .run_ledger import open_ledger
from .prisma import write_prisma_report
from .guides import AOP_LEVELS, write_annotation_guide, write_master_db_guide
from . import guides
from . import ledger_collect
from . import runcontrol
from .runcontrol import RunAborted
from .relevance import run_mean_relevancy_scoring
from .secondary_analysis import convert_network_json_to_excel, analyze_node_relevancy, analyze_edge_relevancy, get_top_network_articles, run_network_overlay_comparison
from .benchmark import run_benchmark, resolve_ground_truth_path, KNOWN_GT_FILENAMES
from .gt_network_validation import run_gt_network_validation
from .validation_report import run_validation_report, run_projection_comparison
from .viz import (
    configure_output,
    load_and_prepare_data, load_full_raw_data, analyze_dispersion,
    plot_cooccurrance_distribution, run_optimization_comparison,
    plot_louvain_community_bars, plot_tsne_louvain_overlap,
    plot_sankey_alluvial, plot_dendrogram, plot_network_graph
)

try:
    from .mesh_stop_words import MESH_STOP_WORDS
except ImportError:
    MESH_STOP_WORDS = set()

def get_version() -> str:
    """The version of the code that is actually running.

    Not the installed distribution's: the packaged build is never pip-installed,
    so metadata there is either absent or - worse - a stale dist-info left in a
    development environment, which reported 3.1.0 as 1.0.0. The module constant
    ships with the code it describes and cannot drift from it. Distribution
    metadata is consulted only if that constant is somehow missing.
    """
    from . import __version__ as pkg_version
    if pkg_version:
        return pkg_version
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

# The figures the user can turn off, in the order they are drawn. Absent from
# the config means on: an existing settings file predates these switches, and
# silently drawing nothing would be the worst possible reading of that.
FIGURE_SWITCHES = [
    ('distribution', 'Figure 1 - Edge weight distribution'),
    ('optimisation', 'Figure 2 - Optimisation trajectory'),
    ('communities',  'Figure 3 - Community composition'),
    ('tsne',         'Figure 4 - t-SNE projection'),
    ('alluvial',     'Figure 5 - AOP alluvial flow'),
    ('dendrogram',   'Figure 6 - Node2Vec dendrogram'),
    ('network',      'Figure 7 - Network graph'),
]


def _figure_on(config, name):
    """Whether this figure is ticked on the Figures tab."""
    value = config.params.get('viz_parameters', {}).get(f'fig_{name}')
    return True if value is None else bool(value)


# Figures 1 and 2 describe the filtering and the optimisation search, so they are
# drawn while that is happening - by the network step, not the figures step.
NETWORK_STEP_FIGURES = ('distribution', 'optimisation')


_FIGURE_OUTPUT_ANNOUNCED = False


def _announce_figure_output(config):
    """Say what figures are being written at, the first time one is written.

    Once per run: Step 3 and Step 4 both draw, and repeating it would be noise.
    """
    global _FIGURE_OUTPUT_ANNOUNCED
    if _FIGURE_OUTPUT_ANNOUNCED:
        return
    _FIGURE_OUTPUT_ANNOUNCED = True
    vp = config.params.get('viz_parameters', {})
    dpi, fmts = configure_output(vp.get('figure_dpi'), vp.get('figure_formats'))
    source = ('your settings' if (vp.get('figure_dpi') or vp.get('figure_formats'))
              else 'built-in defaults - nothing set on the Figures tab')
    print(f"\n    Figures are being written at {dpi} dpi as "
          f"{', '.join(fmts)}  ({source}).")


def _report_skipped_figures(config, step=''):
    """Say what was left out, rather than let it look like a failure.

    A figure that simply does not appear is indistinguishable from one that
    crashed quietly, and the run is long enough that nobody wants to guess.
    """
    off = [label for name, label in FIGURE_SWITCHES if not _figure_on(config, name)]
    if off:
        print(f"\n    {len(off)} figure(s) turned off on the Figures tab, not drawn:")
        for label in off:
            print(f"      - {label}")

    # Running the figures step alone draws four of the six. Without saying so,
    # the two that are missing look like the two that went wrong.
    if step == 'viz':
        elsewhere = [label for name, label in FIGURE_SWITCHES
                     if name in NETWORK_STEP_FIGURES and _figure_on(config, name)]
        if elsewhere:
            print("\n    Drawn by the network step rather than this one, because "
                  "they describe\n    the filtering as it happens:")
            for label in elsewhere:
                print(f"      - {label}")


def _gt_wanted(config):
    """Is the ground-truth analysis in play for this run?

    Reference mode implies it. The intent was always "absent from the config
    means follow reference data", but the settings form writes every key it
    knows about on save, so the key is never absent once the window has been
    used once - and an explicit false then defeated reference mode silently.
    That mattered beyond the benchmark itself: this same answer decides whether
    the subgraph centralities are computed at the network stage, so a reference
    run was quietly producing a network missing six node attributes that the
    published one has.
    """
    if bool(config.get('control_flags', 'use_reference_data')):
        return True
    return bool(config.params.get('benchmark', {}).get('run_ground_truth_analysis'))


# What the bundled reference corpus cannot do, and why. Each of these needs the
# cleaned citation database, which only retrieval produces - and retrieval does
# not run against a corpus that is already built. Skipping them is not a
# limitation of the mode; it IS the mode. The corpus ships with the work already
# done, which is why figures and benchmarking take minutes rather than a night.
REFERENCE_SKIPS = {
    'data_ops': ('Retrieval', 'that corpus is already built - there is nothing to fetch'),
    'network': ('Network construction',
                'that corpus ships with its network already built, filtered and scored'),
    'secondary': ('Secondary analysis',
                  'it reads the cleaned citation database, which only retrieval produces'),
}


def _reference_skip(config, step, requested):
    """Should this step be skipped because the reference corpus is in use?

    Asking for one of them directly is refused outright - the user asked for a
    specific thing that cannot happen. Inside a full run it is skipped and the
    run continues, because the steps that CAN work against the shipped network
    still have real work to do, and killing the whole run over a step that was
    never needed would be worse.
    """
    if not config.use_reference_data or step not in REFERENCE_SKIPS:
        return False
    label, why = REFERENCE_SKIPS[step]
    if requested == step:
        print("\n" + "<" * 30 + ">" * 30)
        print(f"{label} does not run against the bundled reference corpus")
        print("<" * 30 + ">" * 30)
        print(f"  'Use bundled reference data' is ticked, and {why}.")
        print("\n  To do this with your own corpus, untick it on the Search tab.")
        print("  Your own search term, dates and prefix come back as you left them.")
        print("\n  With the reference corpus, run the figures or the benchmark:")
        print("  everything they need is already here.")
        sys.exit(1)
    print(f"\n>>> SKIPPING: {label} - {why}.")
    return True


def _network_term_index(network_path):
    """Every term in the finished network, for checking a typed name against.

    Returns {} when the network cannot be read - the caller then says nothing,
    because a check that cannot run must not produce a warning either way.
    """
    try:
        with open(network_path, 'r', encoding='utf-8') as fh:
            nodes = json.load(fh).get('elements', {}).get('nodes', [])
        return {str(n.get('data', {}).get('id', '')): True for n in nodes if n.get('data')}
    except (OSError, ValueError, AttributeError):
        return {}


def _warn_unknown_term(name, known, what):
    """Say plainly that a typed term is not in the network, and why that matters.

    Nothing fails on a name that does not match - the query returns an empty set
    and the run carries on - so the only symptom is an analysis that produced no
    output, with no indication why. One clear sentence about the cause is worth
    more than a guess at which of the several possible mistakes it was.
    """
    if not known or name in known:
        return
    print(f"\n  [!] {what}: no term named {name!r} in this network.")
    print(f"      Terms must match a network term exactly, including capitals.")
    print(f"      Nothing will be produced for it; the rest of the step runs.")
    print(f"      The full term list is the Nodes tab of "
          f"<prefix>_final_network_nodes_and_edges.xlsx.")


def _build_relevance_db_if_missing(config, include_subgraph_weightings=False):
    """Score the corpus into a relevance database when there isn't one yet.

    The benchmark needs per-article relevance scores. Against the bundled
    reference corpus there were none and there could not be: the database is two
    gigabytes, far too large to ship, so the benchmark simply refused to run and
    the shipped corpus could not be validated by anyone.

    It does not need to be shipped. Everything required to build it is already
    on the machine - the consensus network, which is bundled, and the master
    annotation database, which the user built. So build it.

    This is not free: it scans the whole master database, which is the same pass
    the network step makes and takes a while. It is done once and kept.
    """
    target = config.files['relevance_db']
    if os.path.exists(target) and os.path.getsize(target) > 0:
        return

    network = config.files['consensus_lcc']
    master = config.files['master_db']
    if not os.path.exists(network) or not os.path.exists(master):
        return          # _ensure_prerequisites will say what is missing, and why

    print("\n" + "<" * 30 + ">" * 30)
    print("<<< No relevance database yet - building it >>>")
    print("<" * 30 + ">" * 30)
    print(f"  network : {network}")
    print(f"  master  : {master}")
    print(f"  writing : {target}")
    print("  This scans the master annotation database once. It takes a while,")
    print("  and it is kept - the next benchmark starts straight away.\n")

    # The rescored network goes beside the databases rather than over the input.
    # In reference mode the input IS the shipped corpus, which is read-only and
    # must stay exactly as published.
    scored = Path(config.databases_dir) / f"{config.db_prefix}_scored_for_benchmark.json"

    weightings = [("betweenness_centrality", "MRS_betweenness_centrality"),
                  ("pagerank_centrality", "MRS_pagerank_centrality"),
                  ("eigenvector_centrality", "MRS_eigenvector_centrality")]
    if include_subgraph_weightings:
        weightings += [("betweenness_subgraph_centrality", "MRS_betweenness_subgraph_centrality"),
                       ("pagerank_subgraph_centrality", "MRS_pagerank_subgraph_centrality"),
                       ("eigenvector_subgraph_centrality", "MRS_eigenvector_subgraph_centrality")]

    run_mean_relevancy_scoring(
        input_nodes_file=network,
        output_nodes_file=str(scored),
        master_db_path=master,
        relevance_db_path=target,
        id_key='id',
        weightings=weightings,
        start_date_param=config.get('analysis_parameters', 'context_start_date'),
        end_date_param=config.get('analysis_parameters', 'context_end_date'),
        entrez_email=config.get('credentials', 'entrez_email'),
        entrez_api_key=config.get('credentials', 'entrez_api_key'),
        calculate_full_centrality=config.get('analysis_parameters', 'calculate_full_centrality'))

    print(f"\n  [+] Relevance database built: {target}")


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

def _run_file_check(config, repair=False, deep=False):
    """Report every damaged artefact, optionally delete them. Returns an exit code.

    Deleting is offered because these files are all reproducible: the pipeline
    rebuilds whatever is absent. Keeping a corrupt one is strictly worse than
    having none, because a present-but-wrong file is what the next step will
    happily read.
    """
    from . import integrity

    print("\n" + "<" * 30 + ">" * 30)
    print("CHECKING FILES" + (" (full integrity check)" if deep else ""))
    print("<" * 30 + ">" * 30)
    print(f"  Data     : {config.active_raw_dir}")
    print(f"  Processed: {config.active_source_dir}")
    print(f"  Results  : {config.results_dir}\n")

    artifacts = integrity.scan(config, deep=deep)
    for a in artifacts:
        mark = {integrity.OK: '[+]', integrity.MISSING: '[ ]'}.get(a.status, '[!]')
        note = 'not built yet' if a.status == integrity.MISSING else a.detail
        print(f"  {mark} {a.label:44} {note}")

    print()
    print(integrity.summarise(artifacts))

    suspect = [a for a in integrity.problems(artifacts) if a.status == integrity.SUSPECT]
    if suspect:
        print("\n  These open but do not look right. They are left alone: the")
        print("  pipeline rebuilds what it cannot use, and deleting something")
        print("  merely odd is not a call to make automatically.")
        for a in suspect:
            print(f"      {a.label} - {a.detail}")
        print("  Tools -> Check and repair files in the Workbench can remove them.")

    broken = [a for a in integrity.problems(artifacts) if a.broken]
    if not broken:
        return 0

    if not repair:
        print("\n  Add --repair-files to delete the damaged files listed above.")
        print("  Nothing has been changed.")
        return 1

    costly = [a for a in broken if a.cost == integrity.COST_HOURS]
    if costly:
        print("\n  [!] This will delete files that take HOURS to rebuild:")
        for a in costly:
            print(f"        {a.label} ({a.size_mb:,.0f} MB)")

    removed, freed, failures = integrity.remove(broken)
    print(f"\n  Removed {len(removed)} file(s), {freed/1e6:,.0f} MB freed.")
    for a, err in failures:
        print(f"  [!] Could not remove {a.path}: {err}")
    step = integrity.resume_step(artifacts)
    if step:
        print(f"\n  Now resume from: {integrity.STEP_LABEL.get(step, step)}")
        print(f"    python -m mesh_aop.cli --step {step}"
              + (" --build-database" if step == 'baseline' else ""))
    return 0 if not failures else 1


def _preflight(config, step):
    """Check the inputs a step depends on before it starts, not four hours in.

    Only the master database is examined closely, because it is the one artefact
    whose replacement costs a night rather than a coffee break, and because the
    ways it goes wrong - a sync client serving a placeholder, a build cut short
    by a full disk - all produce a file that looks entirely normal until
    something tries to read from the middle of it.

    Never fatal on its own: a warning here and a clear failure later beats
    refusing to start over a check that might itself be wrong.
    """
    from . import integrity

    needs_master = step in ('all', 'network', 'secondary', 'benchmark')
    if not needs_master or not os.path.exists(config.files['master_db']):
        return True

    print("\n<<< Checking the master annotation database >>>")
    started = time.time()
    status, detail, rows = integrity.check_master_db(config.files['master_db'])
    took = time.time() - started

    if status == integrity.OK:
        print(f"  [+] Master database OK - {detail} ({took:.1f}s)")
        return True

    print("\n" + "<" * 30 + ">" * 30)
    print(f"[!] THE MASTER ANNOTATION DATABASE IS {status.upper()}")
    print("<" * 30 + ">" * 30)
    print(f"    {config.files['master_db']}")
    print(f"    {detail}")
    if status == integrity.SUSPECT:
        print("\n    The file opens, so the run will continue - but the results")
        print("    may be built on incomplete annotations. Worth checking before")
        print("    you rely on them.")
        return True
    print("\n    Scoring reads this database for every article, so continuing")
    print("    would fail later and waste the intervening hours.")
    print("\n    To fix it:")
    print("      Workbench:  Tools -> Check and repair files")
    print("      Terminal :  python -m mesh_aop.cli --step baseline "
          "--build-database --rebuild-corrupt")
    print(f"\n    Instructions are also saved beside the database itself, as")
    print(f"      {guides.MASTER_DB_GUIDE_NAME}")
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
                """Prefer this run's AOP level for a term when it was actually assigned (not 'Unassigned')."""
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

    parser.add_argument('--step', choices=['all', 'baseline', 'process', 'data_ops', 'network', 'secondary', 'viz', 'benchmark'], default='all',
                        help=(
                            "Specify which segment of the pipeline to execute:\n"
                            "  all        : Run the entire pipeline sequentially (Default).\n"
                            "  baseline   : Step 0 only - build the master annotation database, then stop.\n"
                            "               Matches no later stage, so nothing else runs. Pair with --build-database.\n"
                            "  process    : Step 1 - Process MeSH XML and initialize Master Annotations.\n"
                            "  data_ops   : Step 2 - Collect Entrez data and build SQLite PMIDs database.\n"
                            "  network    : Step 3 - Construct network, filter via GLF/SA, and calculate MRS.\n"
                            "  secondary  : Step 3.5 - Execute targeted node/edge queries and metadata hydration.\n"
                            "  viz        : Step 4 - Generate biological plots, joint distributions, and HTML networks.\n"
                            "  benchmark  : Step 5 - Validate ground-truth citations and benchmark ranking performance."
                        ))

    # Default None, NOT 'mesh_config.json'. A bare relative name is resolved
    # against the working directory, so the pipeline opened a different file
    # from the one the Workbench had just written: the window saves to the
    # per-user settings path (paths.config_path()), while this opened
    # ./mesh_config.json wherever the process happened to be started. On any
    # installed copy that file does not exist, so every run silently used
    # factory defaults - an empty search term, no credentials, the default data
    # folder - and reported the search term as unset while the user was looking
    # at it in the form. Passing None lets MeshConfig resolve the same path the
    # Workbench does. An explicit --config still overrides it.
    parser.add_argument('--config', default=None,
                        help="Path to the configuration JSON file. Defaults to this copy's settings file.")

    parser.add_argument('--interactive', action='store_true',
                        help="Launch the interactive command-line wizard to modify parameters before execution.")

    # Forces Step 1 to fetch the descriptor XML and rebuild the stop-word
    # vocabulary from it. For one run only, like the Step 0 switches: the
    # equivalent config setting would refetch on every run afterwards, which is
    # a download nobody asked for repeated indefinitely. This is what the
    # Generate button on the Database screen sends.
    parser.add_argument('--refresh-mesh-support', action='store_true',
                        help="Re-download the MeSH descriptor file and rebuild the stop-word vocabulary from it.")

    # Memory, not speed, is what this controls. Each parser process is a Python
    # interpreter with an XML tree in it, so the database build's footprint is
    # roughly this number times the per-worker peak - which the build now prints
    # as it goes. Left unset it is derived from the RAM the machine reports.
    parser.add_argument('--max-workers', type=int, default=None, metavar='N',
                        help="Parser processes for the database build. Default: chosen from available RAM.")

    # Step 0 switches. These stay on the command line rather than in the config
    # file because they describe one run: a rebuild that persisted would repeat
    # itself, at several hours a time, on every run afterwards. Without the
    # wizard there was no way to reach Step 0 at all.
    parser.add_argument('--build-database', action='store_true',
                        help="Run Step 0 first: download the PubMed baseline and compile the master annotation database.")
    parser.add_argument('--skip-baseline-download', action='store_true',
                        help="With --build-database, compile from archives already on disk instead of downloading them again.")
    parser.add_argument('--with-updates', action='store_true',
                        help="With --build-database, also fetch the daily update files published since the baseline.")
    parser.add_argument('--rebuild-corrupt', action='store_true',
                        help="With --build-database, delete an unreadable master database before rebuilding it.")

    parser.add_argument('--check-files', action='store_true',
                        help="Check every file this project depends on and report anything damaged, then exit.")
    parser.add_argument('--repair-files', action='store_true',
                        help="With --check-files, delete the damaged files so the next run rebuilds them.")
    parser.add_argument('--deep-check', action='store_true',
                        help="With --check-files, run SQLite's full integrity check. Thorough, and slow on a large database.")

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

    if args.check_files:
        sys.exit(_run_file_check(config, repair=args.repair_files, deep=args.deep_check))

    if args.interactive:
        reset_triggered = run_interactive_wizard(config, args.step)
        if reset_triggered:
            config.reset()
            print("\n<<< Pipeline halted. Factory defaults have been restored. Run again to execute. >>>\n")
            sys.exit(0)
        else:
            config.save()

    # Applied after the wizard so an explicit flag wins, and after save() so the
    # switches never reach the config file.
    if args.build_database:
        config.params['_run_baseline_etl'] = True
        config.params['_skip_baseline_download'] = args.skip_baseline_download
        config.params['_run_baseline_updates'] = args.with_updates
        config.params['_delete_corrupt_db'] = args.rebuild_corrupt

    # A project needs a name before it writes anything. The prefix ships empty
    # so a fresh install does not inherit this project's name, which means it
    # has to be checked - without one, every file would be called _pmids.db,
    # _final_network_with_relevance.json and so on, and two projects would
    # collide on the first run.
    if not str(config.prefix or '').strip():
        print("\n" + "<" * 30 + ">" * 30)
        print("[CRITICAL ERROR] No project prefix is set")
        print("<" * 30 + ">" * 30)
        print("  Every file a run produces is named with it, so it is what keeps")
        print("  one analysis apart from another.")
        print("\n  Set 'Project prefix' on the Search tab - a short name for the")
        print("  question you are asking, such as SkinSens or ACD_2026. Letters,")
        print("  digits and underscores.")
        print("\n  Or tick 'Use bundled reference data', which supplies its own.")
        sys.exit(1)

    # What the reference corpus replaced, said out loud. A silent override is
    # the right behaviour and a silent override nobody is told about is not.
    _overrides = getattr(config, 'reference_overrides', [])
    if _overrides:
        print("\n" + "<" * 30 + ">" * 30)
        print("<<< Using the bundled reference corpus >>>")
        print("<" * 30 + ">" * 30)
        print("  These settings describe the published run and are used instead")
        print("  of yours. Your saved settings are NOT changed - untick")
        print("  'Use bundled reference data' and they come straight back.\n")
        for key, was, now in _overrides:
            shown = repr(was) if was not in (None, '') else '(empty)'
            print(f"    {key:<40} {shown}  ->  {now!r}")
        print("\n  Everything else is yours: figure resolution and formats, which")
        print("  figures are drawn, the folders, and your NCBI credentials.")

    total_start = time.time()

    # Counts handed back by the stages that actually execute. A resumed run
    # leaves these empty and the ledger falls back to reading the artefacts.
    build_stats, filter_stats = {}, {}

    # Resolution and file types, set before anything is drawn. This has to
    # happen here rather than in Step 4: Step 3 draws the topological and
    # distribution figures too, and the PRISMA flow is drawn after every step.
    _vp = config.params.get('viz_parameters', {})
    _dpi, _fmts = configure_output(_vp.get('figure_dpi'), _vp.get('figure_formats'))
    # Announced only by the steps that actually draw something. It used to print
    # at startup whatever the run was doing, so building a database or
    # collecting data began by reporting figure settings neither would ever use
    # - and when the settings were absent it reported the module defaults,
    # 300 dpi and jpeg/tif, which is why the line looked hardcoded.
    # Announced by whichever step first draws a figure, not here. On a full run
    # this was the very first line printed - figure settings reported before the
    # pipeline had done anything that could produce a figure. See
    # _announce_figure_output.

    # The rebuild instructions live beside the database, refreshed each run so
    # they always name the current paths. Cheap, and it means the data folder is
    # never just an unexplained pile of gigabytes.
    try:
        write_master_db_guide(config.files['master_db'],
                              workspace_dir=str(config.params.get('directories', {})
                                                .get('etl_workspace_dir', '') or ''))
    except Exception:
        pass

    # Check the expensive input before spending hours on work that depends on it.
    if not _preflight(config, args.step):
        sys.exit(1)

    opt_hist_path = os.path.join(os.path.dirname(config.files['full_network']), f"{config.prefix}_optimization_history.json")
    run_anno_path = os.path.join(config.results_dir, f"{config.prefix}_run_annotations.csv")

    try:
        # <<< STEP 0: Offline Master DB Generation >>>
        if config.params.get('_run_baseline_etl', False):
            runcontrol.checkpoint("before Step 0")
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
                master_db_path=config.files['master_db'],
                workspace_dir=(config.params.get('directories', {})
                               .get('etl_workspace_dir', '').strip() or None),
                max_workers=args.max_workers
            )
            baseline_mgr.run_downloads(
                include_baseline=not config.params.get('_skip_baseline_download', False),
                include_updates=config.params.get('_run_baseline_updates', False)
            )
            baseline_mgr.compile_database()

        # <<< STEP 1: MeSH Support Files Processing >>>
        if args.step in ['all', 'process']:
            runcontrol.checkpoint("before Step 1")
            print("\n>>> STARTING: Step 1 - MeSH Raw Data Processing")

            # The descriptor XML is only read when the support files are being
            # (re)built. When that is the case, guarantee a current desc<year>.xml
            # is present - auto-downloading it from NLM if missing or superseded -
            # so users never have to fetch it by hand. Otherwise extraction is
            # skipped and the path is a harmless placeholder.
            forced = bool(args.refresh_mesh_support)
            needs_support_build = (
                forced
                or config.get('search_parameters', 'update_mesh_support_files')
                or not os.path.exists(config.files['mesh_terms_csv'])
                or not os.path.exists(config.files['mesh_stopwords_py'])
            )
            if needs_support_build:
                xml_file_path = ensure_mesh_descriptor_xml(str(config.active_raw_dir))
            else:
                xml_file_path = os.path.join(config.active_raw_dir, 'desc2025.xml')

            try:
                process_raw_mesh_data(
                    xml_file=xml_file_path,
                    output_csv=config.files['mesh_terms_csv'],
                    output_py=config.files['mesh_stopwords_py'],
                    force_update=forced or config.get('search_parameters', 'update_mesh_support_files')
                )
            except (PermissionError, OSError) as exc:
                # The stop-word list lives inside the package, and an installed
                # copy under Program Files or /opt is read-only. That is not a
                # reason to fail the run: the shipped list is present and in
                # use, and it only goes stale when the MeSH release year moves.
                if not os.path.exists(config.files['mesh_stopwords_py']):
                    raise
                print(f"\n  [!] The MeSH support files could not be rewritten: {exc}")
                print( "      This copy is installed read-only, so the shipped "
                       "stop-word list is used as it is.")
                print( "      That is fine unless you are moving to a newer MeSH "
                       "release; for that, use a portable copy.")

            # Pre-populate Master Annotations Library
            _initialize_master_annotations(config.files['mesh_terms_csv'], config.files['annotations'])

        # <<< STEP 2: Entrez Data Collection >>>
        _skip_retrieval = (args.step in ['all', 'data_ops']
                          and _reference_skip(config, 'data_ops', args.step))
        if args.step in ['all', 'data_ops'] and not _skip_retrieval:
            runcontrol.checkpoint("before Step 2")
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
        _skip_network = (args.step in ['all', 'network']
                        and _reference_skip(config, 'network', args.step))
        if args.step in ['all', 'network'] and not _skip_network:
            runcontrol.checkpoint("before Step 3")
            print("\n>>> STARTING: Step 3 - Network Construction & Statistical Filtering")

            _ensure_prerequisites({"Cleaned Database": config.files['cleaned_db']}, "Step 3 (Network Construction)")

            if not os.path.exists(config.files['full_network']):
                build_stats = run_network_construction(
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
                filter_stats = run_consensus_filtering_and_lcc(
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

            # The subgraph PageRank exists purely to be scored against the
            # whole-corpus one, so it is driven by the ground-truth analysis switch
            # rather than a toggle of its own. It has to be decided here, at the
            # network stage, because it changes which centralities get computed.
            _bench = config.params.get('benchmark', {})
            include_subgraph_weightings = _gt_wanted(config)

            run_step_6 = True
            if os.path.exists(config.files['consensus_lcc']):
                try:
                    with open(config.files['consensus_lcc'], 'r') as f:
                        nodes = json.load(f).get('elements', {}).get('nodes', [])
                        if nodes:
                            data0 = nodes[0].get('data', {})
                            has_communities = data0.get('filtered_louvain_community_id') is not None
                            # Re-run when subgraph centralities are wanted but absent,
                            # otherwise a cached network would silently skip them.
                            needs_subgraph = include_subgraph_weightings and any(
                                data0.get(k) is None for k in (
                                    'pagerank_subgraph_centrality',
                                    'betweenness_subgraph_centrality',
                                    'eigenvector_subgraph_centrality',
                                )
                            )
                            if has_communities and not needs_subgraph:
                                run_step_6 = False
                except Exception:
                    pass

            if run_step_6:
                run_community_detection(
                    network_file_path=config.files['consensus_lcc'],
                    random_seed=config.get('analysis_parameters', 'random_seed'),
                    compute_subgraph_centrality=include_subgraph_weightings
                )

            if not os.path.exists(config.files['final_network']):
                run_mean_relevancy_scoring(
                    input_nodes_file=config.files['consensus_lcc'],
                    output_nodes_file=config.files['final_network'],
                    master_db_path=config.files['master_db'],
                    relevance_db_path=config.files['relevance_db'],
                    id_key='id',
                    # Whole-graph centralities always; the subgraph set as well when
                    # enabled, so centrality SCOPE (whole corpus vs consensus subgraph)
                    # and TYPE (betweenness / PageRank / eigenvector) form a full 3x2
                    # grid rather than varying scope for only some algorithms. The
                    # corpus scan dominates runtime, so extra weightings are near-free.
                    weightings=(
                        [("betweenness_centrality", "MRS_betweenness_centrality"),
                         ("pagerank_centrality", "MRS_pagerank_centrality"),
                         ("eigenvector_centrality", "MRS_eigenvector_centrality")]
                        + ([("betweenness_subgraph_centrality", "MRS_betweenness_subgraph_centrality"),
                            ("pagerank_subgraph_centrality", "MRS_pagerank_subgraph_centrality"),
                            ("eigenvector_subgraph_centrality", "MRS_eigenvector_subgraph_centrality")]
                           if include_subgraph_weightings else [])
                    ),
                    start_date_param=config.get('analysis_parameters', 'context_start_date'),
                    end_date_param=config.get('analysis_parameters', 'context_end_date'),
                    entrez_email=config.get('credentials', 'entrez_email'),
                    entrez_api_key=config.get('credentials', 'entrez_api_key'),
                    calculate_full_centrality=config.get('analysis_parameters', 'calculate_full_centrality')
                )

            # Named for what is in it. "_export.xlsx" told a reader nothing:
            # this is the finished network - every term with its scores, every
            # relation with its weights - as two spreadsheet tabs.
            excel_path = os.path.join(config.secondary_dir,
                                      f"{config.prefix}_final_network_nodes_and_edges.xlsx")
            convert_network_json_to_excel(config.files['final_network'], excel_path)

            # Generate Run-Specific Annotations Template
            _generate_run_annotations(config.files['final_network'], config.files['annotations'], run_anno_path)

            # <<< AOP-INDEPENDENT VISUALIZATIONS >>>
            print("\n<<< Generating Topological & Distribution Figures (AOP-Independent) >>>")
            _, edge_df, _ = load_and_prepare_data(config.files['final_network'], run_anno_path)
            raw_edge_df, raw_all_edges = load_full_raw_data(config.files['full_network'])

            if raw_edge_df is not None and _figure_on(config, 'distribution'):
                _announce_figure_output(config)
                analyze_dispersion(raw_edge_df, str(config.figures_dir), config.prefix)
                plot_cooccurrance_distribution(raw_edge_df, edge_df, str(config.figures_dir), config.prefix)

            if raw_all_edges and os.path.exists(opt_hist_path) and _figure_on(config, 'optimisation'):
                run_optimization_comparison(
                    history_json_path=opt_hist_path,
                    output_dir=str(config.figures_dir),
                    file_prefix=config.prefix
                )

        # <<< STEP 3.5: Secondary Analysis >>>
        sec_params = config.params.get('secondary_analysis', {})

        # Safely determine if we should run secondary analysis based on the execution mode
        should_run_sec = True if args.step == 'secondary' else sec_params.get('export_top_articles', True)

        _skip_secondary = (args.step in ['all', 'secondary']
                           and _reference_skip(config, 'secondary', args.step))

        if args.step in ['all', 'secondary'] and should_run_sec and not _skip_secondary:
            runcontrol.checkpoint("before Secondary Analysis")
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
            sort_metric = sec_params.get('sort_metric', 'F1')
            linear_weight_ars = sec_params.get('linear_weight_ars', 0.5)

            if sec_params.get('export_top_articles', True):
                get_top_network_articles(
                    db_path=config.files['relevance_db'],
                    cleaned_db_path=config.files['cleaned_db'],
                    results_dir=str(config.secondary_dir),
                    file_prefix=config.prefix, limit=export_limit, exclude_reviews=exclude_rev,
                    entrez_email=api_email, entrez_api_key=api_key,
                    sort_metric=sort_metric, linear_weight_ars=linear_weight_ars
                )
            if args.step == 'secondary':
                # The names below have to match a term in the network exactly.
                # A typo, or a heading that simply is not in this corpus, is not
                # an error - the query just returns nothing - so it is checked
                # here and said out loud with the nearest matches offered.
                known = _network_term_index(config.files['final_network'])

                t_nodes = config.get('secondary_analysis', 'target_nodes') or ''
                if t_nodes:
                    for node_name in [n.strip() for n in t_nodes.split(';') if n.strip()]:
                        _warn_unknown_term(node_name, known, 'Target node')
                        analyze_node_relevancy(
                            node_name=node_name, db_path=config.files['relevance_db'],
                            cleaned_db_path=config.files['cleaned_db'],
                            results_dir=str(config.secondary_dir),
                            file_prefix=config.prefix, limit=export_limit, exclude_reviews=exclude_rev,
                            entrez_email=api_email, entrez_api_key=api_key,
                            sort_metric=sort_metric, linear_weight_ars=linear_weight_ars
                        )

                t_edges = config.get('secondary_analysis', 'target_edges') or ''
                if t_edges:
                    for edge_str in [e.strip() for e in t_edges.split(';') if e.strip()]:
                        parts = [p.strip() for p in edge_str.split(' - ')]
                        if len(parts) != 2:
                            print(f"\n  [!] Target edge {edge_str!r} is not two terms "
                                  f"separated by ' - ' (a space, a hyphen, a space).")
                            print( "      Nothing was analysed for it.")
                        if len(parts) == 2:
                            for side in parts:
                                _warn_unknown_term(side, known, 'Target edge term')
                            analyze_edge_relevancy(
                                node1=parts[0], node2=parts[1], db_path=config.files['relevance_db'],
                                cleaned_db_path=config.files['cleaned_db'],
                                results_dir=str(config.secondary_dir),
                                file_prefix=config.prefix, limit=export_limit, exclude_reviews=exclude_rev,
                                entrez_email=api_email, entrez_api_key=api_key,
                                sort_metric=sort_metric, linear_weight_ars=linear_weight_ars
                            )

        # <<< STEP 3.5b: Optional network overlay comparison (default off) >>>
        if (args.step in ['all', 'secondary'] and not _skip_secondary
                and sec_params.get('compare_networks', False)):
            # Look in the active processed folder (data/processed, or reference_processed
            # when reference data is in use); PNG/CSV outputs go to the results section.
            run_network_overlay_comparison(
                network_files=sec_params.get('comparison_networks', ''),
                # Both, in that order: a project built since the networks got
                # their own folder has them there, one built before it has them
                # flat in the processed folder, and a user comparing across two
                # projects may well have one of each.
                search_dir=[str(config.networks_dir), str(config.active_source_dir)],
                results_dir=str(config.secondary_dir),
                file_prefix=config.prefix,
                random_seed=config.get('analysis_parameters', 'random_seed') or 42
            )

        if args.step == 'all':
            if config.params.get('control_flags', {}).get('pause_for_annotation', False):
                anno = os.path.abspath(run_anno_path)
                # Instructions on screen AND in a file beside the one to edit:
                # a paused run can sit for days, and the window that explained
                # what to do is usually long closed by the time anyone returns.
                guide = write_annotation_guide(anno, config)
                print("\n" + "<" * 30 + ">" * 30)
                print("PIPELINE PAUSED - the network needs AOP levels before figures")
                print("<" * 30 + ">" * 30)
                print("\n  The consensus network is built and scored. Each MeSH term in")
                print("  it now needs a biological level assigned, which is a judgement")
                print("  the pipeline cannot make for you.")
                print(f"\n  1. Open this file and fill in the 'aop_level' column:")
                print(f"       {anno}")
                print(f"\n     Allowed levels, in pathway order:")
                print(f"       {', '.join(AOP_LEVELS)}")
                print("     Leave a term as 'Unassigned' if it does not belong to the")
                print("     pathway; it is kept in the network and omitted from the")
                print("     AOP-specific figures. Save as CSV, keeping the semicolons.")
                print(f"\n  2. Then run the figures step:")
                print("       Workbench:  Run -> Step 4 (figures)")
                print("       Terminal :  python -m mesh_aop.cli --step viz")
                if guide:
                    print(f"\n  These instructions are also saved beside the file:")
                    print(f"       {guide}")
                # A marker the Workbench can recognise. Prefixed and on one line
                # so it survives the log reader without being mistaken for prose.
                print(f"\n[PAUSED-FOR-ANNOTATION] {anno}")
                sys.exit(0)
            else:
                print(f"AFK Override Active. Proceeding to visualization with 'Unassigned' AOP levels...")

        # <<< STEP 4: Visualization (AOP-Dependent) >>>
        if args.step in ['all', 'viz']:
            runcontrol.checkpoint("before Step 4")
            print("\n>>> STARTING: Step 4 - Biological Figure Generation")

            # The run annotations are normally written by the network step. When
            # this step is run on its own - which is the whole point of reference
            # mode, where the network already exists - there is nothing to have
            # written them, and the step used to stop dead asking for a file the
            # user had no way to produce. Rebuild it from the stratum dictionary
            # instead; both inputs ship with the reference corpus.
            if not os.path.exists(run_anno_path):
                if (os.path.exists(config.files['final_network'])
                        and os.path.exists(config.files['annotations'])):
                    print("\n    Run annotations absent - rebuilding from the "
                          "stratum dictionary.")
                    _generate_run_annotations(config.files['final_network'],
                                              config.files['annotations'],
                                              run_anno_path)

            _ensure_prerequisites({
                "Final Network JSON": config.files['final_network'],
                "Run-Specific Annotations CSV": run_anno_path
            }, "Step 4 (Biological Figure Generation)")

            is_afk = not config.params.get('control_flags', {}).get('pause_for_annotation', False)
            _sync_run_to_master(run_anno_path, config.files['annotations'], is_afk=is_afk)

            node_df, edge_df, G = load_and_prepare_data(config.files['final_network'], run_anno_path)

            figs = str(config.figures_dir)
            _announce_figure_output(config)
            if _figure_on(config, 'communities'):
                plot_louvain_community_bars(node_df, figs, config.prefix)
            if _figure_on(config, 'tsne'):
                plot_tsne_louvain_overlap(node_df, G, figs, config.prefix)
            if _figure_on(config, 'alluvial'):
                plot_sankey_alluvial(G, node_df, figs, config.prefix)
            if _figure_on(config, 'dendrogram'):
                plot_dendrogram(G, node_df, figs, config.prefix,
                                random_seed=config.get('analysis_parameters', 'random_seed') or 42)
            if _figure_on(config, 'network'):
                plot_network_graph(
                    G, node_df, figs, config.prefix,
                    metric=(config.params.get('viz_parameters', {})
                            .get('network_color_metric', '')),
                    random_seed=config.get('analysis_parameters', 'random_seed') or 42)
            _report_skipped_figures(config, args.step)

        # <<< STEP 5: Ground-Truth Validation & Benchmarking >>>
        gt_enabled = False
        if args.step == 'benchmark':
            runcontrol.checkpoint("before Step 5")
            print("\n>>> STARTING: Step 5 - Ground-Truth Validation & Benchmarking")

            bench_params = config.params.get('benchmark', {})

            # The naive-query baseline needs one MeSH heading. It used to be
            # defaulted to this project's outcome, which meant every other
            # corpus was quietly benchmarked against a heading it might not
            # even contain. Empty is now the default, so: reference mode
            # supplies the reference heading, and an own-data run without one
            # says so and goes on without the baseline rather than inventing it.
            primary_node = (config.get('benchmark', 'primary_node') or '').strip()
            if not primary_node and config.get('control_flags', 'use_reference_data'):
                primary_node = 'Dermatitis, Allergic Contact'
            if not primary_node:
                print("\n  [!] No primary node set, so the naive-query baseline is")
                print("      skipped. Everything else in the benchmark still runs.")
                print("      Set Benchmark -> Primary node to the MeSH heading that")
                print("      names your outcome to enable it.")
            else:
                # A heading that is not in the network gives an empty baseline
                # and no error, which is indistinguishable from a baseline that
                # genuinely found nothing.
                _warn_unknown_term(primary_node,
                                   _network_term_index(config.files['final_network']),
                                   'Primary node')

            # The bundled ground truth describes the reference corpus, so scoring
            # against it only means something when that corpus is in play. Default
            # the analysis on with reference data and off without it; an explicit
            # run_ground_truth_analysis in the config always wins.
            use_reference = bool(config.get('control_flags', 'use_reference_data'))
            gt_enabled = _gt_wanted(config)
            if not gt_enabled:
                print("  [i] Ground-truth analysis is off for this run.")
                print("      Tick 'Run ground-truth analysis' on the Benchmark tab, "
                      "or turn on reference data.")

        if gt_enabled:
            # Own-data runs leave ground_truth_csv empty and auto-detect a file the
            # user dropped in data/raw/ (a recognized name from KNOWN_GT_FILENAMES).
            # The bundled OECD curated set is substituted ONLY when running against
            # the reference corpus, so a user's own file is never silently shadowed.
            # The bundled OECD set, named by where it actually is rather than by
            # a path relative to the working directory. The old spelling -
            # 'data/reference_processed/oecd_ground_truth_curated.xlsx' - was
            # resolved against config.root, which is the directory the process
            # was started in, so on an installed copy it found nothing and the
            # run stopped saying no ground-truth file was configured.
            configured_gt = bench_params.get('ground_truth_csv', '')
            if not configured_gt and use_reference:
                bundled = Path(config.active_source_dir) / 'oecd_ground_truth_curated.xlsx'
                configured_gt = str(bundled) if bundled.exists() else ''

            # Resolve the ground-truth file. A configured name may be a bare filename
            # (looked up in the ACTIVE raw directory -- data/raw/ for own data,
            # data/reference_raw/ under reference data), an absolute path, or a path
            # relative to the project root. Empty triggers auto-detection in raw/.
            gt_path = resolve_ground_truth_path(config.active_raw_dir, configured_gt, root=config.root)
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
            nc_path = resolve_ground_truth_path(config.active_raw_dir, nc_filename, root=config.root) if nc_filename else None

            _build_relevance_db_if_missing(config, include_subgraph_weightings=_gt_wanted(config))
            _ensure_prerequisites({
                "Relevance Database": config.files['relevance_db']
            }, "Step 5 (Benchmarking)")

            print(f"  [+] Using ground-truth file: {os.path.basename(gt_path)}")

            # All --step benchmark artifacts (article-ranking benchmark, ground-truth
            # node/edge validation with its figures, and the validation/projection
            # report) are consolidated under results/benchmark/ so the outputs of this
            # step live in one place rather than scattered across results/.
            bench_dir = config.benchmark_dir
            os.makedirs(bench_dir, exist_ok=True)
            os.makedirs(bench_dir / 'validation', exist_ok=True)

            # Node/edge convergent validation runs first: it takes minutes rather
            # than the benchmark's tens of minutes, and answers a different
            # question (is the network's vocabulary and wiring reproduced?), so a
            # long ranking run never blocks getting the structural result.
            if bench_params.get('run_network_validation', True):
                if os.path.exists(config.files['final_network']):
                    try:
                        run_gt_network_validation(
                            ground_truth_path=gt_path,
                            master_db_path=str(config.files['master_db']),
                            final_network_path=str(config.files['final_network']),
                            output_dir=str(bench_dir),
                            figures_dir=str(bench_dir),
                            file_prefix=config.prefix,
                            weight_key=bench_params.get(
                                'network_validation_weight_key', 'MRS_pagerank_centrality'),
                            min_articles_per_node=bench_params.get('min_articles_per_node', 2),
                            pool_size=bench_params.get('background_pool_size', 50000),
                            random_seed=config.get('analysis_parameters', 'random_seed') or 42
                        )
                    except Exception as e:
                        print(f"\n  [!] WARNING: node/edge validation failed ({e});"
                              f" continuing to the article ranking benchmark.")
                else:
                    print("  [!] Skipping node/edge validation: final network not found.")

            # Consolidated validation report. Evaluates every weighting the pipeline
            # produced across four framings (corpus-wide, within the source query,
            # outside it, and the hybrid reading list) and writes a drill-down
            # workbook plus a narrative report. Independent of run_benchmark, which
            # answers the narrower "how does the primary scorer rank the corpus".
            if bench_params.get('run_validation_report', True):
                if os.path.exists(config.files['final_network']):
                    try:
                        run_validation_report(
                            final_network_file=str(config.files['final_network']),
                            master_db_path=str(config.files['master_db']),
                            pmid_db_path=str(config.files['pmids_db']),
                            ground_truth_file=gt_path,
                            output_dir=str(bench_dir / 'validation'),
                            project_prefix=f"{config.prefix}_",
                            primary_query_term=primary_node,
                            n_boot=bench_params.get('validation_report_n_boot', 2000),
                            random_seed=config.get('analysis_parameters', 'random_seed') or 42,
                            start_date=config.get('analysis_parameters', 'context_start_date'),
                            end_date=config.get('analysis_parameters', 'context_end_date'),
                        )
                    except Exception as e:
                        print(f"\n  [!] WARNING: validation report failed ({e});"
                              f" continuing to the article ranking benchmark.")
                else:
                    print("  [!] Skipping validation report: final network not found.")

            # Article-scoring PROJECTION comparison (normalised vs plain sum vs
            # MRS-weighted vs bipartite vs BM25/baselines). Reuses the validation
            # report's scoring-pool cache; emits a CSV + figure only, no HTML.
            if bench_params.get('run_projection_comparison', True):
                if os.path.exists(config.files['final_network']):
                    try:
                        run_projection_comparison(
                            final_network_file=str(config.files['final_network']),
                            master_db_path=str(config.files['master_db']),
                            pmid_db_path=str(config.files['pmids_db']),
                            ground_truth_file=gt_path,
                            output_dir=str(bench_dir / 'validation'),
                            project_prefix=f"{config.prefix}_",
                            primary_query_term=primary_node,
                            n_boot=bench_params.get(
                                'projection_comparison_n_boot',
                                bench_params.get('validation_report_n_boot', 2000)),
                            random_seed=config.get('analysis_parameters', 'random_seed') or 42,
                            start_date=config.get('analysis_parameters', 'context_start_date'),
                            end_date=config.get('analysis_parameters', 'context_end_date'),
                        )
                    except Exception as e:
                        print(f"\n  [!] WARNING: projection comparison failed ({e});"
                              f" continuing to the article ranking benchmark.")
                else:
                    print("  [!] Skipping projection comparison: final network not found.")

            try:
                run_benchmark(
                    resolved_csv_path=gt_path,
                    relevance_db_path=str(config.files['relevance_db']),
                    output_dir=str(bench_dir),
                    file_prefix=f"{config.prefix}_benchmark",
                    primary_node=primary_node,
                    negative_control_csv=nc_path,
                    random_seed=config.get('analysis_parameters', 'random_seed') or 42,
                    n_boot=bench_params.get('n_boot', 25),
                    n_perm=bench_params.get('n_perm', 25)
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

    except RunAborted as e:
        # Asked to stop, and it reached a checkpoint before being killed. The
        # files on disk are consistent, which is the whole point of stopping
        # here rather than being terminated mid-write.
        print("\n" + "<"*30 + ">"*30)
        print(f"[!] RUN STOPPED - {e}")
        print("<"*30 + ">"*30)
        print("    It stopped at a safe point: no file was left half-written,")
        print("    and completed steps are kept. Running again resumes from the")
        print("    last step that finished.")
        sys.exit(130)
    except Exception as e:
        print("\n" + "<"*30 + ">"*30)
        print(f"[CRITICAL ERROR] Pipeline Execution Failed: {e}")
        print("<"*30 + ">"*30)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\n\n[!] PIPELINE CANCELLED BY USER (Ctrl+C).")
        sys.exit(130)

    total_time = time.time() - total_start

    # <<< RUN LEDGER & PRISMA FLOW REPORT >>>
    # Written last, from the counts the stages returned plus whatever the run
    # left on disk, so a resumed run reports the same totals as a fresh one.
    try:
        ledger = open_ledger(config.results_dir, config.prefix)
        ledger_collect.collect_all(
            ledger, config, args.step, version=get_version(),
            build_stats=build_stats, filter_stats=filter_stats,
            comparison_files=(config.params.get('secondary_analysis', {})
                              .get('comparison_networks', '')
                              if config.params.get('secondary_analysis', {})
                              .get('compare_networks', False) else '')
        )
        ledger.record('run', 'runtime_minutes', f"{total_time/60:.2f}")
        if ledger.save():
            print(f"\n  [+] Run ledger written: {ledger.path}")
        for path in write_prisma_report(ledger, config):
            print(f"  [+] Workflow report written: {path}")
    except Exception as e:
        # A missing report must never turn a completed run into a failed one.
        print(f"\n  [!] Could not write the run ledger / workflow report: {e}")

    print("\n" + "<"*30 + ">"*30)
    print(f"PIPELINE COMPLETED SUCCESSFULLY in {total_time/60:.2f} minutes.")
    print("<"*30 + ">"*30)

if __name__ == "__main__":
    main()