# -*- coding: utf-8 -*-
"""
ledger_collect.py - fills the run ledger from whatever the run left on disk.

The pipeline skips any stage whose output already exists, so a run that is
resumed - or re-run to regenerate figures - executes almost nothing and would
otherwise report almost nothing. Every count here is therefore derived from the
artefacts themselves: the PMID databases, the fetch logs, the subgraph JSONs.
Read a second time, they give the same answer as the first run did.

Where a stage IS executing, the in-flight counts it returns are preferred - they
are exact, and cost nothing - and these functions fill the gaps around them.
Nothing here raises: a ledger is a record of the work, never a gate on it.
"""

import json
import os
import sqlite3

from .run_ledger import open_ledger, ledger_path  # noqa: F401  (re-exported)


def _safe(fn, *args, **kwargs):
    """Run a probe, swallow anything it throws, and return None on failure."""
    try:
        return fn(*args, **kwargs)
    except Exception:
        return None


def _gen_counts(db_path):
    """{generation: n} from a pmids database, or None if it cannot be read."""
    if not db_path or not os.path.exists(str(db_path)):
        return None
    conn = sqlite3.connect(f'file:{db_path}?mode=ro', uri=True)
    try:
        cur = conn.cursor()
        cur.execute("SELECT generation, COUNT(*) FROM pmids_table GROUP BY generation")
        return {str(g): int(n) for g, n in cur.fetchall()}
    finally:
        conn.close()


def _mesh_coverage(db_path):
    """(with_mesh, without_mesh) in a cleaned database, or None."""
    if not db_path or not os.path.exists(str(db_path)):
        return None
    conn = sqlite3.connect(f'file:{db_path}?mode=ro', uri=True)
    try:
        cur = conn.cursor()
        cur.execute("PRAGMA table_info(pmids_table)")
        if 'mesh_terms' not in [row[1] for row in cur.fetchall()]:
            return None
        cur.execute("SELECT COUNT(*) FROM pmids_table WHERE mesh_terms IS NOT NULL AND TRIM(mesh_terms) != ''")
        with_mesh = int(cur.fetchone()[0])
        cur.execute("SELECT COUNT(*) FROM pmids_table")
        return with_mesh, int(cur.fetchone()[0]) - with_mesh
    finally:
        conn.close()


def _line_count(path):
    """Data rows in a TSV log - the header, if any, does not count as a record."""
    if not path or not os.path.exists(str(path)):
        return 0
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        rows = [ln for ln in fh if ln.strip()]
    if rows and not rows[0].split('\t')[0].strip().isdigit():
        rows = rows[1:]
    return len(rows)


def _resolve_pair(metrics, stats, prior, path, lazy=True):
    """Resolve a (nodes, edges) pair from stats, then the ledger, then the file.

    Each half is settled on its own: a stage that returned only one of the two
    should not lose it just because its partner is missing. The file is read at
    most once, and only if something is still unresolved - it can be several
    hundred megabytes.
    """
    out = []
    from_disk = None
    for metric in metrics:
        value = stats.get(metric)
        if value is None:
            value = prior.get(metric)
        if value is None and path:
            if lazy and from_disk is None:
                from_disk = _safe(_graph_size, path) or (None, None)
            if from_disk:
                value = from_disk[metrics.index(metric)]
        out.append(value)
    return out


def _graph_size(path):
    """(nodes, edges) from a Cytoscape JSON, or None."""
    if not path or not os.path.exists(str(path)):
        return None
    with open(path, 'r', encoding='utf-8') as fh:
        elements = json.load(fh).get('elements', {})
    return len(elements.get('nodes', [])), len(elements.get('edges', []))


def _community_count(path):
    """Distinct Louvain communities in the filtered network, or None."""
    if not path or not os.path.exists(str(path)):
        return None
    with open(path, 'r', encoding='utf-8') as fh:
        nodes = json.load(fh).get('elements', {}).get('nodes', [])
    ids = {n.get('data', {}).get('filtered_louvain_community_id') for n in nodes}
    ids.discard(None)
    return len(ids)


# --------------------------------------------------------------------- stages

def record_run(ledger, config, step, version=''):
    """The search itself: what was asked for, and of what."""
    search = config.params.get('search_parameters', {})
    analysis = config.params.get('analysis_parameters', {})
    ledger.clear_stage('run')
    ledger.record('run', 'file_prefix', config.prefix)
    ledger.record('run', 'pipeline_step', step)
    if version:
        ledger.record('run', 'workbench_version', version)
    ledger.record('run', 'search_term', search.get('search_term', ''),
                  'The PubMed query that seeded generation P0')
    ledger.record('run', 'search_start_date', search.get('start_date', ''))
    ledger.record('run', 'search_end_date', search.get('end_date', ''))
    ledger.record('run', 'generations_requested', search.get('generations_n', ''),
                  'Citation hops expanded from the seed set')
    ledger.record('run', 'context_start_date', analysis.get('context_start_date', ''))
    ledger.record('run', 'context_end_date', analysis.get('context_end_date', ''))
    ledger.record('run', 'random_seed', analysis.get('random_seed', ''))
    ledger.record('run', 'reference_data',
                  'on' if config.params.get('control_flags', {}).get('use_reference_data') else 'off')


def record_identification(ledger, config):
    """Records retrieved, by citation generation, before and after cleaning."""
    raw = _safe(_gen_counts, config.files.get('pmids_db'))
    cleaned = _safe(_gen_counts, config.files.get('cleaned_db'))
    if raw is None and cleaned is None:
        return False

    ledger.clear_stage('identification')
    if raw:
        for gen in sorted(raw):
            ledger.record('identification', f'records_retrieved_{gen}', raw[gen],
                          'P0 is the seed search; G1+ are citation hops')
        ledger.record('identification', 'records_retrieved_total', sum(raw.values()))
    if cleaned:
        for gen in sorted(cleaned):
            ledger.record('identification', f'records_retained_{gen}', cleaned[gen])
        total = sum(cleaned.values())
        ledger.record('identification', 'records_retained_total', total)
        if raw:
            ledger.record('identification', 'records_removed_final_generation',
                          sum(raw.values()) - total,
                          'The outermost generation is dropped: its citation links are incomplete')
    return True


def record_screening(ledger, config, stats=None):
    """What annotation retrieval and the term screen let through."""
    stats = stats or {}
    prior = ledger.snapshot_stage('screening')
    ledger.clear_stage('screening')

    failed = _safe(_line_count, config.files.get('failed_mesh')) or 0
    empty = _safe(_line_count, config.files.get('empty_mesh')) or 0
    ledger.record('screening', 'records_failed_mesh_retrieval', failed,
                  'PubMed did not answer for these PMIDs')
    ledger.record('screening', 'records_no_mesh_available', empty,
                  'Retrieved, but carrying no MeSH annotation')

    coverage = _safe(_mesh_coverage, config.files.get('cleaned_db'))
    if coverage and 'articles_with_mesh' not in stats:
        stats = dict(stats, articles_with_mesh=coverage[0], articles_without_mesh=coverage[1])

    notes = {
        'articles_screened': 'Records reaching the term screen',
        'articles_with_mesh': 'Records carrying at least one MeSH heading',
        'articles_without_mesh': 'Records excluded for having no MeSH heading',
        'annotations_total': 'Every heading on every record, subheadings included',
        'annotations_subheadings_removed': 'Qualifier forms (Term/qualifier) dropped in favour of the main heading',
        'stop_words_in_list': 'Size of the MeSH stop-word list applied as a filter',
        'terms_before_stop_words': 'Distinct headings before the stop-word filter',
        'terms_removed_as_stop_words': 'Distinct headings matched by the stop-word list',
        'terms_after_stop_words': 'Distinct headings surviving to network construction',
    }
    for metric in ('articles_screened', 'articles_with_mesh', 'articles_without_mesh',
                   'annotations_total', 'annotations_subheadings_removed',
                   'stop_words_in_list', 'terms_before_stop_words',
                   'terms_removed_as_stop_words', 'terms_after_stop_words'):
        value = stats.get(metric, prior.get(metric))
        if value is not None:
            ledger.record('screening', metric, value, notes.get(metric, ''))
    return True


def record_pruning(ledger, config, build_stats=None, filter_stats=None):
    """The full graph, the two optimisers, and what their consensus discarded."""
    build_stats = build_stats or {}
    filter_stats = filter_stats or {}
    # Preferred in order: counts from a stage that just ran, then what this
    # ledger already recorded, then the artefact itself. The middle step keeps
    # a resumed run from re-reading a several-hundred-megabyte network JSON to
    # recover a number it wrote down last time.
    prior = ledger.snapshot_stage('pruning')
    ledger.clear_stage('pruning')

    n, e = _resolve_pair(('full_network_nodes', 'full_network_edges'),
                         build_stats, prior, config.files.get('full_network'))
    if n is not None:
        ledger.record('pruning', 'full_network_nodes', n, 'Every term the corpus produced')
    if e is not None:
        ledger.record('pruning', 'full_network_edges', e, 'Every co-occurrence the corpus produced')

    for key, path_key in (('glf', 'glf_subgraph'), ('sa', 'sa_subgraph')):
        label = 'Graph Likelihood Filtering' if key == 'glf' else 'Simulated Annealing'
        n, e = _resolve_pair((f'{key}_nodes', f'{key}_edges'), filter_stats, prior,
                             config.files.get(path_key))
        if n is not None:
            ledger.record('pruning', f'{key}_nodes', n, f'{label} optimiser')
        if e is not None:
            ledger.record('pruning', f'{key}_edges', e, f'{label} optimiser')

    notes = {
        'target_edges': 'Edge budget the optimisers were asked to hit',
        'pruning_nodes_considered': 'Nodes in the unfiltered graph',
        'pruning_edges_considered': 'Edges in the unfiltered graph',
        'pruning_nodes_excluded': 'Nodes neither optimiser kept',
        'pruning_edges_excluded': 'Edges neither optimiser kept',
        'consensus_nodes': 'Nodes in the consensus network',
        'consensus_edges': 'Edges both optimisers agreed on',
        'consensus_nodes_excluded': 'Nodes kept by only one optimiser',
        'consensus_edges_excluded': 'Edges kept by only one optimiser',
    }
    for metric, note in notes.items():
        value = filter_stats.get(metric, prior.get(metric))
        if value is not None:
            ledger.record('pruning', metric, value, note)
    return True


def record_included(ledger, config, filter_stats=None):
    """The largest connected component - what the analysis actually reports on."""
    filter_stats = filter_stats or {}
    prior = ledger.snapshot_stage('included')
    ledger.clear_stage('included')

    n, e = _resolve_pair(('lcc_nodes', 'lcc_edges'), filter_stats, prior,
                         config.files.get('consensus_lcc'))
    if n is not None:
        ledger.record('included', 'lcc_nodes', n, 'MeSH terms in the largest connected component')
    if e is not None:
        ledger.record('included', 'lcc_edges', e, 'Co-occurrence relations in the LCC')

    notes = {
        'components_found': 'Separate components in the consensus network',
        'components_excluded': 'Smaller components dropped when the LCC was taken',
        'largest_excluded_component': 'Size of the largest component that was dropped',
        'lcc_nodes_excluded': 'Consensus nodes outside the LCC',
        'lcc_edges_excluded': 'Consensus edges outside the LCC',
    }
    for metric, note in notes.items():
        value = filter_stats.get(metric, prior.get(metric))
        if value is not None:
            ledger.record('included', metric, value, note)

    communities = _safe(_community_count, config.files.get('final_network'))
    if communities is None:
        communities = _safe(_community_count, config.files.get('consensus_lcc'))
    if communities:
        ledger.record('included', 'louvain_communities', communities,
                      'Communities detected within the LCC')

    scored = _safe(_graph_size, config.files.get('final_network'))
    if scored:
        ledger.record('included', 'terms_scored_mrs', scored[0],
                      'Terms carrying a Mean Relevancy Score')
    return True


def record_unipartite(ledger, config, comparison_files=''):
    """The optional overlay comparison against other networks."""
    names = [n.strip() for n in str(comparison_files).split(';') if n.strip()]
    if not names:
        return False
    ledger.clear_stage('unipartite')
    ledger.record('unipartite', 'comparison_networks', len(names),
                  'Networks overlaid against this run')
    for i, name in enumerate(names, 1):
        ledger.record('unipartite', f'comparison_network_{i}', name)
    return True


def collect_all(ledger, config, step, version='', build_stats=None,
                filter_stats=None, screening_stats=None, comparison_files=''):
    """Fill every stage, preferring in-flight counts and falling back to disk."""
    _safe(record_run, ledger, config, step, version)
    _safe(record_identification, ledger, config)
    _safe(record_screening, ledger, config, screening_stats or build_stats)
    _safe(record_pruning, ledger, config, build_stats, filter_stats)
    _safe(record_included, ledger, config, filter_stats)
    _safe(record_unipartite, ledger, config, comparison_files)
    return ledger
