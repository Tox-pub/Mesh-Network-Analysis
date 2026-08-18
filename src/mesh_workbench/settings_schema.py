# -*- coding: utf-8 -*-
"""
schema.py - what the settings form shows, and what it tells the user about it.

One entry per configurable value. The form is generated from this, so adding a
setting here is all it takes to expose it; nothing in the UI hard-codes a field.

Each FIELD carries the help text shown in the description pane when the control
is focused: what it does, what the default is, and - where it matters - what goes
wrong if it is changed carelessly. Those caveats are the ones the pipeline
otherwise only reveals by failing, so they belong in front of the user.

`key` is the dotted path into mesh_config.json.
`kind` is one of: text, int, float, bool, choice, path.
"""

from collections import namedtuple

Field = namedtuple('Field', 'key label kind default what deflt note choices')
Field.__new__.__defaults__ = (None, None)          # note, choices optional


def F(key, label, kind, default, what, deflt, note=None, choices=None):
    return Field(key, label, kind, default, what, deflt, note, choices)


STEPS = [
    ('all',       'all - full pipeline',                    '14 h 20 min'),
    ('process',   'process - MeSH ingest + stop words',     '40 min'),
    ('data_ops',  'data_ops - retrieval & citations',       '6 h'),
    ('network',   'network - co-occurrence + consensus',    '3 h 10 min'),
    ('secondary', 'secondary - top articles & export',      '1 min'),
    ('viz',       'viz - figures',                          '4 min'),
    ('benchmark', 'benchmark - ground truth & validation',  '69 min'),
]

TABS = [
    ('Search', [
        F('search_parameters.search_term', 'Search term', 'text',
          'Dermatitis, Allergic Contact [Mesh]',
          'The PubMed query that defines the starting article set (P0).',
          'Default: the AOP-40 primary heading. Use [Mesh] to search the indexed '
          'heading rather than free text.',
          'Changing this invalidates every downstream result.'),
        F('search_parameters.start_date', 'Start date', 'text', '1950/01/01',
          'Publication window for the initial query.', 'Default: 1950/01/01.',
          'Keep this identical to the context window on the Analysis tab, or the '
          'query set and the scored corpus describe different periods.'),
        F('search_parameters.end_date', 'End date', 'text', '2025/01/01',
          'Publication window for the initial query.', 'Default: 2025/01/01.',
          'Keep this identical to the context window on the Analysis tab.'),
        F('search_parameters.generations_n', 'Citation generations', 'int', 1,
          'How many citation hops to expand from the query result set.',
          'Default: 1. 0 = query only; 1 = adds cited and citing articles; '
          '2 = adds a second hop.',
          'Each generation multiplies retrieval time and disk use. Generation 2 '
          'took this project from 146k to 7.9M articles.'),
        F('control_flags.custom_file_prefix', 'Project prefix', 'text', 'DAC_Mesh',
          'Prefix applied to every output file, so runs do not overwrite each other.',
          'Default: DAC_Mesh.'),
    ]),
    ('Analysis', [
        F('analysis_parameters.random_seed', 'Random seed', 'int', 42,
          'Seeds every stochastic component: sampling, permutations, bootstraps '
          'and tie-breaking in ranking.', 'Default: 42.',
          'Changing it shifts results slightly. Record whatever you use - it is '
          'required to reproduce a run.'),
        F('analysis_parameters.context_start_date', 'Context start', 'text', '1950/01/01',
          'Publication window for scoring and validation.', 'Default: 1950/01/01.',
          'This defines the corpus everywhere - relevance database, validation and '
          'benchmark all use it. Articles outside it are never scored or evaluated.'),
        F('analysis_parameters.context_end_date', 'Context end', 'text', '2025/01/01',
          'Publication window for scoring and validation.', 'Default: 2025/01/01.',
          'This defines the corpus everywhere. Articles outside it are never scored.'),
        F('analysis_parameters.betweenness_k_samples', 'Betweenness samples', 'int', 1000,
          'Node sample size for estimating betweenness on the unfiltered graph.',
          'Default: 1000.',
          'Betweenness on the consensus subgraph is always exact; only the '
          'whole-graph estimate samples. Lower is faster and noisier.'),
        F('analysis_parameters.calculate_full_centrality', 'Whole-graph centrality', 'bool', False,
          'Also compute centrality across the entire unfiltered co-occurrence graph.',
          'Default: off.',
          'Adds hours on a 13,558-node graph. Needed only to compare (full) '
          'against (subgraph) scope.'),
        F('analysis_parameters.eigenvector_max_iter', 'Eigenvector max iterations', 'int', 1000,
          'Iteration cap for the eigenvector centrality solver.', 'Default: 1000.',
          'Raise it only if the solver reports non-convergence.'),
    ]),
    ('Network', [
        F('network_parameters.lambda_val', 'Lambda', 'float', 1.0,
          'Scaling term in the linear weighted additive node model.', 'Default: 1.0.',
          'Rarely changed. Values far from 1 make node weights hard to compare '
          'with published runs.'),
        F('network_parameters.node_weight_factors.centrality', 'Weight: centrality',
          'float', 0.45, 'Contribution of topological centrality to a node weight.',
          'Default: 0.45.', 'The four weights must total 1.00.'),
        F('network_parameters.node_weight_factors.article_rank', 'Weight: article rank',
          'float', 0.15, 'Contribution of normalised article count.',
          'Default: 0.15.', 'The four weights must total 1.00.'),
        F('network_parameters.node_weight_factors.rank_median_cit', 'Weight: median citations',
          'float', 0.20, 'Contribution of normalised median citation count.',
          'Default: 0.20.', 'The four weights must total 1.00.'),
        F('network_parameters.node_weight_factors.rank_total_cit', 'Weight: total citations',
          'float', 0.20, 'Contribution of normalised total citation count.',
          'Default: 0.20.', 'The four weights must total 1.00.'),
    ]),
    ('Consensus', [
        F('simulation_parameters.target_num_edges', 'Target edge count', 'int', 500,
          'Edges each heuristic searches for before the consensus is taken.',
          'Default: 500, applied to GLF and SA alike.',
          'A property of the search, not the result. The consensus intersection is '
          'smaller - 500 and 500 gave 328 edges on the reference run.'),
        F('simulation_parameters.glf_iterations', 'GLF iterations', 'int', 5000000,
          'Metropolis proposals for Graph Likelihood Filtering.',
          'Default: 5,000,000.',
          'Below about 1M the search may not converge, the two heuristics agree '
          'less, and the consensus shrinks. Above 5M gains little for a lot of time.'),
        F('simulation_parameters.sa_iterations', 'SA iterations', 'int', 5000000,
          'Metropolis proposals for simulated annealing.', 'Default: 5,000,000.',
          'Below about 1M the search may not converge and the consensus shrinks.'),
        F('simulation_parameters.sa_temp_start', 'SA start temperature', 'float', 5000.0,
          'Initial temperature for simulated annealing.', 'Default: 5000.',
          'Too low and SA behaves greedily, converging on the same local optimum '
          'as GLF and defeating the point of taking a consensus.'),
        F('simulation_parameters.sa_cooling_rate', 'SA cooling rate', 'float', 0.999995,
          'Geometric factor applied to the temperature each iteration.',
          'Default: 0.999995.',
          'Must be just under 1. Cooling faster freezes the search early; 0.99999 '
          'or lower effectively turns SA into a greedy search.'),
    ]),
    ('Benchmark', [
        F('benchmark.ground_truth_csv', 'Ground truth file', 'path',
          'data/reference_processed/oecd_ground_truth_curated.xlsx',
          'Curated positive set used for validation and benchmarking.',
          'Default: the bundled OECD AOP-40 set (96 PMIDs). Leave empty to '
          'auto-detect a file in data/raw.',
          'Supply your own as .xlsx or .csv with a PMID column. Articles carrying '
          'no network term are unreachable by every method and cap achievable recall.'),
        F('benchmark.primary_node', 'Primary node', 'text', 'Dermatitis, Allergic Contact',
          'MeSH heading used for the naive-query baseline.',
          'Default: the search term primary heading.',
          'Must exist in the final network or the baseline is skipped.'),
        F('benchmark.validation_report_n_boot', 'Bootstrap replicates', 'int', 2000,
          'Resamples used for confidence intervals.', 'Default: 2000.',
          'Sets the resolution of empirical p-values - B = 2000 cannot report '
          'anything below 0.0005. Lower is faster but widens intervals.'),
        F('benchmark.background_pool_size', 'Background pool', 'int', 50000,
          'Random articles sampled to estimate base rates and permutation nulls.',
          'Default: 50000.',
          'Smaller pools make enrichment and the null distributions noisier.'),
        F('benchmark.min_articles_per_node', 'Min articles per node', 'int', 2,
          'Ground-truth terms must appear in this many articles to enter the '
          'reference network.', 'Default: 2.',
          'At 1 the reference network fills with terms mentioned once and the '
          'overlap statistics lose meaning.'),
        F('benchmark.run_network_validation', 'Node / edge validation', 'bool', True,
          'Test network overlap with the ground truth against permutation nulls.',
          'Default: on.'),
        F('benchmark.run_projection_comparison', 'Projection comparison', 'bool', True,
          'Compare the eight article-scoring projections.', 'Default: on.',
          'Adds roughly 20 minutes to the benchmark step.'),
    ]),
    ('Folders', [
        F('directories.output_dir', 'Output folder', 'text', '',
          'Where every result is written - workbooks, figures, reports, logs.',
          'Default: empty, meaning the project\'s own results folder.',
          'Point this somewhere else before a trial run. A step re-run with the '
          'default writes over the existing results, and the outputs a manuscript '
          'was written from are not recoverable from the pipeline alone.'),
        F('directories.input_dir', 'Input folder', 'text', '',
          'Where the pipeline looks for source data.',
          'Default: empty, meaning the project\'s own data folder.',
          'Changing this without moving the databases will make the run fail at '
          'the first step that needs them.'),
        F('control_flags.pause_for_annotation', 'Pause for AOP annotation', 'bool', True,
          'Stop after the network is built so biological strata can be assigned '
          'by hand before the figures are drawn.', 'Default: on.',
          'With this on, a full run stops partway and waits. Turn it off for an '
          'unattended run - strata then stay Uncategorized and the biological '
          'figures lose their meaning.'),
        F('control_flags.use_reference_data', 'Use bundled reference data', 'bool', False,
          'Run against the shipped reference corpus instead of your own retrieval.',
          'Default: off.'),
    ]),
    ('Credentials', [
        F('credentials.entrez_email', 'NCBI e-mail', 'text', '',
          'Identifies you to Entrez, as NCBI requires.', 'Default: none.',
          'Required for retrieval. NCBI may block requests that do not carry it.'),
        F('credentials.entrez_api_key', 'NCBI API key', 'text', '',
          'Raises the Entrez rate limit from 3 to 10 requests per second.',
          'Default: none.',
          'Optional but strongly recommended - retrieval is several times faster '
          'with one. Stored on this machine only.'),
    ]),
]

WEIGHT_KEYS = [f.key for f in TABS[2][1] if f.key.startswith('network_parameters.node_weight_factors')]


def get(cfg, dotted, default=None):
    """Read a dotted path out of the config dict."""
    cur = cfg
    for part in dotted.split('.'):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


def put(cfg, dotted, value):
    """Write a dotted path into the config dict, creating levels as needed."""
    parts = dotted.split('.')
    cur = cfg
    for p in parts[:-1]:
        cur = cur.setdefault(p, {})
    cur[parts[-1]] = value
