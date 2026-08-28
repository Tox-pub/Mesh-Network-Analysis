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
        F('search_parameters.search_term', 'Search term', 'text', '',
          'The PubMed query that defines the starting article set (P0).',
          'Default: empty. There is no sensible query for another research '
          'question, so this has to be filled in before retrieval will run. Use '
          '[Mesh] to search the indexed heading rather than free text.',
          'The bundled reference corpus was built with '
          '"Dermatitis, Allergic Contact [Mesh]" over 1950/01/01 to 2025/01/01. '
          'Enter that to reproduce it. This field, the dates and the generation '
          'count together define one corpus: change any of them and the run '
          'produces a different, self-consistent set of results rather than a '
          'wrong one. Give it its own project prefix and both sets survive '
          'side by side.'),
        F('search_parameters.start_date', 'Start date', 'text', '1950/01/01',
          'Publication window for the initial query. YYYY/MM/DD, YYYY/MM or YYYY.',
          'Default: 1950/01/01, which is earlier than PubMed goes - so it takes '
          'everything. Leave it empty for the same effect.',
          'Keep this identical to the context window on the Analysis tab, or the '
          'query set and the scored corpus describe different periods.'),
        F('search_parameters.end_date', 'End date', 'text', '2025/01/01',
          'Publication window for the initial query. YYYY/MM/DD, YYYY/MM or YYYY.',
          'Default: 2025/01/01. Leave it empty to search every date PubMed has, '
          'or type TODAY for the date the run starts.',
          'Both dates are EDAT - the date the record entered PubMed, not the '
          'date the paper was published. Keep this identical to the context '
          'window on the Analysis tab. A fixed date is what makes a search '
          'repeatable: TODAY and an empty field both move, so the same query '
          'run next month returns a different corpus.'),
        F('search_parameters.generations_n', 'Citation generations', 'int', 1,
          'How many citation hops to expand from the query result set.',
          'Default: 1. 0 = the query result only; 1 = adds the articles it '
          'cites and those citing it; 2 = adds a second hop outward.',
          'Each hop multiplies retrieval time and disk use, and the growth is '
          'steep - as an example, moving from generation 1 to generation 2 took '
          'one project from about 146,000 articles to 7.9 million. Choose this '
          'deliberately: for most projects generation 0 or 1 is the appropriate '
          'choice, and 2 is worth it only when the question genuinely needs the '
          'wider literature.'),
        F('search_parameters.update_mesh_support_files', 'Force-refresh MeSH support files',
          'bool', False,
          'Rebuild the MeSH descriptor and stop-word support files instead of '
          'reusing the ones already on disk.',
          'Default: off - they are built once and kept.',
          'Turn this on after moving to a new MeSH release year. It adds time to '
          'the process step and changes the vocabulary every later step sees.'),
        F('control_flags.custom_file_prefix', 'Project prefix', 'text', 'DAC_Mesh',
          'Prefix applied to every output file, so runs do not overwrite each other.',
          'Default: DAC_Mesh.'),
    ]),
    ('Folders', [
        # Two settings decide everything. The three after them are overrides
        # that predate the two - they are in the original upload, while
        # data_dir arrived later - and are normally left blank. Presenting all
        # five as equals made it look as though results lived inside output.
        F('directories.data_dir', 'Data folder', 'text', '',
          'Everything downloaded and everything built from it. Two subfolders '
          'are created inside: raw (the PubMed archives and the master '
          'annotation database) and processed (networks and score databases).',
          'Default: empty, meaning a private folder under your user profile. '
          'ABOUT 52 GB ENDS UP HERE - roughly 44 GB of downloaded archive and '
          'an 8 GB database built from it.',
          'Choose this before downloading anything, and choose a drive with '
          '60 GB free. The archive can be deleted afterwards from the Data '
          'setup screen, which gets most of it back, but moving ~52 GB later '
          'is far slower than picking the right drive now.'),
        F('directories.results_dir', 'Results folder', 'text', '',
          'Your own outputs: figures, workbooks, the run ledger and the PRISMA '
          'report. A figures subfolder is created inside it.',
          'Default: empty, meaning Documents\\MeSH Workbench for an installed '
          'copy, or the program folder for a portable one. Small - megabytes, '
          'not gigabytes.',
          'Deliberately not inside the data folder: results are your work and '
          'the data folder is a cache that can be deleted and rebuilt. Two '
          'people sharing one installed copy also get their own results here '
          'instead of overwriting each other.'),

        F('_heading.advanced', 'Advanced - normally left blank', 'heading', '',
          'The three settings below override parts of the two above. Leave them '
          'empty unless you have a specific reason, such as splitting the '
          'downloads and the working files across two drives.', '', ''),

        F('directories.input_dir', 'Override: raw data folder', 'text', '',
          'Replaces the raw subfolder of the Data folder - where the PubMed '
          'archive and the master annotation database are kept.',
          'Default: empty, meaning <Data folder>\\raw. Leave it that way unless '
          'the archive has to live on a different drive from everything else.',
          'This is an override, not an extra folder: setting it PINS the raw '
          'location, so later changing the Data folder will not move it. If the '
          'databases already exist, moving them by hand is your job - the run '
          'fails at the first step that needs them.'),
        F('directories.output_dir', 'Override: working files folder', 'text', '',
          'Replaces the processed subfolder of the Data folder - the networks '
          'and score databases the pipeline builds.',
          'Default: empty, meaning <Data folder>\\processed.',
          'A second, older behaviour to know about: if the Results folder above '
          'is blank, results are written here too. That is why this setting can '
          'look as though results live inside it. Set the Results folder '
          'explicitly and the two stay separate.'),
        F('directories.etl_workspace_dir', 'Override: database build workspace',
          'text', '',
          'Scratch space used while the master annotation database is compiled.',
          'Default: empty, meaning the system temp folder.',
          'The workspace holds a full working copy of the database - several GB '
          '- and the system temp folder is on the system drive. If C: is short '
          'of space, point this at a roomier disk or the build can run out '
          'partway through, hours in.'),

        F('_heading.behaviour', 'Run behaviour', 'heading', '', '', '', ''),

        F('control_flags.pause_for_annotation', 'Pause for AOP annotation', 'bool', False,
          'Stop after the network is built so biological strata can be assigned '
          'by hand before the figures are drawn.', 'Default: off.',
          'Off, a full run completes unattended and every term stays '
          'Uncategorized, which drains the biological figures of meaning. Turn '
          'it on to stop partway and assign the strata first - that is the '
          'intended workflow for a real analysis, and it needs you present.'),
        F('control_flags.use_reference_data', 'Use bundled reference data', 'bool', False,
          'Run against the shipped reference corpus instead of your own retrieval.',
          'Default: off - the normal case is analysing your own corpus.',
          'Turn it on to reproduce the published reference run, or to draw '
          'figures before you have retrieved anything. While it is on, the run '
          'scores the shipped network rather than yours.'),
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
    ('Analysis', [
        F('analysis_parameters.random_seed', 'Random seed', 'int', 42,
          'Seeds every stochastic component: sampling, permutations, bootstraps '
          'and tie-breaking in ranking.', 'Default: 42.',
          'Changing it shifts results slightly. Record whatever you use - it is '
          'required to reproduce a run.'),
        F('analysis_parameters.context_start_date', 'Context start', 'text', '1950/01/01',
          'Publication window for scoring and validation. YYYY/MM/DD, YYYY/MM or YYYY.',
          'Default: 1950/01/01. Leave it empty to score every date available.',
          'This defines the corpus everywhere - relevance database, validation and '
          'benchmark all use it. Articles outside it are never scored or evaluated.'),
        F('analysis_parameters.context_end_date', 'Context end', 'text', 'TODAY',
          'Publication window for scoring and validation. YYYY/MM/DD, YYYY/MM or YYYY.',
          'Default: TODAY, resolved to the date the run starts. Leave it empty '
          'to score every date available.',
          'This defines the corpus everywhere - articles outside it are never '
          'scored. Two cautions: TODAY moves, so the same analysis repeated '
          'next month covers a different corpus; and it can reach past your '
          'local data, since the master database only holds what the last '
          'baseline and update download contained. Set a fixed date at or '
          'before that point for a reproducible run.'),
        F('analysis_parameters.betweenness_k_samples', 'Betweenness samples', 'int', 1000,
          'Node sample size for estimating betweenness on the unfiltered graph.',
          'Default: 1000.',
          'Betweenness on the consensus subgraph is always exact; only the '
          'whole-graph estimate samples. Lower is faster and noisier.'),
        F('analysis_parameters.calculate_full_centrality', 'Whole-graph centrality', 'bool', False,
          'Also compute centrality across the entire unfiltered co-occurrence graph.',
          'Default: off, which is much faster.',
          'Centrality on the consensus subgraph is computed either way - this '
          'adds the whole-graph figures on top, and costs hours on a '
          '13,558-node graph. Turn it on only to compare (full) against '
          '(subgraph) scope.'),
        F('analysis_parameters.eigenvector_max_iter', 'Eigenvector max iterations', 'int', 1000,
          'Iteration cap for the eigenvector centrality solver.', 'Default: 1000.',
          'Raise it only if the solver reports non-convergence.'),
        F('analysis_parameters.eigenvector_tol', 'Eigenvector tolerance', 'float', 1.0e-6,
          'Convergence threshold for the eigenvector centrality solver.',
          'Default: 0.000001.',
          'Raising it lets the solver stop on a vector that has not settled; '
          'lowering it can exhaust the iteration cap above and fail the step.'),
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
    ('Secondary', [
        F('secondary_analysis.export_top_articles', 'Export top articles', 'bool', True,
          'Write the ranked article table for the whole network once the run ends.',
          'Default: on.'),
        F('secondary_analysis.export_limit', 'Export limit', 'int', 500,
          'How many articles the export keeps.', 'Default: 500.'),
        F('secondary_analysis.sort_metric', 'Sort metric', 'text', 'F1',
          'How relevance and citation standing are combined into one ranking.',
          "Default: F1. Accepts 'F1' or 'Linear'.",
          'F1 is the harmonic mean, so an article weak on either side drops '
          'sharply. Linear is the weighted average set by the field below.'),
        F('secondary_analysis.linear_weight_ars', 'Linear ARS weight', 'float', 0.5,
          'Share of the Linear metric given to relevance; citations take the rest.',
          'Default: 0.5, and must fall between 0.01 and 0.99.',
          'Read only when the sort metric is Linear. 0.70 means 70% relevance, '
          '30% citations.'),
        F('secondary_analysis.exclude_reviews', 'Exclude reviews', 'bool', True,
          'Drop review articles from the exports.', 'Default: on.',
          'Reviews gather citations across a whole field, so leaving them in '
          'pushes primary studies down a citation-weighted ranking.'),
        F('secondary_analysis.target_nodes', 'Target nodes', 'text', '',
          'Limit the secondary step to named terms.',
          'Default: empty, meaning the whole network.',
          'Semicolon-separated, e.g. "Term A; Term B". Only the secondary step '
          'reads this; a full run ignores it.'),
        F('secondary_analysis.target_edges', 'Target edges', 'text', '',
          'Limit the secondary step to named term pairs.',
          'Default: empty, meaning the whole network.',
          'Written as "NodeA - NodeB; NodeC - NodeD".'),
        F('secondary_analysis.compare_networks', 'Compare multiple networks', 'bool', False,
          'Also compare this run against other saved networks.', 'Default: off.'),
        F('secondary_analysis.comparison_networks', 'Networks to compare', 'text', '',
          'The networks the comparison reads.',
          'Default: empty. Read only when comparison is on.',
          'Comma-separated file names, e.g. "Graph_1.json","Graph_2.graphml". '
          'A bare name resolves to the processed data folder.'),
    ]),
    ('Benchmark', [
        F('benchmark.ground_truth_csv', 'Ground truth file', 'path',
          '',
          'Curated positive set used for validation and benchmarking.',
          'Default: empty, which auto-detects a file you drop in data/raw and '
          'otherwise falls back to the bundled set when reference data is in '
          'use. That bundled example is '
          'data/reference_processed/oecd_ground_truth_curated.xlsx - the OECD '
          'AOP-40 set of 96 PMIDs for the allergic contact dermatitis query.',
          'Supply your own as .xlsx or .csv with a PMID column. Articles carrying '
          'no network term are unreachable by every method and cap achievable recall.'),
        F('benchmark.primary_node', 'Primary node', 'text', '',
          'A single MeSH heading, used as the naive-query baseline the ranking '
          'is measured against.',
          'Default: empty. Enter the one heading that best names your outcome - '
          'for example "Dermatitis, Allergic Contact", which is a self-contained '
          'MeSH heading and the one the bundled reference corpus uses.',
          'It must be a real MeSH heading, spelled as MeSH spells it, and it '
          'must appear in your final network - otherwise the baseline is '
          'skipped and the benchmark reports without it.'),
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
        F('benchmark.run_ground_truth_analysis', 'Run ground-truth analysis', 'bool', False,
          'Score the network against the curated positive set.',
          'Default: off. Turning on "Use bundled reference data" forces it on '
          'regardless, because the published reference run cannot be '
          'reproduced without it.',
          'This also decides whether subgraph PageRank is computed, and the '
          'network step reads it - so set it before building the network, not '
          'afterwards.'),
        F('benchmark.n_boot', 'Benchmark bootstrap resamples', 'int', 25,
          "Resamples for the benchmark's own confidence intervals.",
          'Default: 25.',
          'Separate from "Bootstrap replicates" above, which belongs to the '
          'validation report. n = 200 takes roughly two hours.'),
        F('benchmark.n_perm', 'Permutation iterations', 'int', 25,
          'Null-model iterations behind the node and edge enrichment tests.',
          'Default: 25.',
          'Sets the floor on reportable p-values. The published run used 1000 '
          'node and 500 edge permutations.'),
        F('benchmark.network_validation_weight_key', 'Node validation weight key',
          'text', 'MRS_pagerank_centrality',
          'Node attribute ranked when testing overlap against the ground truth.',
          'Default: MRS_pagerank_centrality.',
          'Written as [MRS_]{betweenness|pagerank|eigenvector}[_subgraph]_centrality. '
          'A name the network does not carry fails the validation step.'),
    ]),
    ('Figures', [
        F('_heading.which_figures', 'Which figures to draw', 'heading', '',
          'Untick anything you do not need. Each is drawn from the finished '
          'network, so leaving one out costs nothing but that figure - and a '
          'few of them are slow.', '', ''),
        F('viz_parameters.fig_distribution', 'Fig 1 - Edge weight distribution',
          'bool', True,
          'Co-occurrence counts before and after consensus filtering, overlaid.',
          'Default: on.',
          'The evidence that filtering removed the long tail rather than the '
          'signal. Drawn by the network step, not the figures step.'),
        F('viz_parameters.fig_optimisation', 'Fig 2 - Optimisation comparison',
          'bool', True,
          'GLF against simulated annealing over the filtering search.',
          'Default: on.',
          'Only has anything to draw if the optimisation history was written, '
          'which the network step does as it runs.'),
        F('viz_parameters.fig_communities', 'Fig 3 - Community composition',
          'bool', True,
          'AOP-level make-up of each Louvain community, as stacked bars.',
          'Default: on.',
          'Needs AOP levels assigned. With everything left Unassigned it draws '
          'one bar and says nothing.'),
        F('viz_parameters.fig_tsne', 'Fig 5 - t-SNE with communities',
          'bool', True,
          't-SNE of the graph distance matrix, coloured by community.',
          'Default: on.',
          'The slowest of the static figures on a large network, and the '
          'projection is stochastic - it is reproducible only because the '
          'random seed on the Analysis tab fixes it.'),
        F('viz_parameters.fig_alluvial', 'Fig 6 - AOP alluvial flow',
          'bool', True,
          'Interactive Sankey from stressors through to adverse outcomes.',
          'Default: on.',
          'Writes both a labelled and an unlabelled HTML version, plus the '
          'connection table behind them. The one figure that shows the pathway '
          'as a pathway.'),
        F('viz_parameters.fig_dendrogram', 'Node2Vec dendrogram',
          'bool', True,
          'Ward clustering of Node2Vec embeddings, leaves coloured by AOP level.',
          'Default: on.',
          'Trains an embedding first, so it is the most expensive figure here. '
          'Reproducible run to run: the walks are seeded and gensim is trained '
          'single-threaded so the result does not depend on thread timing.'),
        F('_heading.figure_output', 'Output format', 'heading', '',
          'How every figure ticked above is written to disk.', '', ''),
        F('viz_parameters.figure_dpi', 'Resolution (dpi)', 'int', 300,
          'Dots per inch for every raster figure the run produces.',
          'Default: 300, the usual journal minimum for a raster figure.',
          '600 is what print production normally asks for, and roughly '
          'quadruples both the time to write each figure and its file size. '
          'Values are clamped to 50-1200.'),
        F('viz_parameters.figure_formats', 'File formats', 'text', 'jpeg,tif',
          'Which image files to write for each figure.',
          'Default: jpeg,tif - a preview to look at and a lossless copy to '
          'submit.',
          'Comma-separated, from jpeg, png, tif, pdf and svg. TIFF is written '
          'with LZW compression because journals ask for lossless and an '
          'uncompressed 600 dpi panel runs to tens of megabytes. Listing more '
          'formats writes every figure more than once.'),
    ]),
]

# Found by key, not by tab position: indexing TABS positionally breaks silently
# the moment a tab is inserted before this one.
WEIGHT_KEYS = [f.key for _, fields in TABS for f in fields
               if f.key.startswith('network_parameters.node_weight_factors')]

# Real settings only. A heading carries no key into the config file, and code
# that walks the schema to save, validate or audit must skip it.
SETTINGS = [(tab, f) for tab, fields in TABS for f in fields if f.kind != 'heading']


# Standing explanations pinned to the foot of a tab, for things that are true of
# the tab as a whole rather than of one control. The description pane below the
# form only shows the field you are on, so anything you need to know BEFORE
# choosing a value has nowhere else to appear.
TAB_NOTES = {
    'Search': (
        'The project prefix decides which run you are working on',
        'Every file a run produces is named with this prefix. That is what makes '
        'a run resumable, and what decides whether you continue one or overwrite '
        'it.\n'
        '\n'
        'KEEP THE SAME PREFIX to pick up where you left off. Each step checks for '
        'its own output first and skips the work if it is already there, so a run '
        'interrupted after retrieval will go straight to the network step, and a '
        'run that reached the figures will only redraw them. This is how you '
        'resume after a pause, an abort, or a machine that was shut down.\n'
        '\n'
        'CHANGE THE PREFIX to start clean, or to keep an earlier analysis intact '
        'while you try something different. Two prefixes never touch each other\'s '
        'files, so both sets of results survive side by side.\n'
        '\n'
        'Re-running a step whose output already exists overwrites it. Before a run '
        'starts, anything that would be replaced is listed for you to confirm - '
        'and if the prefix has not been used for that step, there is nothing to '
        'warn about and you are not asked.'),
}


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
