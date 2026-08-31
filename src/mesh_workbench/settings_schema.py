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


# No time estimates. They were measured on one machine and are dominated by
# how much memory it has and how fast its disk is - the figure shown was
# several times the real runtime on hardware with enough RAM, which is worse
# than saying nothing. "all" excludes the benchmark; that IS worth saying,
# because otherwise it is discovered by looking for results never produced.
STEPS = [
    ('all',       'all - full pipeline (not the benchmark)'),
    ('process',   'process - MeSH ingest + stop words'),
    ('data_ops',  'data_ops - retrieval & citations'),
    ('network',   'network - co-occurrence + consensus'),
    ('secondary', 'secondary - top articles & export'),
    ('viz',       'viz - figures'),
    ('benchmark', 'benchmark - ground truth & validation'),
]

TABS = [
    ('Search', [
        # First on the first tab, deliberately. It decides most of the fields
        # below it, and putting it here means a user sees them grey out the
        # moment it is ticked instead of finding them locked later with no
        # visible cause.
        F('control_flags.use_reference_data',
          'Use bundled reference data (demonstration only)', 'bool', False,
          'FOR DEMONSTRATION ONLY. Analyse the finished reference network that '
          'shipped with the program, so you can see what a completed run looks '
          'like - the figures, the report, the benchmark - in minutes, before '
          'committing to a real retrieval of your own.',
          'Default: off - the normal case is analysing your own corpus.',
          'FOR DEMONSTRATION ONLY. This is not a way to do new research; it is '
          'a way to look at finished output before you spend a day making your '
          'own.\n'
          '\n'
          'A first real run means a long PubMed download and a full rebuild '
          'before you see a single figure. Tick this instead and the program '
          'analyses a network that already exists - the allergic contact '
          'dermatitis network published with this software, included in the '
          'download, already built and scored. Figures, the workflow report '
          'and the benchmark come out in minutes, so you can judge whether the '
          'outputs are what you want and learn where everything lands on disk '
          'first.\n'
          '\n'
          'What you cannot do with it is change the corpus. The articles behind '
          'this network are fixed; the retrieval that produced them is not '
          'repeated and cannot be varied. Findings from it are the published '
          'findings, not yours - untick it and run your own search for that.\n'
          '\n'
          'The network, the AOP strata and the ground-truth set are read from '
          'the program folder. The settings describing that corpus - search '
          'term, both date windows, citation depth, random seed, benchmark '
          'primary node - are set for you and shown greyed out, because '
          'anything else would label the figures with a query that did not '
          'make them. Retrieval will not run: there is nothing to fetch.\n'
          '\n'
          'Yours either way: figure resolution and formats, which figures are '
          'drawn, your folders, your credentials, and your own project. Your '
          'saved settings are not altered - untick and they return. Files it '
          'produces are named Reference_ so they cannot be confused with '
          'yours.\n'
          '\n'
          'The reference networks are copied into your own networks folder the '
          'first time you tick this, and everything downstream reads them from '
          'there. They are yours to open, edit or delete like any other result. '
          'Copies already in that folder are left alone, so your changes '
          'survive later runs; delete one and the pristine original comes back '
          'from the program folder.'),
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
        F('control_flags.custom_file_prefix', 'Project prefix', 'text', '',
          'Prefix applied to every file this project produces, so two projects '
          'never overwrite each other.',
          'Ships empty - name it for your own question. Ticking "Use bundled '
          'reference data" sets it to DAC_Mesh and locks it, because that is '
          'the name the published corpus was built under.',
          'This is what separates one analysis from another. Keep it the same '
          'to resume a run: each step checks for its own output first and skips '
          'work already done. Change it to start clean while leaving an earlier '
          'analysis intact - two prefixes never touch the same files.'),
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
          'Your own outputs: figures, workbooks, the run ledger and the workflow '
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
        F('simulation_parameters.target_num_edges', 'Target edge count (per optimiser)', 'int', 500,
          'How many edges EACH optimiser searches for. Not how many you get.',
          'Default: 500 - so 500 for graph likelihood filtering (GLF) and 500 '
          'for simulated annealing (SA), independently.',
          'The network is built by consensus. Two different optimisers each '
          'pick the set of edges they think best explains the co-occurrence '
          'counts, starting from different places and searching differently, '
          'and only the edges BOTH chose are kept. An edge one optimiser liked '
          'and the other did not is discarded, which is what makes the result a '
          'consensus rather than one heuristic\'s opinion.\n'
          '\n'
          'So the number you get is always smaller than the number you ask for, '
          'and how much smaller tells you how well the two agreed. On the '
          'reference run 500 and 500 agreed on 328 edges. Raising this widens '
          'both searches and usually widens the consensus too, but never to the '
          'number set here.'),
        F('simulation_parameters.glf_iterations',
          'Graph likelihood filtering (GLF) iterations', 'int', 5000000,
          'Metropolis proposals for graph likelihood filtering - GLF, the first '
          'of the two optimisers whose consensus forms the network.',
          'Default: 5,000,000.',
          'Below about 1M the search may not converge, the two heuristics agree '
          'less, and the consensus shrinks. Above 5M gains little for a lot of time.'),
        F('simulation_parameters.sa_iterations',
          'Simulated annealing (SA) iterations', 'int', 5000000,
          'Metropolis proposals for simulated annealing - SA, the second of the '
          'two optimisers whose consensus forms the network.',
          'Default: 5,000,000.',
          'Below about 1M the search may not converge and the consensus shrinks.'),
        F('simulation_parameters.sa_temp_start',
          'Simulated annealing (SA) start temperature', 'float', 5000.0,
          'Initial temperature for simulated annealing.', 'Default: 5000.',
          'Too low and SA behaves greedily, converging on the same local optimum '
          'as GLF and defeating the point of taking a consensus.'),
        F('simulation_parameters.sa_cooling_rate',
          'Simulated annealing (SA) cooling rate', 'float', 0.999995,
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
          'Export the articles behind one or more named terms.',
          'Default: empty, meaning no per-term export.',
          'Semicolon-separated, e.g. "Skin; Haptens".\n'
          '\n'
          'An article is returned when the term is one of the MeSH headings it '
          'was indexed with AND that heading is a node in your network - the '
          'articles that actually put the term there. Matching is on the WHOLE '
          'heading, so "Skin" returns articles indexed under Skin and not those '
          'under Skin Diseases, Skin Absorption or Skin Tests.\n'
          '\n'
          'Because of that, the name has to be spelled exactly as the network '
          'spells it: "Dermatitis, Allergic Contact", not "allergic contact '
          'dermatitis". Take the spellings from the network Excel export or '
          'from the network JSON in Cytoscape. A name that is not in the '
          'network finds nothing, so what you type is checked first and the '
          'nearest matches are offered.\n'
          '\n'
          'Only the secondary step reads this; a full run says it skipped it.'),
        F('secondary_analysis.target_edges', 'Target edges', 'text', '',
          'Export the articles behind one or more named relationships.',
          'Default: empty, meaning no per-edge export.',
          'Written as "NodeA - NodeB; NodeC - NodeD" - a space, a hyphen, a '
          'space between the two headings.\n'
          '\n'
          'An article is returned when BOTH headings are among the MeSH terms '
          'it was indexed with, so these are the articles that put the edge '
          'there rather than either node alone. Both names match on the whole '
          'heading and both must be spelled as the network spells them.'),
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
          'Rank the ground-truth articles eight different ways and compare how '
          'well each recovers them.', 'Default: on.',
          'The question it answers is which way of turning a network into an '
          'article score actually finds the curated papers. Eight are compared, '
          'as a 4x2 grid: four node weightings - raw centrality, MRS, article '
          'count and the adjusted node weight - each applied two ways, summed '
          'over the MeSH terms an article carries and averaged over them. Sum '
          'favours articles indexed with many relevant terms; mean favours '
          'articles that are tightly on-topic.\n'
          '\n'
          'Each is scored by BEDROC, enrichment factor, recall and mean average '
          'precision, with bootstrap intervals, so the comparison is between '
          'ranked lists rather than between single numbers. Adds roughly 20 '
          'minutes to the benchmark step.'),
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
          'choice', 'MRS_pagerank_centrality',
          'Which node score is tested against the ground truth.',
          'Default: MRS_pagerank_centrality.',
          'The validation asks a single question: if the terms in the network '
          'are ranked by this number, do the ones that appear in the curated '
          'ground-truth papers come out on top? It reports the rank '
          'correlation against how often each term appears in those papers, '
          'the AUC separating attested terms from unattested ones, and - '
          'because a term that is simply common everywhere would score well by '
          'accident - the same figures after controlling for corpus-wide '
          'frequency.\n'
          '\n'
          'The choice is a 2x3x2 grid: raw or MRS, betweenness or PageRank or '
          'eigenvector, whole-corpus or subgraph. MRS is the score adjusted by '
          'how strongly the literature supports the term; subgraph means '
          'centrality measured within the consensus network rather than the '
          'whole corpus.\n'
          '\n'
          'The six subgraph options only exist when the ground-truth analysis '
          'was on at the time the network was BUILT. Picking one the network '
          'does not carry fails the validation step.',
          ['MRS_pagerank_centrality',
           'MRS_betweenness_centrality',
           'MRS_eigenvector_centrality',
           'MRS_pagerank_subgraph_centrality',
           'MRS_betweenness_subgraph_centrality',
           'MRS_eigenvector_subgraph_centrality',
           'pagerank_centrality',
           'betweenness_centrality',
           'eigenvector_centrality',
           'pagerank_subgraph_centrality',
           'betweenness_subgraph_centrality',
           'eigenvector_subgraph_centrality']),
    ]),
    ('Figures', [
        F('_heading.which_figures', 'Which figures to draw', 'heading', '',
          'Untick anything you do not need. Each is drawn from the finished '
          'network, so leaving one out costs nothing but that figure - and a '
          'few of them are slow.', '', ''),
        F('viz_parameters.fig_distribution', 'Figure 1 - Edge weight distribution',
          'bool', True,
          'How often each pair of MeSH terms appears in the same article, plotted before and after consensus filtering.',
          'Default: on.',
          'The evidence that filtering removed the long tail rather than the '
          'signal. Drawn by the network step, not the figures step.'),
        F('viz_parameters.fig_optimisation', 'Figure 2 - Optimisation trajectory',
          'bool', True,
          'The two optimisers - generalised local filtering and simulated annealing - traced over the search for the consensus subgraph.',
          'Default: on.',
          'Only has anything to draw if the optimisation history was written, '
          'which the network step does as it runs.'),
        F('viz_parameters.fig_communities', 'Figure 3 - Community composition',
          'bool', True,
          'Stacked bars showing which AOP levels make up each Louvain community, so you can see whether the communities follow the pathway.',
          'Default: on.',
          'Needs AOP levels assigned. With everything left Unassigned it draws '
          'one bar and says nothing.'),
        F('viz_parameters.fig_tsne', 'Figure 4 - t-SNE projection',
          'bool', True,
          'A two-dimensional projection of the graph distance matrix, each node coloured by its Louvain community - do the communities separate?',
          'Default: on.',
          'The slowest of the static figures on a large network, and the '
          'projection is stochastic - it is reproducible only because the '
          'random seed on the Analysis tab fixes it.'),
        F('viz_parameters.fig_alluvial', 'Figure 5 - AOP alluvial flow',
          'bool', True,
          'An interactive Sankey tracing flow from stressors through molecular, cellular, tissue and organ levels to adverse outcomes.',
          'Default: on.',
          'Writes both a labelled and an unlabelled HTML version, plus the '
          'connection table behind them. The one figure that shows the pathway '
          'as a pathway.'),
        F('viz_parameters.fig_dendrogram', 'Figure 6 - Node2Vec dendrogram',
          'bool', True,
          'Ward clustering of Node2Vec embeddings, leaves coloured by AOP level - terms that sit near each other in the learned space cluster together.',
          'Default: on.',
          'Trains an embedding first, so it is the most expensive figure here. '
          'Reproducible run to run: the walks are seeded and gensim is trained '
          'single-threaded so the result does not depend on thread timing.'),
        F('viz_parameters.fig_network', 'Figure 7 - Network graph',
          'bool', True,
          'The consensus network drawn, with every term labelled and the nodes '
          'coloured by one metric.',
          'Default: on.',
          'An overview rather than a working view - enough to see where the '
          'hubs are and how the communities sit against each other. Open the '
          'network JSON in Cytoscape for anything closer. The layout is '
          'seeded from the random seed on the Analysis tab, so the same '
          'network draws the same way every time.'),
        F('viz_parameters.network_color_metric', 'Figure 7 - colour by',
          'choice', 'MRS_pagerank_centrality',
          'Which per-node number decides the colour of each term in Figure 7.',
          'Default: MRS_pagerank_centrality.',
          'Colours are viridis, min-max scaled across the terms drawn, so the '
          'full range of the scale is used however narrow the spread is - dark '
          'purple is the lowest value present, yellow the highest. The colour '
          'bar carries the real numbers. A metric this network does not have - '
          'the subgraph centralities, when the ground-truth analysis was off - '
          'falls back to the first one it does have, and says so.',
          ['MRS_pagerank_centrality',
           'MRS_betweenness_centrality',
           'MRS_eigenvector_centrality',
           'MRS_pagerank_subgraph_centrality',
           'MRS_betweenness_subgraph_centrality',
           'MRS_eigenvector_subgraph_centrality',
           'pagerank_centrality',
           'betweenness_centrality',
           'eigenvector_centrality',
           'pagerank_subgraph_centrality',
           'betweenness_subgraph_centrality',
           'eigenvector_subgraph_centrality',
           'adjusted_node_weight',
           'article_count',
           'degree',
           'clustering_coefficient',
           'major_topic_proportion',
           'article_count_rank_normalized',
           'rank_norm_mean_cit',
           'rank_norm_median_cit',
           'rank_norm_total_cit']),
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

    'Secondary': (
        'Target nodes and edges are NOT part of a full run - start it yourself',
        'Choosing "all - full pipeline" runs the export of top network-wide '
        'articles, but it does NOT run the target nodes and edges you type '
        'below. Each target pulls hundreds of records from NCBI, which would '
        'add a long tail to every run, so they wait until you ask for them.\n'
        '\n'
        'To run them, pick "secondary - targeted queries & exports" from the '
        'step list and press Run. Everything they need is on disk after a full '
        'run, so this can be done at any point afterwards - days later is fine, '
        'as long as the project prefix is the same.\n'
        '\n'
        'WHERE TO FIND THE NAMES TO TYPE. A target has to be a MeSH heading '
        'that is actually in your network, spelled exactly as the network '
        'spells it - "Dermatitis, Allergic Contact", not "allergic contact '
        'dermatitis". Two places give you the real list:\n'
        '\n'
        '  1. THE NETWORK JSON, in your networks folder. This is the easiest '
        'way to see what you are choosing between. Open it in Cytoscape - File '
        '> Import > Network from File - and the whole graph is drawn for you: '
        'click any node or edge and its exact id appears in the table below, '
        'ready to copy. Any tool that reads Cytoscape JSON will do; Gephi and '
        'the networkx library read these files too.\n'
        '\n'
        '  2. THE NETWORK EXCEL EXPORT in your results folder, written by '
        '"Export network to Excel" on this tab. One row per node and per edge, '
        'so you can sort by centrality and copy the names straight out.\n'
        '\n'
        'The "Open network folder" button on the Results screen takes you to '
        'the JSON files.\n'
        '\n'
        'Type several nodes separated by semicolons: Skin; Haptens. An edge is '
        'two headings with a space-hyphen-space between them: Skin - '
        'Dermatitis, Allergic Contact. A name that is not in the network is not '
        'an error - the query simply finds nothing - so the run checks what you '
        'typed first and offers the nearest matches it can see.\n'
        '\n'
        'This step needs the cleaned citation database from your own retrieval, '
        'so it cannot run against the bundled reference corpus, which ships as '
        'a finished network with no article lists behind it.'),

    'Benchmark': (
        'The benchmark is NOT part of a full run - start it yourself',
        'Choosing "all - full pipeline" runs retrieval, the network, the '
        'secondary analysis and the figures. It does NOT run the benchmark. '
        'That is deliberate: the benchmark scores every weighting against every '
        'article in the corpus, which takes a good while on a full-sized one, '
        'and most runs do not need it.\n'
        '\n'
        'To run it, pick "benchmark - ground truth & validation" from the step '
        'list and press Run. Everything it needs is already on disk by then, so '
        'it can be run at any point after a full run finishes - including days '
        'later, as long as the project prefix is the same.\n'
        '\n'
        'WHERE TO PUT THE GROUND-TRUTH FILE. The benchmark scores your network '
        'against a set of articles you already believe are relevant, and it '
        'cannot start without one. Put that file in your RAW DATA folder - the '
        'Folders tab shows exactly where that is on this machine - under one of '
        'these names, and it is picked up on its own:\n'
        '\n'
        '    ground_truth_pmids.csv        ground_truth.csv\n'
        '    ground_truth_pmids.txt        ground_truth.txt\n'
        '    oecd_resolved_citations.csv\n'
        '\n'
        'Any other name works too if you type it into "Ground truth CSV" below; '
        'a bare name is looked for in that same raw folder, and a full path is '
        'used as given. Accepted contents: a CSV or TSV with a PMID column, a '
        'single column of PMIDs, or a plain text file with one PMID per line. '
        'If nothing is found the run stops and prints the folder it searched.\n'
        '\n'
        'The file you used is copied into the results as benchmark/inputs, '
        'under your project prefix, so the numbers stay readable months later '
        'when the original has moved or changed.\n'
        '\n'
        'WHAT THE NETWORK ITSELF LOOKS LIKE. The benchmark reports on the '
        'network JSON in your networks folder. To see the thing being scored - '
        'which nodes and edges survived, and how they connect - open that JSON '
        'in Cytoscape (File > Import > Network from File), or in any tool that '
        'reads Cytoscape JSON. The "Open network folder" button on the Results '
        'screen takes you to the files. The primary node named below has to be '
        'a heading spelled exactly as it appears there.\n'
        '\n'
        'Its results go to a benchmark folder inside your results folder, kept '
        'apart from the figures and the secondary analysis, and grouped by the '
        'question each one answers:\n'
        '\n'
        '  inputs              the ground truth the run was scored against,\n'
        '                      kept under your project prefix so the numbers\n'
        '                      stay readable later\n'
        '  ranking             how well each weighting ranks articles\n'
        '  ranking_validation  every weighting compared, with confidence\n'
        '                      intervals and the projection comparison\n'
        '  network_validation  whether the network\'s own terms and links are\n'
        '                      reproduced by the ground truth\n'
        '\n'
        'One setting here reaches backwards. "Run ground-truth analysis" also '
        'decides whether the subgraph centralities are computed, and those are '
        'written while the NETWORK is built - so it has to be set before the '
        'network step, not before the benchmark. Turning on "Use bundled '
        'reference data" forces it on, because the published reference run '
        'cannot be reproduced without those six node attributes.'),
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
