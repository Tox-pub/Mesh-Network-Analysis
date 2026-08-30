# -*- coding: utf-8 -*-
"""
config_parser.py - configuration loading and path resolution.

Defines the MeshConfig object that every stage of the pipeline shares. It
implements a two-tier configuration system: a set of immutable factory defaults
overlaid by the user's local mesh_config.json, with helpers to read, update,
save, and reset individual parameters.

From the resolved settings it derives every working directory and file path the
pipeline uses (raw inputs, processed databases, network JSONs, results, logs),
switching automatically between the user's own data set and the bundled
reference data set.
"""

import os
import json
import copy
from datetime import date
from pathlib import Path

class MeshConfig:
    """Parses user configuration, applies factory defaults, and generates absolute paths."""

    def __init__(self, config_path: str = None, workspace_root: str = None):
        """
        Initializes the configuration object.

        Args:
            config_path: Path to the user's local mesh_config.json
            workspace_root: The root directory containing the data/ and results/ folders.
                            Defaults to the Current Working Directory if not provided.
        """
        # An installed copy is shared and often read-only, so settings live in
        # the user's own profile; a portable copy keeps them in its folder. See
        # paths.py - both front-ends resolve this the same way.
        from . import paths as _paths
        self.config_path = Path(config_path) if config_path else _paths.config_path()
        self.root = Path(workspace_root) if workspace_root else Path.cwd()

        self._define_factory_defaults()
        self.load_config()
        self.refresh_paths()

    def _define_factory_defaults(self):
        """Immutable baseline settings for the entire pipeline."""
        self.factory_defaults = {
            "control_flags": {
                # Off: the default is a user analysing their own corpus. Turn it
                # on to reproduce the published reference run.
                "use_reference_data": False,
                "custom_file_prefix": "DAC_Mesh",
                # Off, matching the documented behaviour: the run completes
                # unattended with every term left 'Unassigned'. Turn it on to
                # stop after Step 3 and assign AOP levels by hand first.
                "pause_for_annotation": False
            },
            "directories": {
                "input_dir": "",
                "output_dir": "",
                # Scratch space for the database build. Empty means the system
                # temp directory, which is usually the system drive - the
                # workspace holds a full copy of the master database, so a small
                # C: is a real constraint.
                "etl_workspace_dir": "",
                # Downloaded archives and the databases built from them. Empty
                # means the per-user location for an installed copy, or the
                # program folder for a portable one. See paths.py.
                "data_dir": "",
                # The user's own outputs. Kept apart from the databases so an
                # installed copy can give each account its own results.
                "results_dir": ""
            },
            "credentials": {
                "entrez_email": "",
                "entrez_api_key": ""
            },
            "search_parameters": {
                # Blank: a fresh install should not silently inherit this
                # project's query. The reference corpus was built with
                # "Dermatitis, Allergic Contact [Mesh]" over 1950/01/01 to
                # 2025/01/01, which the settings help repeats.
                "search_term": "",
                "start_date": "1950/01/01",
                # A fixed date, not TODAY. A moving end date makes the same
                # search return a different corpus every time it is run, which
                # is the one thing a reproducible pipeline must not do quietly.
                # TODAY still works if it is typed in; it is just not the
                # default, and the help says what it costs.
                "end_date": "2025/01/01",
                "generations_n": 1,
                "update_mesh_support_files": False
            },
            "analysis_parameters": {
                # Hours on a 13,558-node graph. Subgraph centrality is computed
                # either way, so off is the useful default and on is opt-in.
                "calculate_full_centrality": False,
                "random_seed": 42,
                "betweenness_k_samples": 1000,
                "context_start_date": "1950/01/01",
                "context_end_date": "TODAY",
                # Both of these have fallbacks in network.py, so the pipeline
                # ran without them - but the settings window showed its own
                # defaults for keys the file did not contain, which is how a
                # displayed value and an applied value come to differ. Seeded
                # here so there is one answer rather than two.
                "eigenvector_max_iter": 1000,
                "eigenvector_tol": 1.0e-6
            },
            # Figure output. 300 dpi is the usual journal minimum for raster
            # figures and is fast; 600 is what print production asks for and
            # roughly quadruples both time and file size.
            "viz_parameters": {
                "figure_dpi": 300,
                "figure_formats": "jpeg,tif",
                # One switch per figure, all on. Absent from a user's file means
                # on as well - a settings file written before these existed must
                # not be read as "draw nothing".
                "fig_distribution": True,
                "fig_optimisation": True,
                "fig_communities": True,
                "fig_tsne": True,
                "fig_alluvial": True,
                "fig_dendrogram": True,
                "fig_network": True,
                "network_color_metric": "MRS_pagerank_centrality"
            },
            "network_parameters": {
                "lambda_val": 1.0,
                "node_weight_factors": {
                    "centrality": 0.45,
                    "article_rank": 0.15,
                    "rank_median_cit": 0.20,
                    "rank_total_cit": 0.20
                }
            },
            "simulation_parameters": {
                "target_num_edges": 500,
                "glf_iterations": 5000000,
                "sa_iterations": 5000000,
                "sa_temp_start": 5000.0,
                "sa_cooling_rate": 0.999995
            },
            "benchmark": {
                # Empty by default so an own-data run auto-detects a file dropped in
                # data/raw/ (e.g. ground_truth_pmids.csv). The bundled OECD curated
                # set is substituted by the CLI only when 'Use Reference Data' is on,
                # so a user's own file is never silently shadowed.
                "ground_truth_csv": "",
                "negative_control_csv": "",
                # Empty: this is the user's outcome heading, not this
                # project's. Shipping one silently benchmarked every
                # corpus against a heading it may not even contain.
                "primary_node": "",
                "n_boot": 25,
                "n_perm": 25,
                "run_network_validation": True,
                "run_projection_comparison": True,
                "network_validation_weight_key": "MRS_pagerank_centrality",
                "min_articles_per_node": 2,
                "background_pool_size": 50000,
                # Bootstrap depth for the validation report, kept separate from
                # n_boot so the cheap benchmark and the expensive report can be
                # tuned independently.
                "validation_report_n_boot": 2000
                # run_ground_truth_analysis is deliberately NOT seeded here. Its
                # default is not a constant: absent, it follows 'Use reference
                # data', so the bundled corpus benchmarks itself and a user's own
                # corpus does not until they ask. Writing a fixed value would
                # pin it and lose that. The settings window still offers it, and
                # saving from there pins it on purpose.
            },
            # Step 3.5. Declared here so saved values survive a reload: the
            # wizard used to create this section in memory only, and load_config
            # dropped whatever it had written because the section was unknown.
            "secondary_analysis": {
                "export_top_articles": True,
                "export_limit": 500,
                "sort_metric": "F1",
                "linear_weight_ars": 0.5,
                "exclude_reviews": True,
                "target_nodes": "",
                "target_edges": "",
                "compare_networks": False,
                "comparison_networks": ""
            }
        }

    def load_config(self):
        """Loads local JSON and deep-merges it over the factory defaults."""
        self.params = copy.deepcopy(self.factory_defaults)
        if self.config_path.exists():
            try:
                with open(self.config_path, 'r') as f:
                    user_params = json.load(f)
                for section, values in user_params.items():
                    # A leading underscore marks a switch that belongs to one run
                    # only - _run_baseline_etl and friends. Reading those back
                    # would repeat a multi-hour rebuild on every later run, so
                    # they are deliberately not restored.
                    if section.startswith('_'):
                        continue
                    # setdefault, not a membership test: skipping sections absent
                    # from the defaults silently discarded everything the user
                    # had saved under them.
                    if isinstance(values, dict):
                        self.params.setdefault(section, {}).update(values)
                    else:
                        self.params[section] = values
            except Exception as e:
                print(f"\n[!] Failed to parse {self.config_path}. Using factory defaults. ({e})")

    def refresh_paths(self):
        """Triggers path resolution if directories or control flags are updated."""
        self._resolve_dynamic_values()
        self._build_directories()
        self._map_files()

    # HOW THE PUBLISHED REFERENCE CORPUS WAS BUILT.
    #
    # Ticking "Use bundled reference data" analyses that corpus, so these are
    # the settings that describe it. They are applied over whatever is in the
    # form, because the alternative is a run that draws the dermatitis network
    # while its ledger and PRISMA report describe the user's own query - a
    # provenance document stating something false, which is the one thing a
    # PRISMA report exists not to do.
    #
    # Only what defines the CORPUS is here. Figure resolution, file formats,
    # which figures to draw, the folders, the credentials and the worker count
    # are the user's: they describe presentation and environment, and changing
    # them does not change what is being reproduced.
    REFERENCE_RUN = {
        ('search_parameters', 'search_term'): 'Dermatitis, Allergic Contact [Mesh]',
        ('search_parameters', 'start_date'): '1950/01/01',
        ('search_parameters', 'end_date'): '2025/01/01',
        ('search_parameters', 'generations_n'): 1,
        ('analysis_parameters', 'context_start_date'): '1950/01/01',
        ('analysis_parameters', 'context_end_date'): '2025/01/01',
        # Fixes t-SNE, the Node2Vec walks and the Figure 7 layout, so the
        # published figures come out as published rather than merely equivalent.
        ('analysis_parameters', 'random_seed'): 42,
        ('benchmark', 'primary_node'): 'Dermatitis, Allergic Contact',
    }

    # Values the user may type that mean "the day the run starts". Matched
    # case- and space-insensitively because a settings form is typed into by
    # hand: "today" and " TODAY " meant the same thing to everyone except the
    # code, which sent them to PubMed literally.
    _TODAY_WORDS = {'today', 'now'}
    _DYNAMIC_DATES = (('analysis_parameters', 'context_end_date'),
                      ('search_parameters', 'end_date'))

    def _resolve_dynamic_values(self):
        """Handles date interpretation and reference flags.

        Resolution is kept OUT of self.params, in a shadow dict that get() reads
        through. It used to be written back over the stored value, which meant
        the first save() after any load replaced "TODAY" with that day's date -
        permanently. The setting could be chosen once and then silently became a
        frozen date that never moved again, and nothing said so. self.params is
        now exactly what the user typed, and save() persists exactly that.
        """
        self.use_reference_data = self.params['control_flags']['use_reference_data']
        self.prefix = "DAC_Mesh" if self.use_reference_data else self.params['control_flags']['custom_file_prefix']

        # Databases derived from the shipped corpus are named apart from the
        # user's own. Without this they collide outright: the factory default
        # project prefix IS "DAC_Mesh", and the reference corpus forces that
        # same prefix, so a user who had scored their own corpus and then ticked
        # "Use bundled reference data" wrote both to
        # DAC_Mesh_mean_relevancy.db - one overwriting the other, with no
        # warning and no way to tell afterwards which corpus a score came from.
        #
        # The bundled NETWORKS keep their published names; they are read-only
        # and are what the reference corpus is. Only what we derive is renamed.
        self.db_prefix = f"Reference_{self.prefix}" if self.use_reference_data else self.prefix
        today = date.today().strftime("%Y/%m/%d")
        self._resolved = {}

        # The reference corpus first, so a TODAY typed into a date field cannot
        # then move a window that is supposed to be fixed at what was published.
        #
        # These go in the shadow dict, exactly like the TODAY substitution, and
        # for the same reason: self.params stays what the user typed, so save()
        # keeps it and unticking the box brings their own settings straight
        # back. Nothing of theirs is overwritten on disk.
        self.reference_overrides = []
        if self.use_reference_data:
            for (section, key), value in self.REFERENCE_RUN.items():
                was = self.params.get(section, {}).get(key)
                self._resolved[(section, key)] = value
                if str(was if was is not None else '') != str(value):
                    self.reference_overrides.append((f'{section}.{key}', was, value))

        for section, key in self._DYNAMIC_DATES:
            if (section, key) in self._resolved:
                continue                      # the reference window is fixed
            raw = self.params.get(section, {}).get(key)
            if isinstance(raw, str) and raw.strip().lower() in self._TODAY_WORDS:
                self._resolved[(section, key)] = today

    def _annotations_path(self):
        """The AOP stratum dictionary, preferring the user's own copy.

        It lives in the raw data folder, which on a fresh install is empty - so
        the figure step failed on a clean machine even in reference mode, asking
        for a run-annotations file it could not generate without this. The copy
        shipped beside the reference corpus is the fallback, so an install can
        reproduce the reference figures before retrieving anything.
        """
        from . import paths as _paths
        own = self.active_raw_dir / "aop_annotations_master.csv"

        # Reference mode reproduces a published run, so the published stratum
        # dictionary wins over anything the user has annotated for their own
        # corpus. Their file is not touched; it is simply not the one that
        # describes the reference network.
        if self.use_reference_data:
            bundled = _paths.bundled_reference_dir(self.root) / "aop_annotations_master.csv"
            if bundled.exists():
                return bundled

        if own.exists():
            return own
        # The packaged build copies the dictionary in beside the reference
        # corpus; a source checkout keeps it in data/raw, which is not the
        # active raw folder while reference mode is on. Try both, so reference
        # mode draws its figures from a clone as readily as from an installer.
        for fallback in (_paths.bundled_reference_dir(self.root) / "aop_annotations_master.csv",
                         self.root / "data" / "raw" / "aop_annotations_master.csv"):
            if fallback.exists():
                return fallback
        return own

    def get_default_directories(self, use_ref_data: bool):
        """Return the default (input, output) directories for the active data mode.

        Reference mode reads the corpus shipped with the program, which the
        packaging step places beside the packages rather than under data/ - so
        it is resolved by paths.py rather than assumed to sit below the project
        root. Own-data mode writes to the user's data folder, which for an
        installed copy is under their profile and not beside the executable.
        """
        from . import paths as _paths
        if use_ref_data:
            ref_out = _paths.bundled_reference_dir(self.root)
            return str(self.root / 'data' / 'reference_raw'), str(ref_out)
        base_data = Path(self.params.get('directories', {}).get('data_dir', '').strip()
                         or _paths.default_user_data_dir(self.root))
        return str(base_data / 'raw'), str(base_data / 'processed')

    def _build_directories(self):
        """Constructs and ensures the existence of required directories with override support."""
        custom_in = self.params.get('directories', {}).get('input_dir', '').strip()
        custom_out = self.params.get('directories', {}).get('output_dir', '').strip()

        default_in, default_out = self.get_default_directories(self.use_reference_data)

        # THE RAW FOLDER IS ALWAYS THE USER'S, like the databases.
        #
        # It used to follow the reference checkbox to data/reference_raw, which
        # is not shipped, does not exist on an installed copy, and is resolved
        # against the working directory - so with the box ticked, Step 0 would
        # have downloaded a second 44 GB copy of the PubMed baseline into
        # whatever folder the program happened to be started in, and Step 1 a
        # second copy of the MeSH descriptor XML beside it. Nothing in the
        # benchmark path touches these, which is the only reason it went
        # unnoticed.
        #
        # Downloads are downloads whichever corpus is being analysed. Only the
        # networks and the ground-truth files come from the bundle.
        self.active_raw_dir = Path(custom_in) if custom_in else self._user_raw_dir()

        from . import paths as _paths
        custom_results = self.params.get('directories', {}).get('results_dir', '').strip()

        if custom_out:
            self.active_source_dir = Path(custom_out)
        else:
            self.active_source_dir = Path(default_out)

        # Results are the user's own work and are kept apart from the databases:
        # an installed copy puts them in the user's documents, not beside the
        # program, so two people sharing one install do not overwrite each other.
        # output_dir still covers both when set, for configs that predate this.
        if custom_results:
            self.results_dir = Path(custom_results)
        elif custom_out:
            self.results_dir = Path(custom_out)
        else:
            self.results_dir = Path(_paths.default_user_results_dir(self.root))

        # Working files were written flat into one folder, so the network JSONs
        # a user might want to open sat among the SQLite databases the pipeline
        # depends on - which is how someone deletes or edits a database meaning
        # to touch a network. They are separated now.
        self.networks_dir = self.active_source_dir / 'networks'

        # THE DATABASES ARE ALWAYS THE USER'S, never the shipped corpus.
        #
        # "Use bundled reference data" switches the working folder to the corpus
        # inside the program, and that used to drag the databases with it - so
        # with the box ticked the pipeline looked for the master annotation
        # database inside a read-only folder that has never contained one, while
        # the 8 GB the user had actually built sat in their own data folder,
        # invisible. That is why the benchmark could not run against the
        # reference corpus: not a missing file, a wrong folder.
        #
        # Only the networks and the ground-truth files are bundled. Every
        # database is built on this machine, is expensive, and is shared between
        # projects, so all of them live together under the user's data folder
        # whatever the checkbox says.
        self.databases_dir = self._user_processed_dir() / 'databases'

        self.figures_dir = self.results_dir / 'figures'
        # Results likewise: figures, secondary analysis and benchmarking each
        # produce a handful of files, and mixed together in one folder none of
        # them is findable.
        self.secondary_dir = self.results_dir / 'secondary_analysis'
        self.benchmark_dir = self.results_dir / 'benchmark'
        # Logs are machinery, not results. A portable copy keeps everything
        # together; an installed one keeps them in the user's private area.
        self.log_dir = (self.results_dir / 'logs'
                        if _paths.is_portable(self.root) else _paths.log_dir())

        # databases_dir is always the user's, so it is always created. The
        # networks folder is not, when it would land inside the shipped corpus:
        # that is read-only, never written to, and copied into the installers,
        # so creating it there would only ship an empty directory to everyone.
        wanted = [self.active_raw_dir, self.active_source_dir, self.databases_dir,
                  self.figures_dir, self.secondary_dir, self.benchmark_dir,
                  self.log_dir]
        if not self.use_reference_data:
            wanted.insert(2, self.networks_dir)
        for d in wanted:
            try:
                d.mkdir(parents=True, exist_ok=True)
            except OSError:
                # Reference mode reads a corpus inside the installed program
                # folder, which may be read-only. Nothing is written there, so
                # a directory that cannot be created is not an error - the
                # lookup below falls back to the files that are already present.
                pass

    def _user_data_root(self):
        """The user's data folder, whatever the reference checkbox says."""
        from . import paths as _paths
        return Path(self.params.get('directories', {}).get('data_dir', '').strip()
                    or _paths.default_user_data_dir(self.root))

    def _user_raw_dir(self):
        """Where downloads and the master database went before this move."""
        custom = self.params.get('directories', {}).get('input_dir', '').strip()
        return Path(custom) if custom else self._user_data_root() / 'raw'

    def _user_processed_dir(self):
        """Where the project databases were written flat, before this move."""
        custom = self.params.get('directories', {}).get('output_dir', '').strip()
        return Path(custom) if custom else self._user_data_root() / 'processed'

    def _settled_db(self, name, *legacy_dirs):
        """Where a database lives: the databases folder, or where it already is.

        Same rule as _settled, but the fallbacks are given explicitly because a
        database can predate two different moves - the master annotation
        database was written to raw/ for most of this project's life, and the
        project databases were written flat into processed/. Neither is worth
        making anyone rebuild: the master one is a 44 GB download and several
        hours of compilation.
        """
        preferred = self.databases_dir / name
        if preferred.exists():
            return preferred
        for old in legacy_dirs:
            candidate = Path(old) / name
            if candidate.exists():
                return candidate
        return preferred

    def _settled(self, preferred_dir, name):
        """Where a working file lives: the new subfolder, or where it already is.

        The networks and databases used to be written flat into the processed
        folder. An existing project has hours of work sitting there - a cleaned
        citation database, a relevance database - and moving files under someone
        without asking is not something this should do on their behalf. So a
        file already present in the old flat location is used from there; only
        new ones are written to the subfolder. Nothing is orphaned and nothing
        is silently relocated.

        Reference mode benefits from the same rule: the bundled corpus is flat
        and read-only, and is found without special-casing it.
        """
        legacy = self.active_source_dir / name
        if not (preferred_dir / name).exists() and legacy.exists():
            return legacy
        return preferred_dir / name

    def _map_files(self):
        """Maps all explicit file paths required by the pipeline modules."""
        self.files = {
            "annotations": self._annotations_path(),
            # All four databases together, under the user's data folder. The
            # legacy locations are honoured where a file is already sitting in
            # one, so nobody rebuilds anything for this.
            "master_db": self._settled_db("master_mesh_database.db",
                                          self.active_raw_dir,
                                          self._user_raw_dir()),
            "pmids_db": self._settled_db(f"{self.db_prefix}_pmids.db",
                                         self.active_raw_dir,
                                         self._user_raw_dir()),

            "full_network": self._settled(self.networks_dir, f"{self.prefix}_full_network_data.json"),
            "glf_subgraph": self._settled(self.networks_dir, f"{self.prefix}_glf_optimal_subgraph.json"),
            "sa_subgraph": self._settled(self.networks_dir, f"{self.prefix}_sa_optimal_subgraph.json"),
            "consensus_lcc": self._settled(self.networks_dir, f"{self.prefix}_consensus_lcc_network.json"),
            "final_network": self._settled(self.networks_dir, f"{self.prefix}_final_network_with_relevance.json"),

            "cleaned_db": self._settled_db(f"{self.db_prefix}_cleaned_pmids.db",
                                           self.active_source_dir,
                                           self._user_processed_dir()),
            "relevance_db": self._settled_db(f"{self.db_prefix}_mean_relevancy.db",
                                             self.active_source_dir,
                                             self._user_processed_dir()),

            "failed_mesh": self.log_dir / f"{self.prefix}_failed_mesh_fetches.tsv",
            "empty_mesh": self.log_dir / f"{self.prefix}_empty_mesh_pmids.tsv",
            "log_file": self.log_dir / f"{self.prefix}_processing_errors.log",

            "mesh_marc": self.active_raw_dir / "20250301_marc_full2025.bin",
            "mesh_ascii": self.active_raw_dir / "d2025.bin",
            # The user's, like the databases. It is DERIVED - Step 1 writes it
            # from desc2025.xml - so pointing it at the shipped corpus meant a
            # reference run would try to write the MeSH vocabulary into a
            # read-only folder. It is also the same vocabulary either way: it
            # describes the MeSH release, not the corpus being analysed.
            "mesh_terms_csv": self._settled(self._user_processed_dir(), "mesh_terms.csv"),

            # The package's own module, resolved from this file - NOT
            # root/src/mesh_aop/. root is the working directory, so that path
            # pointed at wherever the process happened to be started: it never
            # existed, which made Step 1 rebuild the support files on every
            # single run, and it wrote the regenerated stop-word list to a stray
            # src/mesh_aop/ tree under the launch directory where nothing would
            # ever read it. The list actually in force is the one imported from
            # this package, so that is the file the pipeline must be talking
            # about.
            "mesh_stopwords_py": Path(__file__).resolve().parent / "mesh_stop_words.py"
        }

    # Credentials are read from the environment when the config leaves them blank, so
    # a working setup never requires writing a secret into a tracked file. The config
    # still wins if it is populated, which keeps existing local setups working.
    _CREDENTIAL_ENV = {
        "entrez_email": "MESH_ENTREZ_EMAIL",
        "entrez_api_key": "MESH_ENTREZ_API_KEY",
    }

    def get(self, section: str, key: str):
        """Utility to safely retrieve a parameter.

        Everything that consumes settings comes through here, which is why the
        TODAY substitution lives here rather than in the stored value: callers
        get a real date, the file keeps the word.
        """
        resolved = getattr(self, '_resolved', {})
        if (section, key) in resolved:
            return resolved[(section, key)]
        value = self.params.get(section, {}).get(key)
        value = self._unquote(section, key, value)
        if section == "credentials" and not value:
            return os.environ.get(self._CREDENTIAL_ENV.get(key, ""), "") or value
        return value

    # Fields whose value is matched against a NAME in the data - a MeSH heading,
    # a node, an edge. Quotes typed around one become part of the string and the
    # match then fails silently: the network holds a node called
    #     Dermatitis, Allergic Contact
    # and the lookup asks for
    #     "Dermatitis, Allergic Contact"
    # which is a different string. Nothing errors; the analysis simply comes
    # back empty, and the reason is invisible. Quoting a value with a comma in
    # it is a completely reasonable thing to try, so it is accepted and undone
    # rather than rejected.
    #
    # NOT applied to the search term: PubMed queries use quotes meaningfully,
    # and stripping them there would change what was searched for.
    _UNQUOTE = {
        ('benchmark', 'primary_node'),
        ('secondary_analysis', 'target_nodes'),
        ('secondary_analysis', 'target_edges'),
    }

    @staticmethod
    def _strip_quotes(text):
        """Remove one matched pair of surrounding quotes, if there is one."""
        s = text.strip()
        for q in ('"', "'"):
            if len(s) >= 2 and s.startswith(q) and s.endswith(q):
                return s[1:-1].strip()
        return s

    def _unquote(self, section, key, value):
        """Undo quotes typed around a name that has to match the data exactly."""
        if (section, key) not in self._UNQUOTE or not isinstance(value, str):
            return value
        # These two are semicolon-separated lists, so each entry is handled.
        if key in ('target_nodes', 'target_edges'):
            parts = [self._strip_quotes(p) for p in value.split(';')]
            return ';'.join(p for p in parts if p)
        return self._strip_quotes(value)

    def update(self, section: str, key: str, value):
        """Updates a parameter in memory."""
        if section not in self.params:
            self.params[section] = {}
        self.params[section][key] = value

    def save(self):
        """Writes current memory params back to the JSON file and refreshes paths."""
        # Per-run switches are held out of the file entirely; see load_config.
        persist = {k: v for k, v in self.params.items() if not k.startswith('_')}
        with open(self.config_path, 'w') as f:
            json.dump(persist, f, indent=4)
        print(f"    [+] Configuration saved successfully to {self.config_path.name}")
        self.refresh_paths()

    def reset(self):
        """Erases the local JSON file and reloads from immutable factory defaults."""
        if self.config_path.exists():
            os.remove(self.config_path)
            print(f"    [+] {self.config_path.name} erased.")
        self.load_config()
        self.refresh_paths()
        print("    [+] Factory defaults fully restored.")