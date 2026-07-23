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

    def __init__(self, config_path: str = "mesh_config.json", workspace_root: str = None):
        """
        Initializes the configuration object.

        Args:
            config_path: Path to the user's local mesh_config.json
            workspace_root: The root directory containing the data/ and results/ folders.
                            Defaults to the Current Working Directory if not provided.
        """
        self.config_path = Path(config_path)
        self.root = Path(workspace_root) if workspace_root else Path.cwd()

        self._define_factory_defaults()
        self.load_config()
        self.refresh_paths()

    def _define_factory_defaults(self):
        """Immutable baseline settings for the entire pipeline."""
        self.factory_defaults = {
            "control_flags": {
                "use_reference_data": True,
                "custom_file_prefix": "DAC_Mesh"
            },
            "directories": {
                "input_dir": "",
                "output_dir": ""
            },
            "credentials": {
                "entrez_email": "",
                "entrez_api_key": ""
            },
            "search_parameters": {
                "search_term": "Dermatitis, Allergic Contact [Mesh]",
                "start_date": "1950/01/01",
                "end_date": "2025/01/01",
                "generations_n": 1,
                "update_mesh_support_files": False
            },
            "analysis_parameters": {
                "calculate_full_centrality": True,
                "random_seed": 42,
                "betweenness_k_samples": 1000,
                "context_start_date": "1950/01/01",
                "context_end_date": "TODAY"
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
                "ground_truth_csv": "data/reference_processed/oecd_ground_truth_curated.xlsx",
                "negative_control_csv": "",
                "primary_node": "Dermatitis, Allergic Contact",
                "n_boot": 25,
                "n_perm": 25,
                "run_network_validation": True,
                "network_validation_weight_key": "CRS_pagerank_centrality",
                "min_articles_per_node": 2,
                "background_pool_size": 50000
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
                    if section in self.params and isinstance(values, dict):
                        self.params[section].update(values)
            except Exception as e:
                print(f"\n[!] Failed to parse {self.config_path}. Using factory defaults. ({e})")

    def refresh_paths(self):
        """Triggers path resolution if directories or control flags are updated."""
        self._resolve_dynamic_values()
        self._build_directories()
        self._map_files()

    def _resolve_dynamic_values(self):
        """Handles date interpretation and reference flags."""
        self.use_reference_data = self.params['control_flags']['use_reference_data']
        self.prefix = "DAC_Mesh" if self.use_reference_data else self.params['control_flags']['custom_file_prefix']
        if self.params['analysis_parameters']['context_end_date'] == "TODAY":
            self.params['analysis_parameters']['context_end_date'] = date.today().strftime("%Y/%m/%d")

    def get_default_directories(self, use_ref_data: bool):
        """Return the default (input, output) directories for the active data mode (user vs reference)."""
        base_data = self.root / 'data'
        default_in = base_data / 'reference_raw' if use_ref_data else base_data / 'raw'
        default_out = base_data / 'reference_processed' if use_ref_data else base_data / 'processed'
        return str(default_in), str(default_out)

    def _build_directories(self):
        """Constructs and ensures the existence of required directories with override support."""
        custom_in = self.params.get('directories', {}).get('input_dir', '').strip()
        custom_out = self.params.get('directories', {}).get('output_dir', '').strip()

        default_in, default_out = self.get_default_directories(self.use_reference_data)

        self.active_raw_dir = Path(custom_in) if custom_in else Path(default_in)

        if custom_out:
            self.active_source_dir = Path(custom_out)
            self.results_dir = Path(custom_out)
        else:
            self.active_source_dir = Path(default_out)
            self.results_dir = self.root / 'results'

        self.figures_dir = self.results_dir / 'figures'
        self.log_dir = self.results_dir / 'logs'

        self.active_raw_dir.mkdir(parents=True, exist_ok=True)
        self.active_source_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        self.log_dir.mkdir(parents=True, exist_ok=True)

    def _map_files(self):
        """Maps all explicit file paths required by the pipeline modules."""
        self.files = {
            "annotations": self.active_raw_dir / "aop_annotations_master.csv",
            "master_db": self.active_raw_dir / "master_mesh_database.db",
            "pmids_db": self.active_raw_dir / f"{self.prefix}_pmids.db",

            "full_network": self.active_source_dir / f"{self.prefix}_full_network_data.json",
            "glf_subgraph": self.active_source_dir / f"{self.prefix}_glf_optimal_subgraph.json",
            "sa_subgraph": self.active_source_dir / f"{self.prefix}_sa_optimal_subgraph.json",
            "consensus_lcc": self.active_source_dir / f"{self.prefix}_consensus_lcc_network.json",
            "final_network": self.active_source_dir / f"{self.prefix}_final_network_with_relevance.json",

            "cleaned_db": self.active_source_dir / f"{self.prefix}_cleaned_pmids.db",
            "relevance_db": self.active_source_dir / f"{self.prefix}_contextual_relevance.db",

            "failed_mesh": self.log_dir / f"{self.prefix}_failed_mesh_fetches.tsv",
            "empty_mesh": self.log_dir / f"{self.prefix}_empty_mesh_pmids.tsv",
            "log_file": self.log_dir / f"{self.prefix}_processing_errors.log",

            "mesh_marc": self.active_raw_dir / "20250301_marc_full2025.bin",
            "mesh_ascii": self.active_raw_dir / "d2025.bin",
            "mesh_terms_csv": self.active_source_dir / "mesh_terms.csv",

            "mesh_stopwords_py": self.root / "src" / "mesh_aop" / "mesh_stop_words.py"
        }

    def get(self, section: str, key: str):
        """Utility to safely retrieve a parameter."""
        return self.params.get(section, {}).get(key)

    def update(self, section: str, key: str, value):
        """Updates a parameter in memory."""
        if section not in self.params:
            self.params[section] = {}
        self.params[section][key] = value

    def save(self):
        """Writes current memory params back to the JSON file and refreshes paths."""
        with open(self.config_path, 'w') as f:
            json.dump(self.params, f, indent=4)
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