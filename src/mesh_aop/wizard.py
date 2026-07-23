# -*- coding: utf-8 -*-
"""
wizard.py - interactive configuration wizard (the `--interactive` flag).

Walks the user through the pipeline's settings block by block - control flags,
NCBI credentials, search parameters, analysis and network/simulation parameters,
and secondary-analysis options - showing current values and only prompting for
the blocks relevant to the chosen step.

It also probes the offline master database for presence, schema version,
completeness, and age, and offers to build, resume, or refresh it, then writes
the confirmed settings back to mesh_config.json (or triggers a factory reset).
"""

import os
import time
import sqlite3
from pathlib import Path

def _mask_credential(cred: str, is_email: bool = False) -> str:
    """Masks sensitive API keys and emails for console display."""
    if not cred:
        return ""

    if is_email and "@" in cred:
        parts = cred.split("@")
        return f"{parts[0][:2]}***@{parts[1]}"

    if len(cred) > 4:
        return f"********{cred[-4:]}"

    return "****"

def _prompt_override(label: str, current_value, cast_type=str, help_text="", is_email=False, must_exist=False):
    """
    Displays the current value and loops until a valid input or blank (skip) is provided.
    Enforces type casting and path existence validation.
    """
    if help_text:
        print(f"      {help_text}")

    while True:
        if cast_type == bool:
            curr_str = "True" if current_value else "False"
            user_in = input(f"  {label} (Current: {curr_str}) [y/n/Enter to skip]: ").strip().lower()

            if user_in == '':
                return current_value
            if user_in in ['y', 'yes', 'true', 't', '1']:
                return True
            if user_in in ['n', 'no', 'false', 'f', '0']:
                return False

            print("  [!] Invalid input. Please enter y/yes or n/no.")
            continue

        display_value = current_value
        if "API Key" in label or "Email" in label:
            display_value = _mask_credential(current_value, is_email)

        user_in = input(f"  {label} (Current: {display_value}) [Enter to skip]: ").strip()

        if user_in == '':
            return current_value

        try:
            val = cast_type(user_in)

            if is_email and "@" not in str(val):
                print("  [!] Invalid email format. Must contain '@'.")
                continue

            if must_exist and not os.path.exists(str(val)):
                print(f"  [!] Path does not exist: {val}. Please provide a valid path.")
                continue

            return val

        except ValueError:
            print(f"  [!] Invalid input type. Expected numerical or text value.")

def _ask_block(title: str, preview_dict: dict) -> bool:
    """Displays the current state of a parameter block and asks if the user wants to update it."""
    print(f"\n<<< {title} >>>")

    for k, v in preview_dict.items():
        print(f"  - {k}: {v}")

    ans = input(f"Update {title}? [y/n/Enter to skip]: ").strip().lower()
    return ans in ['y', 'yes']

def _open_master_for_health_check(master_db_path):
    """Open the master DB for the status check and return (connection, columns).

    A read-only open is tried first, but a WAL-mode database cannot be brought
    online read-only on a OneDrive path (the -shm shared-memory file cannot be
    created), so it would report an empty schema and be misread as corrupt. If
    the read-only open finds no schema, fall back to a normal connection - which
    can apply the WAL even on OneDrive - so a healthy database is not flagged
    corrupt purely because of its journal mode.
    """
    conn = sqlite3.connect(f"file:{master_db_path}?mode=ro", uri=True, timeout=120.0)
    columns = [info[1] for info in conn.execute("PRAGMA table_info(master_mesh_annotations)").fetchall()]
    if columns:
        return conn, columns
    conn.close()
    conn = sqlite3.connect(str(master_db_path), timeout=120.0)
    columns = [info[1] for info in conn.execute("PRAGMA table_info(master_mesh_annotations)").fetchall()]
    return conn, columns

def _classify_master_db(master_db_path, downloaded_files):
    """Inspect the master DB and classify it without ever deleting it.

    Returns (status, is_corrupt, is_incomplete, is_locked, pending_count). A
    locked/busy or otherwise unreadable database is reported as locked - never
    corrupt - so the caller never offers to delete a healthy-but-inaccessible DB.
    """
    if not (master_db_path and os.path.exists(master_db_path)):
        return "Missing", False, False, False, 0
    try:
        conn, columns = _open_master_for_health_check(master_db_path)
        try:
            if not columns:
                return "EMPTY / CORRUPTED", True, False, False, 0
            if 'pub_date' not in columns:
                return "OUTDATED SCHEMA", True, False, False, 0
            record_count = conn.execute("SELECT count(*) FROM master_mesh_annotations").fetchone()[0]
            try:
                parsed_count = conn.execute("SELECT count(*) FROM parsed_files").fetchone()[0]
            except sqlite3.OperationalError:
                parsed_count = 0
        finally:
            conn.close()

        if downloaded_files > 0 and parsed_count < downloaded_files:
            return (f"INCOMPLETE ({parsed_count}/{downloaded_files} files parsed, {record_count:,} records)",
                    False, True, False, downloaded_files - parsed_count)
        return f"Valid & Complete ({record_count:,} records)", False, False, False, 0

    except sqlite3.OperationalError as e:
        if any(s in str(e).lower() for s in ("locked", "busy", "in use")):
            return "IN USE / LOCKED", False, False, True, 0
        return "CORRUPTED", True, False, False, 0
    except sqlite3.DatabaseError:
        return "CORRUPTED", True, False, False, 0
    except Exception as e:
        return f"UNREADABLE ({type(e).__name__})", False, False, True, 0

def _attempt_unlock_master_db(master_db_path):
    """Best-effort recovery of a database reported locked/unreadable.

    Safe to call only after other pipeline runs have been stopped. Clears the
    stale -shm shared-memory index, lets SQLite recover any hot journal and
    checkpoint a leftover WAL on a clean open/close, and leaves the file in
    rollback mode. Returns True if the database can afterwards be opened
    read-only. A lock held by a still-running process cannot be broken and
    returns False.
    """
    master_db_path = str(master_db_path)
    shm = master_db_path + "-shm"
    if os.path.exists(shm):
        try:
            os.remove(shm)
        except OSError:
            pass
    try:
        conn = sqlite3.connect(master_db_path, timeout=30.0)
        conn.execute("PRAGMA wal_checkpoint(TRUNCATE);")
        conn.execute("PRAGMA journal_mode=DELETE;")
        conn.commit()
        conn.close()
    except sqlite3.DatabaseError:
        return False
    try:
        c = sqlite3.connect(f"file:{master_db_path}?mode=ro", uri=True, timeout=10.0)
        ok = bool(c.execute("PRAGMA table_info(master_mesh_annotations)").fetchall())
        c.close()
        return ok
    except sqlite3.DatabaseError:
        return False

def run_interactive_wizard(config, step: str) -> bool:
    """Executes block-by-block configuration wizard. Returns True if a factory reset occurred."""
    print("\n" + "<"*30 + ">"*30)
    print("<<< MeSH AOP Pipeline - Interactive Configuration Wizard >>>")
    print("<"*30 + ">"*30)

    # ---------------------------------------------------------
    # SYSTEM OPTIONS: DOUBLE-CONFIRM RESET
    # ---------------------------------------------------------
    print("\n<<< System Options >>>")
    do_reset = input("Do you want to restore Factory Defaults? (Erases API keys & custom settings) [y/n/Enter to skip]: ").strip().lower()

    if do_reset in ['y', 'yes']:
        confirm = input("  [!] Are you SURE you want to erase all settings? This cannot be undone. [y/n]: ").strip().lower()
        if confirm in ['y', 'yes']:
            print("\n  >>> Factory Defaults will be restored. <<<")
            return True
        else:
            print("  >>> Factory Reset aborted. Proceeding with current configuration. <<<")

    print("\n  Press [Enter] on any prompt to keep the (Current) value.")

    params = config.params

    # ---------------------------------------------------------
    # BLOCK 1: Directories and Control Flags (Always Ask)
    # ---------------------------------------------------------
    use_ref = params['control_flags'].get('use_reference_data', True)
    pause_anno = params['control_flags'].get('pause_for_annotation', False)
    def_in, def_out = config.get_default_directories(use_ref)

    curr_in = params['directories'].get('input_dir', '').strip() or def_in
    curr_out = params['directories'].get('output_dir', '').strip() or def_out

    b1_preview = {
        "Use Reference Data": use_ref,
        "Pause for Annotation": pause_anno,
        "Custom Prefix": params['control_flags'].get('custom_file_prefix', ''),
        "Input Directory": curr_in,
        "Output Directory": curr_out
    }

    if _ask_block("Control Flags & Directories", b1_preview):
        params['control_flags']['use_reference_data'] = _prompt_override(
            "Use Reference Data?",
            use_ref,
            bool,
            "If True, skips Entrez downloading and uses curated reference library."
        )

        params['control_flags']['pause_for_annotation'] = _prompt_override(
            "Pause for Annotation (Occurs after Step 3)?",
            pause_anno,
            bool,
            "Set to True if you want to halt the pipeline to manually categorize nodes.\n      Default is False (AFK Mode: assigns 'Unassigned' and finishes pipeline)."
        )

        new_use_ref = params['control_flags']['use_reference_data']
        new_def_in, new_def_out = config.get_default_directories(new_use_ref)
        display_in = params['directories'].get('input_dir', '').strip() or new_def_in
        display_out = params['directories'].get('output_dir', '').strip() or new_def_out

        params['control_flags']['custom_file_prefix'] = _prompt_override(
            "Custom File Prefix",
            params['control_flags'].get('custom_file_prefix', ''),
            str,
            "'DAC_Mesh' is default."
        )

        user_in = _prompt_override(
            "Input Directory",
            display_in,
            str,
            "Must exist.",
            must_exist=True
        )
        params['directories']['input_dir'] = "" if user_in == new_def_in else user_in

        user_out = _prompt_override(
            "Output Directory",
            display_out,
            str,
            "Will be created."
        )
        params['directories']['output_dir'] = "" if user_out == new_def_out else user_out

    config.refresh_paths()

    # ---------------------------------------------------------
    # BLOCK 1.5: Master Database Status, Health & Age Check
    # ---------------------------------------------------------
    master_db_path = config.files.get('master_db')
    print("\n<<< Master Database Status >>>")

    baseline_dir = config.active_raw_dir / "pubmed_baseline"
    updates_dir = config.active_raw_dir / "pubmed_updates"
    downloaded_files = 0

    if baseline_dir.exists():
        downloaded_files += len(list(baseline_dir.glob("*.xml.gz")))
    if updates_dir.exists():
        downloaded_files += len(list(updates_dir.glob("*.xml.gz")))

    db_age_days = 0
    if master_db_path and os.path.exists(master_db_path):
        db_age_days = (time.time() - os.path.getmtime(master_db_path)) / 86400

    db_status, is_corrupt, is_incomplete, is_locked, pending_count = \
        _classify_master_db(master_db_path, downloaded_files)

    if is_locked:
        print(f"  [!] Master DB is {db_status}: {master_db_path.name if master_db_path else 'Unknown'}")
        print("      The database exists and was NOT modified or deleted, but it could not")
        print("      be read. Usually another pipeline run is holding it open, OneDrive is")
        print("      mid-sync, or a stale lock was left by a process that died (common on")
        print("      Linux). First stop every other instance:")
        print("          Get-Process python | Stop-Process -Force      (Windows)")
        print("          pkill -f mesh_aop                              (Linux/macOS)")
        ans = input("  Then attempt to unlock the database now? [y/n/Enter to skip]: ").strip().lower()
        if ans in ['y', 'yes']:
            if _attempt_unlock_master_db(master_db_path):
                print("  [+] Unlock succeeded. Re-checking database status...")
                db_status, is_corrupt, is_incomplete, is_locked, pending_count = \
                    _classify_master_db(master_db_path, downloaded_files)
                print(f"      Status now: {db_status}")
            else:
                print("  [!] Still locked/unreadable - a live process is holding it, or the")
                print("      machine needs a reboot to clear the lock. The DB was NOT modified.")

    if is_locked:
        params['_run_baseline_etl'] = False
        params['_delete_corrupt_db'] = False

    elif db_status == "Missing" or is_corrupt:
        print(f"  [!] Master DB is {db_status}: {master_db_path.name if master_db_path else 'Unknown'}")

        if is_corrupt:
            if "OUTDATED SCHEMA" in db_status:
                print("      The database schema is outdated (missing the pub_date column).")
                print("      It MUST be deleted and rebuilt to support high-speed Contextual Relevance Scoring.")
            else:
                print("      The database file is malformed (likely due to an interrupted build).")
                print("      It must be deleted and rebuilt from your downloaded XML files.")

        ans = input("  Compile the PubMed Baseline Database now? [y/n/Enter to skip]: ").strip().lower()

        if ans in ['y', 'yes']:
            params['_run_baseline_etl'] = True
            params['_delete_corrupt_db'] = is_corrupt
            ans_up = input("  Download Daily Updates as well? [y/n/Enter to skip]: ").strip().lower()
            params['_run_baseline_updates'] = ans_up in ['y', 'yes']
            params['_skip_baseline_download'] = False
        else:
            params['_run_baseline_etl'] = False
            print("\n  [!] CRITICAL WARNING: You chose to skip Master Database compilation.")
            print("      Step 3 (Contextual Relevance Scoring) will fail or return empty results without a valid database.")

    elif is_incomplete:
        print(f"  [!] Master DB is {db_status}")
        print(f"      The parsing process was interrupted. {pending_count} local files still need to be compiled.")
        ans = input("  Resume compiling the database now? [y/n/Enter to skip]: ").strip().lower()

        if ans in ['y', 'yes']:
            params['_run_baseline_etl'] = True
            params['_delete_corrupt_db'] = False
            params['_skip_baseline_download'] = True

            ans_up = input("  Download and include new Daily Updates as well? [y/n/Enter to skip]: ").strip().lower()
            params['_run_baseline_updates'] = ans_up in ['y', 'yes']
        else:
            params['_run_baseline_etl'] = False
            print("\n  [!] CRITICAL WARNING: You chose to run with an incomplete Master Database.")
            print("      Contextual Relevance Scores will be statistically inaccurate until fully compiled.")

    else:
        if db_age_days > 330:
            print(f"  [!] Master DB Found: {db_status}")
            print(f"      WARNING: Database is {int(db_age_days)} days old. A new yearly NLM Baseline is likely available.")
            ans = input("  Rebuild entirely with the newest Yearly Baseline (~40GB)? [y/n/Enter to skip]: ").strip().lower()

            if ans in ['y', 'yes']:
                params['_run_baseline_etl'] = True
                params['_delete_corrupt_db'] = True
                params['_skip_baseline_download'] = False
                ans_up = input("  Download Daily Updates as well? [y/n/Enter to skip]: ").strip().lower()
                params['_run_baseline_updates'] = ans_up in ['y', 'yes']
            else:
                ans_up = input("  Fetch and compile just the new daily updates instead? [y/n/Enter to skip]: ").strip().lower()
                if ans_up in ['y', 'yes']:
                    params['_run_baseline_etl'] = True
                    params['_run_baseline_updates'] = True
                    params['_skip_baseline_download'] = True
                else:
                    params['_run_baseline_etl'] = False
        else:
            print(f"  [+] Master DB Found: {db_status}")
            ans = input("  Fetch and compile new daily updates? [y/n/Enter to skip]: ").strip().lower()

            if ans in ['y', 'yes']:
                params['_run_baseline_etl'] = True
                params['_run_baseline_updates'] = True
                params['_skip_baseline_download'] = True
            else:
                params['_run_baseline_etl'] = False

    # ---------------------------------------------------------
    # BLOCK 2: Credentials
    # ---------------------------------------------------------
    if step not in ['all', 'data_ops', 'network', 'secondary']:
        print("\n<<< NCBI Credentials >>>")
        print(f"  - [Skipped]: Not required for step '{step}'.")
    else:
        b2_preview = {
            "Email": _mask_credential(params['credentials'].get('entrez_email', ''), is_email=True) or "[Empty]",
            "API Key": _mask_credential(params['credentials'].get('entrez_api_key', '')) or "[Empty]"
        }

        if _ask_block("NCBI Credentials", b2_preview):
            params['credentials']['entrez_email'] = _prompt_override(
                "Entrez Email",
                params['credentials'].get('entrez_email', ''),
                str,
                is_email=True
            )
            params['credentials']['entrez_api_key'] = _prompt_override(
                "Entrez API Key",
                params['credentials'].get('entrez_api_key', '')
            )

    # ---------------------------------------------------------
    # BLOCK 3: Search Parameters
    # ---------------------------------------------------------
    if params['control_flags']['use_reference_data']:
        print("\n<<< Search Parameters >>>")
        print("  - [Skipped]: 'Use Reference Data' is True. Using curated reference library.")
    elif step not in ['all', 'data_ops', 'process']:
        print("\n<<< Search Parameters >>>")
        print(f"  - [Skipped]: Not required for step '{step}'.")
    else:
        b3_preview = {
            "Search Term": params['search_parameters'].get('search_term', ''),
            "Start Date": params['search_parameters'].get('start_date', ''),
            "End Date": params['search_parameters'].get('end_date', ''),
            "Generations": params['search_parameters'].get('generations_n', 1),
            "Update MeSH Files": params['search_parameters'].get('update_mesh_support_files', False)
        }

        if _ask_block("Search Parameters", b3_preview):
            params['search_parameters']['search_term'] = _prompt_override(
                "MeSH Search Term",
                b3_preview["Search Term"]
            )
            params['search_parameters']['start_date'] = _prompt_override(
                "Start Date (YYYY/MM/DD)",
                b3_preview["Start Date"]
            )
            params['search_parameters']['end_date'] = _prompt_override(
                "End Date (YYYY/MM/DD)",
                b3_preview["End Date"]
            )
            params['search_parameters']['generations_n'] = _prompt_override(
                "Citation Generations",
                b3_preview["Generations"],
                int,
                "Incoming/outgoing citations per layer. High values exponentially increase network size."
            )
            params['search_parameters']['update_mesh_support_files'] = _prompt_override(
                "Force Update MeSH Support Files?",
                b3_preview["Update MeSH Files"],
                bool
            )

    # ---------------------------------------------------------
    # BLOCK 4: Analysis Parameters
    # ---------------------------------------------------------
    if step not in ['all', 'network', 'viz']:
        print("\n<<< Analysis Parameters >>>")
        print(f"  - [Skipped]: Not required for step '{step}'.")
    else:
        b4_preview = {
            "Calculate Exact Betweenness": params['analysis_parameters'].get('calculate_full_centrality', False),
            "Random Seed": params['analysis_parameters'].get('random_seed', 42),
            "Betweenness K-Samples": params['analysis_parameters'].get('betweenness_k_samples', 1000),
            "Eigenvector Max Iterations": params['analysis_parameters'].get('eigenvector_max_iter', 1000),
            "Eigenvector Tolerance": params['analysis_parameters'].get('eigenvector_tol', 1.0e-6),
            "Context Start Date": params['analysis_parameters'].get('context_start_date', '1950/01/01'),
            "Context End Date": params['analysis_parameters'].get('context_end_date', '2026/05/25')
        }

        if _ask_block("Analysis Parameters", b4_preview):
            params['analysis_parameters']['calculate_full_centrality'] = _prompt_override(
                "Calculate Exact Betweenness?", b4_preview["Calculate Exact Betweenness"], bool
            )
            params['analysis_parameters']['random_seed'] = _prompt_override(
                "Random Seed", b4_preview["Random Seed"], int
            )
            params['analysis_parameters']['betweenness_k_samples'] = _prompt_override(
                "Betweenness K-Samples", b4_preview["Betweenness K-Samples"], int
            )
            params['analysis_parameters']['eigenvector_max_iter'] = _prompt_override(
                "Eigenvector Max Iterations", b4_preview["Eigenvector Max Iterations"], int
            )
            params['analysis_parameters']['eigenvector_tol'] = _prompt_override(
                "Eigenvector Tolerance", b4_preview["Eigenvector Tolerance"], float,
                "E.g., 1.0e-6 for strict precision, or 1.0e-4 for faster processing."
            )
            params['analysis_parameters']['context_start_date'] = _prompt_override(
                "Context Start Date", b4_preview["Context Start Date"], str,
                "Scopes which articles are relevance-scored (ARS). Term specificity (IC) is always computed across the full corpus, independent of this window."
            )
            params['analysis_parameters']['context_end_date'] = _prompt_override(
                "Context End Date", b4_preview["Context End Date"], str
            )

    # ---------------------------------------------------------
    # BLOCK 5: Network Parameters
    # ---------------------------------------------------------
    if step not in ['all', 'network']:
        print("\n<<< Network Parameters >>>")
        print(f"  - [Skipped]: Not required for step '{step}'.")
    else:
        nw_factors = params['network_parameters'].get('node_weight_factors', {})

        b5_preview = {
            "Lambda Value": params['network_parameters'].get('lambda_val', 1.0),
            "Weight: Centrality": nw_factors.get('centrality', 0.45),
            "Weight: Article Rank": nw_factors.get('article_rank', 0.15),
            "Weight: Rank Median Cit": nw_factors.get('rank_median_cit', 0.20),
            "Weight: Rank Total Cit": nw_factors.get('rank_total_cit', 0.20)
        }

        if _ask_block("Network Parameters", b5_preview):
            params['network_parameters']['lambda_val'] = _prompt_override(
                "Lambda Value",
                b5_preview["Lambda Value"],
                float
            )

            print("      Note: Node weight factors below will be automatically normalized out of 1.0 downstream.")

            nw_factors['centrality'] = _prompt_override(
                "Centrality Weight",
                b5_preview["Weight: Centrality"],
                float
            )
            nw_factors['article_rank'] = _prompt_override(
                "Article Rank Weight",
                b5_preview["Weight: Article Rank"],
                float
            )
            nw_factors['rank_median_cit'] = _prompt_override(
                "Rank Median Cit Weight",
                b5_preview["Weight: Rank Median Cit"],
                float
            )
            nw_factors['rank_total_cit'] = _prompt_override(
                "Rank Total Cit Weight",
                b5_preview["Weight: Rank Total Cit"],
                float
            )
            params['network_parameters']['node_weight_factors'] = nw_factors

    # ---------------------------------------------------------
    # BLOCK 6: Simulation Parameters
    # ---------------------------------------------------------
    if step not in ['all', 'network', 'viz']:
        print("\n<<< Simulation Parameters >>>")
        print(f"  - [Skipped]: Not required for step '{step}'.")
    else:
        b6_preview = {
            "Target Edges": params['simulation_parameters'].get('target_num_edges', 500),
            "GLF Iterations": params['simulation_parameters'].get('glf_iterations', 5000000),
            "SA Iterations": params['simulation_parameters'].get('sa_iterations', 5000000),
            "SA Start Temp": params['simulation_parameters'].get('sa_temp_start', 5000.0),
            "SA Cooling Rate": params['simulation_parameters'].get('sa_cooling_rate', 0.999995)
        }

        if _ask_block("Simulation Parameters", b6_preview):
            params['simulation_parameters']['target_num_edges'] = _prompt_override(
                "Target Edges",
                b6_preview["Target Edges"],
                int
            )
            params['simulation_parameters']['glf_iterations'] = _prompt_override(
                "GLF Iterations",
                b6_preview["GLF Iterations"],
                int
            )
            params['simulation_parameters']['sa_iterations'] = _prompt_override(
                "SA Iterations",
                b6_preview["SA Iterations"],
                int
            )
            params['simulation_parameters']['sa_temp_start'] = _prompt_override(
                "SA Start Temp",
                b6_preview["SA Start Temp"],
                float
            )
            params['simulation_parameters']['sa_cooling_rate'] = _prompt_override(
                "SA Cooling Rate",
                b6_preview["SA Cooling Rate"],
                float
            )

# ---------------------------------------------------------
    # BLOCK 7: Secondary Analysis Parameters
    # ---------------------------------------------------------
    if 'secondary_analysis' not in params:
        params['secondary_analysis'] = {}

    if step not in ['all', 'secondary']:
        print("\n<<< Secondary Analysis Parameters >>>")
        print(f"  - [Skipped]: Not required for step '{step}'.")
    else:
        sec_params = params['secondary_analysis']

        if step == 'all':
            b7_preview = {
                "Export Limit": sec_params.get('export_limit', 500),
                "Sort Metric": sec_params.get('sort_metric', 'F1'),
                "Linear ARS Weight": sec_params.get('linear_weight_ars', 0.5),
                "Export Top Network Articles": sec_params.get('export_top_articles', True),
                "Exclude Review Articles": sec_params.get('exclude_reviews', True)
            }
            if _ask_block("Secondary Analysis Parameters (Post-Pipeline Export)", b7_preview):
                sec_params['export_limit'] = _prompt_override(
                    "Export Limit (Number of articles)", b7_preview["Export Limit"], int
                )
                sec_params['sort_metric'] = _prompt_override(
                    "Sort Metric ('Linear' or 'F1')", b7_preview["Sort Metric"], str,
                    "Linear = Balanced average. F1 = Strict harmonic mean (penalizes articles if EITHER score is low)."
                )
                sec_params['linear_weight_ars'] = _prompt_override(
                    "Linear ARS Weight (Must be between 0.01 and 0.99)", b7_preview["Linear ARS Weight"], float,
                    "Citation weight is automatically (1.0 - ARS). Example: 0.70 means 70% ARS, 30% Citations."
                )
                sec_params['export_top_articles'] = _prompt_override(
                    "Export top overall articles across entire network at end of pipeline?",
                    b7_preview["Export Top Network Articles"], bool
                )
                if sec_params['export_top_articles']:
                    sec_params['exclude_reviews'] = _prompt_override(
                        "Exclude Review Articles from exports?",
                        b7_preview["Exclude Review Articles"], bool
                    )
            sec_params['run_secondary_analysis'] = sec_params.get('export_top_articles', True)
            sec_params['target_nodes'] = ''
            sec_params['target_edges'] = ''

        elif step == 'secondary':
            b7_preview = {
                "Export Limit": sec_params.get('export_limit', 500),
                "Sort Metric": sec_params.get('sort_metric', 'F1'),
                "Linear ARS Weight": sec_params.get('linear_weight_ars', 0.5),
                "Exclude Review Articles": sec_params.get('exclude_reviews', True),
                "Export Top Network Articles": sec_params.get('export_top_articles', True),
                "Target Nodes": sec_params.get('target_nodes', ''),
                "Target Edges": sec_params.get('target_edges', '')
            }
            if _ask_block("Secondary Analysis Parameters (Targeted)", b7_preview):
                sec_params['run_secondary_analysis'] = True
                sec_params['export_limit'] = _prompt_override(
                    "Export Limit (Number of articles)", b7_preview["Export Limit"], int
                )
                sec_params['sort_metric'] = _prompt_override(
                    "Sort Metric ('Linear' or 'F1')", b7_preview["Sort Metric"], str,
                    "Linear = Balanced average. F1 = Strict harmonic mean (penalizes articles if EITHER score is low)."
                )
                sec_params['linear_weight_ars'] = _prompt_override(
                    "Linear ARS Weight (Must be between 0.01 and 0.99)", b7_preview["Linear ARS Weight"], float,
                    "Citation weight is automatically (1.0 - ARS). Example: 0.70 means 70% ARS, 30% Citations."
                )
                sec_params['exclude_reviews'] = _prompt_override(
                    "Exclude Review Articles from exports?", b7_preview["Exclude Review Articles"], bool
                )
                sec_params['export_top_articles'] = _prompt_override(
                    "Export top overall articles across entire network?", b7_preview["Export Top Network Articles"], bool
                )
                sec_params['target_nodes'] = _prompt_override(
                    "Analyze Specific Nodes (Semicolon-separated, e.g., Term A; Term B | Empty to skip)",
                    b7_preview["Target Nodes"], str
                )
                sec_params['target_edges'] = _prompt_override(
                    "Analyze Specific Edges (Format: NodeA - NodeB; NodeC - NodeD | Empty to skip)",
                    b7_preview["Target Edges"], str
                )

    # ---------------------------------------------------------
    # BLOCK 8: Ground-Truth Validation & Benchmarking
    # ---------------------------------------------------------
    if 'benchmark' not in params:
        params['benchmark'] = {}

    if step not in ['all', 'benchmark']:
        print("\n<<< Benchmark Parameters >>>")
        print(f"  - [Skipped]: Not required for step '{step}'.")
    else:
        bench_params = params['benchmark']
        # The bundled ground truth describes the reference corpus, so it is only
        # meaningful when that corpus is in play; default the analysis on/off to match.
        _use_ref = bool(params.get('control_flags', {}).get('use_reference_data', False))
        b8_preview = {
            "Run Ground-Truth Analysis": bench_params.get('run_ground_truth_analysis', _use_ref),
            "Ground Truth File": bench_params.get(
                'ground_truth_csv', 'data/reference_processed/oecd_ground_truth_curated.xlsx'
            ),
            "Primary Node": bench_params.get('primary_node', 'Dermatitis, Allergic Contact'),
            "Bootstrap Resamples": bench_params.get('n_boot', 25),
            "Permutations": bench_params.get('n_perm', 25)
        }
        if _ask_block("Ground-Truth Validation & Benchmark Parameters", b8_preview):
            bench_params['run_ground_truth_analysis'] = _prompt_override(
                "Run ground-truth analysis?", b8_preview["Run Ground-Truth Analysis"], bool,
                "Defaults to on when using reference data, off otherwise."
            )
            bench_params['ground_truth_csv'] = _prompt_override(
                "Ground Truth File (.csv/.xlsx; filename or path)",
                b8_preview["Ground Truth File"], str,
                "The curated positives retrieval is scored against. Leave as-is for the bundled OECD AOP40 set."
            )
            bench_params['primary_node'] = _prompt_override(
                "Primary Disease Node", b8_preview["Primary Node"], str,
                "Splits hits into 'under the primary node' vs 'topology-exclusive' (found without it)."
            )
            bench_params['n_boot'] = _prompt_override(
                "Bootstrap Resamples (n_boot)", b8_preview["Bootstrap Resamples"], int,
                "Warning: n = 200 takes approximately 2 hours."
            )
            bench_params['n_perm'] = _prompt_override(
                "Permutation Null Iterations (n_perm)", b8_preview["Permutations"], int
            )

    print("\n<<< Configuration Update Complete >>>\n")
    return False