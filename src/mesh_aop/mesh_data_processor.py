# -*- coding: utf-8 -*-
"""
mesh_data_processor.py - MeSH vocabulary extraction and stop-word generation (pipeline Step 1).

Parses the NLM MeSH descriptor XML (e.g. desc2025.xml) in a single,
memory-efficient streaming pass and produces the support files the rest of the
pipeline depends on.

Specifically it:
  - extracts every descriptor's name, unique identifier, and tree numbers;
  - derives a stop-word list by excluding terms that fall outside the relevant
    MeSH categories (Anatomy, Diseases, Chemicals/Drugs, Phenomena) and writes
    it as an importable Python module;
  - exports a CSV of unique terms, each tagged with a stop-word flag.
"""

import os
import re
import glob
import datetime
import traceback
import xml.etree.ElementTree as ET
from collections import defaultdict
import pandas as pd
import requests
from tqdm import tqdm
from pathlib import Path

from . import vocabulary

# <<< Constants >>>
# The trees and their names live in vocabulary.py, which the settings window
# also reads and which imports nothing heavier than csv. These are the defaults
# only: a run passes its own choice in, so that four trees is what the program
# starts with rather than what it insists on.
CATEGORIES_TO_KEEP = set(vocabulary.DEFAULT_KEPT)
CHECK_TAGS_TO_STOP = set(vocabulary.SEXES)
TOP_LEVEL_CATEGORY_NAMES = dict(vocabulary.TREE_NAMES)

def _get_missing_nlm_file_message(expected_path: str) -> str:
    """Generates a detailed error message for missing NLM XML data."""
    return (
        f"Raw MeSH XML File Not Found at: {expected_path}\n"
        "To resolve this:\n"
        "1. Go to the NLM Data Distribution page: https://nlmpubs.nlm.nih.gov/projects/mesh/2025/xmlmesh/\n"
        "2. Download the 'desc2025.xml' file.\n"
        "3. Place it in your raw data directory and ensure the configuration matches."
    )

# <<< MeSH descriptor XML acquisition (auto-download if missing or superseded) >>>

MESH_XML_URL_TEMPLATE = "https://nlmpubs.nlm.nih.gov/projects/mesh/{year}/xmlmesh/desc{year}.xml"


def _probe_remote_mesh_xml(year: int):
    """Return (url, content_length) if a genuine desc<year>.xml is published at NLM, else None.

    NLM answers requests for unpublished years with HTTP 200 but redirects to an
    HTML 'bad_url' page, so a bare status check is not enough - we confirm the
    response is really the XML file (xml content-type, no bad_url redirect).
    """
    url = MESH_XML_URL_TEMPLATE.format(year=year)
    try:
        resp = requests.get(url, stream=True, timeout=30, allow_redirects=True)
        is_real = (resp.status_code == 200
                   and 'xml' in resp.headers.get('Content-Type', '').lower()
                   and 'bad_url' not in resp.url)
        length = int(resp.headers.get('Content-Length', 0) or 0)
        resp.close()
        return (url, length) if is_real else None
    except requests.RequestException:
        return None


def _newest_remote_mesh_year():
    """Probe NLM for the newest published descriptor year.

    Returns (year, url, content_length), or None when NLM cannot be reached (the
    caller then falls back to a local copy). Checks next year first so a freshly
    released annual file is adopted as soon as it appears.
    """
    this_year = datetime.date.today().year
    for year in range(this_year + 1, this_year - 3, -1):
        found = _probe_remote_mesh_xml(year)
        if found:
            return (year, found[0], found[1])
    return None


def _newest_local_mesh_xml(raw_dir: str):
    """Return the path to the newest desc<year>.xml already present in raw_dir, or None."""
    candidates = []
    for path in glob.glob(os.path.join(raw_dir, "desc*.xml")):
        match = re.search(r'desc(\d{4})\.xml$', os.path.basename(path))
        if match:
            candidates.append((int(match.group(1)), path))
    return max(candidates)[1] if candidates else None


def _download_mesh_xml(url: str, target_path: str, expected_size: int = 0):
    """Stream a large descriptor XML to target_path with a progress bar and atomic rename."""
    tmp_path = target_path + ".part"
    with requests.get(url, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        total = int(resp.headers.get('Content-Length', expected_size) or 0) or None
        with open(tmp_path, 'wb') as handle, tqdm(
                total=total, unit='B', unit_scale=True, unit_divisor=1024,
                desc=os.path.basename(target_path)) as bar:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                if chunk:
                    handle.write(chunk)
                    bar.update(len(chunk))
    os.replace(tmp_path, target_path)


def ensure_mesh_descriptor_xml(raw_dir: str) -> str:
    """Guarantee a current MeSH descriptor XML exists in raw_dir and return its path.

    Downloads the newest annual desc<year>.xml from NLM when it is missing, only
    partially present, or superseded by a newer release - no user action needed.
    When NLM is unreachable it falls back to the newest local copy, and only
    raises when no usable file exists at all.
    """
    os.makedirs(raw_dir, exist_ok=True)
    newest = _newest_remote_mesh_year()

    if newest is None:
        local = _newest_local_mesh_xml(raw_dir)
        if local:
            print(f"    [!] Could not reach NLM to check for MeSH updates; "
                  f"using existing {os.path.basename(local)}.")
            return local
        raise FileNotFoundError(
            "MeSH descriptor XML is missing and NLM could not be reached to download it.\n"
            "Connect to the internet and re-run, or manually place a descYYYY.xml from\n"
            "https://nlmpubs.nlm.nih.gov/projects/mesh/ into: " + raw_dir)

    year, url, size = newest
    target = os.path.join(raw_dir, f"desc{year}.xml")

    if os.path.exists(target) and (size == 0 or os.path.getsize(target) == size):
        print(f"    [+] MeSH descriptor XML is current: {os.path.basename(target)}")
        return target

    action = "Updating to" if _newest_local_mesh_xml(raw_dir) else "Downloading"
    approx = f" (~{size / 1e6:.0f} MB)" if size else ""
    print(f"    [>] {action} newest MeSH descriptor XML desc{year}.xml{approx} from NLM...")
    _download_mesh_xml(url, target, size)
    print(f"    [+] Saved {os.path.basename(target)} to {raw_dir}")
    return target


def write_stopwords_to_python(output_path: str, grouped_stopwords: dict, missing_tns: list) -> bool:
    """Writes the categorized stop words to a properly formatted Python module."""
    cats_display_list = []

    for cat in sorted(list(CATEGORIES_TO_KEEP)):
        name = TOP_LEVEL_CATEGORY_NAMES.get(cat, "Unknown")
        cats_display_list.append(f"'{cat} - {name}'")

    cats_formatted_str = "[" + ", ".join(cats_display_list) + "]"

    try:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        with open(output_path, "w", encoding="utf-8") as f:
            f.write('"""\nCombined MeSH Stop Word List\n')
            f.write('Generated by extracting Tree Numbers from MeSH XML file.\n"""\n\n')
            f.write("MESH_STOP_WORDS = [\n")
            f.write("    # This list contains MeSH terms that generally do NOT fall into the primary\n")
            f.write(f"    # MeSH categories: {cats_formatted_str}.\n")
            f.write("    # It includes:\n")
            f.write("    #   1. Terms explicitly categorized outside of the 'kept' categories.\n")
            f.write("    #   2. Terms for which categories could not be determined.\n\n")
            f.write("    # --- Terms Categorized Outside A, C, D, G ---\n")

            sorted_cats = sorted(grouped_stopwords.keys())

            for cat in sorted_cats:
                terms = sorted(list(grouped_stopwords[cat]))
                cat_desc = TOP_LEVEL_CATEGORY_NAMES.get(cat, "Unknown Category")

                f.write("\n")
                f.write("    #---------------------------------------------------------------------------\n")
                f.write(f"    # MeSH Category (source of these stop words): {cat} - {cat_desc}\n")
                f.write("    #---------------------------------------------------------------------------\n")

                for term in terms:
                    clean_term = term.replace('\\', '\\\\').replace('"', '\\"')
                    f.write(f'    "{clean_term}",\n')

            if missing_tns:
                f.write("\n")
                f.write("    #---------------------------------------------------------------------------\n")
                f.write("    # Terms With Undetermined Categories (No valid Tree Numbers found)\n")
                f.write("    # These are included as stop words as their category could not be confirmed to be A,C,D,G.\n")
                f.write("    #---------------------------------------------------------------------------\n")

                for term in missing_tns:
                    clean_term = term.replace('\\', '\\\\').replace('"', '\\"')
                    f.write(f'    "{clean_term}",\n')

            f.write("]\n")

        print(f"    [+] Successfully generated stop words file: {output_path}")
        return True

    except Exception as e:
        raise RuntimeError(f"Error writing Python stop word file: {e}")

def extract_all_mesh_data_from_xml(xml_file_path: str, output_csv_path: str, output_py_path: str) -> bool:
    """
    Reads the MeSH XML file in a single pass.
    Extracts Headings, UIs, and Tree Numbers, and generates both the CSV and PY files.
    """
    print(f"\n<<< Starting Unified XML Extraction >>>")
    print(f"Reading from: {xml_file_path}")

    if not os.path.exists(xml_file_path):
        raise FileNotFoundError(_get_missing_nlm_file_message(xml_file_path))

    mesh_terms = []
    mesh_term_to_categories = defaultdict(set)
    all_mh_found = set()

    try:
        # Keep a handle on the root so we can drop already-parsed siblings; without
        # this, iterparse retains every emptied record and memory grows unbounded.
        context = ET.iterparse(xml_file_path, events=('start', 'end'))
        _, root = next(context)
        count = 0

        for event, elem in context:
            if event != 'end' or elem.tag != 'DescriptorRecord':
                continue

            count += 1
            if count % 2000 == 0:
                print(f"  Parsed {count:,} descriptor records...", end='\r')

            ui = elem.findtext('.//DescriptorUI')
            name = elem.findtext('.//DescriptorName/String')

            if ui and name:
                mesh_terms.append([name, ui])
                all_mh_found.add(name)

                tree_list = elem.find('.//TreeNumberList')
                if tree_list is not None:
                    for tree_node in tree_list.findall('TreeNumber'):
                        tree_num = tree_node.text
                        if tree_num and tree_num[0].isalpha():
                            top_cat = tree_num[0].upper()
                            mesh_term_to_categories[name].add(top_cat)

            elem.clear()
            root.clear()

        print(f"  Parsed {count:,} total records. Analyzing hierarchies...\n")

        # --- 1. PROCESS STOP WORDS ---
        stop_words_grouped = defaultdict(set)
        terms_missing_tns = sorted(list(all_mh_found - set(mesh_term_to_categories.keys())))
        final_orphans_set = set(terms_missing_tns)

        for mh, cats in mesh_term_to_categories.items():
            if mh in CHECK_TAGS_TO_STOP:
                final_orphans_set.add(mh)
                continue

            if not any(c in CATEGORIES_TO_KEEP for c in cats):
                excluded_cats = sorted(list(cats))
                primary_cat = excluded_cats[0] if excluded_cats else 'Unknown'
                stop_words_grouped[primary_cat].add(mh)

        for tag in CHECK_TAGS_TO_STOP:
            if tag in all_mh_found:
                 final_orphans_set.add(tag)

        all_stop_words = final_orphans_set.copy()
        for s in stop_words_grouped.values():
            all_stop_words.update(s)

        write_stopwords_to_python(output_py_path, stop_words_grouped, sorted(list(final_orphans_set)))

        # --- 2. PROCESS DATAFRAME & TAG CSV ---
        print(f"\n<<< Building DataFrame & Tagging Stop Words >>>")
        df = pd.DataFrame(mesh_terms, columns=["DescriptorName", "DescriptorUI"])
        initial_len = len(df)

        df.drop_duplicates(subset=['DescriptorName'], keep='first', inplace=True)

        stop_words_lower = {s.lower() for s in all_stop_words}
        df['MeSH_stop_term'] = df['DescriptorName'].astype(str).str.lower().isin(stop_words_lower)
        stop_count = df['MeSH_stop_term'].sum()

        # Which trees each descriptor belongs to, so the stop-word set can be
        # recomputed for a different choice of trees without coming back to
        # this 300 MB file. Without it, changing one checkbox would mean
        # re-parsing the XML, and nobody would change it twice.
        df[vocabulary.TREE_COLUMN] = df['DescriptorName'].map(
            lambda n: ''.join(sorted(mesh_term_to_categories.get(n, ()))))

        os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
        df.to_csv(output_csv_path, index=False)

        print(f"    [+] Successfully extracted {len(df):,} unique terms (from {initial_len:,} raw)")
        print(f"    [+] Tagged {stop_count:,} terms as stop words.")
        print(f"    [+] Recorded tree membership, so the stop-word trees can be")
        print(f"        changed later without re-reading the XML.")
        print(f"    [+] Saved CSV to '{output_csv_path}'")

        return True

    except MemoryError:
        raise RuntimeError("Memory Limit Exceeded during XML processing.")
    except Exception as e:
        traceback.print_exc()
        raise RuntimeError(f"Unexpected error during XML extraction: {e}")


def process_raw_mesh_data(xml_file: str, output_csv: str, output_py: str, force_update: bool = False):
    """
    Orchestrates the unified XML extraction and stop word generation process.
    """
    print("--- MeSH Data Processor ---")
    print(f"XML Input:    {xml_file}")
    print(f"CSV Output:   {output_csv}")
    print(f"PY Output:    {output_py}")
    print(f"Force Update: {force_update}")

    print("\n<<< Starting Processing >>>")

    # Because XML generates both files at the exact same time, we check if EITHER is missing
    needs_processing = not os.path.exists(output_csv) or not os.path.exists(output_py) or force_update

    if needs_processing:
        extract_all_mesh_data_from_xml(xml_file, output_csv, output_py)
    else:
        print("    Skipping XML extraction (Both CSV and PY exist, and Force Update is OFF).")

    print("\n<<< Data Processing Complete >>>")