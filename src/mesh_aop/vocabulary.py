# -*- coding: utf-8 -*-
"""
vocabulary.py - which MeSH terms an analysis is allowed to see.

MeSH is organised into sixteen trees. This program was written for adverse
outcome pathways, so it kept four of them - Anatomy, Diseases, Chemicals and
Drugs, Phenomena and Processes - and treated every term in the other twelve as
a stop word. Geographicals, publication types, occupations and named groups are
noise in that question, and dropping them is most of what makes the resulting
network readable.

They are not noise in every question. A study of where research is done needs
Geographicals; one about who does it needs Disciplines and Occupations. So the
trees are a setting, and this module turns that setting into the actual set of
terms a run must ignore.

HOW IT AVOIDS RE-READING THE XML
    The obvious implementation re-parses the 300 MB descriptor file whenever a
    checkbox changes, which would make the setting unusable. Instead the tree
    letters of every descriptor are written into the terms CSV once, when the
    XML is first processed, and the stop-word set for any choice of trees is a
    filter over that CSV - a second or two, no download, no XML.

    A CSV written by an older version has no such column. That case falls back
    to the list shipped in mesh_stop_words.py, which is the same four-tree
    answer, and says so, rather than silently analysing with no stop words at
    all - a run with none would appear to work and produce a network full of
    countries and journal types.

Standard library only, so the settings window can list the trees without
loading pandas.
"""

import csv
import os

# The sixteen top-level MeSH trees, by their letter. Kept here rather than in
# mesh_data_processor because the settings window needs to list them and has no
# reason to import pandas, requests and tqdm to do it.
TREES = (
    ('A', 'Anatomy'),
    ('B', 'Organisms'),
    ('C', 'Diseases'),
    ('D', 'Chemicals and Drugs'),
    ('E', 'Analytical, Diagnostic and Therapeutic Techniques, and Equipment'),
    ('F', 'Psychiatry and Psychology'),
    ('G', 'Phenomena and Processes'),
    ('H', 'Disciplines and Occupations'),
    ('I', 'Anthropology, Education, Sociology, and Social Phenomena'),
    ('J', 'Technology, Industry, Agriculture'),
    ('K', 'Humanities'),
    ('L', 'Information Science'),
    ('M', 'Named Groups'),
    ('N', 'Health Care'),
    ('V', 'Publication Characteristics'),
    ('Z', 'Geographicals'),
)

TREE_NAMES = dict(TREES)

# Kept by default: the four an adverse outcome pathway is built from. Every
# other tree is excluded, which is what the shipped stop-word list encodes.
DEFAULT_KEPT = ('A', 'C', 'D', 'G')
DEFAULT_EXCLUDED = tuple(letter for letter, _name in TREES
                         if letter not in DEFAULT_KEPT)

# Male and Female are check tags. They sit in no tree at all, so no choice of
# trees can exclude them, and they are attached to a large fraction of all
# clinical articles - which makes them two of the highest-degree nodes in any
# network that keeps them, connected to everything and meaning nothing. Stopped
# by default for that reason; kept by anyone actually studying sex differences.
SEXES = ('Male', 'Female')

# The column the tree letters are written into, and its separator.
TREE_COLUMN = 'TreeCategories'
TERM_COLUMN = 'DescriptorName'


def parse_letters(value):
    """A stored tree selection into a tuple of letters, in TREES order.

    Accepts "B;E;F", "B,E,F", "BEF" and a list, because all four turn up: the
    settings file writes the first, a user editing it by hand writes any of
    them. Anything that is not a known tree letter is dropped.
    """
    if not value:
        return ()
    if isinstance(value, (list, tuple, set, frozenset)):
        found = {str(v).strip().upper() for v in value}
    else:
        text = str(value).upper()
        found = {c for c in text if c.isalpha()}
    return tuple(letter for letter, _name in TREES if letter in found)


def format_letters(letters):
    """A tree selection as it is stored in the settings file."""
    return ';'.join(parse_letters(letters))


def parse_terms(value):
    """The hand-entered term list into a tuple, in the order given.

    Semicolon-delimited, because MeSH headings contain commas as a matter of
    course - "Dermatitis, Allergic Contact" is one heading, not two. Commas are
    NOT accepted as an alternative here for exactly that reason: splitting that
    heading on its comma would silently stop two terms that do not exist.
    """
    if not value:
        return ()
    if isinstance(value, (list, tuple, set, frozenset)):
        items = [str(v).strip() for v in value]
    else:
        items = [part.strip() for part in str(value).split(';')]
    out = []
    for term in items:
        if term and term not in out:
            out.append(term)
    return tuple(out)


def _read_terms_csv(path):
    """(term, tree letters) for every descriptor, or None if unusable.

    None means "this CSV cannot answer the question" - it is missing, or it was
    written before the tree column existed. The caller falls back rather than
    treating an empty answer as "nothing is a stop word".
    """
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path, 'r', encoding='utf-8', newline='') as fh:
            reader = csv.DictReader(fh)
            if not reader.fieldnames or TREE_COLUMN not in reader.fieldnames:
                return None
            term_col = (TERM_COLUMN if TERM_COLUMN in reader.fieldnames
                        else reader.fieldnames[0])
            return [(row.get(term_col) or '', row.get(TREE_COLUMN) or '')
                    for row in reader]
    except (OSError, UnicodeError, csv.Error):
        return None


def needs_tree_column(path):
    """True when this terms CSV exists but predates the tree record.

    That file is what makes the tree selection work, so without the column
    every choice on the Stop words tab silently falls back to the built-in
    four-tree list. A missing file is NOT this case - it is rebuilt anyway.
    """
    return bool(path) and os.path.exists(path) and _read_terms_csv(path) is None


def shipped_stop_words():
    """The four-tree list built into the package."""
    try:
        from .mesh_stop_words import MESH_STOP_WORDS
        return list(MESH_STOP_WORDS)
    except Exception:                                              # noqa: BLE001
        return []


# The set in force for this run. Held here rather than bound into each module
# at import time, which is what it used to be: a frozenset built from the
# shipped list the moment network.py was imported could not be changed
# afterwards, so no setting could ever have affected it. configure() is called
# once by the CLI before the network is built; until then active() answers with
# the shipped list, so a notebook driving these modules directly still works.
_ACTIVE = None


def configure(words):
    """Fix the stop words for this run. None restores the shipped list."""
    global _ACTIVE
    _ACTIVE = None if words is None else frozenset(words)
    return _ACTIVE


def active():
    """The stop words in force."""
    global _ACTIVE
    if _ACTIVE is None:
        _ACTIVE = frozenset(shipped_stop_words())
    return _ACTIVE


def resolve(excluded_trees=DEFAULT_EXCLUDED, keep_sexes=False, extra_terms='',
            mesh_csv_path=None):
    """The stop words in force, and a record of where they came from.

    Returns (frozenset, provenance). The provenance is written to the run
    ledger: a network is only interpretable next to the vocabulary it was built
    from, and "which terms were excluded" is not recoverable from the network
    afterwards.
    """
    excluded = parse_letters(excluded_trees)
    extra = parse_terms(extra_terms)
    rows = _read_terms_csv(mesh_csv_path)

    provenance = {
        'excluded_trees': ';'.join(excluded),
        'excluded_tree_names': [f'{c} - {TREE_NAMES[c]}' for c in excluded],
        'kept_trees': ';'.join(c for c, _n in TREES if c not in excluded),
        'sexes': 'kept' if keep_sexes else 'stopped',
        'hand_entered': list(extra),
    }

    if rows is None:
        # No tree column to filter on. The shipped list is the four-tree
        # answer, so it is right whenever the trees are still the default and
        # wrong - loudly, here - when they are not.
        words = set(shipped_stop_words())
        provenance['source'] = 'shipped list (mesh_stop_words.py)'
        if set(excluded) != set(DEFAULT_EXCLUDED):
            provenance['warning'] = (
                'The chosen trees could not be applied: the MeSH terms file '
                'predates the tree column. The built-in four-tree list was '
                'used instead. Turn on "Rebuild the stop-word list" to '
                'regenerate it from the descriptor XML.')
    else:
        words = set()
        for term, letters in rows:
            term = (term or '').strip()
            if not term:
                continue
            cats = {c for c in (letters or '').upper() if c.isalpha()}
            # A term in no tree is a check tag or an orphan: nothing keeps it,
            # so it is stopped, which is what the original list did too.
            if not cats or not (cats - set(excluded)):
                words.add(term)
        provenance['source'] = f'{os.path.basename(mesh_csv_path)} ({len(rows):,} descriptors)'

    for sex in SEXES:
        if keep_sexes:
            words.discard(sex)
        else:
            words.add(sex)
    words.update(extra)

    provenance['total'] = len(words)
    return frozenset(words), provenance
