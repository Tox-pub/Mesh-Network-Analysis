# -*- coding: utf-8 -*-
"""
integrity.py - is what is on disk actually usable, and if not, what now.

A run of this pipeline can take a night. Its intermediate files are large, and
several of them are written over hours: a machine that sleeps, a disk that fills,
a OneDrive client that dehydrates a file mid-write, and what is left behind looks
complete - right name, plausible size - and fails much later, deep inside a step
that had no way to know its input was rubbish.

The cost is asymmetric. The master annotation database is 44 GB downloaded and
several hours compiled; everything else is derived from it and takes minutes.
So the checks here are deliberately unequal: the master database is verified
structurally and recorded in a health file so damage is caught the moment it
appears, while derived artefacts are checked cheaply and simply thrown away when
they are wrong.

Three things this module will not do:

  * Guess. A file is reported CORRUPT only when opening it actually fails, and
    SUSPECT when it is merely implausible. The two are never conflated - one
    warrants deleting hours of work, the other warrants a look.
  * Delete anything on its own. It reports, and something else acts.
  * Touch the user's results. Every artefact here is machinery and reproducible.
"""

import json
import os
import sqlite3
import time
from datetime import datetime

from . import paths as _paths

# ---------------------------------------------------------------- statuses
OK = 'ok'
MISSING = 'missing'
EMPTY = 'empty'            # exists, zero bytes - a write that never started
CORRUPT = 'corrupt'        # opening or parsing it fails outright
SUSPECT = 'suspect'        # opens, but the contents do not add up
ORPHAN = 'orphan'          # a temp/partial file no completed run would leave

# What a damaged artefact costs to replace, which is what should drive whether
# anyone is willing to delete it.
COST_HOURS = 'hours'
COST_MINUTES = 'minutes'
COST_SECONDS = 'seconds'

HEALTH_SUFFIX = '.health.json'

# Pipeline order. Each entry names the step that produces an artefact, so a
# damaged file can be turned into "resume from here" rather than "start again".
STEP_ORDER = ('baseline', 'process', 'data_ops', 'network', 'secondary', 'viz', 'benchmark')
STEP_LABEL = {
    'baseline': 'Step 0 - build the master annotation database',
    'process': 'Step 1 - process the MeSH support files',
    'data_ops': 'Step 2 - collect articles and build the PMID databases',
    'network': 'Step 3 - construct, filter and score the network',
    'secondary': 'Step 3.5 - secondary analysis',
    'viz': 'Step 4 - figures',
    'benchmark': 'Step 5 - benchmark and validation',
}


class Artifact:
    """One file the pipeline produced, and whether it can be trusted."""

    def __init__(self, key, path, kind, step, cost, label, status=OK,
                 detail='', size=0, removable=True):
        self.key, self.path, self.kind = key, str(path), kind
        self.step, self.cost, self.label = step, cost, label
        self.status, self.detail, self.size = status, detail, size
        self.removable = removable

    @property
    def broken(self):
        return self.status in (EMPTY, CORRUPT, ORPHAN)

    @property
    def size_mb(self):
        return self.size / 1e6

    def __repr__(self):
        return f'<Artifact {self.key} {self.status}>'


# ------------------------------------------------------------------ probes

def _stat(path):
    try:
        return os.path.getsize(_paths.long_path(path))
    except OSError:
        return None


def check_sqlite(path, table=None, quick=True):
    """(status, detail, rows). Opens read-only so a live run is never disturbed."""
    size = _stat(path)
    if size is None:
        return MISSING, 'not present', 0
    if size == 0:
        return EMPTY, 'zero bytes - the write never got started', 0
    try:
        conn = sqlite3.connect(f'file:{_paths.long_path(path)}?mode=ro', uri=True, timeout=30.0)
    except sqlite3.Error as exc:
        return CORRUPT, f'cannot be opened: {exc}', 0
    try:
        if quick:
            # quick_check walks the b-tree structure without the full O(n log n)
            # integrity_check, which on an 8 GB database takes long enough that
            # nobody would ever run it.
            result = conn.execute('PRAGMA quick_check(1)').fetchone()
            if result and str(result[0]).lower() != 'ok':
                return CORRUPT, f'SQLite reports: {result[0]}', 0
        if table:
            cols = conn.execute(f'PRAGMA table_info({table})').fetchall()
            if not cols:
                return SUSPECT, f"the '{table}' table is missing", 0
            rows = conn.execute(f'SELECT count(*) FROM {table}').fetchone()[0]
            if rows == 0:
                return SUSPECT, f"'{table}' is present but empty", 0
            return OK, f'{rows:,} rows', rows
        return OK, 'opens cleanly', 0
    except sqlite3.DatabaseError as exc:
        return CORRUPT, f'unreadable: {exc}', 0
    finally:
        try:
            conn.close()
        except Exception:
            pass


def check_json_graph(path):
    """(status, detail, nodes). A network JSON that will not parse is worthless."""
    size = _stat(path)
    if size is None:
        return MISSING, 'not present', 0
    if size == 0:
        return EMPTY, 'zero bytes', 0
    try:
        with open(_paths.long_path(path), 'r', encoding='utf-8') as fh:
            data = json.load(fh)
    except (OSError, ValueError) as exc:
        # A truncated JSON is the classic signature of a run killed mid-write.
        return CORRUPT, f'not valid JSON: {str(exc)[:90]}', 0
    elements = (data or {}).get('elements')
    if not isinstance(elements, dict):
        return SUSPECT, "no 'elements' section - not a network file", 0
    nodes, edges = elements.get('nodes') or [], elements.get('edges') or []
    if not nodes:
        return SUSPECT, 'contains no nodes', 0
    return OK, f'{len(nodes):,} terms, {len(edges):,} relations', len(nodes)


def check_text(path, min_bytes=1):
    size = _stat(path)
    if size is None:
        return MISSING, 'not present', 0
    if size < min_bytes:
        return EMPTY, 'zero bytes', 0
    return OK, f'{size:,} bytes', size


# ------------------------------------------------- the master database health

def health_path(db_path):
    return str(db_path) + HEALTH_SUFFIX


def write_health(db_path, rows=None, note=''):
    """Record what a good master database looked like, right after building one.

    Size and modification time are what make later damage detectable without
    re-scanning eight gigabytes: if either has moved and no build was run, the
    file has been altered by something that had no business altering it -
    a sync client, a backup tool, a half-finished copy.
    """
    try:
        size = os.path.getsize(_paths.long_path(db_path))
    except OSError:
        return None
    if rows is None:
        _, _, rows = check_sqlite(db_path, 'master_mesh_annotations')
    record = {
        'path': str(db_path),
        'size_bytes': size,
        'modified': os.path.getmtime(_paths.long_path(db_path)),
        'rows': rows,
        'verified': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'note': note,
    }
    try:
        with open(_paths.long_path(health_path(db_path)), 'w', encoding='utf-8') as fh:
            json.dump(record, fh, indent=2)
        return record
    except OSError:
        return None


def read_health(db_path):
    try:
        with open(_paths.long_path(health_path(db_path)), 'r', encoding='utf-8') as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def check_master_db(db_path, deep=False):
    """Status of the master database, cheaply, using the health record if there is one.

    `deep=False` is what runs before every pipeline step: open it, structural
    quick_check, confirm the table has rows, and compare against the recorded
    size. That is a second or two even on 8 GB. `deep=True` runs the full
    integrity_check, which is minutes, and is offered only on request.
    """
    status, detail, rows = check_sqlite(db_path, 'master_mesh_annotations', quick=True)
    if status != OK:
        return status, detail, rows

    if deep:
        try:
            conn = sqlite3.connect(f'file:{_paths.long_path(db_path)}?mode=ro', uri=True, timeout=600.0)
            try:
                result = conn.execute('PRAGMA integrity_check').fetchone()
                if result and str(result[0]).lower() != 'ok':
                    return CORRUPT, f'full integrity check failed: {result[0]}', rows
            finally:
                conn.close()
        except sqlite3.DatabaseError as exc:
            return CORRUPT, f'full integrity check could not run: {exc}', rows

    recorded = read_health(db_path)
    if recorded:
        size = _stat(db_path) or 0
        if recorded.get('size_bytes') and size != recorded['size_bytes']:
            delta = size - recorded['size_bytes']
            return (SUSPECT,
                    f"{rows:,} rows, but the file is {abs(delta)/1e6:,.0f} MB "
                    f"{'larger' if delta > 0 else 'smaller'} than when it was built "
                    f"({recorded.get('verified', 'unknown date')}). Something "
                    f"changed it outside a build.", rows)
        if recorded.get('rows') and rows < recorded['rows']:
            return (SUSPECT,
                    f"{rows:,} rows, down from {recorded['rows']:,} at build time", rows)
    return OK, detail + (' (matches its health record)' if recorded else
                         ' (no health record - built before these were kept)'), rows


# ------------------------------------------------------------------ the scan

def _orphans(directories):
    """Temp and part-files no completed run leaves behind."""
    found = []
    patterns = ('.tmp', '.json.tmp', '.part', '.partial', '.crdownload', '.download')
    for d in directories:
        if not d:
            continue
        d = str(d)
        if not os.path.isdir(_paths.long_path(d)):
            continue
        try:
            names = os.listdir(_paths.long_path(d))
        except OSError:
            continue
        for name in names:
            full = os.path.join(d, name)
            if not os.path.isfile(_paths.long_path(full)):
                continue
            lower = name.lower()
            if any(lower.endswith(p) for p in patterns):
                found.append(Artifact(
                    key=name, path=full, kind='temp', step='network',
                    cost=COST_SECONDS, label=f'Leftover temporary file: {name}',
                    status=ORPHAN, size=_stat(full) or 0,
                    detail='written by a run that did not finish'))
    return found


def scan(config, deep=False):
    """Every artefact this configuration depends on, and its condition.

    Ordered by pipeline step, so the first broken thing in the list is also the
    earliest point a rerun would have to start from.
    """
    f = config.files
    out = []

    def add(key, path, kind, step, cost, label, table=None, removable=True):
        if kind == 'sqlite':
            status, detail, _ = (check_master_db(path, deep=deep) if key == 'master_db'
                                 else check_sqlite(path, table))
        elif kind == 'graph':
            status, detail, _ = check_json_graph(path)
        else:
            status, detail, _ = check_text(path)
        out.append(Artifact(key, path, kind, step, cost, label, status, detail,
                            _stat(path) or 0, removable))

    add('master_db', f['master_db'], 'sqlite', 'baseline', COST_HOURS,
        'Master annotation database', 'master_mesh_annotations', removable=True)
    add('mesh_terms_csv', f['mesh_terms_csv'], 'text', 'process', COST_MINUTES,
        'MeSH term table')
    add('pmids_db', f['pmids_db'], 'sqlite', 'data_ops', COST_HOURS,
        'Retrieved PMIDs', 'pmids_table')
    add('cleaned_db', f['cleaned_db'], 'sqlite', 'data_ops', COST_MINUTES,
        'Cleaned PMIDs with MeSH', 'pmids_table')
    add('full_network', f['full_network'], 'graph', 'network', COST_MINUTES,
        'Unfiltered co-occurrence network')
    add('glf_subgraph', f['glf_subgraph'], 'graph', 'network', COST_HOURS,
        'GLF optimiser subgraph')
    add('sa_subgraph', f['sa_subgraph'], 'graph', 'network', COST_HOURS,
        'Simulated annealing subgraph')
    add('consensus_lcc', f['consensus_lcc'], 'graph', 'network', COST_MINUTES,
        'Consensus network (largest component)')
    add('final_network', f['final_network'], 'graph', 'network', COST_MINUTES,
        'Final scored network')
    add('relevance_db', f['relevance_db'], 'sqlite', 'network', COST_MINUTES,
        'Mean relevancy scores', 'article_relevance_scores')

    # _orphans lists one directory at a time and does not recurse, so the
    # networks and databases subfolders have to be named or a half-written file
    # inside them would never be found.
    out.extend(_orphans([config.active_raw_dir, config.active_source_dir,
                         getattr(config, 'networks_dir', None),
                         getattr(config, 'databases_dir', None),
                         config.results_dir,
                         getattr(config, 'figures_dir', None),
                         getattr(config, 'secondary_dir', None),
                         getattr(config, 'benchmark_dir', None),
                         getattr(config, 'benchmark_inputs_dir', None),
                         getattr(config, 'benchmark_ranking_dir', None),
                         getattr(config, 'benchmark_validation_dir', None),
                         getattr(config, 'benchmark_network_dir', None)]))
    return out


# What each step writes, as (config.files key or a results glob, description).
# Used to answer "will running this replace something I already have?" before a
# run starts, rather than after it has already done so.
STEP_OUTPUTS = {
    'baseline': [('master_db', 'Master annotation database')],
    'process': [('mesh_terms_csv', 'MeSH term table')],
    'data_ops': [('pmids_db', 'Retrieved PMIDs'),
                 ('cleaned_db', 'Cleaned PMIDs with MeSH')],
    'network': [('full_network', 'Unfiltered co-occurrence network'),
                ('glf_subgraph', 'GLF optimiser subgraph'),
                ('sa_subgraph', 'Simulated annealing subgraph'),
                ('consensus_lcc', 'Consensus network'),
                ('final_network', 'Final scored network'),
                ('relevance_db', 'Mean relevancy scores')],
}

# Results-folder outputs, matched by glob against the file prefix.
STEP_RESULT_GLOBS = {
    'network': [('{p}_final_network_nodes_and_edges.xlsx', 'Network workbook'),
                ('{p}_run_annotations.csv', 'AOP annotation template')],
    'secondary': [('{p}_Top_Network_Articles.csv', 'Top articles export')],
    'viz': [('figures/{p}_*.jpeg', 'Figures (JPEG)'),
            ('figures/{p}_*.tif', 'Figures (TIFF)'),
            ('figures/{p}_*.png', 'Figures (PNG)'),
            ('figures/{p}_*.svg', 'Figures (SVG)'),
            ('figures/{p}_*.html', 'Interactive figures')],
    'benchmark': [('benchmark/*', 'Benchmark and validation outputs')],
}

# Which steps a given selection actually runs. 'all' is every one of them.
STEP_CHAIN = {
    'all': ('process', 'data_ops', 'network', 'secondary', 'viz', 'benchmark'),
}


def _step_rewrites(config, step, key):
    """Will this step actually rewrite this output, or skip it because it exists?

    Most steps do rewrite: the network stage rebuilds whatever it is asked for.
    Two exceptions, and warning about either is simply untrue.

    The MeSH support files are built once and reused, and only rebuilt when
    explicitly refreshed.

    And with the bundled reference corpus in use, retrieval, network
    construction and the secondary analysis do not run at all - the corpus ships
    with that work already done. Their outputs are read-only files inside the
    program folder. Warning that a run is about to replace them was alarming and
    false in the same breath: nothing can, and nothing will.
    """
    if step == 'process' and key == 'mesh_terms_csv':
        return bool(config.get('search_parameters', 'update_mesh_support_files'))
    if getattr(config, 'use_reference_data', False) and step in REFERENCE_SKIPPED_STEPS:
        return False
    return True


# Kept in step with cli.REFERENCE_SKIPS. Named here rather than imported to
# avoid a cycle: integrity is imported by the CLI, not the other way round.
REFERENCE_SKIPPED_STEPS = ('data_ops', 'network', 'secondary')


def would_overwrite(config, step):
    """What running `step` with this prefix would replace.

    Returns a list of (description, path, size_bytes), empty when nothing is at
    risk - which is the common case for a prefix that has not been used yet, and
    the case in which the user should not be asked anything at all.

    Only counts files that actually exist, and only those the step will really
    rewrite. The MeSH term table is the case that made the difference: Step 1
    rebuilds it only when it is missing or a refresh was asked for, so once it
    exists a later run skips it - and warning that it was about to be replaced
    was simply untrue. Someone who had just generated it from the Database
    screen was then told, on their very next run, that the run would destroy it.
    """
    import glob as _glob

    steps = STEP_CHAIN.get(step, (step,))
    out, seen = [], set()

    for st in steps:
        for key, label in STEP_OUTPUTS.get(st, ()):
            path = str(config.files.get(key, ''))
            if not path or path in seen:
                continue
            if not os.path.exists(_paths.long_path(path)):
                continue
            if not _step_rewrites(config, st, key):
                continue
            seen.add(path)
            out.append((label, path, _stat(path) or 0))

        for pattern, label in STEP_RESULT_GLOBS.get(st, ()):
            spec = os.path.join(str(config.results_dir),
                                pattern.format(p=config.prefix))
            matches = [m for m in _glob.glob(spec) if os.path.isfile(m)]
            matches = [m for m in matches if m not in seen]
            if not matches:
                continue
            seen.update(matches)
            if len(matches) == 1:
                out.append((label, matches[0], _stat(matches[0]) or 0))
            else:
                total = sum(_stat(m) or 0 for m in matches)
                out.append((f'{label} ({len(matches)} files)',
                            os.path.dirname(matches[0]), total))
    return out


def problems(artifacts):
    """Only the ones that need a decision. Missing is not a problem: it is a run
    that has not happened yet, and the pipeline builds what it lacks."""
    return [a for a in artifacts if a.status in (EMPTY, CORRUPT, ORPHAN, SUSPECT)]


def resume_step(artifacts):
    """The earliest step that would have to run again, or None if all is well.

    Anything downstream of a broken artefact is suspect whether or not it opens,
    because it was computed from something that has since proved unreliable -
    so the answer is the earliest damaged step, not a list of them.
    """
    damaged = [a for a in artifacts if a.broken and a.kind != 'temp']
    if not damaged:
        return None
    return min(damaged, key=lambda a: STEP_ORDER.index(a.step)).step


def summarise(artifacts):
    """A short, plain account of the scan, for a dialog or the console."""
    bad = problems(artifacts)
    if not bad:
        healthy = [a for a in artifacts if a.status == OK]
        return (f"Checked {len(artifacts)} files. Everything present opens "
                f"cleanly ({len(healthy)} in place, "
                f"{len(artifacts) - len(healthy)} not yet built).")
    lines = [f"Checked {len(artifacts)} files. {len(bad)} need attention:", '']
    for a in sorted(bad, key=lambda a: STEP_ORDER.index(a.step)):
        lines.append(f"  {a.label}")
        lines.append(f"      {a.status.upper()}: {a.detail}")
        if a.cost == COST_HOURS and a.status != ORPHAN:
            lines.append("      Rebuilding this one takes hours.")
    step = resume_step(artifacts)
    if step:
        lines += ['', f"Remove the damaged files and resume from:",
                  f"  {STEP_LABEL.get(step, step)}"]
    return '\n'.join(lines)


# Files that belong to another file and are meaningless without it. A stale
# write-ahead log left beside a rebuilt database is worse than useless: SQLite
# will try to replay it into the new file.
_SIDECARS = ('-wal', '-shm', '-journal', HEALTH_SUFFIX)


def remove(artifacts):
    """Delete the given artefacts and their sidecars.

    Returns (removed, bytes_freed, failures). Whoever calls this decides what
    goes in the list; this function does not second-guess it, which is what
    lets the GUI offer a merely-suspect file for deletion while the automatic
    repair path sticks to files that are definitely broken.
    """
    removed, freed, failures = [], 0, []
    for a in artifacts:
        if not a.removable:
            continue
        try:
            size = _stat(a.path) or 0
            os.remove(_paths.long_path(a.path))
            removed.append(a)
            freed += size
        except OSError as exc:
            failures.append((a, str(exc)))
            continue
        for suffix in _SIDECARS:
            side = a.path + suffix
            if os.path.exists(_paths.long_path(side)):
                try:
                    freed += _stat(side) or 0
                    os.remove(_paths.long_path(side))
                except OSError:
                    # A leftover sidecar is worth mentioning but not worth
                    # reporting the whole removal as failed over.
                    pass
    return removed, freed, failures
