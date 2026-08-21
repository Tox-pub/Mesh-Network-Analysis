# -*- coding: utf-8 -*-
"""
uninstall.py - find everything this program put on the machine, and remove it.

Deleting the program folder is not enough. The ETL keeps a working copy of the
master database in the system temp directory, which no user would think to look
for and which is the same size as the database itself. An output folder pointed
somewhere else stays behind too, and a pip install leaves a package and its
console scripts in the environment.

This module only reports and removes; it decides nothing. Both front-ends - the
`mesh-uninstall` command and the Workbench's uninstall screen - drive it, so
they cannot disagree about what "clean" means.

Nothing outside the inventory is ever touched, and results are opt-in: work that
took hours to produce is not deleted because someone clicked through a wizard.
"""

import os
import shutil
import sys
import tempfile
from pathlib import Path

# What a category means, in the order a user should think about it.
#   removable : safe to delete, the program put it there
#   regenerate: deleting costs download or compute time, nothing else
#   precious  : the user's own results; never removed unless asked explicitly
APPLICATION = 'application'
DOWNLOADED = 'downloaded'
DERIVED = 'derived'
RESULTS = 'results'
CACHE = 'cache'
CONFIG = 'config'

CATEGORY_LABEL = {
    APPLICATION: 'Application files',
    DOWNLOADED: 'Downloaded source data',
    DERIVED: 'Built databases',
    RESULTS: 'Your results',
    CACHE: 'Temporary working files',
    CONFIG: 'Settings',
}

# Categories a "remove everything the program installed" default should tick.
# RESULTS is deliberately absent.
DEFAULT_SELECTED = (APPLICATION, DOWNLOADED, DERIVED, CACHE, CONFIG)


class Item:
    """One entry in the inventory: what it is, and what deleting it costs.

    `targets` is usually the single path, but a scatter of stray files - the
    orphaned ETL shards - is presented and removed as one entry.
    """

    __slots__ = ('path', 'category', 'label', 'note', 'bytes', 'is_dir',
                 'targets', 'removable')

    def __init__(self, path, category, label, note='', targets=None,
                 removable=True):
        self.path = Path(path)
        self.category = category
        self.label = label
        self.note = note
        self.targets = [Path(t) for t in targets] if targets else [self.path]
        self.is_dir = self.path.is_dir()
        self.bytes = sum(_size(t) for t in self.targets)
        # Reported so the user knows it exists, but never deleted by this tool.
        self.removable = removable

    @property
    def gb(self):
        return self.bytes / 1e9

    def __repr__(self):
        return f'<Item {self.category} {self.path} {self.gb:.2f}GB>'


def _size(path):
    """Bytes used by a file or a whole tree. Unreadable entries count as zero."""
    try:
        if path.is_file():
            return path.stat().st_size
        if not path.is_dir():
            return 0
    except OSError:
        return 0
    total = 0
    for base, _, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(base, f))
            except OSError:
                pass
    return total


def _add(items, path, category, label, note='', removable=True):
    p = Path(path)
    try:
        if p.exists():
            items.append(Item(p, category, label, note, removable=removable))
    except OSError:
        pass


def is_portable(project_dir):
    """True when this is the extracted Windows bundle rather than a pip install.

    The bundle carries its own interpreter beside the app folder; a pip install
    never does. The distinction matters because uninstalling the bundle means
    deleting the folder the running program lives in.
    """
    root = Path(project_dir)
    return (root / 'python' / 'python.exe').exists() and (root / 'app').is_dir()


def _config_dirs(root):
    """The `directories` block from mesh_config.json, or {} if unreadable.

    An installed copy keeps its settings under the user profile, not beside the
    program, so the file is located the same way the application locates it.
    """
    from . import paths as _paths
    cfg_path = Path(_paths.config_path(root))
    if not cfg_path.exists():
        cfg_path = Path(root) / 'mesh_config.json'
    if not cfg_path.exists():
        return {}
    try:
        import json
        with open(cfg_path, encoding='utf-8') as fh:
            return json.load(fh).get('directories') or {}
    except (OSError, ValueError):
        return {}


def _workspace_base(root):
    """Where the ETL was told to do its scratch work, if anywhere."""
    value = (_config_dirs(root).get('etl_workspace_dir') or '').strip()
    return Path(value) if value else None


def _redirected(root):
    """Folders the user pointed outside the project, read straight from the file.

    Deliberately not via MeshConfig: constructing one calls refresh_paths(),
    which creates every directory it resolves. An uninstaller that mkdirs the
    folders it is about to delete - and recreates them on the rescan afterwards
    - is worse than useless.
    """
    out = []
    dirs = _config_dirs(root)
    for key, label in (('output_dir', 'Results folder (redirected)'),
                       ('input_dir', 'Input folder (redirected)')):
        value = (dirs.get(key) or '').strip()
        if value:
            out.append((Path(value), label))
    return out


def inventory(project_dir, config=None):
    """Everything attributable to this program, largest concern first.

    `config` is accepted and ignored; redirected folders are read from
    mesh_config.json directly. Kept in the signature so existing callers work.
    """
    root = Path(project_dir).resolve()
    items = []
    portable = is_portable(root)

    if portable:
        _add(items, root / 'python', APPLICATION, 'Bundled Python runtime',
             'The interpreter this program runs on.')
        _add(items, root / 'app', APPLICATION, 'Application and pipeline code')
    else:
        # Reported, never deleted. Under an editable install this path is the
        # user's own source checkout - deleting it would destroy a working tree
        # and a git repository. Even for an ordinary install, removing
        # site-packages by hand leaves pip's metadata describing a package that
        # is no longer there. `pip uninstall` is the only correct route, so the
        # entry exists to point at it.
        pkg = Path(__file__).resolve().parent
        _add(items, pkg, APPLICATION, 'Installed mesh_aop package',
             'Removed with pip, not by this tool - the command is shown at the end.',
             removable=False)

    _add(items, root / 'data' / 'raw' / 'pubmed_baseline', DOWNLOADED,
         'PubMed baseline archive',
         'Re-downloadable from NLM, but that is tens of GB and several hours.')
    _add(items, root / 'data' / 'raw' / 'master_mesh_database.db', DERIVED,
         'Master annotation database',
         'Compiled from the archive above. Rebuilding takes hours.')
    _add(items, root / 'data' / 'raw' / 'desc2025.xml', DOWNLOADED,
         'MeSH descriptor file')
    _add(items, root / 'data' / 'processed', DERIVED,
         'Processed databases and networks',
         'Everything the pipeline computed for this project.')

    _add(items, root / 'results', RESULTS, 'Results folder',
         'Figures, workbooks and reports. Not removed unless you ask.')

    # An output folder pointed outside the project is invisible to anyone who
    # just deletes the program folder.
    for target, label in _redirected(root):
        target = target.resolve()
        if root in target.parents or target == root:
            continue                      # already covered by the entries above
        # A redirected folder that is itself a checkout belongs to someone else's
        # project. Report it so it is not a surprise, but never offer to delete a
        # working tree.
        checkout = (target / '.git').exists() or (target.parent / '.git').exists()
        _add(items, target, RESULTS, label,
             'Inside a source checkout - reported only, not removed.'
             if checkout else str(target),
             removable=not checkout)

    # The ETL's scratch space. Never cleaned up by the pipeline, and it holds a
    # full copy of the master database. Both the default location and a
    # configured one are checked: a build that ran with the workspace pointed
    # elsewhere would otherwise leave the largest single leftover untouched.
    tmp = Path(tempfile.gettempdir())
    bases, seen = [], set()
    for base in (_workspace_base(root), tmp):
        if base and str(base).lower() not in seen:
            seen.add(str(base).lower())
            bases.append(base)
    for base in bases:
        where = ('the system temp folder' if base == tmp
                 else f'the configured workspace, {base}')
        _add(items, base / 'mesh_etl_workspace', CACHE, 'ETL working folder',
             f'A working copy of the master database, left in {where}.')

    # Shards used to be written directly into temp rather than inside the
    # workspace, so a machine that ran an older build can still have them loose.
    orphans = sorted(tmp.glob('shard_*.db'))
    if orphans:
        items.append(Item(
            tmp, CACHE, f'Orphaned ETL shards ({len(orphans)} files)',
            'Left behind only when a database build was interrupted.',
            targets=orphans))

    # An installed copy keeps per-user state outside the program folder, so
    # Windows removing what Setup wrote leaves all of this behind.
    from . import paths as _paths
    dirs = _config_dirs(root)
    if not _paths.is_portable(root):
        user_data = (dirs.get('data_dir') or '').strip() or _paths.default_data_dir()
        _add(items, Path(user_data), DERIVED, 'Data folder (your profile)',
             'Downloaded archives and the databases built from them.')
        user_results = (dirs.get('results_dir') or '').strip() or _paths.default_results_dir()
        _add(items, Path(user_results), RESULTS, 'Results folder (your profile)',
             'Your figures, workbooks and reports. Not removed unless you ask.')
        _add(items, _paths.log_dir(), CACHE, 'Run logs')
        _add(items, _paths.config_path(root), CONFIG, 'Saved settings',
             'Your search terms, folders and NCBI credentials.')
        # Only this account's copies. Another user of the same machine keeps
        # theirs, which is correct - and worth saying out loud, because it means
        # an uninstall is not necessarily the end of this program's data.
    _add(items, root / 'mesh_config.json', CONFIG, 'Saved settings',
         'Your search terms, folders and NCBI credentials.')

    items.sort(key=lambda i: i.bytes, reverse=True)
    return items


def remove(items, dry_run=False, on_event=None):
    """Delete the given items. Returns (removed, freed_bytes, failures).

    A failure never stops the pass: a file held open by another program should
    not prevent the other 50 GB from being reclaimed.
    """
    say = on_event or (lambda *_: None)
    removed, freed, failures = 0, 0, []

    for item in items:
        # Defence in depth: a front-end should never offer these, but a bug in
        # one must not turn into a deleted source tree.
        if not item.removable:
            say('skip', f'{item.label} - not removable by this tool')
            continue
        ok = True
        for target in item.targets:
            if dry_run:
                continue
            try:
                if target.is_dir():
                    shutil.rmtree(target)
                elif target.exists():
                    target.unlink()
            except OSError as exc:
                ok = False
                failures.append((target, str(exc)))
                say('fail', f'{target}: {exc}')
        if ok:
            removed += 1
            freed += item.bytes
            say('done', f'{item.label} ({item.gb:.2f} GB)')
    return removed, freed, failures


def pip_hint():
    """The command that removes the package itself, where one applies."""
    return f'{Path(sys.executable).name} -m pip uninstall mesh_aop_network'


def summarise(items):
    """Totals per category, for a front-end to display."""
    out = {}
    for it in items:
        agg = out.setdefault(it.category, [0, 0])
        agg[0] += 1
        agg[1] += it.bytes
    return out
