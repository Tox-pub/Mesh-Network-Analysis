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


def _data_root(root):
    """The user's data folder, if one was configured or can be resolved.

    Read from the file rather than through MeshConfig, for the reason given
    below: constructing one creates every directory it resolves, and an
    uninstaller must not create the folders it is about to remove.
    """
    value = (_config_dirs(root).get('data_dir') or '').strip()
    if value:
        return Path(value)
    try:
        from . import paths as _paths
        return Path(_paths.default_user_data_dir(root))
    except Exception:
        return None


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
    # The daily update files are fetched separately from the yearly baseline and
    # into their own folder, so an inventory that named only the baseline left
    # them behind - and they accumulate for as long as the corpus is kept
    # current.
    _add(items, root / 'data' / 'raw' / 'pubmed_updates', DOWNLOADED,
         'PubMed daily update archives',
         'Records published since the baseline snapshot. Re-downloadable, and '
         'they grow with every update run.')
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
    # The data folder is checked too: when the temp directory is RAM-backed -
    # which /tmp is on most current Linux - the build puts its workspace beside
    # the data instead, and a leftover there is the same eight gigabytes.
    beside_data = _data_root(root)
    bases, seen = [], set()
    for base in (_workspace_base(root), tmp,
                 (beside_data / 'etl_workspace') if beside_data else None):
        if base and str(base).lower() not in seen:
            seen.add(str(base).lower())
            bases.append(base)
    for base in bases:
        where = ('the system temp folder' if base == tmp
                 else f'the configured workspace, {base}')
        _add(items, base / 'mesh_etl_workspace', CACHE, 'ETL working folder',
             f'A working copy of the master database, left in {where}.')

    # Shards are written inside the workspace rather than directly into temp,
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

    _add_shortcuts(items, root)

    items.sort(key=lambda i: i.bytes, reverse=True)
    return items


def _shortcut_target(path):
    """Where a .lnk points, or '' if it cannot be read.

    Read through the shell so a shortcut is resolved the way Explorer resolves
    it. Without a target there is no way to tell one program's shortcut from
    another's, and deleting by name alone would take someone else's.
    """
    try:
        import win32com.client                                   # noqa: F401
    except ImportError:
        pass
    else:
        try:
            shell = win32com.client.Dispatch('WScript.Shell')
            return str(shell.CreateShortCut(str(path)).TargetPath or '')
        except Exception:                                        # noqa: BLE001
            return ''
    # No pywin32 in the bundled interpreter, so read the path out of the file.
    # A .lnk holds its target twice, once UTF-16 and once ANSI; matching a
    # drive-rooted path ending in an executable finds either. This is enough to
    # decide ownership and deliberately no more than that.
    import re
    try:
        raw = Path(path).read_bytes()
    except OSError:
        return ''
    pattern = re.compile(r'[A-Za-z]:\\[^<>:"|?*\x00]{2,240}?'
                         r'(?:pythonw\.exe|python\.exe|\.bat)',
                         re.IGNORECASE)
    for encoding in ('utf-16-le', 'latin-1'):
        try:
            text = raw.decode(encoding, errors='ignore')
        except Exception:                                        # noqa: BLE001
            continue
        match = pattern.search(text)
        if match:
            return match.group(0).strip()
    return ''


def shortcut_locations():
    """Every folder a Start-menu or Desktop shortcut could have been put in."""
    if sys.platform != 'win32':
        return []
    out = []
    for var, tail in (
        ('APPDATA', ('Microsoft', 'Windows', 'Start Menu', 'Programs')),
        ('ProgramData', ('Microsoft', 'Windows', 'Start Menu', 'Programs')),
        ('USERPROFILE', ('Desktop',)),
        ('PUBLIC', ('Desktop',)),
    ):
        base = os.environ.get(var)
        if base:
            out.append(Path(base).joinpath(*tail))
    return out


def find_shortcuts(project_dir=None):
    """Shortcuts on this machine that point INTO a MeSH Workbench install.

    Both installers make these, and so does the portable
    'Create desktop shortcut.bat'. Nothing removed them: an interrupted MSI
    uninstall leaves them behind pointing at a folder that may no longer
    exist, and the portable one was never tracked by anything at all.

    Ownership is decided by the target, never by the name, so a shortcut
    belonging to something else is never touched.
    """
    found = []
    if sys.platform != 'win32':
        return found
    root = str(Path(project_dir).resolve()).lower() if project_dir else None
    for folder in shortcut_locations():
        try:
            entries = list(folder.glob('*.lnk'))
        except OSError:
            continue
        for lnk in entries:
            target = _shortcut_target(lnk).lower()
            if not target:
                continue
            ours = ('mesh workbench' in target
                    or (root and target.startswith(root)))
            if ours:
                found.append(lnk)
    return found


def _add_shortcuts(items, root):
    """Report Start-menu and Desktop shortcuts as their own entry."""
    links = find_shortcuts(root)
    if not links:
        return
    where = sorted({('Desktop' if 'desktop' in str(p.parent).lower()
                     else 'Start menu') for p in links})
    items.append(Item(
        links[0], APPLICATION,
        f'Shortcuts ({", ".join(where)})',
        f'{len(links)} shortcut(s) pointing at this program. '
        f'An interrupted uninstall leaves these behind.',
        targets=links))


def _clear_readonly(func, path, _exc):
    """rmtree error handler: drop the read-only bit and try that entry again.

    shutil.rmtree refuses a read-only file with "Access is denied", which reads
    like a permissions problem the user has to solve and is usually just an
    attribute. Files arriving from an archive or a sync client carry it often.
    """
    import stat
    try:
        os.chmod(path, stat.S_IWRITE)
        func(path)
    except OSError:
        raise


def _delete(target, attempts=3):
    """Remove a file or tree, tolerating read-only bits and transient locks.

    OneDrive, the search indexer and antivirus all take brief handles on files
    they have just seen. A single attempt turns that into a permanent failure
    reported to the user; a short retry usually turns it into nothing at all.
    """
    import time
    last = None
    for attempt in range(attempts):
        try:
            if target.is_dir():
                # onexc replaced onerror in 3.12; support both so this works
                # under the bundled interpreter and an older system one.
                if sys.version_info >= (3, 12):
                    shutil.rmtree(target, onexc=lambda f, p, e: _clear_readonly(f, p, e))
                else:
                    shutil.rmtree(target, onerror=lambda f, p, e: _clear_readonly(f, p, e))
            elif target.exists():
                target.unlink()
            return None
        except OSError as exc:
            last = exc
            time.sleep(0.4 * (attempt + 1))
    return last


def _is_running_from(target):
    """True if this interpreter lives inside `target`.

    Windows will not delete a loaded DLL, so an uninstaller running on the
    Python it is trying to remove cannot finish the job from inside - it fails
    on its own tcl or its own python3.dll. Recognising that is what lets it be
    handled rather than reported as a mystery.
    """
    try:
        exe = Path(sys.executable).resolve()
        target = Path(target).resolve()
        return target == exe or target in exe.parents
    except OSError:
        return False


def _safe_to_delete_unattended(folder):
    """Would `rmdir /s /q folder` be a reasonable thing to run with nobody watching?

    The deferred delete runs after this process is gone, so nothing can stop it
    and nothing will report what it did. That deserves more caution than a
    deletion the user is watching: the cost of refusing a legitimate folder is
    one manual delete, and the cost of accepting a wrong one is unbounded.

    So: it must be several levels deep, must not be a home or system directory,
    and must actually look like an installation of this program.
    """
    folder = Path(folder).resolve()
    parts = folder.parts
    if len(parts) < 3:                       # C:\ or C:\something
        return False

    protected = {Path.home().resolve()}
    for var in ('ProgramFiles', 'ProgramW6432', 'ProgramFiles(x86)', 'SystemRoot',
                'WINDIR', 'LOCALAPPDATA', 'APPDATA', 'USERPROFILE', 'TEMP'):
        value = os.environ.get(var)
        if value:
            try:
                protected.add(Path(value).resolve())
            except OSError:
                pass
    for name in ('Documents', 'Desktop', 'Downloads', 'OneDrive'):
        protected.add((Path.home() / name).resolve())
    if folder in protected:
        return False

    # It has to look like ours. Either the folder is named for the program, or
    # it holds the things an installation holds.
    looks_like_ours = (
        'mesh' in folder.name.lower()
        or (folder / 'portable.marker').exists()
        or ((folder / 'app').is_dir() and (folder / 'python').is_dir())
        or (folder.name.lower() == 'python' and (folder.parent / 'app').is_dir())
    )
    return bool(looks_like_ours)


def schedule_delete_after_exit(folder):
    """Have the folder removed once this process has let go of it.

    The only way to delete the interpreter you are running on. A small batch
    file waits for this process ID to disappear, removes the folder, and then
    removes itself; it is launched detached so it outlives us. Windows only -
    every other system allows deleting a file that is open.

    Returns the helper's path, or None if it could not be arranged.
    """
    if os.name != 'nt':
        return None
    folder = Path(folder).resolve()
    if not _safe_to_delete_unattended(folder):
        return None
    folder = str(folder)
    script = Path(tempfile.gettempdir()) / f'mesh_uninstall_{os.getpid()}.bat'
    body = (
        '@echo off\r\n'
        'rem Written by the MeSH Workbench uninstaller. Waits for the program\r\n'
        'rem to close, removes its folder, then removes itself.\r\n'
        f':wait\r\n'
        f'tasklist /FI "PID eq {os.getpid()}" 2>nul | find "{os.getpid()}" >nul\r\n'
        'if not errorlevel 1 (\r\n'
        '  ping -n 2 127.0.0.1 >nul\r\n'
        '  goto wait\r\n'
        ')\r\n'
        f'rmdir /s /q "{folder}" 2>nul\r\n'
        f'del /f /q "%~f0" >nul 2>&1\r\n'
    )
    try:
        script.write_text(body, encoding='utf-8')
        import subprocess
        subprocess.Popen(['cmd', '/c', str(script)],
                         creationflags=0x00000008 | 0x08000000,   # DETACHED, NO_WINDOW
                         close_fds=True)
        return str(script)
    except OSError:
        return None


def remove(items, dry_run=False, on_event=None):
    """Delete the given items. Returns (removed, freed_bytes, failures, deferred).

    A failure never stops the pass: a file held open by another program should
    not prevent the other 50 GB from being reclaimed.

    `deferred` names anything that could not be deleted because this program is
    running from inside it. That is not a failure the user can act on - it is
    arranged to happen the moment the window closes.
    """
    say = on_event or (lambda *_: None)
    removed, freed, failures, deferred = 0, 0, [], []

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
            exc = _delete(target)
            if exc is None:
                continue
            if _is_running_from(target):
                # Expected, and handled: we cannot delete the interpreter we
                # are executing. Hand it to a helper that runs after we exit.
                if schedule_delete_after_exit(target):
                    deferred.append(target)
                    say('deferred', f'{target}: will be removed when the program closes')
                    continue
            ok = False
            failures.append((target, str(exc)))
            say('fail', f'{target}: {exc}')
        if ok:
            removed += 1
            freed += item.bytes
            say('done', f'{item.label} ({item.gb:.2f} GB)')
    return removed, freed, failures, deferred


def package_is_installed():
    """Was this actually installed with pip, or is it just on the path?

    A self-contained bundle unpacks its dependencies into its own interpreter
    but never installs the application: it is put on PYTHONPATH, so there is no
    distribution to uninstall and `pip uninstall` reports it is not installed.
    Telling someone to run it anyway sends them looking for a problem that is
    not there.
    """
    try:
        from importlib.metadata import PackageNotFoundError, distribution
    except ImportError:
        return False
    for name in ('mesh-aop-network', 'mesh_aop_network'):
        try:
            distribution(name)
            return True
        except PackageNotFoundError:
            continue
        except Exception:                                          # noqa: BLE001
            continue
    return False


def bundle_root():
    """The folder to delete for a self-contained copy, or None.

    A bundle is an interpreter, the application and its wheels in one directory
    tree; removing it is removing the program. Recognised by the interpreter
    living inside that tree rather than by a marker, because the Unix bundles
    ship no marker.
    """
    exe = Path(sys.executable).resolve()
    for parent in list(exe.parents)[:4]:
        if (parent / 'python').is_dir() and any(
                (parent / n).exists() for n in ('MeSH Workbench', 'mesh-pipeline',
                                                'MeSH Workbench.bat', 'app')):
            return parent
    return None


def pip_hint():
    """The command that removes the package itself, where one applies.

    The FULL path to the interpreter, not its name. `Path(sys.executable).name`
    gave 'python3.12', which inside a bundle names an interpreter that exists
    only in that bundle - so the command was copied into a shell that had no
    python3.12 on its PATH and failed with "no such file or directory".
    """
    exe = str(Path(sys.executable))
    if ' ' in exe:
        exe = f'"{exe}"'
    return f'{exe} -m pip uninstall mesh_aop_network'


def removal_instructions():
    """How to remove the program itself on THIS copy, as (heading, lines).

    Three different answers, and giving the wrong one wastes real time:
    a pip install is uninstalled with pip, a self-contained bundle is a folder
    to delete, and a Windows install belongs to Add/Remove Programs.
    """
    root = bundle_root()
    if root is not None:
        if sys.platform == 'win32':
            return ('This copy is self-contained - it installed nothing.', [
                'Everything it needs is inside one folder. To finish, delete it:',
                f'    {root}',
                'If it was installed with the .msi or Setup.exe, use',
                'Settings > Apps > MeSH Workbench > Uninstall instead.'])
        return ('This copy is self-contained - it installed nothing.', [
            'Everything it needs is inside one folder, including its own',
            'Python. To finish, delete that folder:',
            f'    rm -rf "{root}"',
            'There is nothing to uninstall with pip: the application was never',
            'installed into an interpreter, only added to its path.'])
    if sys.platform == 'win32':
        return ('Remove the program itself from Windows.', [
            'Settings > Apps > Installed apps > MeSH Workbench > Uninstall.',
            'If you are running a pip install instead, use:',
            f'    {pip_hint()}'])
    if package_is_installed():
        return ('The package itself is still installed.', [
            'To remove it, run:', f'    {pip_hint()}'])
    return ('The program itself was not installed by a package manager.', [
        'It is running from this folder, so deleting the folder removes it:',
        f'    {Path(__file__).resolve().parents[2]}'])


def summarise(items):
    """Totals per category, for a front-end to display."""
    out = {}
    for it in items:
        agg = out.setdefault(it.category, [0, 0])
        agg[0] += 1
        agg[1] += it.bytes
    return out
