# -*- coding: utf-8 -*-
"""
paths.py - where this program keeps things, and who owns them.

Two very different situations, and one function per question so the GUI and the
CLI can never disagree about the answer.

**Portable.** The zip is extracted somewhere and run from there. Everything -
settings, data, results - stays inside that folder. Self-containment is the
entire point: the folder can be moved, copied to a USB stick, or deleted whole.
Detected by a marker file that ships only in the zip.

**Installed.** The program sits in a shared, usually read-only location and may
be used by several people. Writing settings or results beside the executable is
wrong twice over: under Program Files it fails outright, and even where it
succeeds two users would overwrite each other's configuration and outputs. So
per-user state goes under the user's own profile, the way any other Windows
application behaves.

Within the installed case, what goes where is decided by who produced it:

    asked for, and the user's own       results, and the bulk data folder
    machinery, never asked about        settings, logs, caches, ETL scratch

The bulk data - tens of gigabytes of PubMed archives and the databases built
from them - defaults under the user profile but is deliberately askable, because
the profile is on the system drive and the system drive is frequently too small
to hold it.
"""

import os
import sys
from pathlib import Path

APP_DIR_NAME = 'MeSH Workbench'
PORTABLE_MARKER = 'portable.marker'


def program_dir(start=None):
    """The folder the application was launched from."""
    if start:
        return Path(start)
    # src/mesh_aop/paths.py -> up two for a source checkout; in the bundle this
    # is app/mesh_aop/paths.py and the parent of app/ is the program folder.
    return Path(__file__).resolve().parent.parent.parent


def is_portable(program=None):
    """True when this copy owns its folder and should keep everything in it.

    Decided by a marker file rather than by testing whether the folder happens
    to be writable: a writable install directory is not the same thing as a
    portable copy, and guessing would silently change where a user's results go.
    """
    root = Path(program) if program else program_dir()
    if (root / PORTABLE_MARKER).exists():
        return True
    # The bundle keeps the marker beside the launcher, one level above app/.
    return (root.parent / PORTABLE_MARKER).exists()


def user_root():
    """Per-user private state: settings, logs, caches, scratch.

    Never asked about - this is machinery, not the user's work.
    """
    if os.name == 'nt':
        base = os.environ.get('LOCALAPPDATA') or (Path.home() / 'AppData' / 'Local')
    elif sys.platform == 'darwin':
        base = Path.home() / 'Library' / 'Application Support'
    else:
        base = os.environ.get('XDG_DATA_HOME') or (Path.home() / '.local' / 'share')
    return Path(base) / APP_DIR_NAME


def default_results_dir():
    """Where a user's outputs go unless they choose otherwise.

    Somewhere they can find without being told a path: results are documents.
    """
    if os.name == 'nt':
        docs = Path.home() / 'Documents'
        # Redirected Documents (OneDrive) is normal on managed machines.
        onedrive = os.environ.get('OneDrive') or os.environ.get('OneDriveCommercial')
        if onedrive and (Path(onedrive) / 'Documents').is_dir():
            docs = Path(onedrive) / 'Documents'
    else:
        docs = Path.home() / 'Documents'
    return (docs if docs.is_dir() else Path.home()) / APP_DIR_NAME


def default_data_dir():
    """Where the downloaded archives and built databases go by default.

    Under the user profile, so it is private and removable with the account -
    but this is the one bulk location worth asking about, because it is tens of
    gigabytes and the system drive often cannot hold it.
    """
    return user_root() / 'data'


def bundled_reference_dir(program=None):
    """The reference corpus shipped with the program.

    This is part of the program, not user data: it is read and never written.
    The packaging step places it beside the packages in app/, while a source
    checkout keeps it at data/reference_processed - so both are checked. Getting
    this wrong is silent: the run simply reports the reference files missing.
    """
    root = Path(program) if program else program_dir()
    here = Path(__file__).resolve().parent           # .../mesh_aop

    # Walk up from the package itself, because that is the one thing whose
    # location is always known. The three layouts put the corpus in three
    # different places relative to it, and the caller's `program` is often the
    # working directory, which is not where anything is:
    #
    #   source checkout   <repo>/src/mesh_aop        -> <repo>/data/reference_processed
    #   Windows portable  <prog>/app/mesh_aop        -> <prog>/app/reference_processed
    #   Unix bundle       <prog>/app/src/mesh_aop    -> <prog>/reference_processed
    #
    # Only the middle one was found before, and only when the working directory
    # happened to be the program folder. The Unix bundles ship the corpus and
    # could never locate it, so "Use bundled reference data" had nothing to read.
    candidates = [root / 'app' / 'reference_processed',
                  root / 'reference_processed',
                  root / 'data' / 'reference_processed']
    for up in (here.parent, here.parent.parent, here.parent.parent.parent):
        candidates.append(up / 'reference_processed')
        candidates.append(up / 'data' / 'reference_processed')

    for candidate in candidates:
        if candidate.is_dir():
            return candidate
    return root / 'data' / 'reference_processed'


def default_user_data_dir(program=None):
    """Where downloaded archives and built databases go.

    Portable copies keep them in the folder; an installed copy puts them under
    the user profile, private to that account.
    """
    root = Path(program) if program else program_dir()
    return (root / 'data') if is_portable(root) else default_data_dir()


def default_user_results_dir(program=None):
    """Where a user's outputs go."""
    root = Path(program) if program else program_dir()
    return (root / 'results') if is_portable(root) else default_results_dir()


def config_path(program=None):
    """The settings file for this copy."""
    root = Path(program) if program else program_dir()
    if is_portable(root):
        return root / 'mesh_config.json'
    return user_root() / 'mesh_config.json'


def log_dir():
    """Run logs. Machinery: not mixed in with the user's results."""
    return user_root() / 'logs'


def scratch_dir():
    """Default parent for the ETL working folder."""
    return user_root()


def long_path(path):
    """A path Windows will accept even past its 260-character limit.

    Windows APIs reject paths longer than MAX_PATH unless they are absolute,
    fully backslash-separated, and prefixed with \\\\?\\ - which switches the
    kernel to the extended-length form. Without it the failure is a bare
    "No such file or directory" naming a path that plainly exists, which reads
    as a bug in the program rather than a limit of the filesystem.

    This is not hypothetical here: results land under Documents, Documents is
    frequently redirected into OneDrive with a long organisation name in it, and
    the figure filenames are descriptive. Three such segments is all it takes.

    Returns the path unchanged on any other platform, for a UNC path (which
    needs the different \\\\?\\UNC\\ form), or if it is already prefixed.
    """
    p = str(path)
    if os.name != 'nt' or p.startswith('\\\\'):
        return p
    absolute = os.path.abspath(p)
    return '\\\\?\\' + absolute if len(absolute) >= 250 else p


def free_gb(path):
    """Free space on the volume holding `path`, following it up to one that exists.

    A first-run dialog has to report this for a directory that has not been
    created yet, so walk up until something real is found.
    """
    import shutil
    p = Path(path).resolve()
    while not p.exists() and p != p.parent:
        p = p.parent
    try:
        return shutil.disk_usage(str(p)).free / 1e9
    except OSError:
        return None


def describe():
    """Everything resolved at once - for diagnostics and the first-run dialog."""
    portable = is_portable()
    return {
        'portable': portable,
        'program_dir': str(program_dir()),
        'config': str(config_path()),
        'user_root': str(user_root()),
        'default_results': str(default_results_dir()),
        'default_data': str(default_data_dir()),
        'logs': str(log_dir()),
    }


if __name__ == '__main__':
    for k, v in describe().items():
        print(f'  {k:<16} {v}')
