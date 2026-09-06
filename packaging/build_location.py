# -*- coding: utf-8 -*-
"""
build_location.py - the one place a build is allowed to write.

Every build script used to decide this for itself, and they disagreed. The
Windows portable build preferred D: and fell back into the project folder; the
MSI wrote beside whatever tree it was pointed at; the Unix bundles defaulted to
~/Documents. Run them over a few weeks and the output is in four places, three
of them wrong, with two hundred files nobody can account for and two copies of
the same superseded tree on different drives.

So there is one answer now, and it is enforced rather than preferred:

    D:\\mesh_workbench_build

If D: is not attached, the build STOPS. It does not quietly pick somewhere
else. Falling back is what scattered the output in the first place, and the
fallback location was the cloud-synced working copy - so a build nobody asked
for would upload half a gigabyte of reproducible junk.

MESH_BUILD_OUT still overrides, for CI, which has no D: and must write into the
workspace. That is a deliberate answer from the caller rather than a guess.
"""

import os
import sys

#: The only location a build writes to unless it is told otherwise.
BUILD_ROOT = r'D:\mesh_workbench_build'

_ENV = 'MESH_BUILD_OUT'


def _writable(path):
    """Can this actually be written to? Existing is not the same thing.

    D: on another machine may be a read-only optical drive or a disconnected
    network mapping, and a build that discovered that halfway through would
    leave a half-assembled tree behind.
    """
    try:
        os.makedirs(path, exist_ok=True)
        probe = os.path.join(path, '.write_test')
        with open(probe, 'w') as fh:
            fh.write('')
        os.remove(probe)
        return True
    except OSError:
        return False


def resolve(explicit=None, purpose='build'):
    """Where to write, or exit saying why not.

    `explicit` is a --out the user passed; it is honoured without question,
    because naming a path is an answer, not a guess.
    """
    if explicit:
        if not _writable(explicit):
            sys.exit(f'--out is not writable: {explicit}')
        return os.path.abspath(explicit)

    env = os.environ.get(_ENV)
    if env:
        if not _writable(env):
            sys.exit(f'{_ENV} is not writable: {env}')
        return os.path.abspath(env)

    if _writable(BUILD_ROOT):
        return BUILD_ROOT

    sys.exit(
        f'\nThe {purpose} needs {BUILD_ROOT}, and it is not available.\n\n'
        f'  D: is where every artefact for this project lives. Attach it and\n'
        f'  run this again.\n\n'
        f'  Nothing was written. This used to fall back to a folder inside the\n'
        f'  cloud-synced working copy, which is how the output ended up in four\n'
        f'  places at once - so it stops instead.\n\n'
        f'  To build somewhere else on purpose:\n'
        f'      --out <path>            for this run\n'
        f'      set {_ENV}=<path>   for the session (CI uses this)\n')
