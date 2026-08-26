# -*- coding: utf-8 -*-
"""
runcontrol.py - pausing and stopping a run from outside it.

The obvious way to pause a child process is to suspend it: SIGSTOP on POSIX,
NtSuspendProcess on Windows. That was tried first and rejected for two reasons.

It is not reliable. On a managed Windows machine the suspend calls return
STATUS_SUCCESS, SuspendThread reports the expected suspend counts, and the
process carries on regardless - the security layer virtualises them. A pause
button that silently does nothing is worse than no pause button.

It is not safe. A hard suspend lands wherever the process happened to be, which
for this pipeline is often the middle of a multi-gigabyte SQLite transaction or
a half-written JSON. Leave it suspended long enough for a laptop to sleep or a
sync client to notice, and the pause has manufactured exactly the corruption
that integrity.py exists to detect.

So pausing is cooperative. The controller creates a file; the pipeline notices
it at points of its own choosing - between stages, between batches, at the top
of a loop - and waits there, with its transactions closed and its files
consistent. The cost is latency: a pause takes effect at the next checkpoint
rather than instantly, and the caller is told so rather than left wondering.

Stopping is not cooperative. A stop is a stop, and the controller kills the
process; this module only offers a way to notice one and exit tidily if the
checkpoint happens to come first.
"""

import os
import time

ENV_DIR = 'MESH_CONTROL_DIR'
PAUSE_FILE = 'pause.flag'
ABORT_FILE = 'abort.flag'

# How often a waiting process looks again. Short enough that Resume feels
# immediate, long enough to be free.
_POLL_SECONDS = 0.25

_control_dir = None
_announced = False


class RunAborted(RuntimeError):
    """Raised at a checkpoint when the controller has asked the run to stop."""


def set_control_dir(path):
    """Where the flag files live. The GUI sets this through the environment."""
    global _control_dir
    _control_dir = str(path) if path else None
    return _control_dir


def control_dir():
    global _control_dir
    if _control_dir is None:
        _control_dir = os.environ.get(ENV_DIR) or ''
    return _control_dir or None


def _flag(name):
    d = control_dir()
    return os.path.join(d, name) if d else None


def pause_requested():
    p = _flag(PAUSE_FILE)
    return bool(p) and os.path.exists(p)


def abort_requested():
    p = _flag(ABORT_FILE)
    return bool(p) and os.path.exists(p)


def request_pause(on=True):
    """Controller side: ask the run to pause at its next checkpoint."""
    p = _flag(PAUSE_FILE)
    if not p:
        return False
    try:
        if on:
            os.makedirs(os.path.dirname(p), exist_ok=True)
            with open(p, 'w', encoding='utf-8') as fh:
                fh.write(time.strftime('%Y-%m-%d %H:%M:%S'))
        elif os.path.exists(p):
            os.remove(p)
        return True
    except OSError:
        return False


def request_abort():
    p = _flag(ABORT_FILE)
    if not p:
        return False
    try:
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, 'w', encoding='utf-8') as fh:
            fh.write(time.strftime('%Y-%m-%d %H:%M:%S'))
        return True
    except OSError:
        return False


def clear():
    for name in (PAUSE_FILE, ABORT_FILE):
        p = _flag(name)
        if p and os.path.exists(p):
            try:
                os.remove(p)
            except OSError:
                pass


def checkpoint(label=''):
    """A safe point to be paused or stopped at.

    Cheap when nothing is happening - one os.path.exists against a path in the
    user's own temp folder - so it can sit inside a loop without being felt.
    Prints marker lines the Workbench recognises, and readable prose for anyone
    running from a terminal.

    Raises RunAborted if a stop has been requested.
    """
    global _announced
    if abort_requested():
        raise RunAborted(f'stopped by request{" at " + label if label else ""}')
    if not pause_requested():
        return False

    where = f' ({label})' if label else ''
    print(f'\n[RUN-PAUSED]{where}', flush=True)
    print(f'  Paused at a safe point{where}. Nothing is being computed and no '
          f'file is half-written.', flush=True)
    print('  Press Resume in the Workbench, or delete the pause flag:', flush=True)
    print(f'    {_flag(PAUSE_FILE)}', flush=True)
    _announced = True

    waited = 0.0
    while pause_requested():
        if abort_requested():
            print('[RUN-RESUMED]', flush=True)
            raise RunAborted('stopped while paused')
        time.sleep(_POLL_SECONDS)
        waited += _POLL_SECONDS

    print(f'[RUN-RESUMED]', flush=True)
    print(f'  Resumed after {waited/60:.1f} minutes paused.', flush=True)
    return True
