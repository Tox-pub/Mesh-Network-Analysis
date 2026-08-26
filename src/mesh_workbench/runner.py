# -*- coding: utf-8 -*-
"""
runner.py - run a pipeline step in a child process and report progress.

Two decisions worth stating.

SUBPROCESS, NOT A THREAD. The pipeline is long-running, third-party-heavy and
occasionally memory-hungry. In-process it would block the event loop and a hard
failure would take the window with it. As a child process the UI stays live, and
Cancel has something real to terminate.

CARRIAGE RETURNS ARE NOT NEWLINES. tqdm redraws its bar by writing \\r and
overwriting the line. Reading with readline() therefore blocks until a bar
finishes, and progress arrives in one useless burst at the end. This reader takes
raw chunks and splits on both terminators: \\n gives a permanent log line, \\r
gives a transient one that replaces the previous transient. That is what makes a
bar animate instead of freezing.

Progress is inferred from the pipeline's own stdout, so nothing in mesh_aop has
to change to support the UI.
"""

import os
import queue
import re
import subprocess
import sys
import tempfile
import threading
import time

from mesh_aop import runcontrol

# ---- events pushed onto the queue -------------------------------------
LOG = 'log'            # (LOG, text, level)      level: '', 'ok', 'warn', 'err'
TRANSIENT = 'tr'       # (TRANSIENT, text)       redraws in place
PHASE = 'phase'        # (PHASE, name)
PROGRESS = 'prog'      # (PROGRESS, done, total, unit)
DONE = 'done'          # (DONE, returncode, elapsed_seconds)
PAUSED = 'paused'      # (PAUSED, annotation_file_path)
AT_CHECKPOINT = 'ckpt' # (AT_CHECKPOINT, True|False)   reached / left a safe point

# The pipeline prints this when it stops for AOP annotation. A marker line
# rather than pattern-matching the prose around it, so rewording the on-screen
# instructions can never silently break the dialog that depends on them.
_PAUSED_MARK = re.compile(r'^\[PAUSED-FOR-ANNOTATION\]\s*(.+?)\s*$')
# Printed by runcontrol.checkpoint() when the run actually stops and starts
# again, which is what lets the UI distinguish 'pausing' from 'paused'.
_AT_PAUSE = re.compile(r'^\[RUN-PAUSED\]')
_AT_RESUME = re.compile(r'^\[RUN-RESUMED\]')

_PCT_OF = re.compile(r'(\d+)%\s*\(\s*([\d,]+)\s*/\s*([\d,]+)\s*\)')
_TQDM = re.compile(r'(\d+)%\|.*?\|\s*([\d,]+)/([\d,]+)')
_BARE_PCT = re.compile(r'(?<![\d.])(\d{1,3})%')
_PHASE_ARROW = re.compile(r'>>>\s*STARTING:\s*(.+?)\s*$')
# `<<< Title >>>` is a phase header, but the pipeline also prints bare
# `<<<<<<<<>>>>>>>>` rules as separators - require at least one letter inside.
_PHASE_ANGLE = re.compile(r'^<<<\s*(?![<>])(.*[A-Za-z].*?)\s*>>>\s*$')
_PASS = re.compile(r'^\s*(Pass \d+/\d+:.*?)\s*\.\.\.\s*$')
_DONE_OK = re.compile(r'PIPELINE COMPLETED SUCCESSFULLY in ([\d.]+) minutes')
_COUNTED = re.compile(r'^\s*(?:scored pool|Scanned|scanned)\D*([\d,]+)')


def _n(s):
    return int(s.replace(',', ''))


def classify(line):
    """Severity for colouring, from the pipeline's own console conventions."""
    if 'Traceback' in line or 'Error' in line or 'ERROR' in line:
        return 'err'
    if '[!]' in line or 'WARNING' in line or 'Warning' in line:
        return 'warn'
    if '[+]' in line or 'SUCCESS' in line or 'COMPLETED' in line:
        return 'ok'
    return ''


class PipelineRunner:
    """Runs `python -m mesh_aop.cli --step X` and streams its output."""

    def __init__(self, repo_dir, python_exe=None):
        self.repo_dir = repo_dir
        self.python_exe = python_exe or sys.executable
        self.q = queue.Queue()
        self.proc = None
        self._t = None
        self._start = None
        self._cancelled = False
        self._paused = False
        self._pause_asked = False
        self._paused_at = None
        self._paused_total = 0.0
        # A private directory for the pause/stop flags. Under the user's
        # own temp folder, not the results tree: it is machinery, and a
        # stray flag file left among someone's outputs would be baffling.
        self._control = None

    # -- lifecycle ------------------------------------------------------
    def is_running(self):
        return self.proc is not None and self.proc.poll() is None

    def start(self, step, extra=None):
        """Run one pipeline step. `extra` carries per-run switches such as
        --build-database, which are deliberately not stored in the config."""
        if self.is_running():
            raise RuntimeError('a run is already in progress')
        self._cancelled = False
        self._paused = False
        self._pause_asked = False
        self._paused_at = None
        self._paused_total = 0.0
        self._start = time.time()
        self._control = tempfile.mkdtemp(prefix='mesh-run-')
        env = dict(os.environ)
        env[runcontrol.ENV_DIR] = self._control
        src = os.path.join(self.repo_dir, 'src')
        env['PYTHONPATH'] = src + os.pathsep + env.get('PYTHONPATH', '')
        env['PYTHONUNBUFFERED'] = '1'
        env['PYTHONIOENCODING'] = 'utf-8'
        flags = subprocess.CREATE_NEW_PROCESS_GROUP if os.name == 'nt' else 0
        self.proc = subprocess.Popen(
            [self.python_exe, '-u', '-m', 'mesh_aop.cli', '--step', step] + list(extra or []),
            cwd=self.repo_dir, env=env, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
            bufsize=0, creationflags=flags)
        self._t = threading.Thread(target=self._read, daemon=True)
        self._t.start()
        self.q.put((PHASE, 'Starting ' + step))

    # -- pause / resume --------------------------------------------------
    #
    # Cooperative, not a process suspend. Two reasons, both learned the hard way
    # on the machine this was built for.
    #
    # Suspending is not reliable. NtSuspendProcess returned STATUS_SUCCESS and
    # per-thread SuspendThread returned the expected suspend counts, and the
    # child carried on computing regardless - the managed security layer
    # virtualises them. A Pause button that silently does nothing is worse than
    # none at all.
    #
    # Suspending is not safe either. It freezes the process wherever it happens
    # to be, which here is often the middle of a multi-gigabyte SQLite
    # transaction. Leave that suspended long enough for the machine to sleep and
    # the pause has produced exactly the corruption the repair tool exists for.
    #
    # So the pipeline is asked to stop, and stops at a point of its own
    # choosing. The cost is latency, which the UI states rather than hides.

    def can_pause(self):
        return self._control is not None

    def is_paused(self):
        """True once the child has actually reached a checkpoint and stopped.

        Deliberately not "pause was requested": the two are different for as
        long as the current batch takes, and reporting the request as the fact
        would put "Paused" on screen next to a progress bar still moving.
        """
        return self._paused and self.is_running()

    def pause_pending(self):
        return self._pause_asked and not self._paused and self.is_running()

    def pause(self):
        """Ask the run to stop at its next safe point."""
        if not self.is_running() or self._pause_asked:
            return False
        if not self._control_write(runcontrol.PAUSE_FILE):
            return False
        self._pause_asked = True
        self.q.put((LOG, '--- pause requested: the run will stop at its next '
                         'safe point, which may take a few minutes ---', 'warn'))
        return True

    def resume(self):
        """Let a paused run continue."""
        if not self._pause_asked:
            return True
        self._control_remove(runcontrol.PAUSE_FILE)
        self._pause_asked = False
        return True

    def _control_write(self, name):
        if not self._control:
            return False
        try:
            os.makedirs(self._control, exist_ok=True)
            with open(os.path.join(self._control, name), 'w', encoding='utf-8') as fh:
                fh.write(time.strftime('%Y-%m-%d %H:%M:%S'))
            return True
        except OSError as exc:
            self.q.put((LOG, f'--- could not request a pause: {exc} ---', 'err'))
            return False

    def _control_remove(self, name):
        if not self._control:
            return
        try:
            path = os.path.join(self._control, name)
            if os.path.exists(path):
                os.remove(path)
        except OSError:
            pass

    def _mark_paused(self, on):
        """Called when the child reports it has reached, or left, a checkpoint."""
        if on and not self._paused:
            self._paused = True
            self._paused_at = time.time()
        elif not on and self._paused:
            # Time spent paused is not time spent working, and an elapsed clock
            # that counts it makes every duration in the log a lie.
            self._paused_total += time.time() - (self._paused_at or time.time())
            self._paused = False
            self._paused_at = None

    def elapsed(self):
        """Working seconds, with paused time taken out."""
        if not self._start:
            return 0.0
        held = self._paused_total
        if self._paused and self._paused_at:
            held += time.time() - self._paused_at
        return max(0.0, time.time() - self._start - held)

    def cancel(self):
        """Terminate the child and everything it spawned.

        A suspended process cannot be asked to exit, so it is woken first;
        otherwise taskkill's request would sit in a queue that is never read.
        """
        if not self.is_running():
            return
        # Ask it to stop, and release any pause - a run waiting at a
        # checkpoint would otherwise never look up to see the request.
        self._control_write(runcontrol.ABORT_FILE)
        self._control_remove(runcontrol.PAUSE_FILE)
        self._pause_asked = False
        self._cancelled = True
        try:
            if os.name == 'nt':
                subprocess.run(['taskkill', '/F', '/T', '/PID', str(self.proc.pid)],
                               capture_output=True)
            else:
                self.proc.terminate()
        except Exception:
            try:
                self.proc.kill()
            except Exception:
                pass

    # -- reader ---------------------------------------------------------
    def _read(self):
        """Chunked read, then split on real newlines and bar redraws separately.

        os.read returns as soon as anything is available, where a buffered
        read(n) would block until n bytes arrive and stall the log. CRLF is
        normalised first: on Windows every print() ends '\\r\\n', so treating a
        bare '\\r' as a redraw without collapsing CRLF first classifies every
        ordinary line as a transient and the log stays empty.
        """
        fd = self.proc.stdout.fileno()
        buf = ''
        try:
            while True:
                try:
                    raw = os.read(fd, 8192)
                except OSError:
                    break
                if not raw:
                    break
                buf += raw.decode('utf-8', 'replace').replace('\r\n', '\n')
                while '\n' in buf:
                    line, buf = buf.split('\n', 1)
                    self._emit_split(line)
                # a bar redraw arrives without a newline; show the newest frame
                if '\r' in buf:
                    *frames, buf = buf.split('\r')
                    for f in frames:
                        if f.strip():
                            self._emit(f.rstrip(), permanent=False)
        except Exception as exc:                       # pragma: no cover
            self.q.put((LOG, 'reader stopped: %s' % exc, 'err'))
        finally:
            if buf.strip():
                self._emit(buf.rstrip(), permanent=True)
            rc = self.proc.wait() if self.proc else -1
            self._cleanup_control()
            self.q.put((DONE, -1 if self._cancelled else rc, self.elapsed()))

    def _emit_split(self, line):
        """A completed line may still hold earlier bar frames before its text."""
        if '\r' in line:
            *frames, final = line.split('\r')
            for f in frames:
                if f.strip():
                    self._emit(f.rstrip(), permanent=False)
            line = final
        if line.strip():
            self._emit(line.rstrip(), permanent=True)

    def _emit(self, text, permanent):
        # progress first - a bar line carries no other information
        m = _PCT_OF.search(text) or _TQDM.search(text)
        if m:
            self.q.put((PROGRESS, _n(m.group(2)), _n(m.group(3)), 'items'))
            if not permanent:
                self.q.put((TRANSIENT, text))
                return
        elif not permanent:
            m2 = _BARE_PCT.search(text)
            if m2:
                self.q.put((PROGRESS, int(m2.group(1)), 100, '%'))
            self.q.put((TRANSIENT, text))
            return

        pause = _PAUSED_MARK.match(text.strip())
        if pause:
            self.q.put((PAUSED, pause.group(1)))
            return
        if _AT_PAUSE.match(text.strip()):
            self.q.put((AT_CHECKPOINT, True))
            return
        if _AT_RESUME.match(text.strip()):
            self.q.put((AT_CHECKPOINT, False))
            return

        m = _PHASE_ARROW.search(text) or _PHASE_ANGLE.search(text) or _PASS.match(text)
        if m:
            self.q.put((PHASE, m.group(1).strip()))

        self.q.put((LOG, text, classify(text)))

        if _DONE_OK.search(text):
            self.q.put((PROGRESS, 1, 1, 'done'))

    # -- draining -------------------------------------------------------
    def drain(self, limit=400):
        """Pop up to `limit` events. Called from the UI timer, never blocks."""
        out = []
        for _ in range(limit):
            try:
                out.append(self.q.get_nowait())
            except queue.Empty:
                break
        return out

    def _cleanup_control(self):
        """Remove the flag directory once the run it belonged to has ended."""
        if not self._control:
            return
        try:
            import shutil
            shutil.rmtree(self._control, ignore_errors=True)
        finally:
            self._control = None
