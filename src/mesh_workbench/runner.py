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
import threading
import time

# ---- events pushed onto the queue -------------------------------------
LOG = 'log'            # (LOG, text, level)      level: '', 'ok', 'warn', 'err'
TRANSIENT = 'tr'       # (TRANSIENT, text)       redraws in place
PHASE = 'phase'        # (PHASE, name)
PROGRESS = 'prog'      # (PROGRESS, done, total, unit)
DONE = 'done'          # (DONE, returncode, elapsed_seconds)

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

    # -- lifecycle ------------------------------------------------------
    def is_running(self):
        return self.proc is not None and self.proc.poll() is None

    def start(self, step):
        if self.is_running():
            raise RuntimeError('a run is already in progress')
        self._cancelled = False
        self._start = time.time()
        env = dict(os.environ)
        src = os.path.join(self.repo_dir, 'src')
        env['PYTHONPATH'] = src + os.pathsep + env.get('PYTHONPATH', '')
        env['PYTHONUNBUFFERED'] = '1'
        env['PYTHONIOENCODING'] = 'utf-8'
        flags = subprocess.CREATE_NEW_PROCESS_GROUP if os.name == 'nt' else 0
        self.proc = subprocess.Popen(
            [self.python_exe, '-u', '-m', 'mesh_aop.cli', '--step', step],
            cwd=self.repo_dir, env=env, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
            bufsize=0, creationflags=flags)
        self._t = threading.Thread(target=self._read, daemon=True)
        self._t.start()
        self.q.put((PHASE, 'Starting ' + step))

    def cancel(self):
        """Terminate the child and everything it spawned."""
        if not self.is_running():
            return
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
            self.q.put((DONE, -1 if self._cancelled else rc,
                        time.time() - (self._start or time.time())))

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
