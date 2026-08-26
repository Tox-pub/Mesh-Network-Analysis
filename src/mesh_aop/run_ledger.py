# -*- coding: utf-8 -*-
"""
run_ledger.py - the count of what happened at every stage of a run.

The pipeline computes a great many numbers - records retrieved per generation,
terms before and after filtering, what each optimiser kept, what the consensus
discarded - prints them once, and throws them away. Reconstructing them
afterwards means re-running the stage. This keeps them.

One ledger per file prefix, written into the results folder as a semicolon
delimited table:

    stage;metric;value;recorded;note

Semicolon rather than comma because MeSH headings contain commas as a matter of
course ("Dermatitis, Allergic Contact"), and the annotation files already use
this convention. It opens in Excel, reads fine in a text editor, and loads with
one call.

Stages are written at different times - retrieval today, the network next week -
so recording MERGES: a stage re-run overwrites its own rows and leaves the
others alone. Every row carries the timestamp it was written, which is what
makes a mixed-date ledger legible rather than misleading.
"""

import csv
import os
from datetime import datetime

from . import paths as _paths

FIELDS = ('stage', 'metric', 'value', 'recorded', 'note')

# Canonical order, so the file reads top to bottom as the pipeline ran rather
# than in whatever order the stages happened to execute.
STAGE_ORDER = (
    'run',
    'identification',
    'screening',
    'pruning',
    'included',
    'unipartite',
)


class RunLedger:
    """Accumulates counts across stages, keyed on (stage, metric)."""

    def __init__(self, path):
        self.path = str(path)
        self.rows = {}
        self._order = []
        self.load()

    # ------------------------------------------------------------------ io
    def load(self):
        if not os.path.exists(self.path):
            return
        try:
            with open(_paths.long_path(self.path), 'r', encoding='utf-8-sig', newline='') as fh:
                for row in csv.DictReader(fh, delimiter=';'):
                    key = (row.get('stage', ''), row.get('metric', ''))
                    if not key[0]:
                        continue
                    self.rows[key] = row
                    self._order.append(key)
        except (OSError, csv.Error):
            # A malformed ledger must never stop a run; it is a record of the
            # work, not an input to it.
            self.rows, self._order = {}, []

    def save(self):
        os.makedirs(_paths.long_path(os.path.dirname(self.path) or '.'), exist_ok=True)
        try:
            with open(_paths.long_path(self.path), 'w', encoding='utf-8', newline='') as fh:
                w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter=';')
                w.writeheader()
                for key in self._sorted_keys():
                    w.writerow(self.rows[key])
            return True
        except OSError as exc:
            print(f"  [!] Could not write the run ledger: {exc}")
            return False

    def _sorted_keys(self):
        def rank(key):
            stage = key[0]
            idx = STAGE_ORDER.index(stage) if stage in STAGE_ORDER else len(STAGE_ORDER)
            return (idx, self._order.index(key) if key in self._order else 0)
        return sorted(self.rows, key=rank)

    # -------------------------------------------------------------- recording
    def record(self, stage, metric, value, note=''):
        """Set one count. Re-recording replaces the value and the timestamp."""
        key = (stage, metric)
        if key not in self.rows:
            self._order.append(key)
        self.rows[key] = {
            'stage': stage,
            'metric': metric,
            'value': '' if value is None else str(value),
            'recorded': datetime.now().strftime('%Y-%m-%d %H:%M'),
            'note': note,
        }

    def record_many(self, stage, pairs, note=''):
        for metric, value in pairs:
            self.record(stage, metric, value, note)

    def clear_stage(self, stage):
        """Drop a stage's rows before rewriting it.

        A re-run that produces fewer metrics than the last one would otherwise
        leave the surplus behind, and stale rows in a provenance record are
        worse than absent ones.
        """
        for key in [k for k in self.rows if k[0] == stage]:
            del self.rows[key]
            if key in self._order:
                self._order.remove(key)

    # ----------------------------------------------------------------- reading
    def get(self, stage, metric, default=None):
        row = self.rows.get((stage, metric))
        return row['value'] if row else default

    def get_int(self, stage, metric, default=None):
        raw = self.get(stage, metric)
        try:
            return int(float(raw))
        except (TypeError, ValueError):
            return default

    def snapshot_stage(self, stage):
        """{metric: value} for one stage, taken before it is cleared and rewritten.

        A recorded count is a cheaper source than the artefact it came from -
        the unfiltered network is a several-hundred-megabyte JSON, and reading
        it back to learn a number already in this file is pure waste.
        """
        return {metric: row['value'] for (st, metric), row in self.rows.items()
                if st == stage}

    def stage_dates(self):
        """When each stage was last written - what makes a mixed run legible."""
        out = {}
        for (stage, _), row in self.rows.items():
            stamp = row.get('recorded', '')
            if stamp and stamp > out.get(stage, ''):
                out[stage] = stamp
        return out

    def __len__(self):
        return len(self.rows)


def ledger_path(results_dir, prefix):
    """One ledger per file prefix, beside the results it describes."""
    return os.path.join(str(results_dir), f'{prefix}_run_ledger.csv')


def open_ledger(results_dir, prefix):
    return RunLedger(ledger_path(results_dir, prefix))
