# -*- coding: utf-8 -*-
"""
run_all.py - run every check in this folder and report what failed.

    python tests/run_all.py              everything that needs no built data
    python tests/run_all.py --all        everything, including the slow ones

These are standalone scripts rather than pytest cases, deliberately: each one is
readable on its own and prints what it checked, which is what makes a failure
useful six months from now. They exit non-zero when something is wrong, so CI
needs nothing more than this file.

Some need artefacts that only exist after a real run - a built network, the
master annotation database. Those are marked NEEDS_DATA and skipped unless the
files are there, because a check that cannot run is not a check that passed and
should not be reported as one.
"""

import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

# (script, what it covers, artefacts it needs)
SUITES = [
    ('test_paths.py',
     'every path the pipeline resolves lands where it should',
     []),
    ('test_ranking.py',
     'the benchmark ranks positives exactly as sorting the pool would have',
     []),
    ('test_benchmark_layout.py',
     'the benchmark folder groups its inputs and its outputs',
     []),
    ('test_citation_fallback.py',
     'secondary analysis generates the citation data it needs when none was retrieved',
     []),
    ('test_seed_matching.py',
     'a node query returns that heading and not the longer ones containing it',
     []),
    ('test_output_schemas.py',
     'each network JSON carries the attributes and types its stage promises',
     ['data/processed/DAC_Mesh_final_network_with_relevance.json',
      'data/processed/DAC_Mesh_consensus_lcc_network.json']),
    ('test_subgraph_centralities.py',
     'the six subgraph attributes appear exactly when the flag is set',
     ['data/processed/DAC_Mesh_consensus_lcc_network.json']),
    ('test_network_figure.py',
     'Figure 7 draws at every network size and survives the awkward cases',
     ['data/processed/DAC_Mesh_final_network_with_relevance.json',
      'data/raw/aop_annotations_master.csv']),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--all', action='store_true',
                    help='fail rather than skip when a suite needs data that is absent')
    a = ap.parse_args()

    results = []
    for script, what, needs in SUITES:
        absent = [n for n in needs if not os.path.exists(os.path.join(REPO, n))]
        print('\n' + '=' * 70)
        print(f'{script}  -  {what}')
        print('=' * 70)
        if absent and not a.all:
            print(f'  SKIPPED - needs data that is not here: {", ".join(absent)}')
            results.append((script, 'skipped'))
            continue
        started = time.time()
        proc = subprocess.run([sys.executable, os.path.join(HERE, script)], cwd=REPO)
        took = time.time() - started
        state = 'passed' if proc.returncode == 0 else 'FAILED'
        print(f'  -> {state} in {took:,.1f}s')
        results.append((script, state))

    print('\n' + '=' * 70)
    for script, state in results:
        print(f'  {state.upper():<8} {script}')
    failed = [s for s, st in results if st == 'FAILED']
    skipped = [s for s, st in results if st == 'skipped']
    print('=' * 70)
    if skipped:
        print(f'{len(skipped)} skipped for want of built data - '
              f'run with --all to make that a failure instead.')
    print('EVERYTHING PASSED' if not failed else f'{len(failed)} SUITE(S) FAILED')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
