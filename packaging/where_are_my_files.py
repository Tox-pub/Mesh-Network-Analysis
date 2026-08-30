# -*- coding: utf-8 -*-
"""
where_are_my_files.py - what this project expects, and what is actually on disk.

Written because "the databases still have the old prefix" is a question that
cannot be answered by reading code: it depends on what a previous run left
behind, and where. This prints both sides so they can be compared.

    python packaging/where_are_my_files.py
    python packaging/where_are_my_files.py --config /path/to/mesh_config.json
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), 'src'))


def human(n):
    for unit in ('B', 'KB', 'MB', 'GB'):
        if n < 1024 or unit == 'GB':
            return f'{n:,.0f} {unit}' if unit == 'B' else f'{n:,.1f} {unit}'
        n /= 1024


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default=None)
    a = ap.parse_args()

    from mesh_aop.config_parser import MeshConfig
    c = MeshConfig(config_path=a.config)

    print(f'\nproject prefix : {c.prefix}')
    print(f'reference data : {"ON" if c.use_reference_data else "off"}')
    print(f'\nfolders')
    for label, d in (('raw (downloads)', c.active_raw_dir),
                     ('processed', c.active_source_dir),
                     ('networks', c.networks_dir),
                     ('databases', c.databases_dir),
                     ('results', c.results_dir)):
        print(f'  {label:<18} {d}')

    print('\n' + '=' * 78)
    print('WHAT THIS PROJECT EXPECTS')
    print('=' * 78)
    missing = []
    for key in sorted(c.files):
        p = c.files[key]
        there = os.path.exists(p)
        size = human(os.path.getsize(p)) if there and os.path.isfile(p) else ''
        if not there:
            missing.append(key)
        print(f'  {"[+]" if there else "[ ]":<4} {key:<18} {p.name:<46} {size:>10}')

    print('\n' + '=' * 78)
    print('WHAT IS ACTUALLY ON DISK, BY PREFIX')
    print('=' * 78)
    seen = {}
    for folder in {c.active_raw_dir, c.active_source_dir, c.networks_dir,
                   c.databases_dir}:
        if not os.path.isdir(folder):
            continue
        for name in sorted(os.listdir(folder)):
            full = os.path.join(folder, name)
            if not os.path.isfile(full) or not name.endswith(('.db', '.json', '.csv')):
                continue
            # Everything before the first known suffix is the prefix.
            tag = 'no prefix'
            # Longest first: '_pmids.db' is a suffix of '_cleaned_pmids.db', so
            # checking it earlier reads DAC_Mesh_cleaned_pmids.db as belonging
            # to a project called "DAC_Mesh_cleaned".
            for marker in sorted(('_pmids.db', '_cleaned_pmids.db', '_mean_relevancy.db',
                                  '_full_network_data.json', '_consensus_lcc_network.json',
                                  '_final_network_with_relevance.json',
                                  '_glf_optimal_subgraph.json', '_sa_optimal_subgraph.json',
                                  '_optimization_history.json'), key=len, reverse=True):
                if name.endswith(marker):
                    tag = name[:-len(marker)]
                    break
            seen.setdefault(tag, []).append((folder, name, os.path.getsize(full)))

    for tag in sorted(seen):
        mine = (tag == c.prefix)
        note = '  <-- THIS PROJECT' if mine else (
            '  (shared by every project)' if tag == 'no prefix' else
            '  (another project, or a previous prefix)')
        print(f'\n  prefix "{tag}"{note}')
        for folder, name, size in seen[tag]:
            where = os.path.basename(folder) or str(folder)
            print(f'      {where:<12} {name:<48} {human(size):>10}')

    print('\n' + '=' * 78)
    if missing:
        print(f'{len(missing)} of this project\'s files are not built yet:')
        print('   ' + ', '.join(missing))
    else:
        print('Everything this project expects is present.')
    other = [t for t in seen if t not in (c.prefix, 'no prefix')]
    if other:
        print(f'\nFiles under {len(other)} other prefix(es) are also on disk: '
              f'{", ".join(sorted(other))}.')
        print('They belong to earlier runs and are never read by this one. '
              'Delete them\nfrom the Database screen if you no longer want them.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
