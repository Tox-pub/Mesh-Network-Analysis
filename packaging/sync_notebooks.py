# -*- coding: utf-8 -*-
"""
sync_notebooks.py - regenerate src/mesh_aop_notebooks/ from the real modules.

The notebooks are a mirror of the package for people who read or run it in
Colab. Each one is the module's source in a single code cell - that is all they
have ever been - so keeping them current is mechanical, and doing it by hand is
why they went eight releases out of date.

That drift was not cosmetic. The snapshot shipped inside all three installers
still carried the empty-api_key bug that made every retrieval fail with HTTP
400, the master-PMID set that held four gigabytes, and four figures that had
been removed. A second copy of the code in the download, with the bugs still in
it, is worse than no copy.

    python packaging/sync_notebooks.py           rewrite them
    python packaging/sync_notebooks.py --check   fail if any is stale (for CI)
"""

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
PKG = os.path.join(REPO, 'src', 'mesh_aop')
NB = os.path.join(REPO, 'src', 'mesh_aop_notebooks')

# Kept as a plain .py beside the notebooks rather than turned into one. It is
# eight thousand lines of vocabulary with no logic in it, the notebooks import
# it as a module, and a notebook of it would be a quarter-megabyte of JSON
# nobody will ever open.
AS_PLAIN_PY = {'mesh_stop_words'}

METADATA = {
    'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
    'language_info': {'name': 'python'},
}


def source_cell(text):
    """A notebook's single code cell: the module, line by line, as nbformat wants it."""
    lines = text.splitlines(keepends=True)
    return {'cell_type': 'code', 'execution_count': None, 'metadata': {},
            'outputs': [], 'source': lines}


def notebook_for(text):
    return {'cells': [source_cell(text)], 'metadata': METADATA,
            'nbformat': 4, 'nbformat_minor': 4}


def modules():
    """Every module in the package, by name, in a stable order.

    __init__ included: it is where the lazy __getattr__ lives, so it is the
    file that decides what importing mesh_aop actually costs. Leaving it out
    would drop the one module that explains the package's shape.
    """
    return sorted(f[:-3] for f in os.listdir(PKG) if f.endswith('.py'))


def write_if_changed(path, text, check_only):
    """Returns True when the file was (or would be) rewritten."""
    existing = None
    if os.path.exists(path):
        with open(path, 'r', encoding='utf-8') as fh:
            existing = fh.read()
    if existing == text:
        return False
    if not check_only:
        with open(path, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write(text)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true',
                    help='report what is stale and exit non-zero, changing nothing')
    a = ap.parse_args()

    os.makedirs(NB, exist_ok=True)
    changed, names = [], modules()

    for name in names:
        with open(os.path.join(PKG, name + '.py'), 'r', encoding='utf-8') as fh:
            src = fh.read()
        if name in AS_PLAIN_PY:
            target = os.path.join(NB, name + '.py')
            if write_if_changed(target, src, a.check):
                changed.append(name + '.py')
            continue
        target = os.path.join(NB, name + '.ipynb')
        text = json.dumps(notebook_for(src), indent=1, ensure_ascii=False) + '\n'
        if write_if_changed(target, text, a.check):
            changed.append(name + '.ipynb')

    # Anything left over describes a module that no longer exists.
    keep = {n + ('.py' if n in AS_PLAIN_PY else '.ipynb') for n in names}
    allowed = keep | {'environment.yml', 'pyproject.toml', '__pycache__'}
    orphans = [f for f in sorted(os.listdir(NB))
               if f not in allowed and (f.endswith('.ipynb') or f.endswith('.py'))]
    for f in orphans:
        if not a.check:
            os.remove(os.path.join(NB, f))
        changed.append(f + ' (removed - no such module)')

    print(f'{len(names)} modules: {len(names) - len(AS_PLAIN_PY)} notebooks, '
          f'{len(AS_PLAIN_PY)} kept as .py')
    if a.check:
        if changed:
            print('\nSTALE - run: python packaging/sync_notebooks.py')
            for c in changed:
                print('  ', c)
            return 1
        print('every notebook matches its module')
        return 0
    if changed:
        print(f'\nrewrote {len(changed)}:')
        for c in changed:
            print('  ', c)
    else:
        print('already up to date')
    return 0


if __name__ == '__main__':
    sys.exit(main())
