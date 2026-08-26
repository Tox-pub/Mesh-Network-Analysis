# -*- coding: utf-8 -*-
"""
verify_unix_bundle.py - check a Linux/macOS bundle without a Linux or a Mac.

A bundle assembled on Windows cannot be run here, and that is the honest limit
of building it here. But most of the ways it can be wrong are visible from the
outside, and every one of these has actually happened during development:

  * the executable bit lost, because NTFS cannot store it - the launcher and
    the interpreter both fail with "Permission denied";
  * a wheel built for the wrong Python or the wrong architecture, which
    installs on nothing;
  * a dependency with no wheel at all, so the offline install reaches for PyPI
    and fails on a machine with no network;
  * no setuptools, so the one source-only dependency cannot be built;
  * a launcher written with CRLF, which fails at the shebang.

This checks all of them. It is not a substitute for running the thing, and it
does not pretend to be - it is what can be established before that.

    python packaging/verify_unix_bundle.py <bundle.tar.gz>
"""

import argparse
import os
import posixpath
import re
import sys
import tarfile

PY_SERIES = '3.12'
PASS, FAIL = [], []


def check(ok, msg, detail=''):
    (PASS if ok else FAIL).append(msg)
    print(('  [OK] ' if ok else '  [XX] ') + msg + (('  -- ' + detail) if detail and not ok else ''))
    return ok


def verify(path):
    print(f'\n=== {os.path.basename(path)} ({os.path.getsize(path)/1e6:,.0f} MB) ===\n')
    with tarfile.open(path, 'r:gz') as tf:
        members = tf.getmembers()
        by_name = {m.name: m for m in members}
        names = list(by_name)
        root = names[0].split('/')[0]

        def get(rel):
            return by_name.get(posixpath.join(root, rel))

        # ---------------------------------------------------------- layout
        print('-- layout')
        for rel in (f'python/bin/python{PY_SERIES}',
                    f'python/lib/python{PY_SERIES}/tkinter/__init__.py',
                    'app/src/mesh_aop/cli.py',
                    'app/src/mesh_workbench/app.py',
                    'wheels', 'requirements.txt', 'README.txt',
                    'MeSH Workbench', 'mesh-pipeline'):
            check(get(rel) is not None, f'present: {rel}')

        # ------------------------------------------------------ permissions
        print('\n-- permissions (NTFS cannot store these; the builder sets them)')
        for rel in ('MeSH Workbench', 'mesh-pipeline', f'python/bin/python{PY_SERIES}'):
            m = get(rel)
            if m:
                check(m.mode & 0o111, f'executable: {rel}', f'mode {m.mode:o}')
        sos = [m for m in members if m.name.endswith(('.so', '.dylib'))]
        bad_so = [m.name for m in sos if not m.mode & 0o111]
        check(not bad_so, f'all {len(sos)} shared objects are executable',
              f'{len(bad_so)} are not')
        docs = [m for m in members
                if m.name.endswith(('.txt', '.py', '.json')) and m.mode & 0o111]
        check(not docs, 'ordinary files are not marked executable',
              f'{len(docs)} are')

        # ------------------------------------------------------- launchers
        print('\n-- launchers')
        for rel in ('MeSH Workbench', 'mesh-pipeline'):
            m = get(rel)
            if not m:
                continue
            body = tf.extractfile(m).read()
            check(body.startswith(b'#!'), f'{rel}: has a shebang')
            check(b'\r\n' not in body, f'{rel}: LF endings, not CRLF',
                  'CRLF fails at the shebang on Unix')
            check(b'PYTHONPATH' in body, f'{rel}: puts the app on PYTHONPATH')

        launch = tf.extractfile(get('MeSH Workbench')).read().decode()
        check('com.apple.quarantine' in launch,
              'the launcher clears the macOS quarantine flag itself')
        check('--no-index' in launch,
              'the first-run install is offline (--no-index)')

        # ---------------------------------------------------------- wheels
        print('\n-- wheels')
        wheels = [posixpath.basename(n) for n in names
                  if '/wheels/' in n and n.endswith('.whl')]
        sdists = [posixpath.basename(n) for n in names
                  if '/wheels/' in n and n.endswith('.tar.gz')]
        check(bool(wheels), f'{len(wheels)} wheels present')

        # Every wheel must suit this interpreter: either pure Python, or built
        # for cp312 on the target architecture.
        want_abi = 'cp' + PY_SERIES.replace('.', '')
        wrong = []
        arch_tags = set()
        for w in wheels:
            parts = w[:-4].split('-')
            pytag, abi, plat = parts[-3], parts[-2], parts[-1]
            arch_tags.add(plat)
            pure = pytag.startswith('py') and abi == 'none' and plat == 'any'
            if not pure and want_abi not in (pytag, abi):
                wrong.append(w)
        check(not wrong, f'every wheel matches Python {PY_SERIES} or is pure',
              ', '.join(wrong[:3]))

        native = {t for t in arch_tags if t != 'any'}
        expected = ('manylinux' if 'linux' in os.path.basename(path) else 'macosx')
        odd = [t for t in native if expected not in t]
        check(not odd, f'every native wheel targets {expected}', ', '.join(sorted(odd)[:3]))

        # ------------------------------------- can the install run offline?
        print('\n-- offline install')
        reqs = tf.extractfile(get('requirements.txt')).read().decode().split()
        norm = lambda s: re.sub(r'[-_.]+', '-', s).lower()
        have = {norm(w.split('-')[0]) for w in wheels} | {norm(s.rsplit('-', 1)[0]) for s in sdists}
        missing = [r for r in reqs if norm(r) not in have]
        check(not missing, f'all {len(reqs)} requirements are in wheels/',
              'missing: ' + ', '.join(missing))

        if sdists:
            check(any(norm(w.split('-')[0]) == 'setuptools' for w in wheels),
                  'setuptools is bundled, so the source-only dependency can build',
                  'without it the first run fails offline')

        # ------------------------------------------------- reference corpus
        print('\n-- bundled reference data')
        check(get('reference_processed/DAC_Mesh_final_network_with_relevance.json') is not None,
              'the reference network ships, so figures work before any download')
        check(get('reference_processed/aop_annotations_master.csv') is not None,
              'the AOP annotation dictionary ships with it')

    print()
    print('=' * 62)
    print(f'PASSED {len(PASS)}   FAILED {len(FAIL)}')
    if FAIL:
        for f in FAIL:
            print('   FAILED:', f)
    print('\nNot proven here: that it RUNS. This machine cannot execute a Linux\n'
          'or macOS binary. The release workflow installs and launches it on\n'
          'real runners; that is the test this cannot be.')
    return not FAIL


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('bundles', nargs='+')
    a = ap.parse_args()
    sys.exit(0 if all(verify(b) for b in a.bundles) else 1)
