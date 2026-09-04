# -*- coding: utf-8 -*-
"""
verify_windows_bundle.py - check a Windows build before it is published.

The Unix bundles have had a verifier since they existed; the Windows ones have
not, and every fault it would have caught has now happened at least once:

  * the developer's mesh_config.json packaged into the download, carrying their
    NCBI e-mail and API key to anyone who unzipped it;
  * an installer built from a tree that had been run, so bytecode and a data
    folder went out with it;
  * an MSI that harvested none of the application and shipped at 1 MB, because
    WiX does not treat a glob matching no files as an error.

Works on the portable .zip, the assembled folder, or the .msi.

    python packaging/verify_windows_bundle.py <build>\\MeshWorkbench-3.2.0-win64-portable.zip
    python packaging/verify_windows_bundle.py <build>\\MeshWorkbench
    python packaging/verify_windows_bundle.py <build>\\MeSH-Workbench-3.2.0-win64.msi
"""

import argparse
import os
import posixpath
import sys
import zipfile

PASS, FAIL = [], []

# Never in a download. The settings file is the serious one - it holds the
# credentials - but a build tree that has been run also carries results and a
# half-filled data folder, and those are nobody else's business either.
FORBIDDEN_NAMES = ('mesh_config.json', '.env', '.installed')
FORBIDDEN_PREFIXES = ('mesh_config.json.',)
FORBIDDEN_DIRS = ('data/', 'results/', '__pycache__/')

# A finished tree carries these. Their absence means it was assembled wrong.
REQUIRED = (
    'python/python.exe',
    'app/mesh_aop/cli.py',
    'app/mesh_workbench/app.py',
    # Renders the manual for the browser. Listed because losing it is silent:
    # the Help menu would fall back to the in-window reader and say nothing.
    'app/mesh_workbench/mdhtml.py',
    'app/reference_processed/DAC_Mesh_final_network_with_relevance.json',
    'LICENSE',
    'HELP.md',
    # HELP.md links to this in its second paragraph, so a bundle without it
    # ships a manual with a dead link. The portable build did, until now.
    'COMMAND-LINE.md',
    'README.md',
    'MeSH Workbench.bat',
    'Uninstall.bat',
)


def check(ok, msg, detail=''):
    (PASS if ok else FAIL).append(msg)
    print(('  [OK] ' if ok else '  [XX] ') + msg + (f'  -- {detail}' if detail and not ok else ''))
    return ok


def entries(target):
    """(relative posix paths, total uncompressed bytes) for any of the three shapes."""
    if os.path.isdir(target):
        names, total = [], 0
        for base, _, files in os.walk(target):
            for f in files:
                full = os.path.join(base, f)
                names.append(os.path.relpath(full, target).replace('\\', '/'))
                try:
                    total += os.path.getsize(full)
                except OSError:
                    pass
        return names, total

    if target.lower().endswith('.zip'):
        with zipfile.ZipFile(target) as z:
            infos = [i for i in z.infolist() if not i.is_dir()]
        # The zip holds one top-level folder; strip it so paths match the tree.
        root = infos[0].filename.split('/')[0] if infos else ''
        names = [i.filename[len(root) + 1:] if i.filename.startswith(root + '/')
                 else i.filename for i in infos]
        return names, sum(i.file_size for i in infos)

    if target.lower().endswith('.msi'):
        return None, os.path.getsize(target)

    sys.exit(f'not a Windows build: {target}')


def verify(target):
    print(f'\n=== {os.path.basename(target)} ===\n')
    names, total = entries(target)

    if names is None:
        # An MSI is a compound document; reading its file table needs Windows
        # Installer itself. Size is the one thing checkable everywhere, and it
        # is the fault that actually shipped.
        print('-- installer size')
        mb = total / 1e6
        check(mb >= 40, f'the installer is {mb:,.0f} MB, not an empty shell',
              'under 40 MB means it harvested little or none of the application')
        print('\n  An .msi cannot be listed without Windows Installer. Run this '
              'against\n  the portable .zip or the assembled folder to check its '
              'contents.')
        return not FAIL

    lower = [n.lower() for n in names]

    print('-- nothing of the developer in the download')
    leaked = [n for n in names
              if posixpath.basename(n) in FORBIDDEN_NAMES
              or posixpath.basename(n).startswith(FORBIDDEN_PREFIXES)]
    check(not leaked, 'no settings file, .env or install marker is packaged',
          'found: ' + ', '.join(leaked[:4]))

    # Anchored at the top level, except __pycache__ which is unwanted anywhere.
    # Matching 'data/' loosely flags Biopython's substitution_matrices/data and
    # a dozen other legitimate package data folders - a check that cries wolf
    # about a hundred correct files teaches people to ignore it.
    stray_dirs = sorted(
        {d for d in ('data/', 'results/') for n in lower if n.startswith(d)}
        | {'__pycache__/' for n in lower if '__pycache__/' in n})
    check(not stray_dirs, 'no data, results or bytecode folders from a used tree',
          'found at the top level: ' + ', '.join(stray_dirs))

    print('\n-- the application is actually in there')
    for rel in REQUIRED:
        check(rel.lower() in lower, f'present: {rel}')
    mb = total / 1e6
    check(mb >= 150, f'uncompressed tree is {mb:,.0f} MB',
          'under 150 MB means the dependencies are missing')

    print('\n-- the licence and the notices travel with it')
    check('license' in lower, 'LICENSE is included')
    check(any(n.lower() == 'third-party-notices.md' for n in names),
          'THIRD-PARTY-NOTICES.md is included',
          'the bundle redistributes ~47 libraries, four of them copyleft')

    print('\n-- the launcher')
    bats = [n for n in names if n.lower().endswith('.bat')]
    check(bats, f'{len(bats)} .bat launcher(s)')
    check(any('mesh workbench' in b.lower() for b in bats),
          'the one a user double-clicks is there')

    print()
    print('=' * 62)
    print(f'PASSED {len(PASS)}   FAILED {len(FAIL)}')
    for f in FAIL:
        print('   FAILED:', f)
    if not FAIL:
        print('\nNot proven here: that it INSTALLS and RUNS on a machine with no\n'
              'Python. The release workflow does that on a clean runner.')
    return not FAIL


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('targets', nargs='+', help='a .zip, an assembled folder, or an .msi')
    a = ap.parse_args()
    sys.exit(0 if all(verify(t) for t in a.targets) else 1)
