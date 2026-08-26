# -*- coding: utf-8 -*-
"""
build_installer_windows.py - wrap a finished portable build in a Setup.exe.

This does not compile anything. It hands the already-assembled portable tree to
Inno Setup, which packs it into a conventional installer: a wizard with a
directory page, Start-menu and desktop shortcuts, and an Add/Remove Programs
entry. The installed application is byte-for-byte the portable one, so the only
program that ever executes is the PSF-signed python.exe.

Run build_portable_windows.py first; this refuses to run without its output.

    python packaging/build_installer_windows.py
    python packaging/build_installer_windows.py --portable D:\\some\\MeshWorkbench

Setup.exe is itself unsigned. On an unmanaged machine that costs a SmartScreen
prompt; on a managed one it can be refused outright by Defender Exploit Guard,
which blocks executables that have no prevalence and no signature. Only a
code-signing certificate removes that. Install.bat in the portable zip is the
fallback for those machines - it creates no executable at all and installs the
same tree.
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ISS = os.path.join(HERE, 'windows_installer.iss')

# Inno Setup installs per-user as readily as machine-wide, and 7 is not always
# where 6 was, so look in both and in both scopes rather than naming one path.
ISCC_CANDIDATES = [
    os.path.join(os.environ.get('LOCALAPPDATA', ''), 'Programs', f'Inno Setup {v}', 'ISCC.exe')
    for v in (7, 6)
] + [
    os.path.join(os.environ.get('ProgramFiles', ''), f'Inno Setup {v}', 'ISCC.exe')
    for v in (7, 6)
] + [
    os.path.join(os.environ.get('ProgramFiles(x86)', ''), f'Inno Setup {v}', 'ISCC.exe')
    for v in (7, 6)
]


def find_iscc(explicit=None):
    if explicit:
        if os.path.isfile(explicit):
            return explicit
        sys.exit(f'ISCC.exe not found at: {explicit}')
    for path in ISCC_CANDIDATES:
        if path and os.path.isfile(path):
            return path
    sys.exit('Inno Setup not found. Install it from https://jrsoftware.org/isdl.php\n'
             'or pass --iscc with the full path to ISCC.exe.')


def default_portable():
    """The portable tree the sibling script most recently produced."""
    from importlib.util import module_from_spec, spec_from_file_location
    spec = spec_from_file_location('bp', os.path.join(HERE, 'build_portable_windows.py'))
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return os.path.join(mod.BUILD_OUT_DEFAULT, mod.NAME)


def check_not_stale(portable, force=False):
    """Compare the packaged app against the source it should have come from.

    Silent when the tree is current. When it is not, this stops rather than
    warns: a stale installer is indistinguishable from a good one once it has
    been handed to someone, and the only cost of stopping is one rebuild.
    """
    repo = os.path.dirname(HERE)
    src = os.path.join(repo, 'src')
    app = os.path.join(portable, 'app')
    if not os.path.isdir(src) or not os.path.isdir(app):
        return

    def newest(root):
        latest, where = 0.0, None
        for base, dirs, files in os.walk(root):
            dirs[:] = [d for d in dirs if d != '__pycache__']
            for name in files:
                if not name.endswith('.py'):
                    continue
                path = os.path.join(base, name)
                try:
                    mtime = os.path.getmtime(path)
                except OSError:
                    continue
                if mtime > latest:
                    latest, where = mtime, path
        return latest, where

    src_time, src_file = newest(src)
    app_time, _ = newest(app)
    if not src_time or not app_time or app_time >= src_time:
        return

    import datetime
    fmt = lambda t: datetime.datetime.fromtimestamp(t).strftime('%Y-%m-%d %H:%M')
    message = (
        '\n  [!] THIS PORTABLE TREE IS OUT OF DATE\n'
        f'      tree   {portable}\n'
        f'             built from code last changed {fmt(app_time)}\n'
        f'      source {src}\n'
        f'             last changed {fmt(src_time)}'
        f' ({os.path.relpath(src_file, repo)})\n\n'
        '      Packaging this would ship code that is already superseded.\n'
        '      Rebuild the portable tree first:\n\n'
        f'          python packaging/build_portable_windows.py --out <dir>\n\n'
        '      then point this script at it with --portable <dir>/MeshWorkbench.\n'
        '      Pass --force to package it anyway.\n')
    if force:
        print(message)
        print('  --force given: packaging the stale tree regardless.\n')
        return
    sys.exit(message)


# Everything a *run* leaves in a build tree. The zip builder assembles a fresh
# tree each time and never sees any of it; the installers harvest whatever is on
# disk, so a tree someone has actually launched the program from carries extra
# baggage into the release. Both were observed: bytecode added 2,108 files and
# 16 MB to one MSI, and a mesh_config.json written on first run was packaged too.
RUNTIME_FILES = ('mesh_config.json',)
RUNTIME_DIRS = ('data', 'results')


def purge_runtime_artifacts(tree):
    """Strip anything a previous run left behind, so the build is reproducible.

    A shipped mesh_config.json is the more serious of the two: this one held
    only defaults, but a tree the developer had configured would carry their
    search terms - and their NCBI credentials - into a public download.
    """
    removed_dirs = removed_files = 0
    for base, dirs, _ in os.walk(tree):
        for d in list(dirs):
            if d == '__pycache__':
                shutil.rmtree(os.path.join(base, d), ignore_errors=True)
                dirs.remove(d)
                removed_dirs += 1
    for name in RUNTIME_FILES:
        p = os.path.join(tree, name)
        if os.path.isfile(p):
            os.remove(p)
            removed_files += 1
    for name in RUNTIME_DIRS:
        p = os.path.join(tree, name)
        if os.path.isdir(p):
            shutil.rmtree(p, ignore_errors=True)
            removed_files += 1
    if removed_dirs or removed_files:
        print(f'  purged {removed_dirs} __pycache__ dirs, {removed_files} runtime file(s)')
    return removed_dirs + removed_files


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--force', action='store_true',
                    help='package the portable tree even if it is older than the source')
    ap.add_argument('--portable', default=None,
                    help='the assembled portable folder (default: the last build)')
    ap.add_argument('--iscc', default=None, help='full path to ISCC.exe')
    ap.add_argument('--out', default=None,
                    help='where to write Setup.exe (default: beside the portable build)')
    a = ap.parse_args()

    portable = a.portable or default_portable()
    check_not_stale(portable, a.force)
    if not os.path.isdir(portable):
        sys.exit(f'portable build not found: {portable}\n'
                 'Run: python packaging/build_portable_windows.py')
    # Guard against pointing this at a half-built or unrelated folder.
    for needed in ('python\\python.exe', 'app\\launch.py'):
        if not os.path.exists(os.path.join(portable, needed)):
            sys.exit(f'{portable} does not look like a finished build ({needed} missing).')

    purge_runtime_artifacts(portable)
    iscc = find_iscc(a.iscc)
    # Beside the portable build, off the project: Setup.exe is ~100 MB of
    # reproducible output and the working copy is cloud-synced.
    out_dir = a.out or os.path.dirname(os.path.normpath(portable))
    os.makedirs(out_dir, exist_ok=True)
    print(f'Building installer\n  portable : {portable}\n  compiler : {iscc}')

    cmd = [iscc, f'/DPortableDir={portable}', f'/O{out_dir}', ISS]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        # Inno reports the offending line on stdout, not stderr.
        sys.stdout.write(proc.stdout[-4000:])
        sys.stderr.write(proc.stderr[-2000:])
        sys.exit(f'ISCC failed with code {proc.returncode}')

    made = sorted(glob.glob(os.path.join(out_dir, '*.exe')), key=os.path.getmtime)
    if not made:
        sys.exit('ISCC reported success but produced no .exe')
    exe = made[-1]
    print(f'\n  setup    : {exe}  ({os.path.getsize(exe) / 1e6:,.0f} MB)')
    print('\n  Setup.exe is unsigned. On an unmanaged machine that costs a\n'
          '  SmartScreen prompt the user clicks through. On a managed one it can\n'
          '  be refused outright: Defender Exploit Guard rule 01443614 blocks\n'
          '  executables with no prevalence and no signature, and the user sees\n'
          '  0x80070005 with nothing to click. Ship Install.bat from the zip as\n'
          '  the fallback - it creates no executable and installs the same tree.')


if __name__ == '__main__':
    main()
