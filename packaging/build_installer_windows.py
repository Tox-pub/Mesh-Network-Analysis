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

Setup.exe is itself unsigned, so SmartScreen warns once on a machine that has
not seen it before. Only a code-signing certificate removes that, and it never
affects the installed application.
"""

import argparse
import glob
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ISS = os.path.join(HERE, 'windows_installer.iss')
OUT_DIR = os.path.join(HERE, 'dist')

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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--portable', default=None,
                    help='the assembled portable folder (default: the last build)')
    ap.add_argument('--iscc', default=None, help='full path to ISCC.exe')
    a = ap.parse_args()

    portable = a.portable or default_portable()
    if not os.path.isdir(portable):
        sys.exit(f'portable build not found: {portable}\n'
                 'Run: python packaging/build_portable_windows.py')
    # Guard against pointing this at a half-built or unrelated folder.
    for needed in ('python\\python.exe', 'app\\launch.py'):
        if not os.path.exists(os.path.join(portable, needed)):
            sys.exit(f'{portable} does not look like a finished build ({needed} missing).')

    iscc = find_iscc(a.iscc)
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f'Building installer\n  portable : {portable}\n  compiler : {iscc}')

    cmd = [iscc, f'/DPortableDir={portable}', f'/O{OUT_DIR}', ISS]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        # Inno reports the offending line on stdout, not stderr.
        sys.stdout.write(proc.stdout[-4000:])
        sys.stderr.write(proc.stderr[-2000:])
        sys.exit(f'ISCC failed with code {proc.returncode}')

    made = sorted(glob.glob(os.path.join(OUT_DIR, '*.exe')), key=os.path.getmtime)
    if not made:
        sys.exit('ISCC reported success but produced no .exe')
    exe = made[-1]
    print(f'\n  setup    : {exe}  ({os.path.getsize(exe) / 1e6:,.0f} MB)')
    print('\n  Setup.exe is unsigned: SmartScreen will warn once, and the user\n'
          '  clicks More info -> Run anyway. The installed application is the\n'
          '  signed python.exe either way.')


if __name__ == '__main__':
    main()
