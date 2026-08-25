# -*- coding: utf-8 -*-
"""
build_msi_windows.py - wrap a finished portable build in an .msi.

The MSI exists alongside the Setup.exe for a specific reason. An MSI is executed
by msiexec.exe, which is part of Windows and signed by Microsoft, so no new
binary is introduced. A freshly compiled Setup.exe has no signature and no
prevalence, and Defender Exploit Guard refuses it outright on a managed machine
(0x80070005) - which is what happens on the machine this was developed on. The
MSI is the route most likely to work where that rule is enforced.

Requires the WiX Toolset, which is a .NET tool:

    dotnet tool install --global wix
    wix extension add --global WixToolset.UI.wixext

Run build_portable_windows.py first; this refuses to run without its output.

    python packaging/build_msi_windows.py
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WXS = os.path.join(HERE, 'windows_msi.wxs')


def find_wix(explicit=None):
    if explicit:
        if os.path.isfile(explicit):
            return explicit
        sys.exit(f'wix not found at: {explicit}')
    found = shutil.which('wix')
    if found:
        return found
    # dotnet global tools are not always on PATH in a non-login shell.
    candidate = os.path.join(os.path.expanduser('~'), '.dotnet', 'tools', 'wix.exe')
    if os.path.isfile(candidate):
        return candidate
    sys.exit(
        'WiX not found. Install it with:\n'
        '    dotnet tool install --global wix\n'
        '    wix extension add --global WixToolset.UI.wixext\n'
        'or pass --wix with the full path to wix.exe.\n'
        '(.NET SDK is required: https://dotnet.microsoft.com/download)')


def default_portable():
    from importlib.util import module_from_spec, spec_from_file_location
    spec = spec_from_file_location('bp', os.path.join(HERE, 'build_portable_windows.py'))
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return os.path.join(mod.BUILD_OUT_DEFAULT, mod.NAME)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--portable', default=None,
                    help='the assembled portable folder (default: the last build)')
    ap.add_argument('--wix', default=None, help='full path to wix.exe')
    ap.add_argument('--out', default=None,
                    help='where to write the .msi (default: beside the portable build)')
    a = ap.parse_args()

    portable = os.path.normpath(a.portable or default_portable())
    if not os.path.isdir(portable):
        sys.exit(f'portable build not found: {portable}\n'
                 'Run: python packaging/build_portable_windows.py')
    for needed in ('python\\python.exe', 'app\\launch.py'):
        if not os.path.exists(os.path.join(portable, needed)):
            sys.exit(f'{portable} does not look like a finished build ({needed} missing).')

    wix = find_wix(a.wix)
    out_dir = a.out or os.path.dirname(portable)
    os.makedirs(out_dir, exist_ok=True)
    msi = os.path.join(out_dir, 'MeSH-Workbench-3.1.0-win64.msi')

    print(f'Building MSI\n  portable : {portable}\n  compiler : {wix}')
    # -b puts packaging/ on the bind path so license.rtf resolves regardless of
    # the working directory the build was started from.
    cmd = [wix, 'build', WXS,
           '-d', f'PortableDir={portable}',
           '-ext', 'WixToolset.UI.wixext',
           '-arch', 'x64',
           '-b', HERE,
           '-o', msi]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stdout.write(proc.stdout[-4000:])
        sys.stderr.write(proc.stderr[-2000:])
        sys.exit(f'wix failed with code {proc.returncode}')

    if not os.path.exists(msi):
        sys.exit('wix reported success but produced no .msi')
    print(f'\n  msi      : {msi}  ({os.path.getsize(msi) / 1e6:,.0f} MB)')
    print('\n  Installs per user, so msiexec asks for no elevation. Silent install:\n'
          f'    msiexec /i "{os.path.basename(msi)}" /qn')


if __name__ == '__main__':
    main()
