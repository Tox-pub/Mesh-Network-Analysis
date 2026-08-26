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


def dotnet_env():
    """An environment in which wix.exe can find the .NET runtime.

    wix.exe is a .NET apphost: it locates hostfxr.dll through DOTNET_ROOT, or by
    finding dotnet on PATH, or at the machine-wide install. None of those held
    here - the SDK was installed per-user under LocalAppData while only the
    tools directory was added to PATH - so wix aborted with
    "Failed to resolve hostfxr.dll ... 0x80008083", which reads like a WiX
    problem and is not one. Resolving it here makes the build independent of
    whatever the calling shell happens to have inherited.
    """
    env = dict(os.environ)
    if env.get('DOTNET_ROOT') and os.path.isfile(
            os.path.join(env['DOTNET_ROOT'], 'dotnet.exe')):
        return env
    for root in (os.path.join(os.environ.get('LOCALAPPDATA', ''), 'Microsoft', 'dotnet'),
                 os.path.join(os.environ.get('ProgramFiles', ''), 'dotnet'),
                 os.path.join(os.environ.get('ProgramW6432', ''), 'dotnet'),
                 os.path.join(os.path.expanduser('~'), '.dotnet')):
        if root and os.path.isfile(os.path.join(root, 'dotnet.exe')):
            env['DOTNET_ROOT'] = root
            env['PATH'] = root + os.pathsep + env.get('PATH', '')
            print(f'  dotnet   : {root}')
            return env
    return env


def default_portable():
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
    ap.add_argument('--wix', default=None, help='full path to wix.exe')
    ap.add_argument('--out', default=None,
                    help='where to write the .msi (default: beside the portable build)')
    a = ap.parse_args()

    portable = os.path.normpath(a.portable or default_portable())
    check_not_stale(portable, a.force)
    if not os.path.isdir(portable):
        sys.exit(f'portable build not found: {portable}\n'
                 'Run: python packaging/build_portable_windows.py')
    for needed in ('python\\python.exe', 'app\\launch.py'):
        if not os.path.exists(os.path.join(portable, needed)):
            sys.exit(f'{portable} does not look like a finished build ({needed} missing).')

    purge_runtime_artifacts(portable)
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
           '-ext', 'WixToolset.Util.wixext',
           '-arch', 'x64',
           '-b', HERE,
           '-o', msi]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=dotnet_env())
    if proc.returncode != 0:
        sys.stdout.write(proc.stdout[-4000:])
        sys.stderr.write(proc.stderr[-2000:])
        if 'hostfxr' in (proc.stdout + proc.stderr):
            sys.stderr.write(
                '\nwix could not find the .NET runtime. Install the .NET 8 SDK, '
                'or set DOTNET_ROOT to the folder containing dotnet.exe.\n')
        sys.exit(f'wix failed with code {proc.returncode}')

    if not os.path.exists(msi):
        sys.exit('wix reported success but produced no .msi')
    print(f'\n  msi      : {msi}  ({os.path.getsize(msi) / 1e6:,.0f} MB)')
    print('\n  Installs per user, so msiexec asks for no elevation. Silent install:\n'
          f'    msiexec /i "{os.path.basename(msi)}" /qn')


if __name__ == '__main__':
    main()
