# -*- coding: utf-8 -*-
"""
build_portable_windows.py - assemble the download-and-run distribution.

No PyInstaller, no installer, no code signing. The folder carries Python's own
Windows embeddable package, so the only executable a user launches is
python.exe - already signed by the Python Software Foundation. Nothing here is
a freshly compiled unsigned binary, which is what antivirus heuristics and
SmartScreen react to, and what blocked the frozen build.

Layout produced:

    MeshWorkbench/
      python/                embeddable CPython (PSF-signed python.exe)
      python/Lib/site-packages/   the pipeline's dependencies
      app/                   the Workbench and mesh_aop
      MeSH Workbench.bat     the launcher a user double-clicks
      Uninstall.bat          removal, including the temp files outside this tree
      README - Install and First Run.txt

Dependencies are copied from a working virtual environment rather than
pip-installed fresh, so the shipped versions are exactly the ones the published
results were produced with.
"""

import argparse
import os
import shutil
import sys
import time
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
NAME = 'MeshWorkbench'          # zip/folder stem, no spaces
DISPLAY = 'MeSH Workbench'      # what the user sees and clicks
VERSION = '3.0.0'   # tracks the project version in pyproject.toml
PY_VER = '3.12.7'
EMBED_URL = f'https://www.python.org/ftp/python/{PY_VER}/python-{PY_VER}-embed-amd64.zip'
# This script lives in <repo>/packaging, so the project it packages is its parent.
REPO_DEFAULT = os.path.dirname(HERE)
# The environment the release is built from. Dependencies are copied out of it
# rather than installed fresh, so a build ships the versions actually tested.
# Defaults to whichever interpreter is running this script, which is the one the
# builder already chose.
VENV_DEFAULT = os.environ.get('MESH_VENV', sys.prefix)
# The embeddable package deliberately omits tkinter and the Tcl/Tk runtime, so
# the GUI cannot start from it alone. Those files are taken from a full CPython
# install of the same version - which, when this runs inside a venv, is exactly
# what sys.base_prefix points at. Asking Python beats naming one machine's path.
BASE_PY_DEFAULT = os.environ.get('MESH_BASE_PYTHON', sys.base_prefix)


def _default_build_out():
    """Where a release is assembled.

    Off the project by preference: the tree is ~460 MB, is rebuilt from scratch
    every run, and the working copy is cloud-synced, so building in place uploads
    half a gigabyte of reproducible output every time - .gitignore stops git, not
    the sync client.

    D: is used when it is genuinely writable. Testing only that it exists is not
    enough; on another machine D: is as likely to be a read-only optical or
    mapped drive, and the build would fail partway through instead of falling
    back. Override with --out or MESH_BUILD_OUT.
    """
    env = os.environ.get('MESH_BUILD_OUT')
    if env:
        return env
    for candidate in (r'D:\mesh_workbench_build',):
        try:
            os.makedirs(candidate, exist_ok=True)
            probe = os.path.join(candidate, '.write_test')
            with open(probe, 'w') as fh:
                fh.write('')
            os.remove(probe)
            return candidate
        except OSError:
            continue
    return os.path.join(HERE, 'portable')


BUILD_OUT_DEFAULT = _default_build_out()

# Everything the pipeline imports, plus what those need in turn. Copied wholesale
# from the working venv; anything missing shows up immediately as an ImportError
# when the launcher is tested, which is the point of testing it here.
SKIP_DIRS = {'__pycache__', 'pip', 'setuptools', 'wheel', 'pkg_resources',
             '_distutils_hack', 'PyInstaller', '_pyinstaller_hooks_contrib',
             'altgraph', 'pefile', 'pyinstaller_hooks_contrib'}
# Only bytecode. Pruning directories that merely LOOK like test suites breaks
# imports: numpy.testing is public API that scipy touches while loading, so
# dropping "testing" made every scipy-dependent import fail. The couple of
# hundred megabytes saved is not worth guessing which of these are real.
PRUNE = {'__pycache__'}


def fetch_embed(dst, cache):
    """Download the PSF embeddable package (or reuse a cached copy)."""
    os.makedirs(cache, exist_ok=True)
    zpath = os.path.join(cache, os.path.basename(EMBED_URL))
    if not os.path.exists(zpath):
        print(f'  downloading {EMBED_URL}')
        urllib.request.urlretrieve(EMBED_URL, zpath)
    size = os.path.getsize(zpath) / 1e6
    print(f'  embeddable python: {os.path.basename(zpath)} ({size:.0f} MB)')
    os.makedirs(dst, exist_ok=True)
    with zipfile.ZipFile(zpath) as z:
        z.extractall(dst)
    return zpath


def add_tkinter(pydir, base_python):
    """Copy Tk into the embeddable runtime.

    The Windows embeddable package ships without tkinter: no Lib/tkinter, no
    _tkinter.pyd, no Tcl/Tk runtime. For a console tool that is irrelevant, but
    this application IS a Tk window, so the bundle is unusable until these are
    added from a full CPython install of the same version.
    """
    if not os.path.isdir(base_python):
        sys.exit(f'full python install not found (needed for tkinter): {base_python}')
    jobs = [
        (os.path.join(base_python, 'Lib', 'tkinter'),
         os.path.join(pydir, 'Lib', 'tkinter'), 'dir'),
        (os.path.join(base_python, 'tcl'), os.path.join(pydir, 'tcl'), 'dir'),
        (os.path.join(base_python, 'DLLs', '_tkinter.pyd'), pydir, 'file'),
        (os.path.join(base_python, 'DLLs', 'tcl86t.dll'), pydir, 'file'),
        (os.path.join(base_python, 'DLLs', 'tk86t.dll'), pydir, 'file'),
        (os.path.join(base_python, 'DLLs', 'zlib1.dll'), pydir, 'file'),
    ]
    got = 0
    for src, dst, kind in jobs:
        if not os.path.exists(src):
            print(f'  [!] missing {src}')
            continue
        if kind == 'dir':
            shutil.copytree(src, dst, dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns('__pycache__'))
        else:
            shutil.copy2(src, dst)
        got += 1
    print(f'  tkinter runtime: {got}/{len(jobs)} components copied')


def write_pth(pydir):
    """Point the embedded interpreter at our packages.

    The embeddable build ships isolated: `import site` is commented out and
    site-packages is not on the path. Without editing this file the interpreter
    starts but cannot see numpy, and the failure looks like a broken install.
    """
    tag = 'python' + PY_VER.replace('.', '')[:3]
    pth = os.path.join(pydir, tag + '._pth')
    # `Lib` must be listed explicitly. The embeddable build serves the standard
    # library from python312.zip and puts only '.' on the path, so a package
    # added under Lib/ - which is how tkinter arrives - is invisible without it.
    with open(pth, 'w', encoding='utf-8') as fh:
        fh.write(f'{tag}.zip\n'
                 '.\n'
                 'Lib\n'
                 'Lib\\site-packages\n'
                 '..\\app\n'
                 'import site\n')
    print(f'  wrote {os.path.basename(pth)}')


def long(path):
    """Return `path` in the form that lifts the 260-character path limit.

    Some dependencies ship deeply nested files with long names -- gensim's test
    corpora are the worst offenders -- and the copy fails with WinError 3 once
    the destination crosses MAX_PATH. Keeping the output shallow avoids most of
    it, but a caller can point --out anywhere, so do not rely on that.
    The `\\\\?\\` prefix opts the call out of the limit; it requires a fully
    normalised absolute path with no forward slashes, and shutil propagates it
    to everything it joins underneath.
    """
    if os.name != 'nt':
        return path
    p = os.path.abspath(path)
    return p if p.startswith('\\\\?\\') else '\\\\?\\' + p


def tree_size(root):
    """Total size in MB. Walks the extended-length form for the same reason
    the copy does: the tree contains paths past MAX_PATH, and a plain walk
    yields names that then fail to stat."""
    return sum(os.path.getsize(os.path.join(b, f))
               for b, _, fs in os.walk(long(root)) for f in fs) / 1e6


def copy_packages(venv, dst):
    src = os.path.join(venv, 'Lib', 'site-packages')
    if not os.path.isdir(src):
        sys.exit(f'site-packages not found: {src}')
    os.makedirs(dst, exist_ok=True)
    n = kept = 0
    for entry in sorted(os.listdir(src)):
        # .pth files are executed at interpreter start. The ones here belong to
        # setuptools shims and editable installs that this build deliberately
        # leaves behind, so they resolve to nothing and raise on every launch.
        if entry.endswith('.pth'):
            continue
        # .dist-info is kept: several packages read their own version through
        # importlib.metadata at import time, and plotly raises PackageNotFoundError
        # outright without it. All 54 of them together cost under 4 MB.
        if entry in SKIP_DIRS or entry.endswith('.egg-link'):
            continue
        s, d = os.path.join(src, entry), os.path.join(dst, entry)
        if os.path.isdir(s):
            shutil.copytree(long(s), long(d), dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns(*PRUNE))
            kept += 1
        else:
            shutil.copy2(long(s), long(d))
        n += 1
    size = tree_size(dst)
    print(f'  packages: {kept} directories, {n} entries, {size:,.0f} MB')


def copy_app(repo, dst):
    os.makedirs(dst, exist_ok=True)
    skip = shutil.ignore_patterns('__pycache__', '*.pyc', '.ipynb_checkpoints',
                                  '*.ipynb', 'shots', 'dist', 'build', '_staged',
                                  '_spec', '_cache')
    shutil.copytree(os.path.join(repo, 'src', 'mesh_workbench'),
                    os.path.join(dst, 'mesh_workbench'), dirs_exist_ok=True, ignore=skip)
    shutil.copytree(os.path.join(repo, 'src', 'mesh_aop'),
                    os.path.join(dst, 'mesh_aop'), dirs_exist_ok=True, ignore=skip)
    ref = os.path.join(repo, 'data', 'reference_processed')
    if os.path.isdir(ref):
        shutil.copytree(ref, os.path.join(dst, 'reference_processed'),
                        dirs_exist_ok=True, ignore=skip)
    cfg = os.path.join(repo, 'mesh_config.json')
    if os.path.exists(cfg):
        shutil.copy(cfg, os.path.join(dst, 'mesh_config.default.json'))
    with open(os.path.join(dst, 'launch.py'), 'w', encoding='utf-8') as fh:
        fh.write(LAUNCH_PY)
    print('  app + pipeline copied')


LAUNCH_PY = '''# -*- coding: utf-8 -*-
"""Entry point for the portable build."""
import os
import sys

APP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, APP)

# The workbench asks the pipeline to run in a child process; in this build that
# child is the same embedded interpreter, found relative to the app folder.
PY = os.path.join(os.path.dirname(APP), 'python', 'python.exe')

if __name__ == '__main__':
    repo = os.environ.get('MESH_REPO') or os.path.dirname(APP)
    from mesh_workbench.app import main
    main(repo, PY if os.path.exists(PY) else sys.executable)
'''

LAUNCH_BAT = r'''@echo off
REM  MeSH Workbench - portable launcher
REM
REM  Runs the bundled Python, which is the official embeddable build signed by
REM  the Python Software Foundation. Nothing here needs installing and nothing
REM  is written outside this folder.
title MeSH Workbench
cd /d "%~dp0"
if not exist "python\python.exe" (
  echo [X] python\python.exe is missing.
  echo     Extract the whole zip, keeping the folder structure intact.
  echo.
  pause
  exit /b 1
)
start "" "python\pythonw.exe" "app\launch.py" %*
'''

UNINSTALL_BAT = r'''@echo off
REM  Removes the program and everything it downloaded or built - including the
REM  ETL working folder in the system temp directory, which deleting this folder
REM  would otherwise leave behind.
title MeSH Workbench - Uninstall
cd /d "%~dp0"
if not exist "python\python.exe" (
  echo [X] python\python.exe is missing - nothing to run.
  pause
  exit /b 1
)
"python\python.exe" -m mesh_aop.uninstall_cli --project "%~dp0."
echo.
set "FIN="
set /p FIN=Delete the remaining program folder too? [y/N]
if /i not "%FIN%"=="y" goto :done

REM  A folder cannot delete itself while a program inside it is running, so the
REM  removal is handed to a script in TEMP that waits for this one to exit.
> "%TEMP%\mw_cleanup.bat" echo @echo off
>>"%TEMP%\mw_cleanup.bat" echo timeout /t 3 /nobreak ^>nul
>>"%TEMP%\mw_cleanup.bat" echo rmdir /s /q "%~dp0."
>>"%TEMP%\mw_cleanup.bat" echo del "%%~f0"
cd /d "%TEMP%"
start "" /min cmd /c "%TEMP%\mw_cleanup.bat"
echo.
echo The folder will disappear in a few seconds.
timeout /t 2 /nobreak >nul
exit /b 0

:done
echo.
pause
'''

CONSOLE_BAT = r'''@echo off
REM  Same launcher with a console attached, so a startup error stays readable.
title MeSH Workbench (console)
cd /d "%~dp0"
"python\python.exe" "app\launch.py" %*
echo.
echo Exited with code %ERRORLEVEL%
pause
'''

README = r"""MeSH Workbench {ver}
===============================================================

WHAT THIS IS
    A desktop application for building and validating MeSH co-occurrence
    concept networks from PubMed literature.

INSTALLING
    There is no installer. Extract this folder anywhere you can write to -
    Documents, the Desktop, or a USB drive - then open it and double-click

        MeSH Workbench.bat

    Nothing is installed and no administrator rights are needed.

    Windows may take a few seconds to scan the files the first time.

UNINSTALLING
    Run "Uninstall.bat", or open the program and choose Tools -> Uninstall.

    Either one lists what is on disk, with sizes, and lets you pick what goes.
    Your results are kept unless you say otherwise.

    Do not simply delete this folder. Building the database leaves a working
    copy - the same size as the database, several gigabytes - in the Windows
    temp folder, and deleting this folder does not remove it. The uninstaller
    knows about it.

IF THE WINDOW DOES NOT APPEAR
    Run "MeSH Workbench (Troubleshooting).bat" instead. It does the same thing
    but keeps a terminal open, so any error message stays on screen.

BEFORE YOU CAN ANALYSE ANYTHING
    The pipeline works from a local copy of PubMed's annotation data. Nothing
    else runs until that database exists.

    The application opens on Data Setup, which lists what is on disk. Next to
    "Master annotation database" is a Build button; it downloads the PubMed
    baseline and compiles the database for you.

    Expect roughly 44 GB downloaded once and a database of about 8 GB built from
    it, over several hours. The download resumes if interrupted, and the 44 GB
    archive can be deleted afterwards from the same screen.

    Set an output folder on the Folders tab before running a step if you do not
    want results written into the project's own results folder.

REQUIREMENTS
    Windows 10 or 11, 64-bit.
    About 500 MB for this folder, plus room for the data described above.

WHAT IS IN HERE
    MeSH Workbench.bat      the launcher - this is the one to click
    Uninstall.bat           removes the program and its data
    python\                 Python {py}, the official build from python.org
    app\                    the application and the analysis pipeline

    The only program that runs is python.exe, which is signed by the Python
    Software Foundation.

LICENCE AND SOURCE
    https://github.com/Tox-pub/Mesh-Network-Analysis
""".format(ver=VERSION, py=PY_VER)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo', default=os.environ.get('MESH_REPO', REPO_DEFAULT))
    ap.add_argument('--venv', default=VENV_DEFAULT)
    ap.add_argument('--base-python', default=BASE_PY_DEFAULT,
                    help='full CPython install to take the Tk runtime from')
    ap.add_argument('--out', default=os.environ.get('MESH_BUILD_OUT', BUILD_OUT_DEFAULT),
                    help='directory to assemble the release in (default: %(default)s)')
    ap.add_argument('--no-zip', action='store_true')
    a = ap.parse_args()

    out = os.path.join(a.out, NAME)
    if os.path.exists(out):
        shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out, exist_ok=True)
    print(f'Building portable {NAME} {VERSION}')
    t0 = time.time()

    pydir = os.path.join(out, 'python')
    fetch_embed(pydir, os.path.join(a.out, '_cache'))
    add_tkinter(pydir, a.base_python)
    write_pth(pydir)
    copy_packages(a.venv, os.path.join(pydir, 'Lib', 'site-packages'))
    copy_app(a.repo, os.path.join(out, 'app'))

    for fname, body in ((f'{DISPLAY}.bat', LAUNCH_BAT),
                        (f'{DISPLAY} (Troubleshooting).bat', CONSOLE_BAT),
                        ('Uninstall.bat', UNINSTALL_BAT),
                        ('README - Install and First Run.txt', README)):
        with open(os.path.join(out, fname), 'w', encoding='utf-8') as fh:
            fh.write(body)
    print('  launchers + README written')

    total = tree_size(out)
    print(f'\n  folder: {out}\n  size  : {total:,.0f} MB')

    if not a.no_zip:
        zpath = os.path.join(a.out, f'{NAME}-{VERSION}-win64-portable.zip')
        print('  compressing ...')
        with zipfile.ZipFile(zpath, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
            root = long(out)
            for b, _, fs in os.walk(root):
                for f in fs:
                    p = os.path.join(b, f)
                    z.write(p, os.path.join(NAME, os.path.relpath(p, root)))
        print(f'  zip   : {zpath}  ({os.path.getsize(zpath)/1e6:,.0f} MB)')
    print(f'\nDone in {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()
