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
      portable.marker        keeps settings and results in this folder
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
sys.path.insert(0, HERE)
# Where a build may write. Decided in one place for every script in this folder,
# so they cannot disagree about it. See build_location.py.
import build_location                                              # noqa: E402

NAME = 'MeshWorkbench'          # zip/folder stem, no spaces
DISPLAY = 'MeSH Workbench'      # what the user sees and clicks
VERSION = '3.2.10'   # tracks the project version in pyproject.toml
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


# Everything the pipeline imports, plus what those need in turn. Copied wholesale
# from the working venv; anything missing shows up immediately as an ImportError
# when the launcher is tested, which is the point of testing it here.
SKIP_DIRS = {'__pycache__', 'pip', 'setuptools', 'wheel', 'pkg_resources',
             '_distutils_hack', 'PyInstaller', '_pyinstaller_hooks_contrib',
             'altgraph', 'pefile', 'pyinstaller_hooks_contrib'}
# Bytecode, and the shipped test suites - 94 MB of the dependency tree, none of
# which the pipeline imports.
#
# The distinction here is exact and was learned the hard way. "testing" is NOT
# in this set: numpy.testing is public API that scipy reaches for while loading,
# and dropping it broke every scipy-dependent import. "tests" and "test" are the
# suites themselves, along with their fixture corpora - gensim alone ships a
# Wikipedia extract and word2vec models under test/test_data. Removing a
# directory whose name merely resembles a test suite is what caused the earlier
# breakage, so the names are matched exactly rather than by substring.
PRUNE = {'__pycache__', 'tests', 'test'}


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
    # The AOP stratum dictionary travels with the reference corpus. Without it
    # the figure step cannot build its run-annotations file, so a fresh install
    # could not draw a single figure even in reference mode.
    anno = os.path.join(repo, 'data', 'raw', 'aop_annotations_master.csv')
    if os.path.exists(anno):
        os.makedirs(os.path.join(dst, 'reference_processed'), exist_ok=True)
        shutil.copy2(anno, os.path.join(dst, 'reference_processed',
                                        'aop_annotations_master.csv'))
    ref = os.path.join(repo, 'data', 'reference_processed')
    if os.path.isdir(ref):
        shutil.copytree(ref, os.path.join(dst, 'reference_processed'),
                        dirs_exist_ok=True, ignore=skip)
    # The documentation. The Help menu opens these, and without them every entry
    # in it reported the file missing - the bundle shipped no markdown at all.
    # They go one level above app/, beside the launchers, where a user browsing
    # the folder will also find them.
    docs_dst = os.path.dirname(dst)
    for name in ('LICENSE', 'THIRD-PARTY-NOTICES.md', 'CITATION.cff',
                 'HELP.md', 'INSTALL.md', 'README.md', 'COMMAND-LINE.md'):
        src_doc = os.path.join(repo, name)
        if os.path.exists(src_doc):
            shutil.copy2(src_doc, os.path.join(docs_dst, name))
    prov = os.path.join(repo, 'results', 'reference_figures', 'PROVENANCE.md')
    if os.path.exists(prov):
        os.makedirs(os.path.join(dst, 'reference_figures'), exist_ok=True)
        shutil.copy2(prov, os.path.join(dst, 'reference_figures', 'PROVENANCE.md'))

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

# The three .bat launchers are not built here any more - they live as real
# files in packaging/launchers/ so they can be read and reviewed in the
# repository, and are copied verbatim by main().

MARKER_TEXT = """This file marks the folder as a portable copy.

Settings, data and results are kept here rather than under your user profile.
Delete it and this copy behaves as an installed one: settings and results move
to your own profile, and each user of the machine gets their own.
"""

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

    For a Desktop icon, run "Create desktop shortcut.bat" once. It points at
    this folder rather than copying anything, so keep the folder where it is.

    To install it properly instead - Start menu entry, desktop icon, and an
    entry in Settings > Apps - run "Install.bat". That copies the program into
    your user profile and gives each account its own settings and results.
    Nothing needs administrator rights.

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
    Install.bat             installs it properly, if you would rather
    Create desktop shortcut.bat   puts an icon on your Desktop
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
    ap.add_argument('--out', default=None,
                    help=f'directory to assemble the release in '
                         f'(default: {build_location.BUILD_ROOT})')
    ap.add_argument('--no-zip', action='store_true')
    a = ap.parse_args()
    a.out = build_location.resolve(a.out, purpose='Windows portable build')

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

    # Launchers are copied from packaging/launchers rather than written from
    # strings in here, so what a user double-clicks is a file anyone can read and
    # review in the repository. Batch files are copied verbatim to preserve their
    # CRLF line endings, which cmd.exe requires.
    launchers = os.path.join(HERE, 'launchers')
    copied = 0
    for fname in sorted(os.listdir(launchers)):
        if fname.lower().endswith('.bat'):
            shutil.copy2(os.path.join(launchers, fname), os.path.join(out, fname))
            copied += 1
    with open(os.path.join(out, 'README - Install and First Run.txt'),
              'w', encoding='utf-8') as fh:
        fh.write(README)
    # Declares this copy self-contained: settings, data and results stay in the
    # folder rather than going to the user profile. The installer excludes it,
    # which is what makes an installed copy per-user. Decided by a marker rather
    # than by testing writability - a writable folder is not the same thing as a
    # portable copy, and guessing would silently relocate someone's results.
    with open(os.path.join(out, 'portable.marker'), 'w', encoding='utf-8') as fh:
        fh.write(MARKER_TEXT)
    print(f'  launchers ({copied}) + README written')

    total = tree_size(out)
    print(f'\n  folder: {out}\n  size  : {total:,.0f} MB')

    # A tree that has been launched writes mesh_config.json into itself - it is
    # a portable copy, so that is where its settings live. Zipping the tree
    # afterwards puts that file in the download. It was empty this time; a tree
    # the developer had configured would have carried their NCBI e-mail and API
    # key to every person who extracted it. The installers already purge this;
    # the zip did not.
    strays = 0
    for name in ('mesh_config.json', '.installed', 'mesh_config.json.bak'):
        p = os.path.join(out, name)
        if os.path.isfile(long(p)):
            os.remove(long(p))
            strays += 1
    for name in ('data', 'results'):
        p = os.path.join(out, name)
        if os.path.isdir(long(p)):
            shutil.rmtree(long(p), ignore_errors=True)
            strays += 1
    if strays:
        print(f'  purged {strays} file(s) a previous run left in the tree')

    if not a.no_zip:
        zpath = os.path.join(a.out, f'{NAME}-{VERSION}-win64-portable.zip')
        print('  compressing ...')
        # Deflate, at the highest level. LZMA would take roughly another quarter
        # off, but a zip compressed that way cannot be opened by Windows Explorer
        # or PowerShell's Expand-Archive - both refuse it with "unsupported
        # compression method". An archive the user cannot open is not a saving.
        with zipfile.ZipFile(zpath, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as z:
            root = long(out)
            for b, _, fs in os.walk(root):
                for f in fs:
                    p = os.path.join(b, f)
                    z.write(p, os.path.join(NAME, os.path.relpath(p, root)))
        print(f'  zip   : {zpath}  ({os.path.getsize(zpath)/1e6:,.0f} MB)')
    print(f'\nDone in {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()
