# -*- coding: utf-8 -*-
"""
build.py - produce a distributable MeSH AOP Workbench for the host platform.

PyInstaller does not cross-compile: a Windows .exe must be built on Windows, a
.app on macOS, an ELF binary on Linux. This script therefore builds for whatever
it is run on, and the CI workflow (.github/workflows/build.yml) runs it on all
three so a release carries one artifact per platform - the same arrangement any
Python desktop tool ships under.

Stage 1 (here) is the application bundle. Stage 2 is the platform installer:

    Windows   Inno Setup    -> MeshWorkbench-<ver>-win64-setup.exe
    macOS     hdiutil       -> MeshWorkbench-<ver>-macos.dmg
    Linux     tar.gz        -> MeshWorkbench-<ver>-linux-x86_64.tar.gz

onedir, not onefile: onefile unpacks to a temp directory on every launch, which
for a bundle this size adds seconds to startup and confuses antivirus. onedir
starts immediately and is what most Python desktop apps ship.
"""

import argparse
import os
import platform
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
NAME = 'MeshWorkbench'
VERSION = '3.0.0'   # tracks the project version in pyproject.toml
# This script lives in <repo>/packaging, so the project it packages is its parent.
REPO_DEFAULT = os.path.dirname(HERE)

# Packages PyInstaller's static analysis misses. gensim, sklearn and statsmodels
# all import submodules dynamically; matplotlib needs its backend explicitly.
HIDDEN = [
    'gensim', 'gensim.models', 'gensim.models.word2vec',
    'sklearn.manifold', 'sklearn.utils._typedefs', 'sklearn.neighbors._partition_nodes',
    'statsmodels.tsa.statespace._filters', 'statsmodels.tsa.statespace._smoothers',
    # scipy.special._cdflib was removed in scipy 1.15+; listing it only produces
    # a spurious ERROR line in the build log.
    'scipy._lib.messagestream',
    'matplotlib.backends.backend_agg', 'matplotlib.backends.backend_tkagg',
    'community', 'community.community_louvain',
    'Bio', 'Bio.Entrez', 'openpyxl', 'pymarc', 'plotly', 'kaleido',
    'seaborn', 'adjustText', 'tqdm', 'requests',
]
COLLECT_DATA = ['matplotlib', 'plotly', 'seaborn', 'statsmodels', 'sklearn', 'gensim']
EXCLUDE = ['PyQt5', 'PyQt6', 'PySide2', 'PySide6', 'IPython', 'jupyter',
           'notebook', 'pytest', 'sphinx']


def sh(cmd, **kw):
    print('  $ ' + ' '.join(str(c) for c in cmd[:6]) + (' …' if len(cmd) > 6 else ''))
    return subprocess.run(cmd, check=True, **kw)


def stage_pipeline(repo):
    """Copy the pipeline package and its shipped reference data next to the app."""
    src = os.path.join(repo, 'src', 'mesh_aop')
    if not os.path.isdir(src):
        sys.exit(f'pipeline package not found: {src}')
    dst = os.path.join(HERE, '_staged', 'mesh_aop')
    # ignore_errors: the staging area sits on a FAT32 volume where a leftover
    # .ipynb_checkpoints directory can refuse rmdir. Failing the whole build over
    # a stale scratch folder helps nobody, and copytree overwrites what matters.
    if os.path.exists(os.path.dirname(dst)):
        shutil.rmtree(os.path.dirname(dst), ignore_errors=True)
    skip = shutil.ignore_patterns('__pycache__', '*.pyc', '.ipynb_checkpoints',
                                  '*.ipynb', '.git')
    shutil.copytree(src, dst, ignore=skip, dirs_exist_ok=True)
    ref_src = os.path.join(repo, 'data', 'reference_processed')
    ref_dst = os.path.join(HERE, '_staged', 'reference_processed')
    if os.path.isdir(ref_src):
        shutil.copytree(ref_src, ref_dst, ignore=skip, dirs_exist_ok=True)
    cfg = os.path.join(repo, 'mesh_config.json')
    if os.path.exists(cfg):
        shutil.copy(cfg, os.path.join(HERE, '_staged', 'mesh_config.default.json'))
    n = sum(len(f) for _, _, f in os.walk(os.path.dirname(dst)))
    print(f'  staged {n} pipeline files')
    return os.path.join(HERE, '_staged')


def build_bundle(staged, clean=True, console=False):
    """console=True keeps a terminal attached.

    A --windowed build has nowhere to write a traceback: if Python raises while
    importing, the process exits quietly with no window, no message and no entry
    in the Application event log, which makes a startup failure impossible to
    diagnose from the outside. The console variant exists so that failure has a
    voice, and is what to reach for whenever the windowed build "does nothing".
    """
    sep = ';' if os.name == 'nt' else ':'
    name = NAME + ('-debug' if console else '')
    cmd = [sys.executable, '-m', 'PyInstaller', '--noconfirm',
           '--console' if console else '--windowed',
           '--name', name, '--distpath', os.path.join(HERE, 'dist'),
           '--workpath', os.path.join(HERE, 'build'),
           '--specpath', os.path.join(HERE, '_spec'),
           '--paths', HERE, '--paths', staged]
    if clean:
        cmd.append('--clean')
    for h in HIDDEN:
        cmd += ['--hidden-import', h]
    for c in COLLECT_DATA:
        cmd += ['--collect-data', c]
    for e in EXCLUDE:
        cmd += ['--exclude-module', e]
    cmd += ['--add-data', f'{os.path.join(staged, "mesh_aop")}{sep}mesh_aop']
    ref = os.path.join(staged, 'reference_processed')
    if os.path.isdir(ref):
        cmd += ['--add-data', f'{ref}{sep}reference_processed']
    dflt = os.path.join(staged, 'mesh_config.default.json')
    if os.path.exists(dflt):
        cmd += ['--add-data', f'{dflt}{sep}.']
    cmd.append(os.path.join(HERE, 'app_entry.py'))
    t0 = time.time()
    sh(cmd)
    out = os.path.join(HERE, 'dist', name)
    size = sum(os.path.getsize(os.path.join(b, f))
               for b, _, fs in os.walk(out) for f in fs) / 1e6
    print(f'\n  bundle: {out}\n  size  : {size:,.0f} MB   built in {time.time()-t0:.0f}s')
    return out


def package(bundle):
    """Wrap the bundle in whatever the host platform's users expect."""
    system = platform.system()
    out = os.path.join(HERE, 'dist')
    if system == 'Windows':
        iss = os.path.join(HERE, 'installer.iss')
        for exe in (r'C:\Program Files (x86)\Inno Setup 6\ISCC.exe',
                    r'C:\Program Files\Inno Setup 6\ISCC.exe'):
            if os.path.exists(exe):
                sh([exe, iss]); return
        print('\n  [i] Inno Setup not found - the bundle is complete and runnable,\n'
              '      but no setup.exe was produced. Install Inno Setup 6 and rerun\n'
              '      with --package, or ship dist\\%s as a zip.' % NAME)
    elif system == 'Darwin':
        app = os.path.join(out, NAME + '.app')
        dmg = os.path.join(out, f'{NAME}-{VERSION}-macos.dmg')
        if os.path.isdir(app):
            sh(['hdiutil', 'create', '-volname', NAME, '-srcfolder', app,
                '-ov', '-format', 'UDZO', dmg])
    else:
        tar = os.path.join(out, f'{NAME}-{VERSION}-linux-x86_64.tar.gz')
        sh(['tar', '-czf', tar, '-C', out, NAME])
        print(f'  {tar}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo', default=os.environ.get('MESH_REPO', REPO_DEFAULT))
    ap.add_argument('--package', action='store_true', help='also build the installer')
    ap.add_argument('--no-clean', action='store_true')
    ap.add_argument('--console', action='store_true',
                    help='build the diagnostic variant that keeps a terminal')
    a = ap.parse_args()
    print(f'MeSH AOP Workbench {VERSION} - building for {platform.system()} '
          f'{platform.machine()}')
    print(f'  python {sys.version.split()[0]}')
    staged = stage_pipeline(a.repo)
    bundle = build_bundle(staged, clean=not a.no_clean, console=a.console)
    if a.package:
        package(bundle)
    print('\nDone.')


if __name__ == '__main__':
    main()
