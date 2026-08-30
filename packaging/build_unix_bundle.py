# -*- coding: utf-8 -*-
"""
build_unix_bundle.py - a self-contained Linux or macOS package, built anywhere.

This can be run on Windows. That is worth stating plainly, because the received
wisdom is that you cannot build for one operating system on another, and for
compiled software that is true: PyInstaller, Nuitka and appimagetool all need to
run the target's own toolchain, and none of them can cross-compile.

None of that applies here, because this application has nothing to compile.
Every piece of it is one of two things:

  * pure Python - the application itself, which is the same bytes everywhere;
  * a prebuilt wheel - numpy, scipy, scikit-learn and the rest, for which PyPI
    already publishes compiled manylinux and macOS wheels. `pip download
    --platform` fetches them from any machine, because it is a download, not a
    build.

The interpreter is the one piece that is genuinely platform-specific, and it
does not have to be built either: python-build-standalone publishes
redistributable CPython for every target as a tarball. Statically linked,
relocatable, and - checked, because the window depends on it - carrying
_tkinter and its Tcl/Tk runtime.

So the whole package is an assembly job, and assembly is platform-agnostic.

WHAT THIS CANNOT DO is test the result. A bundle built here has never been run
on the system it targets. The CI workflow does that on real runners; this script
is for producing something today, on one machine, without waiting for a tag.

    python packaging/build_unix_bundle.py --target linux
    python packaging/build_unix_bundle.py --target macos-arm --out D:/builds
"""

import argparse
import io
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
NAME = 'MeSH-Workbench'
PY_SERIES = '3.12'

# Each target names three things: the CPython triple published by
# python-build-standalone, the pip platform tag its wheels are built for, and
# what a human calls it.
# wheel_tags are handed to pip ALL AT ONCE, not tried in turn - see
# fetch_wheels. Each project picks its own minimum OS version, so a real
# dependency set spans several tags and no single one describes it.
TARGETS = {
    'linux': dict(
        triple='x86_64-unknown-linux-gnu',
        wheel_tags=['manylinux_2_17_x86_64', 'manylinux2014_x86_64',
                    'manylinux_2_28_x86_64', 'manylinux_2_5_x86_64',
                    'manylinux1_x86_64'],
        label='linux-x86_64', pretty='Linux (Intel/AMD 64-bit)'),
    'linux-arm': dict(
        triple='aarch64-unknown-linux-gnu',
        wheel_tags=['manylinux_2_17_aarch64', 'manylinux2014_aarch64',
                    'manylinux_2_28_aarch64'],
        label='linux-arm64', pretty='Linux (ARM 64-bit)'),
    'macos-arm': dict(
        triple='aarch64-apple-darwin',
        wheel_tags=['macosx_11_0_arm64', 'macosx_12_0_arm64',
                    'macosx_13_0_arm64', 'macosx_14_0_arm64',
                    'macosx_10_9_universal2', 'macosx_10_13_universal2'],
        label='macos-arm64', pretty='macOS (Apple silicon)'),
    'macos-intel': dict(
        triple='x86_64-apple-darwin',
        wheel_tags=['macosx_10_13_x86_64', 'macosx_11_0_x86_64',
                    'macosx_12_0_x86_64', 'macosx_13_0_x86_64',
                    'macosx_10_9_universal2', 'macosx_10_13_universal2'],
        label='macos-x86_64', pretty='macOS (Intel)'),
}

PBS_API = 'https://api.github.com/repos/astral-sh/python-build-standalone/releases/latest'


def version():
    with open(os.path.join(REPO, 'pyproject.toml'), encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('version'):
                return line.split('"')[1]
    return '0.0.0'


def dependencies():
    """The runtime dependencies from pyproject, WITH their version pins.

    The pins were being stripped, so pip took whatever was newest on PyPI on
    the day of the build. That is how a build that worked on Monday failed on
    Tuesday with no change to this file: some project published a release, and
    whether it carried a wheel for this platform was nobody's decision.

    pyproject already states the tested ranges. Honouring them makes the bundle
    contain the versions the results were produced with, and makes two builds a
    week apart agree - which is the same guarantee the Windows build gets by
    copying out of a working virtual environment.

    python-louvain is separated out: it publishes no wheel, only an sdist, and
    pip refuses to mix --platform with a source download. It is pure Python, so
    the sdist is correct for every target anyway.
    """
    deps, inside = [], False
    with open(os.path.join(REPO, 'pyproject.toml'), encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('dependencies'):
                inside = True
                continue
            if inside:
                if line.strip().startswith(']'):
                    break
                spec = line.strip().rstrip(',').strip('"').strip("'").strip()
                if spec and not spec.startswith('#'):
                    deps.append(spec)

    def name_of(spec):
        for sep in ('>', '<', '=', '!', '~', '['):
            spec = spec.split(sep)[0]
        return spec.strip().lower()

    pure = [d for d in deps if name_of(d) == 'python-louvain']
    return [d for d in deps if d not in pure], pure


def _api_headers():
    """Headers for the GitHub API, authenticated when a token is available.

    Unauthenticated GitHub API calls are limited to 60 an hour PER IP, and CI
    runners share outbound addresses - macOS ones heavily. That limit is not
    reached by this build; it is reached by everybody else's build on the same
    address, and then this one gets 403 rate limit exceeded before it has done
    anything. Which is exactly how a green pipeline starts failing on one
    platform with no change to the code.
    
    A token raises the limit to 1000 an hour for the repository. Every GitHub
    Actions workflow is issued one automatically; nothing has to be created.
    Locally there is usually no token and none is needed - one developer does
    not approach 60 an hour.
    """
    headers = {'User-Agent': 'mesh-workbench-build'}
    token = (os.environ.get('GITHUB_TOKEN') or os.environ.get('GH_TOKEN') or '').strip()
    if token:
        headers['Authorization'] = f'Bearer {token}'
    return headers


def _api_json(url, attempts=4):
    """GET some JSON, waiting out a rate limit rather than dying on it.

    Even authenticated, a busy moment can return 403 or 429. The window is
    short and the build is long, so waiting is strictly better than failing a
    twenty-minute job at its first step.
    """
    last = None
    for attempt in range(attempts):
        req = urllib.request.Request(url, headers=_api_headers())
        try:
            return json.load(urllib.request.urlopen(req, timeout=60))
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code not in (403, 429):
                raise
            # GitHub says when the window resets; honour it when it is soon.
            reset = exc.headers.get('X-RateLimit-Reset')
            wait = 15 * (attempt + 1)
            if reset:
                try:
                    wait = max(1, min(120, int(reset) - int(time.time()) + 2))
                except ValueError:
                    pass
            remaining = exc.headers.get('X-RateLimit-Remaining', '?')
            print(f'    GitHub API rate limit (403/{exc.code}, remaining={remaining}); '
                  f'retrying in {wait}s  [{attempt + 1}/{attempts}]')
            time.sleep(wait)
    authed = 'Authorization' in _api_headers()
    sys.exit(
        f'\nGitHub API refused the request after {attempts} attempts: {last}\n'
        f'  authenticated : {authed}\n'
        + ('' if authed else
           '  Set GITHUB_TOKEN. Unauthenticated calls are limited to 60 an hour\n'
           '  per IP address, which CI runners share.\n'))


def fetch_python(triple, dest, stripped=True):
    """Download and unpack a redistributable CPython for the target."""
    print(f'  interpreter: looking up python-build-standalone')
    release = _api_json(PBS_API)
    suffix = 'install_only_stripped' if stripped else 'install_only'
    wanted = None
    for asset in release['assets']:
        n = asset['name']
        if (n.startswith(f'cpython-{PY_SERIES}.') and triple in n
                and n.endswith(f'-{suffix}.tar.gz')):
            wanted = asset
            break
    if not wanted:
        sys.exit(f'no {suffix} build for {triple} in release {release["tag_name"]}')

    print(f'  interpreter: {wanted["name"]}  ({wanted["size"]/1e6:.0f} MB)')
    # Cached: the interpreter is tens of megabytes and does not change between
    # targets rebuilt in the same session.
    cache_dir = os.path.join(HERE, '_cache')
    os.makedirs(cache_dir, exist_ok=True)
    cached = os.path.join(cache_dir, wanted['name'])
    if os.path.exists(cached) and os.path.getsize(cached) == wanted['size']:
        print('    (from cache)')
        blob = open(cached, 'rb').read()
    else:
        # The asset itself comes from a redirect to a CDN rather than from the
        # API, so it is not rate limited the same way - but a transient failure
        # here still costs the whole build, and this download is 30-45 MB.
        blob = None
        for attempt in range(4):
            try:
                req = urllib.request.Request(wanted['browser_download_url'],
                                             headers=_api_headers())
                blob = urllib.request.urlopen(req, timeout=180).read()
                break
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
                if attempt == 3:
                    sys.exit(f'\ncould not download {wanted["name"]}: {exc}')
                wait = 10 * (attempt + 1)
                print(f'    download failed ({exc}); retrying in {wait}s '
                      f'[{attempt + 2}/4]')
                time.sleep(wait)
        with open(cached, 'wb') as fh:
            fh.write(blob)
    with tarfile.open(fileobj=io.BytesIO(blob), mode='r:gz') as tf:
        # The archive contains a single top-level "python/" directory.
        tf.extractall(dest, filter='data')
    root = os.path.join(dest, 'python')
    if not os.path.isdir(root):
        sys.exit('the interpreter archive did not contain python/')

    # The window will not open without these, and a bundle that fails at that
    # point has already been shipped. Check now, on this machine.
    lib = os.path.join(root, 'lib', f'python{PY_SERIES}')
    checks = {
        'tkinter': os.path.isdir(os.path.join(lib, 'tkinter')),
        '_tkinter': any(f.startswith('_tkinter') for f in
                        os.listdir(os.path.join(lib, 'lib-dynload'))
                        if os.path.isdir(os.path.join(lib, 'lib-dynload'))),
        'sqlite3': os.path.isdir(os.path.join(lib, 'sqlite3')),
        'Tcl/Tk runtime': any(d.startswith(('tcl', 'tk')) for d in
                              os.listdir(os.path.join(root, 'lib'))),
    }
    for what, ok in checks.items():
        print(f'    {what:16} {"present" if ok else "MISSING"}')
    if not all(checks.values()):
        sys.exit('the interpreter is missing something the application needs')
    return root


def _run_pip(cmd, what):
    """Run pip and, if it fails, print what it said before giving up.

    check=True raises CalledProcessError, whose message is the exit code and
    nothing else - and with --quiet the output was thrown away too. That
    combination is why a macOS build failure arrived in CI as a bare
    "Process completed with exit code 1", with nothing to act on.
    """
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode == 0:
        return proc
    tail = lambda text: '\n'.join('    ' + line
                                  for line in text.strip().splitlines()[-20:])
    sys.exit('could not download %s\n\n  command:\n    %s\n\n  pip said:\n%s\n%s'
             % (what, ' '.join(cmd), tail(proc.stdout), tail(proc.stderr)))


def describe_host():
    """What this build is running on.

    First thing in the log, because the failure that cost two CI rounds
    differed from a working build only by the pip version.
    """
    import platform
    try:
        pip_v = subprocess.run([sys.executable, '-m', 'pip', '--version'],
                               capture_output=True, text=True).stdout.strip()
    except Exception:
        pip_v = '(pip not reachable)'
    print(f'  host     : {platform.system()} {platform.machine()}')
    print(f'  python   : {sys.version.split()[0]}')
    print(f'  {pip_v}')


def fetch_wheels(tags, dest):
    """Download every dependency as a wheel built for the target platform."""
    os.makedirs(dest, exist_ok=True)
    wheeled, pure = dependencies()
    base = [sys.executable, '-m', 'pip', 'download', '--quiet',
            '--only-binary=:all:', '--python-version', PY_SERIES, '--dest', dest]

    # Every tag in ONE call. pip accepts --platform repeatedly and takes a
    # wheel matching any of them, which is the only way to satisfy a set whose
    # members target different minimum OS versions - and they do: the macOS
    # wheels here span macosx_10_13_universal2, macosx_11_0_arm64 and
    # macosx_12_0_arm64 at once. Trying the tags in turn asks each one to
    # satisfy every package alone, which no single tag can, and it failed on
    # the runner while passing on the developer's older, laxer pip.
    cmd = base + [arg for tag in tags for arg in ('--platform', tag)] + wheeled
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        # Say what pip said. Exiting with only a summary is what made this
        # fail invisibly in CI, with nothing in the log to act on.
        tail = lambda text: '\n'.join('    ' + line
                                      for line in text.strip().splitlines()[-15:])
        sys.exit('could not resolve wheels for: %s\n\n  pip said:\n%s\n%s'
                 % (', '.join(tags), tail(proc.stdout), tail(proc.stderr)))
    print(f'  wheels: {len(os.listdir(dest))} across {len(tags)} platform tags')

    # The pure-Python stragglers, which have no platform to speak of.
    if pure:
        _run_pip([sys.executable, '-m', 'pip', 'download', '--no-deps',
                  '--dest', dest] + pure,
                 'the source-only dependencies (%s)' % ', '.join(pure))
        print(f'  wheels: + {len(pure)} pure-Python source distribution(s)')

        # An sdist has to be BUILT, and the standalone interpreter ships pip
        # but not setuptools - Python 3.12 dropped it from the default install.
        # Without these the first run fails offline with "No module named
        # setuptools", or silently reaches for PyPI and defeats the point.
        _run_pip([sys.executable, '-m', 'pip', 'download', '--no-deps',
                  '--dest', dest, 'setuptools', 'wheel'],
                 'setuptools and wheel')
        print('  wheels: + setuptools and wheel, to build that sdist offline')
    return dest


LAUNCH_SH = """#!/bin/sh
# MeSH Workbench - launcher.
#
# Two jobs before the window opens.
#
# 1. macOS quarantine. Anything downloaded through a browser is tagged
#    com.apple.quarantine, and Gatekeeper then refuses to run the unsigned
#    interpreter inside this folder - "cannot be opened because the developer
#    cannot be verified", with no obvious way forward. Clearing the tag needs
#    no password, no Apple account and no money; it is one command, and this
#    script runs it on the folder it lives in. THIS script is not blocked
#    itself, because a shell script is read by /bin/sh, which Apple signs.
#
# 2. First-run install. The libraries are unpacked from wheels/ into the
#    bundled interpreter, offline. The application itself is not installed at
#    all - it is put on PYTHONPATH, so there is no build step and nothing to
#    go wrong on a machine with no compiler.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
PY="$HERE/python/bin/python3.12"
STAMP="$HERE/.installed"

if [ "$(uname -s)" = "Darwin" ] && command -v xattr >/dev/null 2>&1; then
    if xattr -p com.apple.quarantine "$PY" >/dev/null 2>&1; then
        echo "macOS has quarantined this download. Clearing that flag..."
        if xattr -dr com.apple.quarantine "$HERE" 2>/dev/null; then
            echo "  Done - no password or Apple account needed."
        else
            echo "  Could not clear it automatically. Run this once, by hand:" >&2
            echo "" >&2
            echo "      xattr -dr com.apple.quarantine \\"$HERE\\"" >&2
            echo "" >&2
            exit 1
        fi
    fi
fi

if [ ! -f "$STAMP" ]; then
    echo "First run: unpacking the libraries (about a minute, no network needed)..."
    "$PY" -m pip install --quiet --no-index --find-links "$HERE/wheels" \
        --no-warn-script-location -r "$HERE/requirements.txt"
    touch "$STAMP"
    echo "Done."
fi

# The application runs from source on the path: pure Python, nothing to build.
PYTHONPATH="$HERE/app/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONPATH
exec "$PY" -m mesh_workbench "$@"
"""

UNINSTALL_SH = """#!/bin/sh
# Remove what this program put OUTSIDE this folder.
#
# The bundle is self-contained, but the things it produces are not: settings
# and downloaded data go under ~/.local/share/MeSH Workbench, and results go to
# ~/Documents/MeSH Workbench, so that a database survives replacing the program
# and two people on one machine do not share a results folder.
#
# That means deleting this folder is NOT a complete uninstall - it can leave
# tens of gigabytes behind. This lists what exists, asks, and removes it.
# Afterwards, delete this folder and nothing of the program remains.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
PY="$HERE/python/bin/python3.12"

if [ ! -f "$HERE/.installed" ]; then
    echo "First run: unpacking the libraries (no network needed)..."
    "$PY" -m pip install --quiet --no-index --find-links "$HERE/wheels" \\
        --no-warn-script-location -r "$HERE/requirements.txt"
    touch "$HERE/.installed"
fi

PYTHONPATH="$HERE/app/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONPATH
"$PY" -m mesh_aop.uninstall_cli "$@"

echo
echo "That covered everything outside this folder."
echo "To finish, delete the folder itself:"
echo
echo "    rm -rf \\"$HERE\\""
echo
"""

PIPELINE_SH = """#!/bin/sh
# The pipeline, without the window - for a machine with no desktop.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
PY="$HERE/python/bin/python3.12"

if [ "$(uname -s)" = "Darwin" ] && command -v xattr >/dev/null 2>&1; then
    xattr -dr com.apple.quarantine "$HERE" 2>/dev/null || true
fi

if [ ! -f "$HERE/.installed" ]; then
    echo "First run: unpacking the libraries (no network needed)..."
    "$PY" -m pip install --quiet --no-index --find-links "$HERE/wheels" \
        --no-warn-script-location -r "$HERE/requirements.txt"
    touch "$HERE/.installed"
fi

PYTHONPATH="$HERE/app/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONPATH
exec "$PY" -m mesh_aop.cli "$@"
"""


def repack(target, out_dir):
    """Re-archive a staging folder that is already assembled.

    Only the permission bookkeeping and the compression change, so this is the
    fast path when the packaging is being corrected rather than the contents.
    """
    spec = TARGETS[target]
    stem = f'{NAME}-{version()}-{spec["label"]}'
    staging = os.path.join(out_dir, stem)
    if not os.path.isdir(staging):
        sys.exit(f'nothing to repack: {staging} does not exist')
    archive = os.path.join(out_dir, stem + '.tar.gz')
    if os.path.exists(archive):
        os.remove(archive)
    print(f'\nRe-packing {spec["pretty"]}')
    t0 = time.time()
    with tarfile.open(archive, 'w:gz', compresslevel=6) as tf:
        tf.add(staging, arcname=stem, filter=_unix_permissions)
    print(f'  archive: {archive}  ({os.path.getsize(archive)/1e6:,.0f} MB)'
          f'  in {time.time()-t0:.0f}s')
    return archive


def _unix_permissions(info):
    """Restore the executable bit, which Windows cannot store.

    NTFS has no execute permission, so everything unpacked and re-packed here
    comes out 0666 - and on the far side `./MeSH Workbench` fails with
    "Permission denied", as does the interpreter it would have called. tar
    carries the mode, so the mode is reconstructed from what the file is.

    This is the one genuine difficulty in building a Unix package on Windows,
    and it is a bookkeeping problem rather than a compilation one.
    """
    name = info.name
    if info.isdir():
        info.mode = 0o755
    elif ('/python/bin/' in name
            or name.endswith(('.so', '.dylib'))
            or '.so.' in os.path.basename(name)
            or os.path.basename(name) in ('MeSH Workbench', 'mesh-pipeline', 'mesh-uninstall')):
        info.mode = 0o755
    else:
        info.mode = 0o644
    # Files built on one machine should not carry that machine's account.
    info.uid = info.gid = 0
    info.uname = info.gname = 'root'
    return info


def build(target, out_dir, stripped=True):
    spec = TARGETS[target]
    ver = version()
    stem = f'{NAME}-{ver}-{spec["label"]}'
    staging = os.path.join(out_dir, stem)
    shutil.rmtree(staging, ignore_errors=True)
    os.makedirs(staging, exist_ok=True)

    print(f'\nBuilding {spec["pretty"]}  ({stem})')
    t0 = time.time()

    fetch_python(spec['triple'], staging, stripped=stripped)
    fetch_wheels(spec['wheel_tags'], os.path.join(staging, 'wheels'))

    # The application source, installed at first run from this directory.
    app = os.path.join(staging, 'app')
    os.makedirs(app, exist_ok=True)
    # mesh_config.json is a developer's own settings file - it is gitignored
    # for that reason, and one was being copied into the bundle carrying this
    # project's search term. Empty this time; a configured tree would have
    # shipped its NCBI e-mail and API key to everyone who downloaded it.
    skip = shutil.ignore_patterns('__pycache__', '*.pyc', '.ipynb_checkpoints',
                                  'mesh_config.json', 'mesh_config.json.bak*',
                                  '.env', '*.pyo')
    shutil.copytree(os.path.join(REPO, 'src'), os.path.join(app, 'src'),
                    dirs_exist_ok=True, ignore=skip)

    # Our own licence travels with the program. Redistributing CPython, Tcl/Tk
    # and forty compiled wheels without stating what any of it is licensed
    # under is not a defensible thing to publish.
    for f in ('LICENSE', 'THIRD-PARTY-NOTICES.md', 'CITATION.cff',
                 'pyproject.toml', 'README.md', 'HELP.md', 'INSTALL.md'):
        src = os.path.join(REPO, f)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(app, f))

    # The reference corpus, so the program can draw its figures before anyone
    # has downloaded 44 GB of PubMed.
    ref = os.path.join(REPO, 'data', 'reference_processed')
    if os.path.isdir(ref):
        shutil.copytree(ref, os.path.join(staging, 'reference_processed'),
                        dirs_exist_ok=True, ignore=skip)
    anno = os.path.join(REPO, 'data', 'raw', 'aop_annotations_master.csv')
    if os.path.exists(anno):
        shutil.copy2(anno, os.path.join(staging, 'reference_processed',
                                        'aop_annotations_master.csv'))

    # A requirements file naming exactly what is in wheels/, so the offline
    # install resolves to the bundled files and nothing else.
    wheeled, pure = dependencies()
    with open(os.path.join(staging, 'requirements.txt'), 'w',
              encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(wheeled + pure) + '\n')

    for fname, body in (('MeSH Workbench', LAUNCH_SH),
                        ('mesh-pipeline', PIPELINE_SH),
                        ('mesh-uninstall', UNINSTALL_SH)):
        path = os.path.join(staging, fname)
        with open(path, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write(body)
        os.chmod(path, 0o755)

    is_mac = target.startswith('macos')
    gatekeeper = '''
BEFORE THE FIRST RUN, ON macOS
    macOS tags everything downloaded through a browser as quarantined, and
    refuses to run unsigned programs carrying that tag. This program is
    unsigned - signing requires a paid Apple Developer ID - so you will meet
    that.

    Clearing the tag needs no password, no Apple account and nothing to pay.
    The launcher does it for you: open Terminal, change to this folder, and

        ./"MeSH Workbench"

    The launcher itself is a shell script, which macOS reads with its own
    signed /bin/sh, so it is not blocked and can clear the flag from
    everything else here.

    Start it from Terminal, not by double-clicking in Finder - Finder will
    either open this in a text editor or refuse it outright. After the first
    run, either works.

    To clear the flag by hand instead:

        xattr -dr com.apple.quarantine "<this folder>"

    Nothing here is hidden from you: that command removes one extended
    attribute from these files and touches nothing else on your Mac.
''' if is_mac else ''

    with open(os.path.join(staging, 'README.txt'), 'w',
              encoding='utf-8', newline='\n') as fh:
        fh.write(f'''MeSH Workbench {ver} - {spec["pretty"]}
{gatekeeper}
TO RUN IT

    ./"MeSH Workbench"

That is all. This folder carries its own Python, so nothing needs to be
installed first - no system Python, no tkinter package, no administrator
rights. The first run unpacks the libraries from wheels/ and takes about a
minute; it needs no network.

Without a desktop, the pipeline runs on its own:

    ./mesh-pipeline --step all

TO CHECK IT WORKS WITHOUT DOWNLOADING 44 GB

    Turn on "Use bundled reference data" on the Folders tab of Settings, then
    run Step 4 - Figures. The reference corpus is in this folder, so that
    draws every figure and the PRISMA report from data already on disk.

WHERE YOUR FILES GO

    The program is self-contained, but what it produces is not - so that a
    database survives replacing the program, and two people on one machine do
    not share a results folder. On the first run you are asked where these
    should live. Answer with anything you like; the defaults are:

        ~/.local/share/MeSH Workbench/data      downloads and databases
                                                THIS IS THE BIG ONE, ~52 GB
        ~/Documents/MeSH Workbench              your results and figures
        ~/.local/share/MeSH Workbench           settings and logs

    The Folders tab shows the paths actually in use at any time.

TO UNINSTALL

    ./mesh-uninstall

    That lists everything the program put outside this folder, tells you how
    much space each takes, and asks before removing anything. Then delete this
    folder and nothing remains.

    Deleting the folder ALONE is not enough - it can leave tens of gigabytes
    of downloaded data behind.

The folder itself can be moved or copied anywhere, including onto a USB stick.
''')

    archive = os.path.join(out_dir, stem + '.tar.gz')
    if os.path.exists(archive):
        os.remove(archive)
    print('  packing ...')
    with tarfile.open(archive, 'w:gz', compresslevel=6) as tf:
        tf.add(staging, arcname=stem, filter=_unix_permissions)
    size = os.path.getsize(archive) / 1e6
    print(f'  archive: {archive}  ({size:,.0f} MB)  in {time.time()-t0:.0f}s')
    return archive


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--target', choices=sorted(TARGETS), default='linux')
    ap.add_argument('--all', action='store_true', help='build every target')
    ap.add_argument('--out', default=os.path.join(os.path.expanduser('~'),
                                                  'Documents', 'mesh_workbench_build'))
    ap.add_argument('--full-python', action='store_true',
                    help='use the unstripped interpreter (larger, with symbols)')
    ap.add_argument('--repack', action='store_true',
                    help='re-pack an already-assembled staging folder, no downloads')
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    describe_host()
    targets = sorted(TARGETS) if a.all else [a.target]
    if a.repack:
        made = [repack(t, a.out) for t in targets]
    else:
        made = [build(t, a.out, stripped=not a.full_python) for t in targets]

    print('\nBuilt:')
    for m in made:
        print(f'  {os.path.basename(m):<48} {os.path.getsize(m)/1e6:,.0f} MB')
    print('\nNot tested: these were assembled on this machine and have not been'
          '\nrun on the system they target. The release workflow does that on'
          '\nreal runners.')


if __name__ == '__main__':
    main()
