"""Every build script writes to one place, and refuses rather than guessing.

The artefacts ended up in four folders across two drives because each script
decided this for itself: the Windows portable build preferred D: and fell back
into packaging/ - inside the cloud-synced working copy - while the Unix bundles
defaulted to ~/Documents and the MSI wrote beside whatever tree it was given.

The fallback is the part that mattered. A build run without D: attached did not
fail; it quietly wrote half a gigabyte somewhere else, and nobody noticed until
there were two hundred files nobody could account for.
"""
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, 'packaging')

import build_location                                              # noqa: E402

FAILS = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + ('' if ok or not extra
                                                    else f'   -- {extra}'))
    if not ok:
        FAILS.append(msg)


print('=== 1. there is one location, and it is D: ===')
ck(build_location.BUILD_ROOT == r'D:\mesh_workbench_build',
   f'BUILD_ROOT is {build_location.BUILD_ROOT}')

print('\n=== 2. an explicit --out is honoured without argument ===')
tmp = tempfile.mkdtemp()
ck(build_location.resolve(tmp) == os.path.abspath(tmp),
   'a path passed in is used as given')

print('\n=== 3. MESH_BUILD_OUT is honoured - CI has no D: ===')
old = os.environ.get('MESH_BUILD_OUT')
os.environ['MESH_BUILD_OUT'] = tmp
try:
    ck(build_location.resolve() == os.path.abspath(tmp),
       'the environment override wins over the default')
finally:
    if old is None:
        os.environ.pop('MESH_BUILD_OUT', None)
    else:
        os.environ['MESH_BUILD_OUT'] = old

print('\n=== 4. with no D: and no override, it EXITS - never falls back ===')
# Run in a subprocess with a BUILD_ROOT that cannot exist, because the whole
# point is that resolve() calls sys.exit rather than returning a second choice.
probe = os.path.join(tmp, 'probe.py')
with open(probe, 'w', encoding='utf-8') as fh:
    fh.write(
        'import sys\n'
        f'sys.path.insert(0, {os.path.abspath("packaging")!r})\n'
        'import os\n'
        'os.environ.pop("MESH_BUILD_OUT", None)\n'
        'import build_location\n'
        'build_location.BUILD_ROOT = "\\\\\\\\?\\\\Q:\\\\definitely_not_a_drive"\n'
        'print("RESOLVED:", build_location.resolve())\n')
r = subprocess.run([sys.executable, probe], capture_output=True, text=True)
out = (r.stdout or '') + (r.stderr or '')
ck(r.returncode != 0, f'it exits non-zero: {r.returncode}')
ck('RESOLVED:' not in out, 'and returns no path at all', out[:200])
ck('not available' in out or 'needs' in out, 'saying what is missing', out[:200])
ck('Nothing was written' in out, 'and that nothing was written', out[:200])

print('\n=== 5. the fallback that scattered the output is gone ===')
src = open('packaging/build_portable_windows.py', encoding='utf-8').read()
ck('_default_build_out' not in src,
   'build_portable_windows.py no longer has its own answer')
ck("os.path.join(HERE, 'portable')" not in src,
   'and cannot land inside the cloud-synced project folder')
ck('build_location.resolve' in src, 'it asks build_location instead')

unix = open('packaging/build_unix_bundle.py', encoding='utf-8').read()
ck("'Documents', 'mesh_workbench_build'" not in unix,
   'build_unix_bundle.py no longer defaults to ~/Documents')
ck('build_location.resolve' in unix, 'it asks build_location too')

msi = open('packaging/build_msi_windows.py', encoding='utf-8').read()
ck('BUILD_OUT_DEFAULT' not in msi,
   'build_msi_windows.py no longer reads a constant that was deleted')
ck('build_location' in msi, 'it asks build_location as well')

print('\n=== 6. all three scripts still parse and show the same default ===')
for script in ('build_portable_windows.py', 'build_unix_bundle.py',
               'build_msi_windows.py'):
    r = subprocess.run([sys.executable, os.path.join('packaging', script), '--help'],
                       capture_output=True, text=True)
    ok = r.returncode == 0
    ck(ok, f'{script:<28} --help works', (r.stderr or '')[:160])
    if ok:
        ck('mesh_workbench_build' in r.stdout,
           f'{script:<28} names the shared location in its help')

import shutil                                                       # noqa: E402
shutil.rmtree(tmp, ignore_errors=True)

print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('every build script writes to one place, or stops')
