"""Telling someone how to remove the program has to be right for THEIR copy.

On Linux the uninstaller ended with:

    The package itself is still installed. To remove it:
        python3.12 -m pip uninstall mesh_aop_network

Both halves were wrong. `python3.12` came from Path(sys.executable).name, and
inside a bundle that names an interpreter which exists only in the bundle - so
it was pasted into a shell that had no such command. And the bundle never
pip-installs the application at all; it puts it on PYTHONPATH, so there was no
distribution to uninstall even if the command had run.

Nothing here touches the filesystem: sys.executable is pointed at a fake tree
of each shape and the advice is read back.
"""
import os
import sys
import tempfile

sys.path.insert(0, 'src')

from mesh_aop import uninstall as U                                # noqa: E402

FAILS = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + (f'   -- {extra}' if extra else ''))
    if not ok:
        FAILS.append(msg)


def fake_bundle(box, launcher, exe_rel):
    """A bundle laid out the way the packaging scripts lay one out."""
    root = os.path.join(box, 'MeSH-Workbench-3.2.0')
    os.makedirs(os.path.join(root, 'app'), exist_ok=True)
    exe = os.path.join(root, exe_rel)
    os.makedirs(os.path.dirname(exe), exist_ok=True)
    open(exe, 'w').close()
    open(os.path.join(root, launcher), 'w').close()
    return root, exe


box = tempfile.mkdtemp()
real_exe, real_platform = sys.executable, sys.platform

print('=== 1. the command names the interpreter that actually exists ===')
try:
    root, exe = fake_bundle(box, 'MeSH Workbench', os.path.join('python', 'bin', 'python3.12'))
    sys.executable = exe
    hint = U.pip_hint()
    ck('python3.12 -m pip' != hint.split(os.sep)[-1][:len('python3.12 -m pip')]
       or os.sep in hint,
       f'the hint is a path, not a bare name')
    ck(exe in hint or f'"{exe}"' in hint,
       f'and it is THIS interpreter: {hint}')
    ck(not hint.startswith('python3.12'),
       'it no longer starts with a bare python3.12', hint)

    spaced_box = os.path.join(box, 'has space')
    os.makedirs(spaced_box, exist_ok=True)
    root2, exe2 = fake_bundle(spaced_box, 'MeSH Workbench',
                              os.path.join('python', 'bin', 'python3.12'))
    sys.executable = exe2
    ck(U.pip_hint().startswith('"'), f'a path with a space is quoted: {U.pip_hint()}')
finally:
    sys.executable = real_exe

print('\n=== 2. a self-contained bundle is a folder to delete, not a pip command ===')
try:
    root, exe = fake_bundle(box, 'mesh-pipeline', os.path.join('python', 'bin', 'python3.12'))
    sys.executable = exe
    ck(U.bundle_root() is not None, f'the bundle is recognised: {U.bundle_root()}')
    ck(str(U.bundle_root()) == root, 'and it is the folder holding the launcher',
       f'{U.bundle_root()} != {root}')

    sys.platform = 'linux'
    heading, lines = U.removal_instructions()
    body = '\n'.join(lines)
    ck('self-contained' in heading, f'heading: {heading}')
    ck('rm -rf' in body, 'it says to delete the folder')
    ck(root in body, 'and names the folder')
    ck('pip uninstall' not in body,
       'and does NOT tell the user to pip uninstall anything', body)
    ck('never' in body or 'nothing to uninstall' in body.lower(),
       'and says why there is nothing to uninstall')
finally:
    sys.executable, sys.platform = real_exe, real_platform

print('\n=== 3. Windows says Add/Remove, not a shell command ===')
try:
    root, exe = fake_bundle(box, 'MeSH Workbench.bat', os.path.join('python', 'python.exe'))
    sys.executable = exe
    sys.platform = 'win32'
    heading, lines = U.removal_instructions()
    body = '\n'.join(lines)
    ck(root in body, 'a portable Windows copy names its folder')
    ck('Settings' in body and 'Uninstall' in body,
       'and points at Add/Remove for an installed one', body)
finally:
    sys.executable, sys.platform = real_exe, real_platform

print('\n=== 4. a real pip install still gets the pip command ===')
# this interpreter, running from a source checkout, is not a bundle
ck(U.bundle_root() is None, f'a source checkout is not a bundle: {U.bundle_root()}')
heading, lines = U.removal_instructions()
body = '\n'.join(lines)
if U.package_is_installed():
    ck('pip uninstall' in body, 'an installed package is removed with pip')
    ck(sys.executable in body or f'"{sys.executable}"' in body,
       'using this interpreter by full path')
else:
    ck('pip' not in body or 'not installed' in heading.lower()
       or 'not installed by a package manager' in heading.lower(),
       f'an uninstalled package is not told to pip uninstall: {heading}')
ck(bool(heading) and bool(lines), 'there is always something to say')

print('\n=== 5. both front ends use it ===')
import inspect                                                     # noqa: E402
from mesh_aop import uninstall_cli                                 # noqa: E402
cli_src = inspect.getsource(uninstall_cli.main)
ck('removal_instructions' in cli_src, 'the console uninstaller uses it')
ck('pip_hint()' not in cli_src, 'and no longer prints the bare hint')
app_src = open('src/mesh_workbench/app.py', encoding='utf-8').read()
ck('removal_instructions' in app_src, 'the application dialog uses it')
ck('U.pip_hint()' not in app_src, 'and no longer prints the bare hint')

print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('all uninstall-advice checks passed')
