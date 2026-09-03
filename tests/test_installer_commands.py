"""Every command line an installer builds, split the way Windows splits it.

An MSI uninstall hung at "1 minute remaining" indefinitely. [INSTALLFOLDER]
always ends in a backslash, so &quot;[INSTALLFOLDER]&quot; emitted

    --project "C:\\...\\MeSH Workbench\\"

and Windows reads \\" as an ESCAPED QUOTE. The string never closed, --yes and
--keep-data were swallowed into the --project value, and the uninstaller
reached its "Type REMOVE to continue" prompt with no console to answer it.
It waited forever, and msiexec waited on it.

Two defences, both checked here. The arguments have to parse - verified with
Windows' own CommandLineToArgvW where available - and a prompt with no console
has to decline instead of blocking, so a future mistake in the first cannot
hang an uninstall again.
"""
import argparse
import os
import re
import sys

sys.path.insert(0, 'src')

FAILS = []


def ck(ok, msg, extra=''):
    print(('  [OK] ' if ok else '  [XX] ') + msg + (f'   -- {extra}' if extra else ''))
    if not ok:
        FAILS.append(msg)


def split(cmd):
    """Split a command line exactly as Windows does."""
    if sys.platform == 'win32':
        import ctypes
        f = ctypes.windll.shell32.CommandLineToArgvW
        f.restype = ctypes.POINTER(ctypes.c_wchar_p)
        f.argtypes = [ctypes.c_wchar_p, ctypes.POINTER(ctypes.c_int)]
        n = ctypes.c_int()
        p = f(cmd, ctypes.byref(n))
        return [p[i] for i in range(n.value)]
    # Same rule, for running the suite on a CI runner that is not Windows.
    out, cur, q, i = [], '', False, 0
    while i < len(cmd):
        c = cmd[i]
        if c == '\\' and i + 1 < len(cmd) and cmd[i + 1] == '"':
            cur += '"'
            i += 2
            continue
        if c == '"':
            q = not q
        elif c == ' ' and not q:
            if cur:
                out.append(cur)
            cur = ''
        else:
            cur += c
        i += 1
    if cur:
        out.append(cur)
    return out


def parsed(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument('--project')
    ap.add_argument('--yes', action='store_true')
    ap.add_argument('--keep-data', action='store_true')
    ap.add_argument('--results', action='store_true')
    a, _ = ap.parse_known_args(argv)
    return a


FOLDER = r'C:\Users\someone\AppData\Local\Programs\MeSH Workbench'

print('=== 1. the MSI, whose property always ends in a backslash ===')
wxs = open('packaging/windows_msi.wxs', encoding='utf-8').read()
vals = re.findall(r'Value="([^"]*uninstall_cli[^"]*)"', wxs)
ck(len(vals) == 2, f'both uninstall command lines found: {len(vals)}')
for v in vals:
    cmd = (v.replace('&quot;', '"')
            .replace('[INSTALLFOLDER]', FOLDER + '\\'))
    argv = split(cmd)
    a = parsed(argv[3:])
    tail = ' '.join(argv[4:])
    ck(a.yes, f'--yes survives: {tail}', f'argv={argv}')
    ck(a.project is not None and '"' not in a.project,
       f'--project is a clean path: {a.project!r}')
    ck(os.path.normpath(a.project).rstrip('\\') == FOLDER,
       f'and resolves to the install folder: {os.path.normpath(a.project)!r}')
    if '--keep-data' in v:
        ck(a.keep_data, '--keep-data survives too')

print('\n=== 2. the same line WITHOUT the fix must fail (the test can detect it) ===')
broken = ('"' + FOLDER + '\\python\\python.exe" -m mesh_aop.uninstall_cli '
          '--project "' + FOLDER + '\\" --yes --keep-data')
b = parsed(split(broken)[3:])
ck(not b.yes, 'the unfixed form loses --yes, as it did in the field')
ck('"' in (b.project or ''), 'and swallows the rest of the line into --project')

print('\n=== 3. the portable .bat launchers ===')
for name in ('Uninstall.bat', 'Install.bat'):
    body = open(os.path.join('packaging/launchers', name), encoding='utf-8').read()
    for line in body.splitlines():
        if 'uninstall_cli' not in line:
            continue
        # %~dp0 always ends in a backslash; %TARGET% may.
        cmd = (line.split('uninstall_cli', 1)[1]
                   .replace('%~dp0', FOLDER + '\\')
                   .replace('%TARGET%', FOLDER))
        a = parsed(split('x ' + cmd)[1:])
        ck(a.project is not None and '"' not in a.project,
           f'{name}: --project is clean: {a.project!r}')
        ck(os.path.normpath(a.project).rstrip('\\') == FOLDER,
           f'{name}: resolves to the folder', f'{os.path.normpath(a.project)!r}')
        # and again with a user-supplied trailing backslash
        cmd2 = (line.split('uninstall_cli', 1)[1]
                    .replace('%~dp0', FOLDER + '\\')
                    .replace('%TARGET%', FOLDER + '\\'))
        a2 = parsed(split('x ' + cmd2)[1:])
        ck(a2.project is not None and '"' not in a2.project,
           f'{name}: still clean when the path already ends in a backslash',
           f'{a2.project!r}')

print('\n=== 4. Inno Setup ===')
iss = open('packaging/windows_installer.iss', encoding='utf-8').read()
for m in re.findall(r"ExpandConstant\('([^']*uninstall_cli[^']*)'\)", iss):
    cmd = m.replace('{app}', FOLDER)          # Inno's {app} carries no trailing sep
    a = parsed(split('x ' + cmd)[1:])
    ck(a.yes, f'--yes survives: {m}')
    ck(a.project is not None and '"' not in a.project,
       f'--project is clean: {a.project!r}')

print('\n=== 5. the safety net: no console means decline, never block ===')
from mesh_aop import uninstall_cli                                 # noqa: E402


class NoConsole:
    def isatty(self):
        return False


import subprocess                                                   # noqa: E402
import tempfile                                                     # noqa: E402
import time                                                         # noqa: E402
from mesh_aop.console import has_console                            # noqa: E402

real = sys.stdin
try:
    sys.stdin = NoConsole()
    ck(has_console() is False, 'a redirected stdin is not a console')
    sys.stdin = None
    ck(has_console() is False, 'nor is a missing one')
finally:
    sys.stdin = real

# The case that actually mattered. On Windows the NUL device is a CHARACTER
# device, so isatty() calls it a terminal - and NUL is exactly what an
# installer hands a child process. A test covering only pipes would have
# passed while the bug shipped, which is what happened.
env = dict(os.environ, PYTHONPATH='src')
probe = ('import sys; sys.path.insert(0, "src");'
         'from mesh_aop.console import has_console;'
         'print(sys.stdin.isatty(), has_console())')
out = subprocess.run([sys.executable, '-c', probe], capture_output=True,
                     text=True, stdin=subprocess.DEVNULL, env=env).stdout.split()
if sys.platform == 'win32':
    ck(out[:1] == ['True'], f'isatty() still calls NUL a terminal: {out}')
ck(out[1:] == ['False'], f'has_console() does not: {out}')

# End to end: the exact invocation that hung for twenty minutes must return.
box = tempfile.mkdtemp()
t0 = time.time()
r = subprocess.run([sys.executable, '-m', 'mesh_aop.uninstall_cli', '--project', box],
                   capture_output=True, text=True, timeout=90,
                   stdin=subprocess.DEVNULL, env=env)
took = time.time() - t0
ck(took < 30, f'the MSI invocation returns instead of hanging ({took:.1f}s)')
ck(r.returncode == 1, f'and removes nothing: exit {r.returncode}')
ck('console' in r.stdout.lower(), 'saying why', r.stdout.strip()[-160:])
ck('Type REMOVE' not in r.stdout, 'without ever printing the prompt')

src = open('src/mesh_aop/uninstall_cli.py', encoding='utf-8').read()
# Comments mention input() when explaining why it is guarded, so compare on
# code only - otherwise the explanation reads as the thing being explained.
code = '\n'.join(ln for ln in src.splitlines() if not ln.strip().startswith('#'))
start = code.index('if not a.yes:')
i_console = code.index('_has_console()', start)
i_input = code.index('input(', start)
ck(i_console < i_input, 'the console check comes BEFORE the prompt',
   f'console at {i_console}, input at {i_input}')
ck('return 1' in code[i_console:i_input], 'and returns rather than asking')

print('\n=== 6. nothing else in the pipeline can block on a prompt ===')
for mod in ('cli.py', 'uninstall_cli.py', 'check_env.py'):
    body = open(os.path.join('src/mesh_aop', mod), encoding='utf-8').read()
    hits = []
    for n, line in enumerate(body.splitlines(), 1):
        if not re.search(r'(?<![\w_.])input\s*\(', line):
            continue
        if line.strip().startswith(('#', '*', '"""')) or 'input()' in line:
            continue
        # Far enough past the call to reach the except clause that guards it.
        window = '\n'.join(body.splitlines()[max(0, n - 14):n + 8])
        guarded = ('_has_console' in window or 'isatty' in window
                   or 'EOFError' in window or '_console_answer' in window)
        if not guarded:
            hits.append(f'{mod}:{n}')
    ck(not hits, f'{mod}: every prompt is guarded', f'unguarded at {hits}')

print()
if FAILS:
    print(f'FAILED ({len(FAILS)}):')
    for f in FAILS:
        print('   -', f)
    sys.exit(1)
print('all installer command-line checks passed')
