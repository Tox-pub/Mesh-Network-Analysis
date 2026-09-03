# -*- coding: utf-8 -*-
"""
console.py - is there a person who can answer a question?

One job, because getting it wrong hangs things. A prompt printed where nobody
can reply does not fail: it waits, and whatever launched it waits too. That is
how an MSI uninstall sat at "1 minute remaining" indefinitely, and how a paused
run died the moment it asked about annotations.

sys.stdin.isatty() is NOT the test on Windows. The NUL device is a character
device, so isatty() answers True for it - and NUL is exactly what an installer
or a GUI hands a child process. The only reliable question is whether the
standard input handle is a real console, which GetConsoleMode answers and
nothing else does.
"""

import sys

STD_INPUT_HANDLE = -10
INVALID_HANDLE_VALUE = -1


def has_console():
    """True only when a person could actually type an answer."""
    try:
        if sys.stdin is None:
            return False
        if not sys.stdin.isatty():
            return False                     # pipes and files, on any platform
    except (AttributeError, ValueError):
        return False                         # closed, or replaced with a stub

    if sys.platform != 'win32':
        return True

    try:
        import ctypes
        kernel32 = ctypes.windll.kernel32
        handle = kernel32.GetStdHandle(STD_INPUT_HANDLE)
        if handle in (0, INVALID_HANDLE_VALUE, ctypes.c_void_p(-1).value, None):
            return False
        mode = ctypes.c_uint()
        # Succeeds for a console handle and fails for everything else,
        # including NUL, which isatty() cannot tell apart from a terminal.
        return bool(kernel32.GetConsoleMode(handle, ctypes.byref(mode)))
    except Exception:                                          # noqa: BLE001
        return False


def ask(prompt, default=''):
    """input(), or the default when nothing can answer.

    Prints the question either way, so a log shows both what was asked and what
    was assumed when nobody was there to be asked.
    """
    if not has_console():
        print(prompt)
        print(f"      No console is attached, so this cannot be answered here - "
              f"continuing as '{default or 'skip'}'.")
        return default
    try:
        return input(prompt).strip().lower()
    except (EOFError, KeyboardInterrupt):
        print(f"\n      No answer - continuing as '{default or 'skip'}'.")
        return default
