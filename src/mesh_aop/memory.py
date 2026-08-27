# -*- coding: utf-8 -*-
"""
memory.py - how much memory a stage actually used, reported by the stage itself.

Written because "the process is holding 16 GB" is not a measurement anyone can
act on. Task Manager shows one number, at one moment, for one process, and on
Windows it counts pages the operating system is caching on the program's behalf
and will hand back the instant anything else wants them. That number cannot
distinguish a leak from a large file being read, and the database build reads a
44 GB one.

So the pipeline measures itself: current and peak for this process, the same for
the worker processes it spawns, and the total across both. Peak is the figure
that matters - it is what the machine had to find, and it survives the fall
after a stage finishes.

Standard library only. psutil would do this in one line and is the obvious
choice for a program that already depends on it; adding it here would mean
another wheel in three bundles for one diagnostic.
"""

import os
import sys

_MB = 1024.0 * 1024.0


def _windows():
    """Current and peak working set for this process, in bytes."""
    import ctypes
    from ctypes import wintypes

    class _Counters(ctypes.Structure):
        _fields_ = [('cb', wintypes.DWORD),
                    ('PageFaultCount', wintypes.DWORD),
                    ('PeakWorkingSetSize', ctypes.c_size_t),
                    ('WorkingSetSize', ctypes.c_size_t),
                    ('QuotaPeakPagedPoolUsage', ctypes.c_size_t),
                    ('QuotaPagedPoolUsage', ctypes.c_size_t),
                    ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t),
                    ('QuotaNonPagedPoolUsage', ctypes.c_size_t),
                    ('PagefileUsage', ctypes.c_size_t),
                    ('PeakPagefileUsage', ctypes.c_size_t)]

    counters = _Counters()
    counters.cb = ctypes.sizeof(counters)

    # argtypes and restype are not optional here. GetCurrentProcess returns the
    # pseudo-handle -1, and without a declared signature ctypes passes it as a
    # 32-bit int on a 64-bit build - the call then fails, silently, and every
    # reading comes back empty.
    kernel32 = ctypes.windll.kernel32
    kernel32.GetCurrentProcess.restype = wintypes.HANDLE
    handle = kernel32.GetCurrentProcess()

    # The entry point moved: current Windows exports it from kernel32 as
    # K32GetProcessMemoryInfo, older ones only from psapi.dll.
    for get in (lambda: kernel32.K32GetProcessMemoryInfo,
                lambda: ctypes.WinDLL('psapi').GetProcessMemoryInfo):
        try:
            fn = get()
            fn.argtypes = [wintypes.HANDLE, ctypes.POINTER(_Counters), wintypes.DWORD]
            fn.restype = wintypes.BOOL
        except (AttributeError, OSError):
            continue
        try:
            if fn(handle, ctypes.byref(counters), counters.cb):
                return counters.WorkingSetSize, counters.PeakWorkingSetSize
        except (AttributeError, OSError, ctypes.ArgumentError):
            continue
    return None, None


def _linux():
    """Current and peak resident set from /proc, in bytes."""
    current = peak = None
    try:
        with open('/proc/self/status') as fh:
            for line in fh:
                if line.startswith('VmRSS:'):
                    current = int(line.split()[1]) * 1024
                elif line.startswith('VmHWM:'):      # high-water mark
                    peak = int(line.split()[1]) * 1024
    except OSError:
        pass
    return current, peak


def _rusage(who):
    """Peak resident set from getrusage, in bytes.

    ru_maxrss is kilobytes on Linux and bytes on macOS - a genuine difference
    between the platforms, not a mistake in either.
    """
    try:
        import resource
    except ImportError:
        return None
    try:
        raw = resource.getrusage(who).ru_maxrss
    except (ValueError, OSError):
        return None
    return raw if sys.platform == 'darwin' else raw * 1024


def usage():
    """Return {current, peak, children_peak, total_peak} in bytes; None if unknown.

    children_peak is only available on Unix, where the kernel tracks it. On
    Windows a spawned worker's peak dies with the worker, so the pool's own
    accounting has to supply it - see report_pool().
    """
    if os.name == 'nt':
        current, peak = _windows()
        children = None
    else:
        current, peak = _linux()
        if peak is None:
            peak = _rusage(getattr(__import__('resource'), 'RUSAGE_SELF', 0))
        children = _rusage(getattr(__import__('resource'), 'RUSAGE_CHILDREN', -1))

    total = None
    if peak is not None:
        total = peak + (children or 0)
    return {'current': current, 'peak': peak,
            'children_peak': children, 'total_peak': total}


def fmt(n):
    """Bytes as a figure a person reads, or a dash when it is not known."""
    if n is None:
        return '   -   '
    return f'{n / _MB:,.0f} MB' if n < 4 * 1024 * _MB else f'{n / (1024 * _MB):,.2f} GB'


def report(label, stream=None):
    """Print one line of memory accounting for a finished stage.

    Called after the expensive stages rather than on a timer: the point is the
    high-water mark, and that is already recorded by the time the stage ends.
    """
    stream = stream or sys.stdout
    u = usage()
    if u['peak'] is None and u['current'] is None:
        return u
    parts = [f"  Memory after {label}:",
             f"now {fmt(u['current'])}",
             f"peak {fmt(u['peak'])}"]
    if u['children_peak']:
        parts.append(f"workers peak {fmt(u['children_peak'])}")
        parts.append(f"total peak {fmt(u['total_peak'])}")
    print('   '.join(parts), file=stream)
    return u


def worker_budget(cores, per_worker_mb=700, reserve_mb=2048):
    """How many parser processes this machine can afford.

    The database build spawned one worker per core less one. On a workstation
    that is dozens of Python interpreters, each parsing its own XML file into
    its own batch list, and the memory is the sum of them - which is the
    likeliest explanation for a build that appeared to hold sixteen gigabytes.
    Cores are not the scarce resource here; the work is dominated by gzip
    decompression and by the single SQLite writer they all merge into, so
    capping this costs far less throughput than it sounds like.

    Returns at least 1 - a machine too small by this arithmetic still has to
    run, just slowly.
    """
    total = total_ram()
    if not total:
        return max(1, min(cores, 8))
    usable = max(0, total / _MB - reserve_mb)
    return max(1, min(cores, int(usable // per_worker_mb) or 1))


def total_ram():
    """Physical memory in bytes, or None where it cannot be established."""
    if os.name == 'nt':
        try:
            import ctypes
            from ctypes import wintypes

            class _Status(ctypes.Structure):
                _fields_ = [('dwLength', wintypes.DWORD),
                            ('dwMemoryLoad', wintypes.DWORD),
                            ('ullTotalPhys', ctypes.c_ulonglong),
                            ('ullAvailPhys', ctypes.c_ulonglong),
                            ('ullTotalPageFile', ctypes.c_ulonglong),
                            ('ullAvailPageFile', ctypes.c_ulonglong),
                            ('ullTotalVirtual', ctypes.c_ulonglong),
                            ('ullAvailVirtual', ctypes.c_ulonglong),
                            ('ullAvailExtendedVirtual', ctypes.c_ulonglong)]

            status = _Status()
            status.dwLength = ctypes.sizeof(status)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
                return status.ullTotalPhys
        except Exception:
            return None
        return None
    try:
        return os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
    except (ValueError, OSError, AttributeError):
        return None
