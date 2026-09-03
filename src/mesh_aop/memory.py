# -*- coding: utf-8 -*-
"""
memory.py - how much memory this machine has, so the build can size itself to it.

One question, asked once: how many parser processes can the database build
afford to start. Each is a Python interpreter with an XML parser in it, and the
build must not start one per core regardless of whether the machine has room.

This does not monitor or report anything at run time. It answers that question
and gets out of the way.
"""

import os

_MB = 1024.0 * 1024.0


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


def fmt(n):
    """Bytes as a figure a person reads."""
    if n is None:
        return 'unknown'
    return f'{n / _MB:,.0f} MB' if n < 4 * 1024 * _MB else f'{n / (1024 * _MB):,.1f} GB'


def worker_budget(cores, per_worker_mb=700, reserve_mb=2048):
    """How many parser processes this machine can afford.

    The database build spawned one worker per core less one. On a workstation
    that is a dozen or more Python interpreters, each parsing its own XML file,
    and the memory cost is the sum of them. Cores are not the scarce resource
    here - the work is gzip decompression feeding a single SQLite writer that
    they all merge into one at a time - so capping this costs far less
    throughput than it sounds like.

    Returns at least 1: a machine too small by this arithmetic still has to run,
    just slowly.
    """
    total = total_ram()
    if not total:
        return max(1, min(cores, 8))
    usable = max(0, total / _MB - reserve_mb)
    return max(1, min(cores, int(usable // per_worker_mb) or 1))
