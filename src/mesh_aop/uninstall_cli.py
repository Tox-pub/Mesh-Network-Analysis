# -*- coding: utf-8 -*-
"""
uninstall_cli.py - `mesh-uninstall`, the terminal front-end to uninstall.py.

Identical behaviour on Windows, macOS and Linux. It shows what is on the machine
before asking anything, defaults to leaving results alone, and refuses to guess:
with nothing selected it reports and exits.
"""

import argparse
import os
import sys

from . import uninstall as U


def _fmt(n_bytes):
    if n_bytes >= 1e9:
        return f'{n_bytes / 1e9:,.2f} GB'
    if n_bytes >= 1e6:
        return f'{n_bytes / 1e6:,.1f} MB'
    return f'{n_bytes / 1e3:,.0f} KB'


def _print_inventory(items):
    if not items:
        print('\nNothing belonging to this program was found.\n')
        return
    width = max(len(i.label) for i in items)
    current = None
    for it in sorted(items, key=lambda i: (i.category, -i.bytes)):
        if it.category != current:
            current = it.category
            print(f'\n  {U.CATEGORY_LABEL[current]}')
        mark = '' if it.removable else '  [kept]'
        print(f'    {it.label:<{width}}  {_fmt(it.bytes):>12}{mark}   {it.path}')
        if it.note:
            print(f'    {"":<{width}}  {"":>12}   {it.note}')


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog='mesh-uninstall',
        description='Remove MeSH Workbench and the data it downloaded or built.')
    ap.add_argument('--project', default=None,
                    help='Project directory (default: the current directory).')
    ap.add_argument('--list', action='store_true',
                    help='Show what is on disk and exit, changing nothing.')
    ap.add_argument('--all', action='store_true',
                    help='Select everything found, including your results.')
    ap.add_argument('--keep-data', action='store_true',
                    help='Keep the downloaded archive and the built databases.')
    ap.add_argument('--results', action='store_true',
                    help='Also remove the results folder. Off by default.')
    ap.add_argument('--dry-run', action='store_true',
                    help='Report what would be removed without removing it.')
    ap.add_argument('--yes', action='store_true',
                    help='Do not ask for confirmation. Intended for scripts.')
    a = ap.parse_args(argv)

    root = os.path.abspath(a.project or os.getcwd())

    portable = U.is_portable(root)
    items = U.inventory(root)

    print(f'\nMeSH Workbench - uninstall\n  project: {root}')
    _print_inventory(items)
    total = sum(i.bytes for i in items)
    print(f'\n  {len(items)} item(s), {_fmt(total)} in total.')

    if a.list or not items:
        if items:
            print('\nNothing was changed. Re-run without --list to remove.\n')
        return 0

    # Entries flagged not-removable are reported but never selected, whatever
    # the options say - --all included.
    candidates = [i for i in items if i.removable]
    if a.all:
        chosen = list(candidates)
    else:
        keep = set()
        if a.keep_data:
            keep |= {U.DOWNLOADED, U.DERIVED}
        if not a.results:
            keep.add(U.RESULTS)
        chosen = [i for i in candidates if i.category not in keep]

    if not chosen:
        print('\nEverything found was excluded by the options given. '
              'Nothing to do.\n')
        return 0

    freed = sum(i.bytes for i in chosen)
    skipped = [i for i in items if i not in chosen]
    print(f'\n  Will remove {len(chosen)} item(s), reclaiming {_fmt(freed)}.')
    if skipped:
        print(f'  Keeping {len(skipped)}: ' +
              ', '.join(sorted({U.CATEGORY_LABEL[i.category] for i in skipped})))

    if a.dry_run:
        print('\n  --dry-run: nothing was removed.\n')
        return 0

    if not a.yes:
        print('\n  This cannot be undone.')
        try:
            if input('  Type REMOVE to continue: ').strip() != 'REMOVE':
                print('  Cancelled. Nothing was removed.\n')
                return 1
        except (EOFError, KeyboardInterrupt):
            print('\n  Cancelled. Nothing was removed.\n')
            return 1

    print()
    removed, freed, failures, deferred = U.remove(
        chosen, on_event=lambda kind, msg: print(
            f'  {"removed" if kind == "done" else "FAILED "}  {msg}'))

    print(f'\n  Removed {removed} of {len(chosen)}, reclaimed {_fmt(freed)}.')
    if deferred:
        print('\n  In use by this program itself - scheduled for removal as soon')
        print('  as it exits:')
        for path in deferred[:10]:
            print(f'    {path}')
    if failures:
        print(f'\n  {len(failures)} item(s) could not be removed. Close anything '
              'holding them open - a\n  file manager, an editor, or a sync client '
              'mid-upload - and run this again:')
        for path, err in failures[:10]:
            print(f'    {path}\n        {err}')

    if portable:
        print('\n  To finish, delete this folder.\n')
    else:
        heading, lines = U.removal_instructions()
        print(f'\n  {heading}')
        for line in lines:
            print(f'  {line}')
        print()
    return 0 if not failures else 1


if __name__ == '__main__':
    sys.exit(main())
