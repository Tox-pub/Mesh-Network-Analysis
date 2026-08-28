# -*- coding: utf-8 -*-
"""
make_third_party_notices.py - list what the installers redistribute, and under what.

Every bundle carries around fifty third-party libraries as prebuilt wheels. Most
are permissive and need only attribution; four are not, and one of those -
gensim, LGPL-2.1 - carries a real obligation. Shipping them without saying so is
a licence breach, not an oversight.

The list is taken from the wheels a bundle actually contains, not from the
development environment, because those differ: the environment here holds
build-time tools that are never redistributed.

    python packaging/make_third_party_notices.py <bundle.tar.gz>
    python packaging/make_third_party_notices.py <bundle.tar.gz> --check
"""

import argparse
import json
import posixpath
import re
import subprocess
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
OUT = REPO / 'THIRD-PARTY-NOTICES.md'

# Licences that oblige more than attribution, and what the obligation is.
COPYLEFT = {
    'LGPL': ('LGPL', 'Weak copyleft. The library is redistributed unmodified and '
                     'is dynamically imported, never linked into a derived work. '
                     'Its source and licence are at the URL below.'),
    'MPL':  ('MPL-2.0', 'File-level copyleft. The files are redistributed '
                        'unmodified; source for them is at the URL below.'),
    'GPL':  ('GPL', 'Strong copyleft. If this is here it needs deciding on, not '
                    'documenting.'),
}


def shipped_distributions(bundle):
    """The distribution names a bundle actually redistributes."""
    with tarfile.open(bundle) as tf:
        files = [posixpath.basename(n) for n in tf.getnames()
                 if '/wheels/' in n and (n.endswith('.whl') or n.endswith('.tar.gz'))]
    return sorted({re.split(r'-\d', f, 1)[0].replace('_', '-').lower() for f in files})


def environment_licences():
    """name -> {License, URL, Author, Version} for everything installed here."""
    raw = subprocess.run([sys.executable, '-m', 'piplicenses', '--format=json',
                          '--with-urls', '--with-authors'],
                         capture_output=True, text=True)
    if raw.returncode != 0:
        sys.exit('pip-licenses is needed:\n'
                 f'    {sys.executable} -m pip install pip-licenses')
    return {e['Name'].replace('_', '-').lower(): e for e in json.loads(raw.stdout)}


def classify(licence):
    """Copyleft or not, matched as a licence identifier rather than as text.

    Several projects put their entire licence into the metadata field, and a
    plain substring search finds "MPL" inside "IMPLIED" - which is in the
    warranty disclaimer of every MIT licence ever written. That flagged three
    MIT packages as weak copyleft. Word boundaries, and only the first line
    where the whole text is present.
    """
    head = licence.strip().splitlines()[0] if licence.strip() else ''
    subject = head if len(head) <= 60 and not head.startswith('#') else licence[:200]
    for key, (label, note) in COPYLEFT.items():
        if re.search(rf'\b(?:L?{key}|{key})\b', subject.upper()):
            return label, note
    return None, None


def tidy(licence):
    """Some projects put their whole licence text in the metadata field."""
    first = licence.strip().splitlines()[0].strip() if licence.strip() else ''
    if len(first) > 60 or first.startswith('#'):
        for known in ('MIT', 'BSD', 'Apache', 'MPL', 'LGPL', 'GPL', 'PSF'):
            if known.lower() in licence.lower():
                return known + ' (full text in the package)'
        return 'see the package'
    return first or 'not declared'


def build(bundle):
    shipped = shipped_distributions(bundle)
    known = environment_licences()

    rows, notable, missing = [], [], []
    for name in shipped:
        e = known.get(name)
        if not e:
            missing.append(name)
            continue
        lic = tidy(e['License'])
        label, note = classify(e['License'])
        rows.append((name, e.get('Version', ''), lic, e.get('URL', '') or ''))
        if label:
            notable.append((name, label, note, e.get('URL', '') or ''))

    lines = [
        '# Third-party notices',
        '',
        'MeSH Workbench is distributed as a self-contained bundle: the installer '
        'carries its own Python interpreter and every library it needs, as '
        'prebuilt wheels, unmodified from their published releases.',
        '',
        'This file lists what those are. It is generated from the wheels a built '
        'bundle actually contains — not from a development environment — by',
        '',
        '```bash',
        'python packaging/make_third_party_notices.py <bundle.tar.gz>',
        '```',
        '',
        'MeSH Workbench itself is MIT licensed; see [LICENSE](LICENSE).',
        '',
        '---',
        '',
        '## Licences that ask for more than attribution',
        '',
    ]
    if notable:
        for name, label, note, url in sorted(notable):
            lines += [f'### {name} — {label}', '', note, '']
            if url:
                lines += [f'Source and licence: <{url}>', '']
    else:
        lines += ['None: everything redistributed is permissively licensed.', '']

    lines += ['---', '',
              f'## Everything redistributed ({len(rows)} packages)', '',
              '| Package | Version | Licence | Project |',
              '| :--- | :--- | :--- | :--- |']
    for name, ver, lic, url in rows:
        link = f'<{url}>' if url and url.lower() != 'unknown' else ''
        lines.append(f'| {name} | {ver} | {lic} | {link} |')
    lines.append('')

    if missing:
        lines += ['', '> Not resolvable from the environment this was generated in: '
                  + ', '.join(missing) + '.', '']

    lines += ['---', '',
              'The bundled Python interpreter is CPython, redistributed under the '
              'Python Software Foundation License, from the '
              '[python-build-standalone](https://github.com/astral-sh/python-build-standalone) '
              'project on Linux and macOS and from the '
              '[python.org embeddable package](https://www.python.org/downloads/windows/) '
              'on Windows.', '']
    return '\n'.join(lines) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('bundle', help='a built .tar.gz bundle to read the wheel list from')
    ap.add_argument('--check', action='store_true',
                    help='fail if the file on disk is out of date, changing nothing')
    a = ap.parse_args()

    text = build(a.bundle)
    current = OUT.read_text(encoding='utf-8') if OUT.exists() else None
    if current == text:
        print(f'{OUT.name} is up to date')
        return 0
    if a.check:
        print(f'{OUT.name} is STALE - run: python packaging/make_third_party_notices.py <bundle>')
        return 1
    OUT.write_text(text, encoding='utf-8', newline='\n')
    print(f'wrote {OUT}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
