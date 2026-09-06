# Pre-release checklist

What has to be true before a build is put in front of the public. Written
against 3.2.10; the same list applies to every release after it.

Status as of commit `31f02db` on `foxglove`:

| | |
| :--- | :--- |
| **Done** | 21 items |
| **Blocking** — must be true before publishing | 6 items |
| **Should** — publish without them and it shows | 5 items |
| **Later** — real, not urgent | 4 items |

Updated after the first Linux test run. Nine defects came back from it; all
nine are fixed in `232a61a`, `efb6aa6` and `31f02db`, and the sections below
record which were real and what they turned out to be. One of them — D0 —
made data collection impossible for every user on a default install, so
nothing could have been published before it was found.

---

## A. The download must not contain anything of yours

| | Item | Status |
| :-- | :--- | :--- |
| A1 | `mesh_config.json` kept out of both bundles | **done** — `b87f231` |
| A2 | Verifier fails the bundle if it comes back | **done** — Unix only |
| A3 | Same check for the Windows zip and MSI | **blocking** |
| A4 | No `data/`, `results/`, `.env`, `*.bak` in the tree | **done** (Windows purge, Unix ignore) |
| A5 | Verified on a rebuilt bundle | **done** — Linux 36/36, Windows tree clean |

**A3 is the gap.** `build_portable_windows.py` now purges the settings file,
but nothing inspects the finished `.zip` or `.msi` to confirm it worked — which
is exactly the shape of the bug that shipped a 345 KB installer as a green
build. Until there is a `verify_windows_bundle.py`, the manual check is:

```bash
python -c "import zipfile,sys; z=zipfile.ZipFile(sys.argv[1]); print([n for n in z.namelist() if 'mesh_config' in n or n.endswith('.env')] or 'clean')" MeSH-Workbench-3.2.10-win64-portable.zip
```

This matters more than it looks. Your NCBI e-mail and API key live in that
file. A build made from a configured tree publishes them to anyone who
downloads it, and an API key in a public artifact is a key you have to rotate.

---

## B. Legal and attribution

| | Item | Status |
| :-- | :--- | :--- |
| B1 | `LICENSE` exists (MIT) | **done** |
| B2 | `LICENSE` travels inside both bundles | **done** — `b87f231` |
| B3 | `THIRD-PARTY-NOTICES.md` | **blocking** |
| B4 | `CITATION.cff` | **should** |
| B5 | Copyright line says what you want it to say | **check** — currently `Copyright (c) 2026 Tox-pub` |

**B3.** Each bundle redistributes ~40 third-party libraries as binaries. Most
are permissive (BSD/MIT/Apache) and need only an attribution file. One is not:
**gensim is LGPL-2.1**, and redistributing it obliges you to say so and to
point at its source. The wheels are unmodified upstream, which makes compliance
a document rather than an engineering problem — but the document has to exist.
Generate the inventory from the build environment:

```bash
python -m pip install pip-licenses && pip-licenses --format=markdown --with-urls --with-license-file > THIRD-PARTY-NOTICES.md
```

**B4.** This is a research tool aimed at people who publish. Without a
`CITATION.cff` at the repo root, GitHub shows no "Cite this repository" button
and every citation of it will be formatted differently. It is a twelve-line
file and it is the cheapest thing on this list.

---

## C. It has to actually run, on a machine that is not this one

| | Item | Status |
| :-- | :--- | :--- |
| C1 | Linux: extract and run on a real Linux box | **blocking** — you, flash drive |
| C2 | macOS: extract and run on a real Mac past Gatekeeper | **blocking** — no Mac available yet |
| C3 | Windows MSI on a clean machine with no Python | **blocking** |
| C4 | CI green on all three at the current commit | **blocking** — not re-run since the MSI fix |
| C5 | One full end-to-end run from an installed copy | **blocking** |
| C6 | CI: Linux extract, offline install, imports, Tk window | **done** |
| C7 | CI: MSI size guard (fails under 40 MB) | **done, unproven** |
| C8 | Linux font rendering | **fixed** — `efb6aa6`; confirm on the machine |
| C9 | Retrieval works without an NCBI API key | **done** — proven live |

**C4.** The Windows job went green while producing a 345 KB installer wrapped
around an empty tree. The guard for that is written and committed but has never
run. Push `foxglove` and start the workflow by hand from the Actions tab (no
`publish` tick) — that builds and tests all three without releasing anything.

**C5** is the one that is easy to skip and shouldn't be. Everything tested so
far proves the *package* is sound: the interpreter starts, the imports resolve,
a window opens. Nothing has yet proved the *application* works when installed
— the first-run dialog, the folder settings, the pause-for-annotation popup,
abort, the integrity scan, the ledger, the PRISMA report. Do it from the
installed copy, not the dev tree, using bundled reference data so it takes a
minute instead of a weekend:

> Settings → Folders → tick "Use bundled reference data" → step `viz` → Run.
> Expect eight figures and a PRISMA flow report in the results folder.

**C2.** No Mac, no test. The bundle carries the workaround (the launcher clears
`com.apple.quarantine` itself, no password and nothing to pay) and CI builds and
runs it on a macOS runner — but a CI runner has Gatekeeper's quarantine flag
absent, because the file never arrived through a browser. The one thing that
cannot be tested without a Mac is the exact case a real user hits. Either
borrow a Mac for ten minutes, or publish the macOS build labelled as untested
on a real download. Do not publish it silently as if it had been.

---

## D. Things the public will hit that we haven't handled

| | Item | Status |
| :-- | :--- | :--- |
| D0 | Data collection returned HTTP 400 without an API key | **done** — `232a61a` |
| D1 | NCBI request rate is illegal without an API key | **done** — `232a61a` |
| D2 | 52 GB data-folder warning | **done** |
| D3 | Folders tab explains data vs results | **done** |
| D4 | Search term / citation generation help text | **done** |
| D5 | Getting an NCBI API key is documented | **done** — the run prints the link |
| D6 | What to do when a run breaks | **done** — integrity scan, repair from checkpoint |

**D0 was the one that mattered.** `Entrez.api_key` was assigned
unconditionally, so a user without a key sent NCBI `&api_key=` — an empty
parameter, which it answers with 400 Bad Request. A blank key is the default,
so on a fresh install Step 2 could not run at all. Reproduced live against
PubMed, fixed, and re-run: the same query now retrieves records. The failure
surfaced inside the search, which is why the query and then the date range
looked like the cause; neither was.

**D1** was real too: the delay was a flat `0.11` s (~9/second) whether or not a
key was set, against NCBI's unregistered limit of 3/second — which it enforces
by blocking the IP, not throttling it. The delay now follows the key.

---

## E. The repo itself

| | Item | Status |
| :-- | :--- | :--- |
| E1 | `__pycache__` untracked and ignored | **done** |
| E2 | Stray probe directories removed | **done** |
| E3 | README no longer claims deleting the folder uninstalls it | **done** — `b87f231` |
| E4 | Uninstall documented for all three systems | **done** |
| E5 | No test suite under version control | **should** |
| E6 | `CHANGELOG.md` | **should** |
| E7 | `CONTRIBUTING.md`, issue templates | **later** |
| E8 | Screenshot in the README | **should** |

**E5.** The ~146 checks across 12 suites that have been run during development
live in a scratch directory, not in the repo. Nothing in CI runs a single
assertion about pipeline behaviour — the workflow proves the package installs
and imports, and stops there. A public project that anyone might send a patch
to needs the tests in the tree, if only so a change that breaks the ledger
fails visibly. This does not block a release; it blocks the second one.

**E8.** People decide whether to download a 90 MB scientific tool from a
screenshot. There isn't one.

---

## F. Version bookkeeping — do this at the moment of release

The version is declared in **five** places and nothing keeps them in step:

| File | Line |
| :--- | :--- |
| `pyproject.toml` | `version = "3.2.10"` |
| `src/mesh_aop/__init__.py` | `__version__ = "3.2.10"` |
| `packaging/build_portable_windows.py` | `VERSION = '3.2.10'` |
| `packaging/build_msi_windows.py` | the output filename |
| `packaging/windows_msi.wxs` | `<?define AppVersion = "3.2.10" ?>` |

A mismatch here produces an installer whose Add/Remove Programs entry disagrees
with the application's own About box, and on Windows a wrong `AppVersion`
breaks in-place upgrades — the new MSI installs alongside the old one instead
of replacing it. Grep before tagging:

```bash
grep -rn "3\.2\.0" pyproject.toml src/mesh_aop/__init__.py packaging/
```

Worth doing once: have the builders read the version out of `pyproject.toml`
the way `build_unix_bundle.py` already does, and delete four of the five.

---

## G. Publishing — the last step, and a deliberate one

1. Everything above in section A–D that says **blocking** is resolved.
2. `foxglove` merged to `main`, or the release cut from `foxglove` knowingly.
3. Version bumped in all five places, committed.
4. Tag pushed. **A tag builds and tests. It does not publish** — that was
   changed after tagging released 3.2.0 without anyone deciding to.
5. Download the artifacts from the run and install each one on a real machine.
   This is the last point at which a bad build costs nothing.
6. **Ask before publishing.** Actions → Release → Run workflow → tick
   `publish`. This creates the public GitHub Release. It is not to be run
   without explicit approval, per standing instruction.

Release notes should lead with the download table — one file per system — and
carry the macOS quarantine instructions, because that is the only step a user
cannot guess.
