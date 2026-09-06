# Packaging

Build tooling for distributing **MeSH Workbench**. Nothing here is imported by
the application or the pipeline — these scripts only assemble releases.

| File | Purpose |
| --- | --- |
| `build_portable_windows.py` | Produces the download-and-run Windows zip. Run this first. |
| `build_msi_windows.py` | Wraps that zip's tree in an `.msi`. |
| `windows_msi.wxs` | The WiX source it compiles. |
| `build_unix_bundle.py` | The self-contained Linux and macOS tarballs. |
| `verify_windows_bundle.py`, `verify_unix_bundle.py` | Check a build before it is published. |
| `launchers/` | The `.bat` files a user double-clicks, copied verbatim into the zip. |
| `make_icon.py` | Redraws the application icon into `src/mesh_workbench/assets/`. |

### There is no Setup.exe

There was, built with Inno Setup, and it was dropped. It installed the same tree
as the MSI, from the same portable build, to the same place — the only real
difference being that a freshly compiled `Setup.exe` is an unsigned binary and
the MSI is executed by `msiexec`, which Microsoft signs. On a managed machine
that difference decides whether the install runs at all, and it never favoured
Setup.exe. Two installers that do the same thing is one more artefact to build,
verify, hash, upload and support, for no case the other does not cover.

Anyone who installed an older version with it removes it the same way as the
MSI: **Settings → Apps → MeSH Workbench → Uninstall**.

## Building the Windows release

```bash
python packaging/build_portable_windows.py
```

Produces `MeshWorkbench-<version>-win64-portable.zip`, roughly 150 MB. The user
extracts it and double-clicks **`MeSH Workbench.bat`**.

### Everything is built in one place

`D:\mesh_workbench_build`. Every script here resolves it through
`build_location.py`, so they cannot disagree.

Deliberately **outside the project**: the tree is ~460 MB, is rebuilt from
scratch every run, and the working copy is cloud-synced, so keeping it here
would upload half a gigabyte of reproducible output on every build.
`.gitignore` stops git, not the sync client.

**If D: is not attached, the build stops.** It does not pick somewhere else.
Each script used to decide this for itself — the portable build fell back into
`packaging/portable`, the Unix bundles defaulted to `~/Documents`, and the MSI
wrote beside whatever tree it was handed. The result was four folders across
two drives holding several superseded copies of the same artefacts, and a build
run without D: quietly writing half a gigabyte into the synced project folder.

Override deliberately when you mean to: `--out <path>` for one run, or
`MESH_BUILD_OUT` for a session. CI sets the latter, because a runner has no D:.

Options: `--out` (build directory, or set `MESH_BUILD_OUT`), `--repo` (project
directory), `--venv` (environment to take dependencies from), `--base-python`
(full CPython install to take the Tk runtime from), `--no-zip`.

### What it assembles

```
MeshWorkbench/
  MeSH Workbench.bat                    launcher
  MeSH Workbench (Troubleshooting).bat  same, with a console attached
  README - Install and First Run.txt
  python/                               embeddable CPython + Tk
  app/                                  mesh_workbench, mesh_aop, reference data
```

## Building the installer

```bash
python packaging/build_portable_windows.py   # first - assembles the tree
python packaging/build_msi_windows.py        # then  - wraps it
```

Produces `MeSH-Workbench-<version>-win64.msi`: a directory page, Start-menu and
desktop shortcuts, and an Add/Remove Programs entry.

It installs **the same tree the zip contains**, so the only program that ever
executes is still the PSF-signed `python.exe`. Shortcuts point at
`pythonw.exe app\launch.py` rather than at a `.bat`, so nothing flashes a
console. The `.bat` files and the portable README are excluded from the install -
they would be telling an installed copy to use launchers it does not have.

Installs per-user, so there is no administrator prompt and nothing for a managed
machine to refuse.

Uninstalling through Windows removes what the installer wrote, then offers to
run the application's own uninstaller for the downloaded data and the
temp-folder workspace - which Windows would otherwise leave behind, at tens of
gigabytes. Results are kept.

The MSI takes a while on "Computing space requirements": it carries one
component per file, and there are nearly ten thousand. That is `CostFinalize`
doing its job, not a hang.

## Publishing a release

Four artefacts, all built from the same source:

| Asset | Size | For |
| --- | --- | --- |
| `MeSH-Workbench-<version>-win64.msi` | ~94 MB | Windows. `msiexec` performs the install, and that is signed by Microsoft. |
| `MeshWorkbench-<version>-win64-portable.zip` | ~117 MB | No install at all, or where the installer is refused. Contains `Install.bat`. |
| `MeshWorkbench-<version>-linux-x86_64.tar.gz` | ~400 MB | Linux, self-contained: its own CPython, Tk and wheels. |
| `MeshWorkbench-<version>-macos-<arch>.tar.gz` | ~400 MB | macOS, same shape. **Unvalidated** - assembled on Windows and never run on a Mac. |

### Repository or release?

**The repository holds source only.** Build scripts, launchers, the WiX source,
the icon generator — everything needed to *produce* an artefact, and none of the
artefacts themselves. `packaging/portable/` and `packaging/dist/` are gitignored
precisely so a build cannot land in a commit.

**The release holds the artefacts.** GitHub rejects any file over 100 MB
committed to a repository; release assets are capped at 2 GB each. Three of the
four exceed the repository limit, and the one that would fit belongs in a
release too, because a binary committed to git history stays there forever and
every clone pays for it.

Releases are free, including the bandwidth to serve them, and they do not count
against Git LFS quotas.

### Steps

1. Confirm the version agrees everywhere. Ten files carry it, and a mismatch
   means an asset whose name disagrees with the tag it hangs under.

2. Build them all, onto an internal disk. Building from removable media took
   roughly three times as long here, and the MSI failed mid-read with "the
   volume for a file has been externally altered".

   ```
   python packaging/build_portable_windows.py --out <build dir>
   python packaging/build_msi_windows.py      --portable <build dir>\MeshWorkbench --out <build dir>
   python packaging/build_unix_bundle.py      --platform linux --out <build dir>
   python packaging/build_unix_bundle.py      --platform macos --out <build dir>
   ```

   Then verify, before anything is hashed or uploaded:

   ```
   python packaging/verify_windows_bundle.py <build dir>\MeshWorkbench-3.2.0-win64-portable.zip
   python packaging/verify_unix_bundle.py    <build dir>\MeshWorkbench-3.2.0-linux-x86_64.tar.gz
   ```

3. Hash them all, so a download can be verified:

   ```
   certutil -hashfile "<build dir>\MeshWorkbench-3.2.0-win64-portable.zip" SHA256
   ```

4. Tag the exact commit the artefacts were built from, and push the tag:

   ```
   git tag -a v3.2.0 -m "MeSH Workbench 3.2.0"
   git push origin v3.2.0
   ```

5. On GitHub: **Releases** → **Draft a new release** → choose the tag → attach
   the four files → paste the checksums into the notes → publish.

   With the `gh` CLI it is one command instead:

   ```
   gh release create v3.2.0 --title "MeSH Workbench 3.2.0" --notes-file notes.md *.msi *.zip *.tar.gz
   ```

6. Point the project README's download link at the new release.

A tag is not required to draft a release — GitHub will create one — but tagging
first is worth the extra step, since it records exactly which commit produced
the binaries.

### Why this rather than a frozen executable

The only program a user runs is `python.exe` from the official embeddable
distribution, **signed by the Python Software Foundation**. That side-steps the
whole trust problem: no code-signing certificate to buy, no SmartScreen
reputation to earn, and none of the antivirus false positives PyInstaller
bundles routinely attract. On a managed Windows device a freshly built unsigned
binary can be refused outright, while a signed interpreter runs normally.

Dependencies are copied from a working virtual environment rather than installed
fresh, so a release ships the exact versions the published results were produced
with.

### Three things that will break it if changed carelessly

1. **`Lib` must stay on the path in `python3xx._pth`.** The embeddable build
   serves the standard library from a zip and puts only `.` on `sys.path`, so
   anything added under `Lib/` — which is how tkinter arrives — is invisible
   without it.
2. **Do not prune directories that look like test suites.** `numpy.testing` is
   public API that scipy imports while loading; removing it breaks every
   scipy-dependent import.
3. **Keep `.dist-info`.** Several packages, plotly among them, read their own
   version through `importlib.metadata` and raise without it.

The embeddable package also ships **without tkinter**, so `add_tkinter()` copies
`Lib/tkinter`, `_tkinter.pyd`, the Tcl/Tk DLLs and the `tcl/` runtime from a full
CPython install of the same version.


## macOS and Linux

```bash
python packaging/build_unix_bundle.py --platform linux
python packaging/build_unix_bundle.py --platform macos
```

Self-contained tarballs carrying their own CPython, Tk and every wheel, so the
first run installs offline (`--no-index`) and needs nothing from the machine.
The Linux bundle also drops a `.desktop` entry into the applications menu.

`pip install mesh_aop_network` remains the lighter route on either platform for
anyone who already has Python.

**The macOS bundle is unvalidated.** It is assembled on Windows and has never
been run on a Mac. `verify_unix_bundle.py` checks what can be checked from the
outside — the executable bit NTFS cannot store, wheels built for the wrong
Python or architecture, a launcher written with CRLF — but that is not the same
as running it, and the README should not imply otherwise until it has been.
