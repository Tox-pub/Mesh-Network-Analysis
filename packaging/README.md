# Packaging

Build tooling for distributing **MeSH Workbench**. Nothing here is imported by
the application or the pipeline — these scripts only assemble releases.

| File | Purpose |
| --- | --- |
| `build_portable_windows.py` | Produces the download-and-run Windows zip. Run this first. |
| `build_installer_windows.py` | Wraps that zip's tree in a conventional `Setup.exe`. |
| `windows_installer.iss` | The Inno Setup script it compiles. |
| `launchers/` | The `.bat` files a user double-clicks, copied verbatim into the zip. |
| `make_icon.py` | Redraws the application icon into `src/mesh_workbench/assets/`. |

## Building the Windows release

```bash
python packaging/build_portable_windows.py
```

Produces `MeshWorkbench-<version>-win64-portable.zip`, roughly 150 MB. The user
extracts it and double-clicks **`MeSH Workbench.bat`**.

Output goes to `D:\mesh_workbench_build` by default, deliberately **outside the
project**: the tree is ~460 MB, is rebuilt from scratch every run, and the
working copy is cloud-synced, so keeping it here would upload half a gigabyte of
reproducible output on every build. `.gitignore` stops git, not the sync client.

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
python packaging/build_portable_windows.py     # first - assembles the tree
python packaging/build_installer_windows.py    # then  - wraps it
```

Produces `packaging/dist/MeSH-Workbench-<version>-win64-setup.exe`: a wizard with
a directory page, Start-menu and desktop shortcuts, and an Add/Remove Programs
entry.

It installs **the same tree the zip contains**, so the only program that ever
executes is still the PSF-signed `python.exe`. Shortcuts point at
`pythonw.exe app\launch.py` rather than at a `.bat`, so nothing flashes a
console. The `.bat` files and the portable README are excluded from the install -
they would be telling an installed copy to use launchers it does not have.

Defaults to a per-user install, so there is no administrator prompt and nothing
for a managed machine to refuse; the wizard still offers all-users to anyone with
the rights.

Uninstalling through Windows removes what Setup wrote, then offers to run the
application's own uninstaller for the downloaded data and the temp-folder
workspace - which Windows would otherwise leave behind, at tens of gigabytes.
Results are kept.

**Setup.exe is unsigned.** SmartScreen warns once per machine and the user clicks
*More info -> Run anyway*. Only a code-signing certificate removes that, and it
applies to the installer alone, never to the installed application. This is why
the portable zip stays the fallback for locked-down machines.

## Publishing a release

Three artefacts, all built from the same tree:

| Asset | Size | For |
| --- | --- | --- |
| `MeSH-Workbench-<version>-win64.msi` | ~94 MB | Managed machines. `msiexec` performs the install, and that is signed by Microsoft. |
| `MeSH-Workbench-<version>-win64-setup.exe` | ~68 MB | Ordinary machines. Conventional wizard. |
| `MeshWorkbench-<version>-win64-portable.zip` | ~117 MB | No install at all, or where both installers are refused. Contains `Install.bat`. |

### Repository or release?

**The repository holds source only.** Build scripts, launchers, the WiX and Inno
sources, the icon generator — everything needed to *produce* an artefact, and
none of the artefacts themselves. `packaging/portable/` and `packaging/dist/`
are gitignored precisely so a build cannot land in a commit.

**The release holds the artefacts.** GitHub rejects any file over 100 MB
committed to a repository; release assets are capped at 2 GB each. The zip
exceeds the repository limit, and even the two that would fit belong in a
release, because a binary committed to git history stays there forever and
every clone pays for it.

Releases are free, including the bandwidth to serve them, and they do not count
against Git LFS quotas.

### Steps

1. Confirm the version agrees everywhere. Ten files carry it, and a mismatch
   means an asset whose name disagrees with the tag it hangs under.

2. Build all three, onto an internal disk. Building from removable media took
   roughly three times as long here, and the MSI failed mid-read with "the
   volume for a file has been externally altered".

   ```
   python packaging/build_portable_windows.py   --out <build dir>
   python packaging/build_installer_windows.py  --portable <build dir>\MeshWorkbench --out <build dir>
   python packaging/build_msi_windows.py        --portable <build dir>\MeshWorkbench --out <build dir>
   ```

3. Hash all three, so a download can be verified:

   ```
   certutil -hashfile "<build dir>\MeshWorkbench-3.1.0-win64-portable.zip" SHA256
   ```

4. Tag the exact commit the artefacts were built from, and push the tag:

   ```
   git tag -a v3.1.0 -m "MeSH Workbench 3.1.0"
   git push origin v3.1.0
   ```

5. On GitHub: **Releases** → **Draft a new release** → choose the tag → attach
   the three files → paste the checksums into the notes → publish.

   With the `gh` CLI it is one command instead:

   ```
   gh release create v3.1.0 --title "MeSH Workbench 3.1.0" --notes-file notes.md *.msi *.exe *.zip
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

Not built here. The embeddable-Python approach is Windows-specific and has no
equivalent on other platforms, which do not need one: both ship Python, so

```bash
pip install mesh_aop_network
mesh-workbench
```

is the expected route. A `.command` (macOS) or `.desktop` (Linux) launcher can
wrap that if a double-clickable file is wanted.
