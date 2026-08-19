# Packaging

Build tooling for distributing **MeSH Workbench**. Nothing here is imported by
the application or the pipeline — these scripts only assemble releases.

| File | Purpose |
| --- | --- |
| `build_portable_windows.py` | **The supported route.** Produces the download-and-run Windows zip. |
| `build_frozen.py` | PyInstaller bundle. Kept for reference; see the caveat below. |
| `windows_installer.iss` | Inno Setup script, pairs with `build_frozen.py`. |

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

## The frozen build

`build_frozen.py` produces a PyInstaller bundle and `windows_installer.iss`
wraps it in a conventional Setup wizard. Both work, but the resulting executable
is unsigned and was refused on the managed Windows machine this was developed on,
with no error surfaced — a `--windowed` build has nowhere to write a traceback.
Build with `--console` when diagnosing that.

Use this route only if a signing certificate is available. Otherwise prefer the
portable zip.

## macOS and Linux

Not built here. The embeddable-Python approach is Windows-specific and has no
equivalent on other platforms, which do not need one: both ship Python, so

```bash
pip install mesh_aop_network
mesh-workbench
```

is the expected route. A `.command` (macOS) or `.desktop` (Linux) launcher can
wrap that if a double-clickable file is wanted.
