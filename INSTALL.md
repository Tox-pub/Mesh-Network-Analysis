# Installing MeSH Workbench

**One file per system.** Download the one with your system's name on it from
[Releases](https://github.com/Tox-pub/Mesh-Network-Analysis/releases), and
ignore the other two.

| System | File | What to do |
| :--- | :--- | :--- |
| **Windows** | `MeSH-Workbench-3.1.0-windows.msi` | Double-click it. |
| **Linux** | `MeSH-Workbench-3.1.0-linux.tar.gz` | Extract, then `./packaging/install.sh` |
| **macOS** | `MeSH-Workbench-3.1.0-macos.tar.gz` | Extract, then `./packaging/install.sh` |

Everything installs for the current user only - no administrator rights, and
nothing written outside your own profile.

- [Windows](#windows)
- [Linux and macOS](#linux-and-macos)
- [If a download is blocked](#if-a-download-is-blocked)
- [Where things are kept](#where-things-are-kept)
- [Updating](#updating)
- [Removing it](#removing-it)

---

## Windows

Double-click `MeSH-Workbench-3.1.0-windows.msi`. It asks where to install,
offers a desktop and Start-menu shortcut, and appears in Add/Remove Programs
afterwards.

It carries its own Python, so nothing needs installing first, and it is
executed by `msiexec.exe`, which is part of Windows and signed by Microsoft.
That matters on a managed machine: a freshly built `.exe` installer has no
signature and no prevalence, and Defender Exploit Guard refuses it outright
with `0x80070005` and nothing to click. An MSI introduces no new binary.

Silent install, for deploying to several machines:

```
msiexec /i "MeSH-Workbench-3.1.0-windows.msi" /qn
```

---

## Linux and macOS

**Requirements:** Python 3.11 to 3.13, including tkinter. On most Linux
distributions tkinter is a separate package and is missing by default - and its
absence is only discovered when the window fails to open, so install it now:

```
sudo apt install python3.12 python3.12-venv python3-tk     # Debian, Ubuntu
sudo dnf install python3 python3-tkinter                   # Fedora
brew install python@3.12 python-tk                         # macOS
```

Then:

```
tar -xzf MeSH-Workbench-3.1.0-linux.tar.gz      # or -macos
cd MeSH-Workbench-3.1.0
./packaging/install.sh
```

The script checks for a supported Python and for tkinter **before** changing
anything, creates a private virtual environment, installs the package into it,
and adds a menu entry (Linux) or an application bundle in `~/Applications`
(macOS). `--prefix DIR` puts the environment somewhere else.

Resolving the dependencies needs a network connection and takes a few minutes
the first time - scipy, scikit-learn, gensim, statsmodels and matplotlib are
all sizeable.

To run it without a desktop, the pipeline works headless:

```
~/.local/share/mesh-workbench/venv/bin/mesh-pipeline --step all
```

### Checking the install without downloading 44 GB

Turn on **Use bundled reference data** on the Folders tab, then run
**Step 4 - Figures**. The reference corpus ships with the program, so this
draws all eight figures and the PRISMA report from data already on disk. It is
the quickest honest end-to-end test.

---

## If a download is blocked

Two things can stop a fresh download from running, and they look similar but are
not the same.

**A warning you can dismiss.** *"Windows protected your PC"* means SmartScreen
has not seen this file before. Click **More info → Run anyway**. Expected for
any installer without a purchased code-signing certificate.

**A refusal you cannot dismiss.** *"Blocked an operation that is not allowed by
your IT administrator"*, or error `0x80070005`, means a security policy —
usually Defender Exploit Guard — has refused the file because it is newly
compiled and unsigned. There is nothing to click through.

If that happens the `.msi` is the route that works: it is executed by `msiexec.exe`, which is part of Windows and signed by Microsoft, so no unsigned binary is introduced.
Neither creates a new executable, so neither trips that rule.

**A file that arrived over the network** may also be flagged. Right-click it →
**Properties** → tick **Unblock** → OK. This is the most common reason a copied
download appears to do nothing at all.

---

## Where things are kept

An installed copy keeps each user's settings and results separate, the way any
other program does.

| | Installed | Portable |
| --- | --- | --- |
| Program | `%LOCALAPPDATA%\Programs\MeSH Workbench` | the folder you extracted |
| Settings | `%LOCALAPPDATA%\MeSH Workbench` | the same folder |
| Results | `Documents\MeSH Workbench` | the same folder |
| Downloaded data | your choice, set on the Folders tab | the same folder |

A portable copy is self-contained because it carries a file named
`portable.marker`. Delete that file and it behaves as an installed copy.

**The downloaded data is the large part** — roughly 44 GB for the PubMed
archive and 8 GB for the database built from it. Set the data folder on the
**Folders** tab before starting a download, and put it on a drive with room.

---

## Updating

Install the new version over the old one. The `.msi` replaces it in place, and
your settings and results are kept. For a portable copy, extract the new zip and
delete the old folder.

---

## Removing it

**Installed:** Settings → Apps → MeSH Workbench → Uninstall.

By default this keeps the PubMed data you downloaded, since it is large and
slow to fetch again, and always clears the temporary working copy the database
build leaves in your Windows temp folder. To remove the downloaded data too:

```
msiexec /x "MeSH-Workbench-3.1.0-win64.msi" REMOVEDATA=1
```

**Portable:** run `Uninstall.bat`, then delete the folder.

**macOS and Linux:** `./packaging/install.sh --uninstall`.

Your results are never removed unless you explicitly ask.

To see exactly what is on disk before deciding, without changing anything:

```
mesh-uninstall --list
```

That lists every file the program downloaded, built or installed, with sizes —
including the pieces that live outside the program folder and would otherwise
be left behind.
