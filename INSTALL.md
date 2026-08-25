# Installing MeSH Workbench

Three ways on Windows, one on macOS and Linux. All install for the current user
only — no administrator rights, and nothing written outside your own profile.

- [Windows](#windows)
- [macOS and Linux](#macos-and-linux)
- [If a download is blocked](#if-a-download-is-blocked)
- [Where things are kept](#where-things-are-kept)
- [Updating](#updating)
- [Removing it](#removing-it)

---

## Windows

Download from [Releases](https://github.com/Tox-pub/Mesh-Network-Analysis/releases).
Each download carries its own Python, so nothing needs installing first.

### The installer (`.msi`) — start here

Double-click it. The wizard asks where to install, and whether you want a Start
menu entry and a desktop shortcut. It appears in **Settings → Apps** afterwards
like any other program.

This is the one to try first on a **work or managed computer**. The install is
performed by `msiexec`, which is part of Windows and signed by Microsoft, so it
does not introduce any new program for security policy to object to.

For unattended deployment:

```
msiexec /i "MeSH-Workbench-3.1.0-win64.msi" /qn
```

### The Setup wizard (`.exe`)

A conventional installer, and fine on a personal machine. It is not signed, so
Windows shows *"Windows protected your PC"* the first time — click **More info**
then **Run anyway**.

On some managed machines this one is refused outright rather than merely warned
about; see [below](#if-a-download-is-blocked).

### The portable zip

1. Right-click the zip → **Extract All** → choose any folder you can write to.
2. Open the extracted folder.
3. Double-click **`MeSH Workbench.bat`**.

Nothing is installed and nothing is written outside that folder — settings,
data and results all stay inside it. Move it, copy it to a USB stick, or delete
it; it leaves nothing behind.

To turn it into a proper installation instead, run **`Install.bat`** from inside
the folder. That copies it into your profile, adds shortcuts, and registers it
in Settings → Apps, using only tools Windows already has. Add `/S` to install
silently.

If the window does not appear, run **`MeSH Workbench (Troubleshooting).bat`** —
it keeps a console open so any error stays on screen.

---

## macOS and Linux

There is no packaged download; install from source.

**Requirements:** Python 3.11 to 3.13, including tkinter. On most Linux
distributions tkinter is a separate package and is missing by default:

```
sudo apt install python3.12 python3.12-venv python3-tk     # Debian, Ubuntu
sudo dnf install python3 python3-tkinter                   # Fedora
brew install python@3.12 python-tk                         # macOS
```

Then:

```
git clone https://github.com/Tox-pub/Mesh-Network-Analysis.git
cd Mesh-Network-Analysis
./packaging/install.sh
```

It checks for a supported Python and for tkinter before doing anything, creates
a private virtual environment, installs the package into it, and adds a menu
entry (Linux) or an application bundle in `~/Applications` (macOS).

`./packaging/install.sh --prefix DIR` puts the environment somewhere else.

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

If that happens, use the **`.msi`**, or the portable zip with `Install.bat`.
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
