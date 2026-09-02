# Installing MeSH Workbench

**One file per system.** Download the one with your system's name on it from
[Releases](https://github.com/Tox-pub/Mesh-Network-Analysis/releases), and
ignore the other two.

| System | File | What to do |
| :--- | :--- | :--- |
| **Windows** | `MeSH-Workbench-3.2.0-windows.msi` | Double-click it. |
| **Linux** | `MeSH-Workbench-3.2.0-linux-x86_64.tar.gz` | Extract, then `./"MeSH Workbench"` — it adds itself to your applications menu |
| **macOS** | `MeSH-Workbench-3.2.0-macos-arm64.tar.gz` | Extract, then `./"MeSH Workbench"` |

**Every one carries its own Python.** There is nothing to install first - no
system Python, no `python3-tk`, no administrator rights, and nothing written
outside your own profile.

- [Windows](#windows)
- [Linux](#linux)
- [macOS](#macos)
- [If a download is blocked](#if-a-download-is-blocked)
- [Where things are kept](#where-things-are-kept)
- [Updating](#updating)
- [Removing it](#removing-it)

---

## Windows

Double-click `MeSH-Workbench-3.2.0-windows.msi`. It asks where to install,
offers a desktop and Start-menu shortcut, and appears in Add/Remove Programs
afterwards.

It carries its own Python, so nothing needs installing first, and it is
executed by `msiexec.exe`, which is part of Windows and signed by Microsoft.
That matters on a managed machine: a freshly built `.exe` installer has no
signature and no prevalence, and Defender Exploit Guard refuses it outright
with `0x80070005` and nothing to click. An MSI introduces no new binary.

Silent install, for deploying to several machines:

```
msiexec /i "MeSH-Workbench-3.2.0-windows.msi" /qn
```

---

## Linux

Download `MeSH-Workbench-3.2.0-linux-x86_64.tar.gz`, then:

```
tar -xzf MeSH-Workbench-3.2.0-linux-x86_64.tar.gz
cd MeSH-Workbench-3.2.0-linux-x86_64
./"MeSH Workbench"
```

That is all. The folder carries its own Python, its own Tk, and every library
already compiled - **you do not need a system Python, and you do not need
`python3-tk`**. The first run unpacks the libraries from `wheels/`, takes about
a minute, and needs no network.

**It appears in your applications menu after that first launch.** A `.desktop`
entry has to name an absolute path, and that is not known until you unpack the
folder somewhere, so it cannot be shipped ready-made - the launcher writes it
the first time it runs, to `~/.local/share/applications/mesh-workbench.desktop`.
No root, nothing outside your own home directory, and it is rewritten if you
later move the folder. `./mesh-uninstall` removes it again. If your desktop does
not pick it up immediately, log out and back in.

Without a desktop:

```
./mesh-pipeline --step all
```

---

## macOS

**Yes, this runs on a Mac without paying Apple anything.**

Download `MeSH-Workbench-3.2.0-macos-arm64.tar.gz` (Apple silicon), extract it,
then **open Terminal**, change to the extracted folder and run:

```
./"MeSH Workbench"
```

### Why Terminal, and not a double-click

macOS tags everything downloaded through a browser with a `com.apple.quarantine`
attribute, and Gatekeeper refuses to run unsigned programs carrying it. This
program is unsigned, because signing requires a paid Apple Developer ID. So the
tag has to be cleared once.

Clearing it needs **no password, no Apple account and no payment**. The launcher
does it for you: it is a shell script, which macOS reads with its own signed
`/bin/sh`, so it is not itself blocked and can clear the flag from everything
else in the folder. You will see it say so on the first run.

Start it from Terminal rather than Finder for that first run - Finder will
either open the script in a text editor or refuse it. Afterwards, either works.

### Doing it by hand instead

If you would rather not have a script change attributes on your behalf:

```
xattr -dr com.apple.quarantine "/path/to/MeSH-Workbench-3.2.0-macos-arm64"
```

That removes one extended attribute from the files in that folder and touches
nothing else on your Mac. `xattr` is part of macOS.

### If macOS still refuses

On macOS Sequoia and later, the override moved: **System Settings → Privacy &
Security**, scroll to the bottom, and press **Open Anyway** next to the blocked
item. That button appears only after the first refusal.

### Intel Macs

The published build is for Apple silicon. An Intel build is one command away -
`python packaging/build_unix_bundle.py --target macos-intel` - and can be added
to a release on request.

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
msiexec /x "MeSH-Workbench-3.2.0-win64.msi" REMOVEDATA=1
```

**Windows portable zip:** run `Uninstall.bat`, then delete the folder.

**Linux and macOS (the `.tar.gz` bundle)** — two steps, and the second one is
just deleting a folder:

```
cd MeSH-Workbench-3.2.0-linux-x86_64
```

```
./mesh-uninstall
```

That removes what the program put **outside** its own folder: settings and
downloaded data under `~/.local/share/MeSH Workbench` (`~/Library/Application
Support` on macOS), the applications-menu entry, and — if you ask it to — your
results. It lists everything with sizes and asks before removing anything.

Then delete the folder itself, which is the program:

```
rm -rf MeSH-Workbench-3.2.0-linux-x86_64
```

**There is nothing to `pip uninstall`.** The bundle carries its own Python and
puts the application on that interpreter's path; it never installs it. Earlier
versions ended the uninstall by suggesting `python3.12 -m pip uninstall
mesh_aop_network`, which failed with *no such file or directory* because
`python3.12` exists only inside the bundle — ignore that instruction if you meet
it, and delete the folder instead.

**Installed from source with pip:** `python -m pip uninstall mesh_aop_network`,
after running `mesh-uninstall` to clear the data. Use the same interpreter you
installed it into.

Your results are never removed unless you explicitly ask.

To see exactly what is on disk before deciding, without changing anything:

```
mesh-uninstall --list
```

That lists every file the program downloaded, built or installed, with sizes —
including the pieces that live outside the program folder and would otherwise
be left behind.
