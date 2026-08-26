# MeSH Workbench

Builds and validates **MeSH co-occurrence concept networks** from the PubMed
literature — connecting chemical stressors to adverse outcomes through
biological intermediates, in the shape of an Adverse Outcome Pathway.

A desktop application for Windows, and a command-line pipeline for any platform.

---

## Get it

**One file per system**, from
[Releases](https://github.com/Tox-pub/Mesh-Network-Analysis/releases).

| System | Download | What to do |
| :--- | :--- | :--- |
| Windows | `MeSH-Workbench-<version>-windows.msi` | Double-click it. |
| Linux | `MeSH-Workbench-<version>-linux-x86_64.tar.gz` | Extract, then `./"MeSH Workbench"` |
| macOS | `MeSH-Workbench-<version>-macos-arm64.tar.gz` | Extract, then `./"MeSH Workbench"` |

**Each one carries its own Python.** Nothing needs installing first - no system
Python, no `python3-tk`, no administrator rights, nothing written outside your
own profile. Delete the folder to uninstall.

On Windows the only program that ever executes is `python.exe`, signed by the
Python Software Foundation, and the installer is run by `msiexec.exe`, which is
part of Windows.

On macOS the download is unsigned - signing needs a paid Apple Developer ID -
so macOS quarantines it. The launcher clears that itself on first run, with no
password and nothing to pay; start it from Terminal that first time. See
[INSTALL.md](INSTALL.md#macos).

Full instructions, including what to do if a download is blocked:
**[INSTALL.md](INSTALL.md)**

---

## First run

The application opens on **Data Setup**. The pipeline works from a local copy of
PubMed's annotation data, so nothing can be analysed until that exists — press
**Build** and it downloads and compiles it for you.

Budget roughly **44 GB downloaded once** and about **8 GB** kept afterwards. The
download resumes if interrupted, and the archive can be deleted once the
database is built. You choose which drive holds it.

Once the database exists, everyday analysis runs offline.

---

## What it does

1. **Process** the MeSH vocabulary and build the stop-word set.
2. **Retrieve** an article cohort from a PubMed query, and expand it by citation.
3. **Build** the co-occurrence network, filter it to a consensus subgraph
   (GLF and simulated annealing), and score every term by mean relevancy.
4. **Validate** against a curated ground truth, with permutation nulls and
   bootstrap confidence intervals.
5. **Draw** the result as Adverse Outcome Pathway flows.

Each stage can be run on its own, from the application or the command line:

```
mesh-pipeline --step network
mesh-workbench
mesh-uninstall --list
```

A curated reference corpus ships with the program — the OECD AOP-40 allergic
contact dermatitis set — so the figures can be reproduced before retrieving
anything of your own.

---

## Documentation

| | |
| --- | --- |
| **[INSTALL.md](INSTALL.md)** | Installing, updating and removing, on all three platforms. |
| **[HELP.md](HELP.md)** | How the pipeline works, every setting, the outputs, and troubleshooting. |
| **[packaging/README.md](packaging/README.md)** | Building the releases. For maintainers. |

---

## Citing this work

The bundled ground truth is derived from the OECD Adverse Outcome Pathway
programme's case study on skin sensitisation (AOP 40). See
[HELP.md](HELP.md#citation) for the full citation and how to cite this software.

## Licence

MIT — see [LICENSE](LICENSE).
