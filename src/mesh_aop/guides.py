# -*- coding: utf-8 -*-
"""
guides.py - the written instructions the program leaves next to its files.

Two moments in this pipeline leave someone holding a file and no context: the
run that pauses for AOP annotation, and the master database that took a night to
build. Both are places where the explanation is given once, on a console that is
closed by the time it matters - the annotation pause can sit for days, and nobody
reads a build log a month later when the database will not open.

So the instructions are written to disk beside the file they describe, in plain
text, named so they are the obvious thing to open. This module holds that text
and nothing else: no imports from the rest of the package, so it stays cheap and
can be read from the GUI, the CLI, or a bare interpreter.
"""

import os

# The biological strata, in pathway order. Duplicated in viz.AOP_ORDER for
# plotting; this copy exists so instructions can be produced without importing
# matplotlib, which the annotation pause has no other reason to load.
AOP_LEVELS = ('Stressor', 'Molecular', 'Cellular', 'Tissue', 'Organ',
              'Adverse Outcome', 'Uncategorized')

ANNOTATION_GUIDE_NAME = 'HOW TO ANNOTATE - read me.txt'
MASTER_DB_GUIDE_NAME = 'HOW TO REBUILD THIS DATABASE - read me.txt'


def _wrap(path):
    """The path, and a hint when it is long enough to be worth one."""
    return path


def annotation_guide_text(annotation_path, prefix='', results_dir=''):
    """What to do with a run-annotations file, for the file's own folder."""
    levels = '\n'.join(f'      {lvl}' for lvl in AOP_LEVELS)
    return f"""\
================================================================================
THE RUN IS PAUSED, WAITING FOR YOU
MeSH Workbench - AOP level annotation
================================================================================

Nothing has gone wrong. The pipeline built and scored the network, then stopped
on purpose, because the next step needs a judgement it cannot make.

WHAT IT NEEDS
    Every MeSH term in the final network has to be placed on the adverse outcome
    pathway - is it the stressor, an early molecular event, a cellular one, and
    so on through to the adverse outcome itself. That mapping is the biology,
    and it is yours to supply.

THE FILE TO EDIT
    {annotation_path}

    It is a semicolon-delimited CSV. Open it in Excel, LibreOffice, or a text
    editor. There is one row per MeSH term. Fill in the 'aop_level' column.

    Keep the semicolons. MeSH headings contain commas as a matter of course
    ("Dermatitis, Allergic Contact"), which is exactly why this file does not
    use them as separators. If your spreadsheet offers to save as "CSV
    (comma delimited)", say no and keep the current format.

THE ALLOWED LEVELS, IN PATHWAY ORDER
{levels}

    Spelling and capitalisation have to match. Anything unrecognised is read as
    Uncategorized.

    Use Uncategorized deliberately, not as a shrug: a term marked that way stays
    in the network and in every topological figure, and is simply left out of the
    AOP-specific ones. It is the right answer for a term that genuinely does not
    sit on the pathway.

WHEN YOU ARE DONE
    Save the file, then run the figures step:

        In the Workbench:   Run  ->  Step 4 - Figures
        In a terminal:      python -m mesh_aop.cli --step viz

    Nothing before that step is recomputed, so this is quick. You can edit the
    file and re-run the figures as often as you like.

IF YOU WOULD RATHER NOT ANNOTATE AT ALL
    Turn off "Pause for AOP annotation" on the Folders tab of Settings. The run
    then completes unattended with every term left Uncategorized. The network,
    the scores, the benchmark and the topological figures are all still valid -
    only the AOP-specific figures lose their meaning, because there is no
    pathway left in them to show.

================================================================================
Written automatically when the run paused. Safe to delete.
================================================================================
"""


def master_db_guide_text(db_path, data_dir='', workspace_dir=''):
    """How the master annotation database is built, kept beside the database."""
    where = os.path.dirname(db_path) or data_dir
    return f"""\
================================================================================
THE MASTER ANNOTATION DATABASE
MeSH Workbench - what this file is, and how to rebuild it
================================================================================

    {db_path}

WHAT IT IS
    A local copy of PubMed's MeSH annotations: for every article, its publication
    date and its subject headings. Everything else in this program reads it. It
    is built once and then reused by every run.

WHY YOU SHOULD NOT DELETE IT CASUALLY
    Building it means downloading roughly 44 GB from the NLM and compiling it
    into a database of about 8 GB. That takes several hours and, on a slow
    connection, most of a day. It is by far the longest thing this program does.

    The 44 GB of downloaded archives, on the other hand, ARE safe to delete once
    the database exists. The Data Setup screen has a button for exactly that.

HOW TO REBUILD IT
    In the Workbench:
        Database -> Data setup, then press Build next to "Master annotation
        database". If a damaged database is already there, tick the option to
        delete it first.

    In a terminal:
        python -m mesh_aop.cli --step baseline --build-database

        Add --skip-baseline-download to compile from archives already on disk,
        --with-updates to include the daily files published since the baseline,
        and --rebuild-corrupt to delete an unreadable database before starting.

IF IT WILL NOT OPEN
    A database interrupted mid-build - a power cut, a full disk, a machine that
    slept - can be left unreadable. The program checks it before every run and
    will tell you plainly if it is damaged rather than failing halfway through
    a long step.

    Run Tools -> Check and repair files in the Workbench. It reports what is
    damaged, offers to remove it, and tells you which step to resume from. Only
    the master database has to be rebuilt from scratch; everything else is
    derived and can be regenerated in minutes.

    A full disk is the usual cause. The build needs room for the archives, the
    database, and a working copy of the database at the same time. If this drive
    is short of space, set a different workspace under Settings -> Folders ->
    "Database build workspace"{f' (currently {workspace_dir})' if workspace_dir else ''}.

WHAT ELSE LIVES HERE
    {where}

    Downloaded PubMed archives, the per-project PMID databases, and this file.
    Everything here can be rebuilt; nothing here is your own work. Your results
    are kept somewhere else entirely, so clearing this folder never destroys an
    analysis - only the time it would take to download it all again.

================================================================================
Written automatically by MeSH Workbench. Safe to delete; it will come back.
================================================================================
"""


def _write(path, text):
    """Write a guide, returning its path, or None if it could not be written.

    A guide is a courtesy. Failing to write one must never interrupt whatever
    the program was actually doing.
    """
    try:
        os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(text)
        return path
    except OSError:
        return None


def write_annotation_guide(annotation_path, config=None):
    """Drop the annotation instructions beside the file to be annotated."""
    folder = os.path.dirname(os.path.abspath(annotation_path))
    return _write(os.path.join(folder, ANNOTATION_GUIDE_NAME),
                  annotation_guide_text(os.path.abspath(annotation_path)))


def write_master_db_guide(db_path, workspace_dir=''):
    """Drop the rebuild instructions beside the master database.

    Written next to where the database belongs, whether or not it exists yet -
    someone looking at an empty data folder and wondering what goes in it is
    exactly the person this is for.
    """
    folder = os.path.dirname(os.path.abspath(db_path))
    return _write(os.path.join(folder, MASTER_DB_GUIDE_NAME),
                  master_db_guide_text(os.path.abspath(db_path),
                                       workspace_dir=workspace_dir))
