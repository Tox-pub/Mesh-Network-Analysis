# -*- coding: utf-8 -*-
"""
strata.py - the scheme a project groups its MeSH terms by.

The pipeline used to hard-code the seven levels of an adverse outcome pathway.
That is one scheme among many. The same network is worth stratifying by organ
system, by exposure route, by study design, by anything a question needs, and a
term's stratum is a judgement the program cannot make either way. So the strata
are now whatever the annotation file contains; the AOP seven are the default,
an example rather than a rule.

Two things still have to be decided in one place, and they are decided here.

WHAT THE COLUMN IS CALLED
    Files are written with 'stratum'. Files written before this change say
    'aop_level' and are still read, so a project annotated last month opens
    without being re-annotated.

WHAT ORDER THE STRATA GO IN
    An alluvial flow needs one, and no order can be inferred from names: only
    the person who chose them knows that Molecular comes before Cellular.
    The order is a setting. Anything found in the annotation file but missing
    from that setting is appended rather than dropped, because a stratum
    someone typed is a stratum they meant - the caller is told, so it can say
    so, but nothing is silently discarded.

Standard library only. The annotation pause writes its instructions from here
without loading pandas or matplotlib.
"""

__all__ = ['COLUMN', 'LEGACY_COLUMN', 'UNASSIGNED', 'PLACEHOLDER',
           'DEFAULT_ORDER', 'DEFAULT_ORDER_TEXT', 'parse_order',
           'resolve_order', 'column_of', 'normalise', 'is_unplaced']

# Written to every new file. The old name is still read - see column_of().
COLUMN = 'stratum'
LEGACY_COLUMN = 'aop_level'

# Two words for a term with no stratum, and they do NOT mean the same thing to
# the person doing the annotating:
#
#   PLACEHOLDER  nobody has looked at this term yet - every row starts here
#   UNASSIGNED   somebody looked, and it belongs to none of their groups
#
# The distinction is what lets a long file be worked through over several
# sittings: it separates "not reached" from "finished with". It is preserved in
# the file, and the count of each is reported.
#
# The FIGURES cannot use it. A stratum is a group to draw, and neither of these
# is one, so both are drawn the same way - kept in the network and in every
# topological figure, left out of the ones that show the scheme. normalise() is
# where that collapse happens, and it happens at plotting time only.
#
# They were already collapsed before this was written down, but by accident:
# the template wrote PLACEHOLDER, the figures' categorical only listed
# UNASSIGNED, so every untouched row fell out as NaN and was filled back in.
# Right answer, wrong reason, and nothing said so.
UNASSIGNED = 'Uncategorized'
PLACEHOLDER = 'Unassigned'

# The adverse outcome pathway, which is where this program started and is still
# what most users want. Shipped as the default so an AOP project needs no
# setting at all; replaced wholesale by anyone using a different scheme.
DEFAULT_ORDER = ('Stressor', 'Molecular', 'Cellular', 'Tissue', 'Organ',
                 'Adverse Outcome', UNASSIGNED)

# Semicolons, to match the annotation file. MeSH headings contain commas
# ("Dermatitis, Allergic Contact") and a stratum name may well too.
DEFAULT_ORDER_TEXT = '; '.join(DEFAULT_ORDER)


def parse_order(text):
    """A settings string into a list of strata, in the order given.

    Semicolon-delimited, tolerant of a comma-delimited list because that is
    what someone types when they have not read the label. Blanks and repeats
    are dropped: a duplicate would make the flow diagram ambiguous.
    """
    if not text:
        return []
    raw = text.split(';') if ';' in text else text.split(',')
    out = []
    for name in (s.strip() for s in raw):
        if name and name not in out:
            out.append(name)
    return out


def resolve_order(configured, present):
    """The strata to plot, in order, and the ones the setting did not mention.

    `configured` is the setting, `present` is what the annotation file actually
    contains. The result keeps the configured order, drops configured strata
    nobody used, and appends the rest so that mislabelling a stratum shows up
    as a stratum at the end of the figure instead of a term that vanished.

    Returns (order, appended). `appended` is what the caller should report.
    """
    present = list(dict.fromkeys(present))          # de-duplicate, keep order
    wanted = parse_order(configured) if isinstance(configured, str) \
        else list(configured or [])

    order = [s for s in wanted if s in present]
    appended = [s for s in present if s not in wanted]
    order.extend(appended)

    # Unassigned always sits last: it is the absence of a stratum, not one of
    # them, and in a flow diagram it belongs nowhere else.
    if UNASSIGNED in order:
        order = [s for s in order if s != UNASSIGNED] + [UNASSIGNED]
        appended = [s for s in appended if s != UNASSIGNED]
    return order, appended


def is_unplaced(value):
    """True when this cell means "nobody has placed this term".

    Blank, missing, the template's placeholder, or the word the figures use.
    Compared case-insensitively: a user who types 'unassigned' in a spreadsheet
    means the same thing as one who leaves the capital in place.
    """
    if value is None:
        return True
    text = str(value).strip()
    if not text or text.lower() == 'nan':
        return True
    return text.lower() in (PLACEHOLDER.lower(), UNASSIGNED.lower())


def normalise(value):
    """One cell of the column, as the rest of the program should see it.

    Everything that means "unplaced" becomes UNASSIGNED. Everything else is
    passed through with its surrounding space removed and nothing else changed:
    the user's spelling of their own stratum is not ours to correct.
    """
    return UNASSIGNED if is_unplaced(value) else str(value).strip()


def column_of(columns):
    """Which of the two spellings this file uses, or None.

    The new name wins when a file somehow carries both, which happens if a
    user merges an old annotation file into a new one by hand.
    """
    names = {str(c).lower(): c for c in columns}
    return names.get(COLUMN) or names.get(LEGACY_COLUMN)
