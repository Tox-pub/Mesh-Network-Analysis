# -*- coding: utf-8 -*-
"""
prisma.py - a PRISMA-style flow report for a completed run.

PRISMA 2020 describes a flow of records through a systematic review:
identification, screening, eligibility, inclusion, with the excluded counts and
their reasons set beside each step. This pipeline has the same shape - a seeded
PubMed search expanded over citation generations, a MeSH annotation screen, two
optimisers pruning the co-occurrence graph, and a largest connected component
that everything downstream is reported on - so it is drawn the same way.

The diagram is a rendering of the run ledger, not a second set of calculations.
Whatever the ledger recorded is what appears here; a stage the run never reached
is drawn as pending rather than as zero, because "not run" and "found nothing"
are different claims and a reader cannot tell them apart from a blank box.

Two files come out of it: a figure (in the formats the run is configured for)
and a plain-text summary for the results folder, which is what usually gets
pasted into a methods section.
"""

import os
from datetime import datetime

from . import paths as _paths

# The five stages, in flow order, as (ledger stage, printed heading).
STAGES = (
    ('identification', 'Identification'),
    ('screening', 'Screening'),
    ('pruning', 'Pruning'),
    ('included', 'Included'),
    ('unipartite', 'Comparison (optional)'),
)


def _fmt(n):
    """Thousands separators, and a dash where the ledger has nothing."""
    if n is None:
        return '-'
    try:
        return f"{int(n):,}"
    except (TypeError, ValueError):
        return str(n)


def _generation_rows(ledger, kind):
    """(label, count) per generation, seed first, in citation-hop order."""
    rows = []
    for (stage, metric), row in ledger.rows.items():
        if stage != 'identification' or not metric.startswith(f'records_{kind}_'):
            continue
        gen = metric[len(f'records_{kind}_'):]
        if gen == 'total':
            continue
        try:
            rows.append((gen, int(row['value'])))
        except (TypeError, ValueError):
            continue
    # P0 is the seed search and sorts before the numbered hops.
    return sorted(rows, key=lambda r: (0, 0) if r[0].upper() == 'P0'
                  else (1, int(r[0][1:]) if r[0][1:].isdigit() else 99))


def build_flow(ledger):
    """The flow as data: a list of (heading, [kept lines], [excluded lines]).

    Kept separate from the drawing so the text report and the figure can never
    disagree about what the run did.
    """
    g = ledger.get_int
    flow = []

    # ---------------------------------------------------------- Identification
    kept, excluded = [], []
    for gen, n in _generation_rows(ledger, 'retrieved'):
        label = 'Seed search (P0)' if gen.upper() == 'P0' else f'Citation hop {gen}'
        kept.append(f"{label}: {_fmt(n)}")
    total = g('identification', 'records_retrieved_total')
    if total is not None:
        kept.append(f"Records identified in total: {_fmt(total)}")
    dropped = g('identification', 'records_removed_final_generation')
    if dropped:
        excluded.append(f"Outermost generation removed: {_fmt(dropped)}")
        excluded.append("  (its citation links are incomplete)")
    retained = g('identification', 'records_retained_total')
    if retained is not None:
        kept.append(f"Records carried forward: {_fmt(retained)}")
    flow.append(('Identification', kept, excluded))

    # ---------------------------------------------------------------- Screening
    kept, excluded = [], []
    screened = g('screening', 'articles_screened')
    if screened is not None:
        kept.append(f"Records screened for MeSH: {_fmt(screened)}")
    with_mesh = g('screening', 'articles_with_mesh')
    if with_mesh is not None:
        kept.append(f"Records with MeSH annotation: {_fmt(with_mesh)}")
    for metric, label in (('records_failed_mesh_retrieval', 'MeSH retrieval failed'),
                          ('records_no_mesh_available', 'No MeSH available'),
                          ('articles_without_mesh', 'Records without MeSH')):
        n = g('screening', metric)
        if n:
            excluded.append(f"{label}: {_fmt(n)}")
    total_anno = g('screening', 'annotations_total')
    if total_anno is not None:
        kept.append(f"MeSH annotations retrieved: {_fmt(total_anno)}")
    for metric, label in (('annotations_subheadings_removed', 'Subheading forms removed'),
                          ('terms_removed_as_stop_words', 'Terms removed as stop words')):
        n = g('screening', metric)
        if n:
            excluded.append(f"{label}: {_fmt(n)}")
    stops = g('screening', 'stop_words_in_list')
    if stops:
        excluded.append(f"  (stop-word list: {_fmt(stops)} terms)")
    after = g('screening', 'terms_after_stop_words')
    if after is not None:
        kept.append(f"Distinct MeSH terms retained: {_fmt(after)}")
    flow.append(('Screening', kept, excluded))

    # ------------------------------------------------------------------ Pruning
    kept, excluded = [], []
    fn, fe = g('pruning', 'full_network_nodes'), g('pruning', 'full_network_edges')
    if fn is not None:
        kept.append(f"Unfiltered network: {_fmt(fn)} terms, {_fmt(fe)} relations")
    target = g('pruning', 'target_edges')
    if target is not None:
        kept.append(f"Edge budget set for the optimisers: {_fmt(target)}")
    for key, label in (('glf', 'Graph Likelihood Filtering'), ('sa', 'Simulated Annealing')):
        n, e = g('pruning', f'{key}_nodes'), g('pruning', f'{key}_edges')
        if n is not None:
            kept.append(f"{label}: {_fmt(n)} terms, {_fmt(e)} relations")
    cn, ce = g('pruning', 'consensus_nodes'), g('pruning', 'consensus_edges')
    if cn is not None:
        kept.append(f"Consensus of both optimisers: {_fmt(cn)} terms, {_fmt(ce)} relations")
    for metric, label in (('pruning_nodes_excluded', 'Terms kept by neither optimiser'),
                          ('pruning_edges_excluded', 'Relations kept by neither optimiser'),
                          ('consensus_nodes_excluded', 'Terms kept by only one optimiser'),
                          ('consensus_edges_excluded', 'Relations kept by only one optimiser')):
        n = g('pruning', metric)
        if n:
            excluded.append(f"{label}: {_fmt(n)}")
    flow.append(('Pruning', kept, excluded))

    # ----------------------------------------------------------------- Included
    kept, excluded = [], []
    ln, le = g('included', 'lcc_nodes'), g('included', 'lcc_edges')
    if ln is not None:
        kept.append(f"Largest connected component: {_fmt(ln)} terms, {_fmt(le)} relations")
    for metric, label in (('louvain_communities', 'Communities detected'),
                          ('terms_scored_mrs', 'Terms scored (MRS)')):
        n = g('included', metric)
        if n is not None:
            kept.append(f"{label}: {_fmt(n)}")
    smaller = g('included', 'components_excluded')
    if smaller:
        excluded.append(f"Smaller networks excluded: {_fmt(smaller)}")
        biggest = g('included', 'largest_excluded_component')
        if biggest:
            excluded.append(f"  (largest of them: {_fmt(biggest)} terms)")
    for metric, label in (('lcc_nodes_excluded', 'Terms outside the LCC'),
                          ('lcc_edges_excluded', 'Relations outside the LCC')):
        n = g('included', metric)
        if n:
            excluded.append(f"{label}: {_fmt(n)}")
    flow.append(('Included', kept, excluded))

    # -------------------------------------------------------------- Comparison
    n = g('unipartite', 'comparison_networks')
    if n:
        kept = [f"Networks overlaid for comparison: {_fmt(n)}"]
        for i in range(1, n + 1):
            name = ledger.get('unipartite', f'comparison_network_{i}')
            if name:
                kept.append(f"  {name}")
        flow.append(('Comparison (optional)', kept, []))

    return [(head, k, e) for head, k, e in flow if k or e]


# ------------------------------------------------------------------- text form

def render_text(ledger, flow=None):
    """The flow as a plain-text report - what gets pasted into a methods section."""
    flow = flow if flow is not None else build_flow(ledger)
    dates = ledger.stage_dates()
    width = 78
    out = ['=' * width,
           'PRISMA-STYLE FLOW REPORT',
           'MeSH Workbench - MeSH co-occurrence network pipeline',
           '=' * width, '']

    term = ledger.get('run', 'search_term', '') or '(not recorded)'
    out.append(f"Search term      : {term}")
    start = ledger.get('run', 'search_start_date', '')
    end = ledger.get('run', 'search_end_date', '')
    if start or end:
        out.append(f"Search window    : {start or '?'} to {end or '?'}")
    gens = ledger.get('run', 'generations_requested', '')
    if gens:
        out.append(f"Citation hops    : {gens} generation(s) expanded from the seed set")
    ctx_s = ledger.get('run', 'context_start_date', '')
    ctx_e = ledger.get('run', 'context_end_date', '')
    if ctx_s or ctx_e:
        out.append(f"Relevance window : {ctx_s or '?'} to {ctx_e or '?'}")
    seed = ledger.get('run', 'random_seed', '')
    if seed:
        out.append(f"Random seed      : {seed}")
    ver = ledger.get('run', 'workbench_version', '')
    if ver:
        out.append(f"Workbench version: {ver}")
    out.append(f"Report generated : {datetime.now().strftime('%Y-%m-%d %H:%M')}")

    # Stages written on different days are legible only if the dates are shown.
    if len(set(dates.values())) > 1:
        out += ['', 'Stages were recorded at different times:']
        for stage, heading in STAGES:
            if stage in dates:
                out.append(f"  {heading:<22} {dates[stage]}")
    out.append('')

    for heading, kept, excluded in flow:
        out.append('-' * width)
        out.append(heading.upper())
        out.append('-' * width)
        for line in kept:
            out.append(f"  {line}")
        if excluded:
            out.append('')
            out.append('  Excluded at this step:')
            for line in excluded:
                out.append(f"    - {line}" if not line.startswith('  ') else f"    {line.strip()}")
        out.append('')

    out.append('=' * width)
    out.append('Counts are read from the run ledger written alongside this report.')
    out.append('A dash means the pipeline had not reached that stage when it was written.')
    out.append('=' * width)
    return '\n'.join(out)


# ----------------------------------------------------------------- figure form

def render_figure(ledger, output_dir, file_prefix, flow=None):
    """Draw the flow diagram. Returns the paths written, or [] if it could not."""
    flow = flow if flow is not None else build_flow(ledger)
    if not flow:
        return []

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    from .viz import _FIGURE_DPI, _FIGURE_FORMATS, _FORMAT_KWARGS

    # Height follows content: a run with many citation hops has a taller
    # Identification box, and fixed geometry would either clip it or leave the
    # short stages swimming in white space. One axis unit is one inch here, so
    # the per-line figures are the rendered line heights (font size x
    # linespacing, in points, over 72) with a little padding on top and bottom.
    kept_line, excl_line, pad = 9.0 * 1.55 / 72, 8.5 * 1.55 / 72, 0.45
    heights = [max(1.0,
                   len(k) * kept_line + pad,
                   (len(e) + 1) * excl_line + pad)  # +1 for the "Excluded" heading
               for _, k, e in flow]
    gap = 0.55
    fig_h = sum(heights) + gap * (len(flow) - 1) + 1.4
    fig = plt.figure(figsize=(13.5, fig_h))
    # The axes has to fill the canvas exactly, or a y unit is not an inch and
    # the height arithmetic above silently under-reserves - text then runs out
    # through the bottom of its box.
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, fig_h)
    ax.axis('off')

    stage_colour = '#1f4e79'
    kept_face = '#eaf1f8'
    excl_face = '#f7ece9'
    excl_edge = '#a4442f'

    ax.text(5.0, fig_h - 0.35, 'PRISMA-style flow of records through the pipeline',
            ha='center', va='center', fontsize=15, fontweight='bold', color=stage_colour)
    term = ledger.get('run', 'search_term', '')
    subtitle = f"Search: {term}" if term else 'Search term not recorded'
    gens = ledger.get('run', 'generations_requested', '')
    if gens:
        subtitle += f"   |   {gens} citation generation(s)"
    ax.text(5.0, fig_h - 0.75, subtitle, ha='center', va='center',
            fontsize=9.5, color='#444444')

    top = fig_h - 1.15
    centres = []
    for (heading, kept, excluded), h in zip(flow, heights):
        bottom = top - h
        mid = (top + bottom) / 2.0
        centres.append((mid, top, bottom))

        # Stage label, rotated down the left margin as in the published layout.
        # Rotated text is bounded by the box HEIGHT, so a short stage gets a
        # smaller label rather than one that runs past the bar.
        label_size = min(10.5, max(6.5, (h - 0.15) * 72 / (0.62 * len(heading))))
        ax.text(0.30, mid, heading.upper(), ha='center', va='center',
                fontsize=label_size, fontweight='bold', color='white', rotation=90)
        ax.add_patch(FancyBboxPatch((0.05, bottom), 0.5, h,
                                    boxstyle='round,pad=0.01,rounding_size=0.05',
                                    facecolor=stage_colour, edgecolor='none', zorder=0))

        ax.add_patch(FancyBboxPatch((0.75, bottom), 5.15, h,
                                    boxstyle='round,pad=0.01,rounding_size=0.06',
                                    facecolor=kept_face, edgecolor=stage_colour,
                                    linewidth=1.1))
        ax.text(0.95, top - 0.22, '\n'.join(kept) or '(no counts recorded)',
                ha='left', va='top', fontsize=9, color='#111111', linespacing=1.55)

        if excluded:
            ax.add_patch(FancyBboxPatch((6.35, bottom + 0.05), 3.55, h - 0.1,
                                        boxstyle='round,pad=0.01,rounding_size=0.06',
                                        facecolor=excl_face, edgecolor=excl_edge,
                                        linewidth=1.0, linestyle='--'))
            ax.text(6.55, top - 0.27, 'Excluded\n' + '\n'.join(excluded),
                    ha='left', va='top', fontsize=8.5, color='#5a2116', linespacing=1.55)
            ax.add_patch(FancyArrowPatch((5.90, mid), (6.35, mid),
                                         arrowstyle='-|>', mutation_scale=13,
                                         color=excl_edge, linewidth=1.1, linestyle='--'))
        top = bottom - gap

    for (mid, _, bottom), (next_mid, next_top, _) in zip(centres, centres[1:]):
        ax.add_patch(FancyArrowPatch((3.3, bottom), (3.3, next_top),
                                     arrowstyle='-|>', mutation_scale=15,
                                     color=stage_colour, linewidth=1.4))

    stamp = ledger.get('run', 'workbench_version', '')
    ax.text(9.95, 0.12,
            f"MeSH Workbench{' v' + stamp if stamp else ''}  -  "
            f"{datetime.now().strftime('%Y-%m-%d')}",
            ha='right', va='bottom', fontsize=7.5, color='#888888')

    os.makedirs(_paths.long_path(output_dir), exist_ok=True)
    written = []
    # SVG always, on top of whatever raster formats are configured. This
    # diagram is almost entirely numbers someone will want to quote, and text
    # in an SVG can be selected and copied out of any browser; text in a JPEG
    # has to be retyped. The .txt report and the ledger carry the same figures
    # for anyone who would rather have them as text or as a table.
    formats = list(_FIGURE_FORMATS)
    if 'svg' not in formats:
        formats.append('svg')
    for ext in formats:
        path = os.path.join(output_dir, f"{file_prefix}_prisma_flow.{ext}")
        try:
            # matplotlib's default is to write SVG text as outlined paths,
            # which looks identical and cannot be selected - exactly what the
            # SVG is here to avoid. 'none' keeps the glyphs as <text>.
            with matplotlib.rc_context({'svg.fonttype': 'none'}):
                fig.savefig(_paths.long_path(path), dpi=_FIGURE_DPI, bbox_inches='tight',
                            facecolor='white', **_FORMAT_KWARGS.get(ext, {}))
            written.append(path)
        except Exception as e:
            print(f"    [!] Error saving the PRISMA figure as .{ext}: {e}")
    plt.close(fig)
    return written


def write_prisma_report(ledger, config):
    """Write both forms of the report. Returns the paths written."""
    if not len(ledger):
        return []
    flow = build_flow(ledger)
    if not flow:
        return []

    paths = []
    text_path = os.path.join(str(config.results_dir), f"{config.prefix}_prisma_flow_report.txt")
    try:
        os.makedirs(_paths.long_path(str(config.results_dir)), exist_ok=True)
        with open(_paths.long_path(text_path), 'w', encoding='utf-8') as fh:
            fh.write(render_text(ledger, flow))
        paths.append(text_path)
    except OSError as e:
        print(f"    [!] Could not write the PRISMA text report: {e}")

    try:
        paths += render_figure(ledger, str(config.figures_dir), config.prefix, flow)
    except Exception as e:
        print(f"    [!] Could not draw the PRISMA figure: {e}")
    return paths
