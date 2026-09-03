# -*- coding: utf-8 -*-
"""
citation.py - how to cite this program, and what it is built on.

One place, read by the Help menu and available to the CLI, so the citation a
user copies and the one in CITATION.cff cannot drift apart.

The method references are here rather than buried in a docstring because they
are the ones a reader of the resulting paper needs: the edge filter and the
ground-truth set are both other people's work, and a methods section that
describes them without citing them is incomplete.
"""

from importlib.metadata import PackageNotFoundError, version as _pkg_version

TITLE = ('MeSH Workbench: MeSH co-occurrence concept networks for '
         'Adverse Outcome Pathways')
AUTHOR_FAMILY = 'Sax'
AUTHOR_GIVEN = 'Jakob'
AFFILIATION = 'Karolinska Institutet'
YEAR = 2026
REPO = 'https://github.com/Tox-pub/Mesh-Network-Analysis'


def package_version(default='3.2.0'):
    """The installed version, or the fallback.

    `version()` does not only raise when the package is absent - on an editable
    install whose metadata carries no Version field it returns None, and that
    went straight into a citation reading "Version None". Anything falsy is
    treated as not found.
    """
    for name in ('mesh_aop_network', 'mesh-aop-network'):
        try:
            found = _pkg_version(name)
        except (PackageNotFoundError, Exception):       # noqa: B014
            continue
        if found:
            return str(found)
    return default


# Work this program implements rather than invents. GLF is Dianati's Global
# Likelihood Filter - the name is his, and the codebase spelled it three
# different wrong ways ("Graph Likelihood Filtering", "Global-Local
# Filtering", "global likelihood filtering") before this was checked against
# the paper.
REFERENCES = [
    ('Edge filtering (GLF)',
     'Dianati, N. (2016). Unwinding the hairball graph: Pruning algorithms '
     'for weighted complex networks. Physical Review E, 93, 012304. '
     'https://doi.org/10.1103/PhysRevE.93.012304',
     'The Global Likelihood Filter this pipeline runs as one of its two '
     'optimisers. The same paper introduces the Marginal Likelihood Filter '
     '(MLF), used here for per-edge p-values.'),
    ('Ground-truth set (reference corpus)',
     'OECD (2014). The Adverse Outcome Pathway for Skin Sensitisation '
     'Initiated by Covalent Binding to Proteins. OECD Series on Testing and '
     'Assessment No. 168. OECD Publishing, Paris. '
     'https://doi.org/10.1787/9789264221444-en',
     'AOP 40. The curated positive set the bundled reference corpus is '
     'validated against comes from the references of this document.'),
    ('Community detection',
     'Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008). '
     'Fast unfolding of communities in large networks. Journal of Statistical '
     'Mechanics: Theory and Experiment, 2008(10), P10008.',
     'The Louvain method, used to partition the consensus network.'),
    ('Early-recognition metric',
     'Truchon, J.-F., & Bayly, C. I. (2007). Evaluating virtual screening '
     'methods: Good and bad metrics for the "early recognition" problem. '
     'Journal of Chemical Information and Modeling, 47(2), 488-508.',
     'BEDROC, the benchmark\'s headline ranking metric.'),
    ('Article citation counts',
     'Hutchins, B. I., Baker, K. L., Davis, M. T., et al. (2019). The NIH '
     'Open Citation Collection: A public access, broad coverage resource. '
     'PLOS Biology, 17(10), e3000385.',
     'The iCite API, which supplies the citation side of the Article Impact '
     'Score.'),
]


def citation_text(version=None):
    """A single-line citation, for a bibliography."""
    v = version or package_version()
    return (f'{AUTHOR_FAMILY}, {AUTHOR_GIVEN[0]}. ({YEAR}). {TITLE} '
            f'(Version {v}) [Computer software]. {REPO}')


def bibtex(version=None):
    v = version or package_version()
    return (f'@software{{sax{YEAR}meshworkbench,\n'
            f'  author       = {{{AUTHOR_FAMILY}, {AUTHOR_GIVEN}}},\n'
            f'  title        = {{{TITLE}}},\n'
            f'  year         = {{{YEAR}}},\n'
            f'  version      = {{{v}}},\n'
            f'  url          = {{{REPO}}},\n'
            f'  organization = {{{AFFILIATION}}}\n'
            f'}}')


def citation_markdown(version=None):
    """The whole Cite panel, as Markdown."""
    v = version or package_version()
    out = [
        '# Cite this Program',
        '',
        'If you use this software, or a network it produced, in published '
        'work, please cite it.',
        '',
        '## Citation',
        '',
        '```',
        citation_text(v),
        '```',
        '',
        '## BibTeX',
        '',
        '```',
        bibtex(v),
        '```',
        '',
        '## Also cite the methods this rests on',
        '',
        'These are other people\'s work, implemented here. A methods section '
        'that describes the filtering or the ground truth without citing them '
        'is incomplete.',
        '',
    ]
    for what, ref, why in REFERENCES:
        out += [f'**{what}**', '', f'{ref}', '', f'{why}', '']
    out += ['---', '',
            'A machine-readable version ships with the program as '
            'CITATION.cff.', '']
    return '\n'.join(out)
