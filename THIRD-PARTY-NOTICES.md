# Third-party notices

MeSH Workbench is distributed as a self-contained bundle: the installer carries its own Python interpreter and every library it needs, as prebuilt wheels, unmodified from their published releases.

This file lists what those are. It is generated from the wheels a built bundle actually contains — not from a development environment — by

```bash
python packaging/make_third_party_notices.py <bundle.tar.gz>
```

MeSH Workbench itself is MIT licensed; see [LICENSE](LICENSE).

---

## Licences that ask for more than attribution

### certifi — MPL-2.0

File-level copyleft. The files are redistributed unmodified; source for them is at the URL below.

Source and licence: <https://github.com/certifi/python-certifi>

### gensim — LGPL

Weak copyleft. The library is redistributed unmodified and is dynamically imported, never linked into a derived work. Its source and licence are at the URL below.

Source and licence: <https://radimrehurek.com/gensim/>

### orjson — MPL-2.0

File-level copyleft. The files are redistributed unmodified; source for them is at the URL below.

Source and licence: <https://github.com/ijl/orjson>

### tqdm — MPL-2.0

File-level copyleft. The files are redistributed unmodified; source for them is at the URL below.

Source and licence: <https://tqdm.github.io>

---

## Everything redistributed (45 packages)

| Package | Version | Licence | Project |
| :--- | :--- | :--- | :--- |
| adjusttext | 1.4.0 | MIT | <https://github.com/Phlya/adjustText> |
| biopython | 1.87 | LicenseRef-Biopython-License-Agreement | <https://biopython.org/> |
| certifi | 2026.6.17 | Mozilla Public License 2.0 (MPL 2.0) | <https://github.com/certifi/python-certifi> |
| charset-normalizer | 3.4.7 | MIT | <https://github.com/jawah/charset_normalizer/blob/master/CHANGELOG.md> |
| choreographer | 1.3.0 | MIT (full text in the package) | <https://github.com/plotly/choreographer> |
| colorama | 0.4.6 | BSD License | <https://github.com/tartley/colorama> |
| contourpy | 1.3.3 | BSD License | <https://github.com/contourpy/contourpy> |
| cycler | 0.12.1 | BSD License | <https://matplotlib.org/cycler/> |
| et-xmlfile | 2.0.0 | MIT License | <https://foss.heptapod.net/openpyxl/et_xmlfile> |
| fonttools | 4.63.0 | MIT | <http://github.com/fonttools/fonttools> |
| gensim | 4.4.0 | LGPL-2.1-only | <https://radimrehurek.com/gensim/> |
| idna | 3.18 | BSD-3-Clause | <https://github.com/kjd/idna> |
| joblib | 1.5.3 | BSD-3-Clause | <https://joblib.readthedocs.io> |
| kaleido | 1.3.0 | The MIT License (MIT) | <https://github.com/plotly/kaleido> |
| kiwisolver | 1.5.0 | BSD License | <https://github.com/nucleic/kiwi> |
| logistro | 2.0.1 | MIT License | <https://github.com/geopozo/logistro> |
| matplotlib | 3.11.0 | Python Software Foundation License | <https://matplotlib.org> |
| narwhals | 2.22.1 | MIT | <https://github.com/narwhals-dev/narwhals> |
| networkx | 2.8.8 | BSD License | <https://networkx.org/> |
| numpy | 1.26.4 | BSD License | <https://numpy.org> |
| openpyxl | 3.1.5 | MIT License | <https://openpyxl.readthedocs.io> |
| orjson | 3.11.9 | MPL-2.0 AND (Apache-2.0 OR MIT) | <https://github.com/ijl/orjson> |
| packaging | 26.2 | Apache-2.0 OR BSD-2-Clause | <https://github.com/pypa/packaging> |
| pandas | 3.0.3 | BSD License | <https://pandas.pydata.org> |
| patsy | 1.0.2 | BSD License | <https://github.com/pydata/patsy> |
| pillow | 12.2.0 | MIT-CMU | <https://python-pillow.github.io> |
| platformdirs | 4.10.0 | MIT | <https://github.com/tox-dev/platformdirs> |
| plotly | 6.8.0 | MIT | <https://plotly.com/python/> |
| pymarc | 5.3.1 | MIT (full text in the package) |  |
| pyparsing | 3.3.2 | MIT | <https://github.com/pyparsing/pyparsing/> |
| python-dateutil | 2.9.0.post0 | Apache Software License; BSD License | <https://github.com/dateutil/dateutil> |
| python-louvain | 0.16 | BSD License | <https://github.com/taynaud/python-louvain> |
| requests | 2.34.2 | Apache Software License | <https://github.com/psf/requests> |
| scikit-learn | 1.9.0 | BSD-3-Clause | <https://scikit-learn.org> |
| scipy | 1.17.1 | BSD License | <https://scipy.org/> |
| seaborn | 0.13.2 | BSD License | <https://github.com/mwaskom/seaborn> |
| simplejson | 4.1.1 | MIT OR AFL-2.1 | <https://github.com/simplejson/simplejson> |
| six | 1.17.0 | MIT License | <https://github.com/benjaminp/six> |
| smart-open | 7.7.0 | MIT License | <https://github.com/piskvorky/smart_open> |
| statsmodels | 0.14.6 | BSD License | <https://www.statsmodels.org/> |
| threadpoolctl | 3.6.0 | BSD License | <https://github.com/joblib/threadpoolctl> |
| tqdm | 4.68.3 | MPL-2.0 AND MIT | <https://tqdm.github.io> |
| tzdata | 2026.2 | Apache-2.0 | <https://github.com/python/tzdata> |
| urllib3 | 2.7.0 | MIT | <https://github.com/urllib3/urllib3/blob/main/CHANGES.rst> |
| wrapt | 2.2.2 | BSD-2-Clause | <https://github.com/GrahamDumpleton/wrapt> |


> Not resolvable from the environment this was generated in: setuptools, wheel.

---

The bundled Python interpreter is CPython, redistributed under the Python Software Foundation License, from the [python-build-standalone](https://github.com/astral-sh/python-build-standalone) project on Linux and macOS and from the [python.org embeddable package](https://www.python.org/downloads/windows/) on Windows.

