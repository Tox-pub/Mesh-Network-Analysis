# Reference figures — provenance

These are the figures from the run the paper reports, not output regenerated
later. They accompany the reference corpus in `data/reference_processed/` and
describe the same network: **173 MeSH terms, 328 co-occurrence relations**, one
citation generation (G1) expanded from `"Dermatitis, Allergic Contact"[MeSH]`.

Anyone can reproduce them by turning on **Use bundled reference data** and
running the pipeline; these files are what that run produced when it was done
for publication, kept so a reader can compare rather than take it on trust.

## Where they came from

| Source | What it contributed |
| :--- | :--- |
| `fig_rework/Consolidated/Figures/G1` | The main figures — communities, joint plot, t-SNE, dumbbell, scatter panels, dendrogram, alluvial, and the ground-truth and Fig3 validation panels |
| `fig_rework/Consolidated/Figures/Supplementary` | Edge-weight and dispersion distributions, the GLF/SA trajectory, the full-scorer enrichment supplement, NPMI edge validation |
| `fig_rework/Consolidated/HTML_Reports/G1` | The interactive alluvial flows and the validation report |
| `fig_rework/Consolidated/Tables/G1` | The alluvial connections table underlying Figure 6 |

Copied unchanged, one file for one file. Nothing here was re-rendered, resized,
or recompressed on the way in.

## Which stage produces each

Running with reference data on reproduces these in `results/figures/`:

- **Step 3 (network)** — `distribution_plot`, `EdgeWeightDistribution_Comparison`,
  `GLF_SA_Trajectory`
- **Step 4 (figures)** — `Louvain_Community_Composition`, `Joint_Plot_Final`,
  `T-SNE_Final`, `Alluvial_Clean` / `Alluvial_Labeled` (+ the connections table),
  `MRS_Dumbbell_Plot_Top20`, `Panel_C_Scatter_Betweenness`,
  `Panel_D_Scatter_PageRank`, `Dendrogram_Node2Vec_AOP`
- **Step 5 (benchmark & validation)** — `benchmark_enrichment`, the `GT_*` panels,
  `Fig3A`–`Fig3D`, `FigS_full_scorer_enrichment`,
  `Network_edge_validation_NPMI`, `validation_report.html`

A run also writes `*_prisma_flow.*` and `*_run_ledger.csv`, which post-date this
figure set and so have no counterpart here.

## Formats

JPEG or PNG for reading, TIFF (LZW) for submission, HTML where the figure is
interactive. A live run writes whatever the **Figures** tab specifies — 300 dpi
JPEG + TIFF by default — so regenerated files may differ in format from these
without differing in content.
