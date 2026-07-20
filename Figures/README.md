# Figures

Code that renders the manuscript's main and supplemental figures. Two
orchestrator R Markdown documents drive per-panel render scripts and
composition:

- **`make_main_figures.Rmd`** — Figure 1 (isoform discovery + NMD response) and
  Figure 2 (productive-output response).
- **`make_supplemental_figures.Rmd`** — the supplements whose render code lives in
  this repository (ER-stress pathway across cell types, normalization sanity,
  reproducibility battery, SRSF isopair structures, cell-type-specific NMD,
  NMD machinery).

Figures 3–5 (Isopair splice-event attribution and the deep-learning model) are
produced in separate repositories and are not built here. Several S1/S2
supplements (SR↔LR correlation, isoform length, pairwise expression, SQANTI3
categories, mashr sharing, proportion-vs-expression) are rendered in the Isopair
repository and are noted as out of scope inside the supplemental Rmd.

## Pipeline

Each figure is built in two stages:

1. **Per-panel render scripts** (`make_*.R`) read the analysis data bundle and
   write panel PNGs into `fig_panels/`.
2. **Stitching scripts** (`fig_panels/*.py`, matplotlib) compose the panels into
   the multipanel figure.

Run order is handled by the Rmds; see each Rmd's per-figure table for the
panel→script→input mapping.

## Data

These scripts read the analysis data bundle (`nmd_fig_data/`: DGEList objects,
mashr posterior/lfsr CSVs, and example-gene annotation). That bundle is **not**
included here — it contains primary-cell-derived data and is available through
the study data deposit (see the top-level `README`). Place or symlink the bundle
in the working directory before running, so the scripts' relative `nmd_fig_data/`
path resolves.

`make_pathway_allct.R` uses `pathview`; if the KEGG `hsa04141` template is not
cached in `fig_panels/`, it is fetched from KEGG at run time (needs internet).

## Requirements

R (data.table, dplyr, tidyr, tibble, ggplot2, scales, patchwork, magick, edgeR,
pathview, Isopair, org.Hs.eg.db, AnnotationDbi) and Python 3 with matplotlib.
