# Cross-species GO-transfer divergence-tolerance benchmark — code & processed data

Code and processed result tables underlying the manuscript *"How far can you borrow? Divergence limits
of cross-species Gene Ontology enrichment in non-model organisms."*

## Contents
- `code/` — analysis pipeline (S0–S16: data download, proteome prep, focal truth, orthology
  mapping, GO transfer, enrichment, set-level metrics, ID50, WGD, robustness, orthogonal validation,
  cross-clade, guide table) plus the figure/table generators (`make_figures.R`, `export_tables.R`).
- `data/` — processed result tables:
  - `figdata_setlevel.tsv` — per-list set-level metrics (the master table).
  - `figdata_ic.tsv`, `figdata_refmeta.tsv` — IC-stratified metrics; reference metadata.
  - `crossclade_ci.tsv`, `id50.tsv`, `category.tsv`, `guide_table_unified.tsv`,
    `nonreciprocity.tsv`, `wgd_fungi.tsv`, `string_conservation.tsv`, `expr_conservation.tsv`,
    `regression_clade.txt`, `master_table.tsv`, `eggnog_compare.tsv`, etc.
- `RESULTS_SUMMARY.md` — manuscript-ready results narrative with key numbers.
- `DATA_SOURCES.md` — public databases and versions used (raw data are not redistributed here).

## Reproducing
Raw inputs are obtained from public resources (see `DATA_SOURCES.md`) by `code/S0_download.*`.
The pipeline runs per track via `code/run_track.sh <config>` and `code/S13_run_all.sh`; figures and
tables are produced by `code/make_figures.R` and `code/export_tables.R` from the `data/` tables.
Two pinned conda environments are used (Python: diamond, goatools, pandas, scipy, statsmodels,
biopython; R 4.4: clusterProfiler, GOSemSim, topGO, GO.db, ggplot2).

## Note on raw data
All raw inputs derive from public databases and are reproducible via the scripts; they are not
redistributed here to avoid licensing/size issues. See `DATA_SOURCES.md`.

## License
Code under the MIT License (`LICENSE`); processed data tables under CC-BY-4.0.

## Citation
[Manuscript citation — to add upon acceptance.]
