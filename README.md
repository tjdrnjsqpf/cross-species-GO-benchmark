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
- `code/revision/`, `data/revision/` — analyses added during peer review (revision):
  per-reference ortholog-identity distributions (`R1C1_identity_iqr.py` -> `identity_iqr.tsv`),
  per-list mapping-rate/identity by functional category with covariate-robustness OLS
  (`R1C7_category_mapping.py`, `r1c7_ols_robustness.py`), gene-list inventory
  (`R1C4_merge_lists_meta.py`, `make_tableS4_genelists.py`), bootstrap CIs for ID50/floor
  (`make_fig1b_ci.py`) and for the reliability guide with random-term nulls
  (`make_table1_ci_null.py`), the reference-panel overview figure
  (`make_figS1_panel_overview.py`), and clade-boundary-annotated tolerance panels
  (`make_fig1a_boundaries.R`). Scripts run from `code/revision/` and read/write `data/`
  and `data/revision/`; the three `R1C*` scripts run against the full pipeline output
  (`results/<track>/...`) like the numbered pipeline stages.
- `supplementary/Supplementary_Tables.xlsx` — the assembled Supplementary Tables S1-S12 exactly as
  cited in the manuscript (one sheet per table; the README sheet inside maps each table).
- `RESULTS_SUMMARY.md` — manuscript-ready results narrative with key numbers.

## Where each manuscript table/figure lives
| Manuscript item | Ready-made | Underlying data in this repo |
|---|---|---|
| Table S1 (reference metadata + scientific names) | xlsx sheet S1_refmeta | `data/revision/tableS1_refmeta_with_names.tsv` |
| Table S2 (full per-list set-level results) | — (too large for xlsx) | `data/figdata_setlevel.tsv` |
| Table S3 (per-category ortholog coverage) | xlsx sheet S3_category_mapping | `data/revision/category_mapping.tsv` |
| Table S4 (gene-list inventory) | xlsx sheet S4_gene_lists | `data/revision/tableS4_*.tsv` |
| Table S5 (STRING check) | xlsx sheet S5_STRING | `data/string_conservation.tsv` |
| Table S6 (ID50/floor + fit stats + CIs) | xlsx sheet S6_ID50_floor | `data/id50.tsv` + `data/revision/fig1b_ci.tsv` |
| Table S7 (null Wang baselines) | xlsx sheet S7_null_wang | `data/null_wang.tsv` |
| Table S8 (orthology-method means) | xlsx sheet S8_eggnog_means | `data/eggnog_compare.tsv` |
| Table S9 (expression conservation) | xlsx sheet S9_expression | `data/expr_conservation.tsv` |
| Table S10 (random-orthology control) | xlsx sheet S10_random_orthology | `data/random_orthology_<track>.tsv` |
| Table S11 (size x IC robustness) | xlsx sheet S11_size_x_IC | derived from `data/figdata_setlevel.tsv` (size, list_ic) |
| Table S12 (reliability guide + CIs + null) | xlsx sheet S12_reliability_full | `data/guide_table_unified.tsv` (+ `code/revision/make_table1_ci_null.py`) |
| Table 1 (main text, with CIs + null rows) | `data/revision/table1_extended.tsv` | same script |
| Fig S1 (reference-panel overview) | `data/revision/FigS1_panel_overview.pdf` | `code/revision/make_figS1_panel_overview.py` |
| Fig 1A boundary panels / Fig 1B CIs | `data/revision/Fig1*_boundaries.pdf`, `Fig1B_ID50_floor_CI.pdf` | `code/revision/` |
- `DATA_SOURCES.md` — public databases and versions used (raw data are not redistributed here).

## Reproducing
**Paths.** Scripts resolve the project root automatically (the parent of `code/`), so run them from
inside `code/` or set `export GOTX_ROOT=/path/to/this/repo`. External tools (`diamond`, `emapper.py`,
`Rscript`, `python3`) are expected on `$PATH` (e.g. via the conda environments below).

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
