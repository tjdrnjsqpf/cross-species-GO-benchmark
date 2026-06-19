# Data sources

All raw inputs derive from the public resources below and are reproducible via `code/S0_download.*`.
Record the exact version/release and access date used (fill in [ ]) for the Methods/Data-availability statement.

| Resource | Use | Version / release | Access date |
|---|---|---|---|
| UniProt Reference Proteomes | reference & focal proteomes (one protein/gene) | [release] | [date] |
| GO Annotation Database (GOA) / GAF | GO annotations (focal truth = experimental codes) | [release] | [date] |
| Gene Ontology (go-basic.obo) | ontology structure, true-path propagation, IC | [release] | [date] |
| QuickGO | taxon GO downloads for refs lacking a dedicated GAF | [release] | [date] |
| eggNOG / eggNOG-mapper | orthologous-group transfer | [version] | [date] |
| Bgee | cross-species expression conservation (orthogonal validation) | [release] | [date] |
| EBI Expression Atlas | real differentially expressed gene (DEG) sets | [release] | [date] |
| TimeTree | divergence times (My) | [version] | [date] |
| STRING | network-conservation secondary check (Supplementary) | [version] | [date] |

Focal/reference species and accessions are listed in `data/figdata_refmeta.tsv`.
