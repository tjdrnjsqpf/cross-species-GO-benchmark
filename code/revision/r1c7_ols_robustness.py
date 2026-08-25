#!/usr/bin/env python3
"""R1-C7 local integration — does the conserved-vs-specific effect survive
controlling for orthology coverage (frac_mapped) and gene-level identity?

Merge server granular table (category_mapping.tsv) with the Wang scores
(category.tsv: BP, RBH, all-evidence) and compare:
  M0 (manuscript): wang ~ pident_z * category + clade      (reproduce baseline)
  M1 (+coverage) : M0 + frac_mapped + med_pident_genes_z
Outputs: r1c7_ols_robustness.txt
"""
import os
import pandas as pd
import statsmodels.formula.api as smf

HERE = os.path.dirname(os.path.realpath(__file__))
ROOT = os.path.join(HERE, "../..")
TABD = os.path.join(HERE, "../../data/revision")

cm = pd.read_csv(f"{ROOT}/data/revision/category_mapping.tsv", sep="\t")
cat = pd.read_csv(f"{ROOT}/data/category.tsv", sep="\t")

cm_bp = cm[cm.namespace == "biological_process"][
    ["track", "ref", "list_id", "category", "n_genes", "frac_mapped", "med_pident_genes"]]
d = cat.merge(cm_bp, on=["track", "ref", "list_id"], suffixes=("", "_srv"))
lines = [f"merged rows: {len(d)} (category.tsv {len(cat)}, server BP {len(cm_bp)})"]
# category consistency check between the two files
mism = (d.category != d.category_srv).sum()
lines.append(f"category mismatch between files: {mism}")
d = d[d.category.isin(["conserved", "specific", "other"])].dropna(subset=["wang_bma", "median_pident", "frac_mapped"])

# standardize
for c, z in [("median_pident", "pident_z"), ("med_pident_genes", "gene_pident_z")]:
    d[z] = (d[c] - d[c].mean()) / d[c].std()

d["category"] = pd.Categorical(d.category, categories=["specific", "conserved", "other"])  # ref level = specific

MODELS = [
    ("M0  wang ~ pident_z*category + clade (manuscript baseline)",
     "wang_bma ~ pident_z * C(category) + C(clade)"),
    ("M1  M0 + frac_mapped",
     "wang_bma ~ pident_z * C(category) + C(clade) + frac_mapped"),
    ("M2  M0 + frac_mapped + gene-level identity",
     "wang_bma ~ pident_z * C(category) + C(clade) + frac_mapped + gene_pident_z"),
]
for name, f in MODELS:
    m = smf.ols(f, data=d).fit()
    lines.append(f"\n=== {name} ===  (n={int(m.nobs)}, R2={m.rsquared:.3f})")
    for term in m.params.index:
        if "category" in term or term in ("frac_mapped", "gene_pident_z", "pident_z"):
            lines.append(f"  {term:45s} {m.params[term]:+.4f}  p={m.pvalues[term]:.2e}")

out = "\n".join(lines)
open(os.path.join(TABD, "r1c7_ols_robustness.txt"), "w").write(out + "\n")
print(out)
