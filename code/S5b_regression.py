#!/usr/bin/env python3
"""S5b — regression: which factor dominates set-level F1?
 f1 ~ median_pident + ref_richness + method + aspect + evset
Outputs results/metrics/regression_summary.txt
"""
import os
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf

ROOT = "/var2/lsg/Claude_Code/Cross-species-GeneOntology"
import os as _os
OUT = _os.environ.get("GOTX_OUT", ROOT + "/results")
m = pd.read_csv(f"{OUT}/metrics/setlevel_metrics.tsv", sep="\t")
d = m[m.ic_bin == "all"].copy()
d = d.dropna(subset=["f1"])
# standardize continuous predictors for comparable coefficients
for c in ["median_pident", "ref_richness"]:
    d[c+"_z"] = (d[c]-d[c].mean())/d[c].std()

out = []
print(f"{'outcome':10s} {'distance(+id)':>14s} {'richness':>10s} {'R2':>6s}")
for outcome in ["f1", "recall", "precision"]:
    dd = d.dropna(subset=[outcome])
    mod = smf.ols(f"{outcome} ~ median_pident_z + ref_richness_z + C(method) + C(aspect) + C(evset)", data=dd).fit()
    out.append(f"\n################## OUTCOME = {outcome}  (n={len(dd)}, R2={mod.rsquared:.3f}) ##################")
    out.append(mod.summary().as_text())
    out.append(f"# standardized: distance(median_pident_z)={mod.params['median_pident_z']:+.3f}  richness(ref_richness_z)={mod.params['ref_richness_z']:+.3f}")
    print(f"{outcome:10s} {mod.params['median_pident_z']:>+14.3f} {mod.params['ref_richness_z']:>+10.3f} {mod.rsquared:>6.3f}")

with open(f"{OUT}/metrics/regression_summary.txt","w") as f:
    f.write("\n".join(out))
print("\n[S5b] median_pident_z>0 => closer(higher %id) raises the metric.")
print("[S5b] richness effect: should be +recall, ~0 or -precision (the richness tradeoff).")
print("[S5b] wrote results/metrics/regression_summary.txt")
