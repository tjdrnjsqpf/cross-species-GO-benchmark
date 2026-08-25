#!/usr/bin/env python3
"""Revision R1-10 — Table 1 extended with 95% bootstrap CIs and the random-term null.

Cell = mean per-reference Wang (transferred vs truth; RBH, all-evidence, synthetic
lists) within each focal x distance-band x aspect; CI = percentile bootstrap over
the references in the cell (10,000 resamples; cells with n=1 get no CI).
Null = expected Wang between the truth set and a size-matched random term set
(null_truth_rand_mean from the pipeline's S20 baseline), per focal x aspect.

Inputs : DATA_DEPOSIT/data/crossclade_ci.tsv, null_wang.tsv, guide_table_unified.tsv
Outputs: table1_extended.tsv, table1_extended.md (+ verification vs guide table)
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.realpath(__file__))
DATA = os.path.join(HERE, "../../data")
TABD = os.path.join(HERE, "../../data/revision")
RNG = np.random.default_rng(20260825)
NBOOT = 10000

null = pd.read_csv(f"{DATA}/null_wang.tsv", sep="\t")
guide = pd.read_csv(f"{DATA}/guide_table_unified.tsv", sep="\t")

# per-list rows, matching S15 exactly: RBH, all-evidence, SYN_ lists (small synthetic)
sl = pd.read_csv(f"{DATA}/figdata_setlevel.tsv", sep="\t")
cc = sl[(sl.method == "rbh") & (sl.evset == "all") &
        (sl.list_id.str.startswith("SYN_"))].dropna(subset=["wang_bma"]).copy()

def band(pid):
    return "near" if pid > 70 else ("mid" if pid >= 55 else "far")

cc["band"] = cc.median_pident.map(band)

ORDER = ["mammal", "fish", "insect", "plant_rice", "plant_arabidopsis", "fungi"]
LABEL = {"mammal": "Mammal", "fish": "Fish", "insect": "Insect",
         "plant_rice": "Plant-R", "plant_arabidopsis": "Plant-A", "fungi": "Fungi"}
BANDS = ["near", "mid", "far"]
ASPECTS = ["BP", "MF", "CC"]

rows = []
for t in ORDER:
    for b in BANDS:
        row = dict(Focal=LABEL[t], Distance=b)
        for asp in ASPECTS:
            g = cc[(cc.track == t) & (cc.band == b) & (cc.aspect == asp)]
            if len(g) == 0:
                row[asp] = ""
                continue
            mu = g.wang_bma.mean()                     # pooled-list mean, as in S15
            refs = g.ref.unique()
            if len(refs) > 1:                          # cluster bootstrap over references
                per_ref = {r: v.wang_bma.values for r, v in g.groupby("ref")}
                boots = []
                for _ in range(NBOOT):
                    pick = RNG.choice(refs, size=len(refs), replace=True)
                    boots.append(np.concatenate([per_ref[r] for r in pick]).mean())
                lo, hi = np.percentile(boots, [2.5, 97.5])
                row[asp] = f"{mu:.2f} [{lo:.2f}-{hi:.2f}]"
            else:
                row[asp] = f"{mu:.2f} [n=1]"
            row[f"_n_{asp}"] = len(refs)
        row["n"] = row.pop("_n_BP", "")
        rows.append(row)
    # null row for this focal
    nr = dict(Focal=LABEL[t], Distance="random-term null", n="")
    for asp in ASPECTS:
        nv = null[(null.track == t) & (null.aspect == asp)].null_truth_rand_mean
        nr[asp] = f"{float(nv.iloc[0]):.2f}" if len(nv) else ""
    rows.append(nr)

tbl = pd.DataFrame(rows)[["Focal", "Distance", "n", "BP", "MF", "CC"]]
tbl.to_csv(os.path.join(TABD, "table1_extended.tsv"), sep="\t", index=False)

# markdown render
md = ["| Focal | Distance | n | BP | MF | CC |", "|---|---|---|---|---|---|"]
for _, r in tbl.iterrows():
    md.append(f"| {r.Focal} | {r.Distance} | {r.n} | {r.BP} | {r.MF} | {r.CC} |")
open(os.path.join(TABD, "table1_extended.md"), "w").write("\n".join(md) + "\n")
print("\n".join(md))

# ---- verification against the published guide table (rbh rows) ----
g = guide[guide.method == "rbh"].copy()
g["track"] = g.clade.map({"mammal": "mammal", "fish": "fish", "insect": "insect",
                          "plant_rice": "plant_rice", "plant_arabidopsis": "plant_arabidopsis",
                          "fungi": "fungi"})
g["band"] = g.distance.str.split(" ").str[0]
bad = 0
for _, r in g.iterrows():
    sub = cc[(cc.track == r.track) & (cc.band == r.band) & (cc.aspect == r.aspect)]
    if len(sub) == 0:
        continue
    if abs(sub.wang_bma.mean() - r.wang) > 0.006 or sub.ref.nunique() != r.n_ref:
        bad += 1
        print(f"[mismatch] {r.track} {r.band} {r.aspect}: guide {r.wang} (n={r.n_ref}) vs recomputed {sub.wang_bma.mean():.3f} (n={sub.ref.nunique()})")
print(f"\n[check] guide-table reproduction: {len(g)-bad}/{len(g)} cells match (means & n)")
