#!/usr/bin/env python3
"""S9 — conserved (metabolic/housekeeping) vs taxon-specific (developmental/immune) tolerance.
Tests Primmer's core caveat: conserved processes transfer over-optimistically; taxon-specific
processes collapse sooner. BP only (where the contrast is meaningful). Uses existing results.
Outputs: results/crossclade/{category.tsv, fig_conserved_vs_specific.png, category_stats.txt}
"""
import os
import pandas as pd, numpy as np
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = "/var2/lsg/Claude_Code/Cross-species-GeneOntology"
TRACKS = {"fish":"fish","mammal":"mammal","plant_rice":"plant","plant_arabidopsis":"plant",
          "fungi":"fungi","insect":"insect"}
OUTD = f"{ROOT}/results/crossclade"; os.makedirs(OUTD, exist_ok=True)
cat_color = {"conserved":"#2c7fb8","specific":"#de2d26","other":"#969696"}

rows = []
for track, clade in TRACKS.items():
    base = f"{ROOT}/results/{track}"
    if not os.path.exists(f"{base}/metrics/semantic_sim.tsv"): continue
    sem = pd.read_csv(f"{base}/metrics/semantic_sim.tsv", sep="\t")
    meta = pd.read_csv(f"{base}/enrichment/lists_meta.tsv", sep="\t")[["list_id","category","namespace"]]
    m = pd.read_csv(f"{base}/metrics/setlevel_metrics.tsv", sep="\t")
    pid = m[m.ic_bin=="all"].groupby("ref").median_pident.first()
    rec = m[(m.ic_bin=="all")&(m.method=="rbh")&(m.evset=="all")][["list_id","aspect","ref","recall"]]
    sem = sem[sem.list_id.str.startswith("SYN_")].copy()
    sem["ref"]=sem.label.str.split(".").str[0]; sem["method"]=sem.label.str.split(".").str[1]; sem["evset"]=sem.label.str.split(".").str[2]
    sem = sem[(sem.method=="rbh")&(sem.evset=="all")&(sem.aspect=="BP")]
    sem = sem.merge(meta, on="list_id", how="left")
    sem = sem.merge(rec, on=["list_id","aspect","ref"], how="left")
    sem["track"]=track; sem["clade"]=clade; sem["median_pident"]=sem.ref.map(pid)
    rows.append(sem)
df = pd.concat(rows, ignore_index=True)
df.to_csv(f"{OUTD}/category.tsv", sep="\t", index=False)

out = []
# overall: conserved vs specific mean (distance-matched mid range 30-65% id)
mid = df[(df.median_pident>=30)&(df.median_pident<=68)]
out.append("=== Wang (BP, RBH, all-evidence) by category, distance-matched 30-68% id ===")
piv = mid[mid.category.isin(["conserved","specific"])].groupby(["clade","category"]).wang_bma.mean().unstack()
piv["conserved_minus_specific"] = piv.get("conserved") - piv.get("specific")
out.append(piv.round(3).to_string())

# regression: does conserved decay slower with distance? interaction category x distance
d = df[df.category.isin(["conserved","specific"])].dropna(subset=["wang_bma","median_pident"]).copy()
d["pid_z"] = (d.median_pident - d.median_pident.mean())/d.median_pident.std()
mod = smf.ols("wang_bma ~ pid_z * C(category, Treatment('specific')) + C(clade)", data=d).fit()
out.append("\n=== regression: wang ~ pid_z * category(ref=specific) + clade ===")
main_key = [k for k in mod.params.index if k.startswith("C(category") and "T.conserved" in k and "pid_z" not in k]
if main_key:
    k = main_key[0]
    out.append(f"  conserved main effect = {mod.params[k]:+.3f} (p={mod.pvalues[k]:.3g})  [>0 => conserved higher overall]")
ix = [k for k in mod.params.index if 'pid_z:' in k]
for k in ix:
    out.append(f"  interaction {k} = {mod.params[k]:+.3f} (p={mod.pvalues[k]:.3g})  [>0 => conserved decays slower]")
out.append(f"  R2={mod.rsquared:.3f}, n={len(d)}")

# ---- figure: wang vs %identity by category, facet per track ----
tracks_present = [t for t in TRACKS if t in set(df.track)]
agg = df.groupby(["track","category","ref","median_pident"]).wang_bma.mean().reset_index()
ncol=3; nrow=int(np.ceil(len(tracks_present)/ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(5*ncol, 4*nrow), sharey=True, squeeze=False)
for i,track in enumerate(tracks_present):
    ax=axes[i//ncol][i%ncol]
    for cat in ["conserved","specific","other"]:
        t=agg[(agg.track==track)&(agg.category==cat)].sort_values("median_pident",ascending=False)
        if t.empty: continue
        ax.plot(t.median_pident,t.wang_bma,marker="o",color=cat_color[cat],label=cat,alpha=.8,ms=4)
    ax.set_title(track); ax.invert_xaxis(); ax.grid(alpha=.3); ax.set_xlabel("% identity")
axes[0][0].set_ylabel("Wang sim (BP)"); axes[0][0].legend(fontsize=8)
fig.suptitle("Conserved (metabolic) vs taxon-specific (developmental/immune) GO-transfer tolerance — BP, RBH")
fig.tight_layout(); fig.savefig(f"{OUTD}/fig_conserved_vs_specific.svg"); fig.savefig(f"{OUTD}/fig_conserved_vs_specific.pdf"); plt.close(fig)

open(f"{OUTD}/category_stats.txt","w").write("\n".join(out))
print("\n".join(out)); print("\n[S9] wrote -> ", OUTD)
