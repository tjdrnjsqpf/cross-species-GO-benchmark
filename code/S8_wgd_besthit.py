#!/usr/bin/env python3
"""S8 — re-test the plant-WGD penalty in the BEST-HIT method (vs RBH).
Hypothesis: WGD paralogs in plants make best-hit mis-assign orthologs more than in
fish/mammal; RBH filters this. So the best-hit GO-transfer penalty (and steepness)
should be larger for plant.
Three lines of evidence, all from existing on-disk results (no pipeline re-run):
 (1) genomic paralog signature: non-reciprocity & many-to-one collapse of best-hit maps
 (2) GO-transfer method penalty (RBH - besthit) for wang / precision / recall, by clade
 (3) best-hit-only clade x distance interaction (is plant steeper when paralogs unfiltered?)
Outputs: results/crossclade/{wgd_besthit.txt, nonreciprocity.tsv, fig_besthit_vs_rbh.png,
         fig_method_penalty.png}
"""
import os, glob
import pandas as pd, numpy as np
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = "/var2/lsg/Claude_Code/Cross-species-GeneOntology"
TRACKS = {"fish":("zebrafish","fish"), "mammal":("mouse","mammal"),
          "plant_rice":("rice","plant"), "plant_arabidopsis":("arabidopsis","plant")}
OUTD = f"{ROOT}/results/crossclade"
clade_color = {"fish":"#1b9e77","mammal":"#d95f02","plant":"#7570b3"}
out = []

# ---------- (1) genomic paralog signature from mapping files ----------
nr = []
for track,(focal,clade) in TRACKS.items():
    for bh in glob.glob(f"{ROOT}/results/{track}/mapping/{focal}__*.besthit.tsv"):
        ref = os.path.basename(bh).split("__")[1].split(".")[0]
        rbhf = bh.replace(".besthit.",".rbh.")
        b = pd.read_csv(bh, sep="\t"); r = pd.read_csv(rbhf, sep="\t")
        n_bh, n_rbh = len(b), len(r)
        if n_bh == 0: continue
        frac_nonrecip = 1 - n_rbh/n_bh                       # best-hits failing reciprocity
        many_to_one  = 1 - b.ref_acc.nunique()/n_bh          # focal paralogs collapsing on one ref
        nr.append(dict(track=track, clade=clade, ref=ref, n_bh=n_bh, n_rbh=n_rbh,
                       frac_nonrecip=frac_nonrecip, many_to_one=many_to_one,
                       median_pident=b.pident.median()))
nr = pd.DataFrame(nr); nr.to_csv(f"{OUTD}/nonreciprocity.tsv", sep="\t", index=False)
out.append("=== (1) genomic paralog signature (best-hit maps), mean by clade ===")
out.append(nr.groupby("clade")[["frac_nonrecip","many_to_one"]].mean().round(3).to_string())
# controlled for distance: restrict to mid-range refs (45-70% id) so distance is comparable
mid = nr[(nr.median_pident>=45)&(nr.median_pident<=72)]
out.append("\n  same, restricted to 45-72% identity refs (distance-matched):")
out.append(mid.groupby("clade")[["frac_nonrecip","many_to_one"]].mean().round(3).to_string())

# ---------- load GO metrics for all tracks ----------
sem_all, set_all = [], []
for track,(focal,clade) in TRACKS.items():
    s = pd.read_csv(f"{ROOT}/results/{track}/metrics/semantic_sim.tsv", sep="\t")
    m = pd.read_csv(f"{ROOT}/results/{track}/metrics/setlevel_metrics.tsv", sep="\t")
    pid = m[m.ic_bin=="all"].groupby("ref").median_pident.first()
    s = s[s.list_id.str.startswith("SYN_")].copy()
    s["ref"]=s.label.str.split(".").str[0]; s["method"]=s.label.str.split(".").str[1]; s["evset"]=s.label.str.split(".").str[2]
    s=s[s.evset=="all"]; s["track"]=track; s["clade"]=clade; s["median_pident"]=s.ref.map(pid)
    sem_all.append(s)
    mm=m[(m.ic_bin=="all")&(m.evset=="all")].copy(); mm["track"]=track; mm["clade"]=clade
    set_all.append(mm)
sem=pd.concat(sem_all,ignore_index=True); setm=pd.concat(set_all,ignore_index=True)

# ---------- (2) method penalty (RBH - besthit) by clade ----------
out.append("\n=== (2) GO-transfer method penalty = RBH - besthit (positive = besthit worse) ===")
# wang from semantic
w = sem.groupby(["clade","ref","method"]).wang_bma.mean().unstack("method")
w["penalty"]=w["rbh"]-w["besthit"]
# precision/recall from setlevel
pr = setm.groupby(["clade","ref","method"])[["precision","recall"]].mean().reset_index()
prp = pr.pivot_table(index=["clade","ref"], columns="method", values=["precision","recall"])
prp[("precision","penalty")]=prp[("precision","rbh")]-prp[("precision","besthit")]
prp[("recall","penalty")]=prp[("recall","rbh")]-prp[("recall","besthit")]
out.append("  Wang penalty (mean by clade):")
out.append(w.groupby("clade").penalty.mean().round(3).to_string())
out.append("  precision penalty (mean by clade):")
out.append(prp[("precision","penalty")].groupby("clade").mean().round(3).to_string())
out.append("  recall penalty (mean by clade):")
out.append(prp[("recall","penalty")].groupby("clade").mean().round(3).to_string())

# ---------- (3) best-hit-only clade x distance interaction ----------
out.append("\n=== (3) clade x distance interaction, BEST-HIT only (cf. RBH was plant ns) ===")
for meth in ["besthit","rbh"]:
    d=sem[sem.method==meth].dropna(subset=["wang_bma","median_pident"]).copy()
    d["pid_z"]=(d.median_pident-d.median_pident.mean())/d.median_pident.std()
    mod=smf.ols("wang_bma ~ pid_z * C(clade, Treatment('fish')) + C(aspect)", data=d).fit()
    inter={k:(v,mod.pvalues[k]) for k,v in mod.params.items() if "pid_z:" in k}
    line=", ".join(f"{k.split('T.')[-1].rstrip(']')}: slope_diff={v:+.3f}(p={p:.3f})" for k,(v,p) in inter.items())
    out.append(f"  [{meth}] {line}")

# ---------- figures ----------
# besthit vs rbh curves per clade (BP, wang)
agg=sem.groupby(["clade","track","method","ref","aspect","median_pident"]).wang_bma.mean().reset_index()
bp=agg[agg.aspect=="BP"]
fig,axes=plt.subplots(1,2,figsize=(12,4.8),sharey=True)
for ax,meth in zip(axes,["besthit","rbh"]):
    for track,(focal,clade) in TRACKS.items():
        t=bp[(bp.track==track)&(bp.method==meth)].sort_values("median_pident",ascending=False)
        if t.empty: continue
        ax.plot(t.median_pident,t.wang_bma,marker="o",color=clade_color[clade],alpha=.8,label=track,ms=5)
    ax.set_title(f"method = {meth}"); ax.set_xlabel("median ortholog % identity"); ax.invert_xaxis(); ax.grid(alpha=.3)
axes[0].set_ylabel("Wang semantic similarity (BP)"); axes[1].legend(fontsize=8)
fig.suptitle("WGD re-test: best-hit vs RBH tolerance by clade (does plant suffer more without reciprocity?)")
fig.tight_layout(); fig.savefig(f"{OUTD}/fig_besthit_vs_rbh.svg"); fig.savefig(f"{OUTD}/fig_besthit_vs_rbh.pdf"); plt.close(fig)

# method penalty bars
fig,axes=plt.subplots(1,3,figsize=(12,4))
pen_w=w.groupby("clade").penalty.mean(); pen_p=prp[("precision","penalty")].groupby("clade").mean(); pen_r=prp[("recall","penalty")].groupby("clade").mean()
for ax,(pen,lab) in zip(axes,[(pen_w,"Wang"),(pen_p,"precision"),(pen_r,"recall")]):
    cl=["fish","mammal","plant"]; ax.bar(cl,[pen.get(c,np.nan) for c in cl],color=[clade_color[c] for c in cl])
    ax.set_title(f"{lab} penalty (RBH - besthit)"); ax.axhline(0,color="k",lw=.6); ax.grid(alpha=.3,axis="y")
fig.suptitle("Best-hit penalty by clade (higher = best-hit hurts more; WGD predicts plant largest)")
fig.tight_layout(); fig.savefig(f"{OUTD}/fig_method_penalty.svg"); fig.savefig(f"{OUTD}/fig_method_penalty.pdf"); plt.close(fig)

open(f"{OUTD}/wgd_besthit.txt","w").write("\n".join(out))
print("\n".join(out))
print("\n[S8] wrote wgd_besthit.txt + figs ->", OUTD)
