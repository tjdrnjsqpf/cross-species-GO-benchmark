#!/usr/bin/env python3
"""S13d — compare 3 transfer methods (besthit / rbh / eggnog) on the eval subset.
Q: does eggNOG-OG-curated orthology beat reciprocity (rbh) and loose besthit for GO transfer?
Uses Wang semantic similarity + set-level recall/precision on the SAME synthetic lists & refs.
Inputs: results/eggnog_eval/<track>/metrics/{semantic_sim,setlevel_metrics}.tsv
Outputs: results/crossclade/{eggnog_compare.tsv, eggnog_compare.txt, fig_eggnog_methods.svg/.pdf}
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT="/var2/lsg/Claude_Code/Cross-species-GeneOntology"
OUTD=f"{ROOT}/results/crossclade"
TRACKS=["fish","mammal","plant_rice","plant_arabidopsis","fungi","insect"]
rows=[]
for t in TRACKS:
    base=f"{ROOT}/results/eggnog_eval/{t}/metrics"
    sp=f"{base}/semantic_sim.tsv"; mp=f"{base}/setlevel_metrics.tsv"
    if not os.path.exists(sp):
        print(f"[S13d] {t}: no eval metrics yet, skip"); continue
    sem=pd.read_csv(sp,sep="\t")
    sem["ref"]=sem.label.str.split(".").str[0]; sem["method"]=sem.label.str.split(".").str[1]; sem["evset"]=sem.label.str.split(".").str[2]
    sem=sem[(sem.evset=="all")&(sem.list_id.str.startswith("SYN_"))]
    m=pd.read_csv(mp,sep="\t")
    pid=m[m.ic_bin=="all"].groupby("ref").median_pident.first()
    mm=m[(m.ic_bin=="all")&(m.evset=="all")]
    for (ref,method),g in sem.groupby(["ref","method"]):
        sub=mm[(mm.ref==ref)&(mm.method==method)]
        rows.append(dict(track=t,ref=ref,method=method,median_pident=pid.get(ref,np.nan),
                         wang=g.wang_bma.mean(), recall=sub.recall.mean(), precision=sub.precision.mean()))
df=pd.DataFrame(rows)
if df.empty:
    print("[S13d] no eval data; run S13a->S13c first"); raise SystemExit
df.to_csv(f"{OUTD}/eggnog_compare.tsv",sep="\t",index=False)

out=["=== eggNOG vs besthit vs rbh — mean over eval subset (synthetic, all-evidence) ==="]
agg=df.groupby("method")[["wang","recall","precision"]].mean().round(3)
out.append(agg.to_string())
# paired delta vs besthit at matched track/ref
piv=df.pivot_table(index=["track","ref","median_pident"],columns="method",values=["wang","recall","precision"])
out.append("\n=== paired means (method - besthit) at matched ref ===")
for met in ["wang","recall","precision"]:
    sub=piv[met]
    if {"besthit","rbh","eggnog"}.issubset(sub.columns):
        out.append(f"  {met}: rbh-besthit={ (sub['rbh']-sub['besthit']).mean():+.3f}   "
                   f"eggnog-besthit={ (sub['eggnog']-sub['besthit']).mean():+.3f}   "
                   f"eggnog-rbh={ (sub['eggnog']-sub['rbh']).mean():+.3f}")
out.append("\n=== per-track Wang by method (BP+MF+CC pooled) ===")
out.append(df.pivot_table(index="track",columns="method",values="wang").round(3).to_string())
open(f"{OUTD}/eggnog_compare.txt","w").write("\n".join(out)); print("\n".join(out))

# figure: precision-recall by method (eggnog should raise precision if it filters paralogs)
fig,ax=plt.subplots(1,2,figsize=(12,4.8))
col={"besthit":"#7570b3","rbh":"#1b9e77","eggnog":"#d95f02"}
for meth,g in df.groupby("method"):
    ax[0].scatter(g.recall,g.precision,s=30,alpha=.6,color=col.get(meth,"k"),label=meth)
ax[0].set_xlabel("recall");ax[0].set_ylabel("precision");ax[0].set_title("set-level P-R by orthology method");ax[0].legend();ax[0].grid(alpha=.3)
mb=df.groupby("method")[["recall","precision","wang"]].mean()
mb.plot.bar(ax=ax[1]);ax[1].set_title("mean recall/precision/Wang by method");ax[1].grid(alpha=.3,axis="y")
fig.tight_layout();fig.savefig(f"{OUTD}/fig_eggnog_methods.svg");fig.savefig(f"{OUTD}/fig_eggnog_methods.pdf")
print(f"\n[S13d] wrote -> {OUTD}/eggnog_compare.tsv/.txt + fig_eggnog_methods.svg/.pdf")
