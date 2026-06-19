#!/usr/bin/env python3
"""S11 — fungi-internal WGD test (addendum track D core): post-WGD vs pre-WGD reference
species, SAME focal (S. cerevisiae, itself post-WGD), distance-matched. Independent test of
the plant-WGD hypothesis within fungi. Predicts: post-WGD refs have more paralog ambiguity
(non-reciprocity), benefit more from RBH, and transfer slightly worse than pre-WGD at equal %id.
Outputs: results/crossclade/{wgd_fungi.txt, fig_wgd_fungi.png}
"""
import os, glob, yaml
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT="/var2/lsg/Claude_Code/Cross-species-GeneOntology"; OUTD=f"{ROOT}/results/crossclade"
cfg=yaml.safe_load(open(f"{ROOT}/config/track_fungi.yaml"))
wgd={s:v.get("wgd","none") for s,v in cfg["species"].items()}
focal=cfg["focal"]; base=f"{ROOT}/results/fungi"
out=[]

# genomic paralog signature from mapping files
nr=[]
for bh in glob.glob(f"{base}/mapping/{focal}__*.besthit.tsv"):
    ref=os.path.basename(bh).split("__")[1].split(".")[0]
    b=pd.read_csv(bh,sep="\t"); r=pd.read_csv(bh.replace(".besthit.",".rbh."),sep="\t")
    if len(b)==0: continue
    nr.append(dict(ref=ref,wgd=wgd.get(ref,"none"),n_bh=len(b),n_rbh=len(r),
                   frac_nonrecip=1-len(r)/len(b), many_to_one=1-b.ref_acc.nunique()/len(b),
                   median_pident=b.pident.median()))
nr=pd.DataFrame(nr)

# GO metrics
sem=pd.read_csv(f"{base}/metrics/semantic_sim.tsv",sep="\t")
m=pd.read_csv(f"{base}/metrics/setlevel_metrics.tsv",sep="\t")
sem=sem[sem.list_id.str.startswith("SYN_")].copy()
sem["ref"]=sem.label.str.split(".").str[0]; sem["method"]=sem.label.str.split(".").str[1]; sem["evset"]=sem.label.str.split(".").str[2]
wang=sem[(sem.evset=="all")&(sem.aspect=="BP")].groupby(["ref","method"]).wang_bma.mean().unstack("method")
pr=m[(m.ic_bin=="all")&(m.evset=="all")].groupby(["ref","method"])[["precision","recall"]].mean()

tab=nr.set_index("ref").join(wang, how="left")
tab["wgd"]=tab.index.map(lambda r: wgd.get(r,"none"))
tab.to_csv(f"{OUTD}/wgd_fungi.txt".replace(".txt",".tsv"),sep="\t")

# controlled contrast: only Saccharomycetaceae with clean post/pre label, similar distance
grp=nr[nr.wgd.isin(["post","pre"])]
out.append("=== fungi-internal: post-WGD vs pre-WGD reference species (focal S.cerevisiae) ===")
out.append(grp.groupby("wgd")[["median_pident","frac_nonrecip","many_to_one"]].mean().round(3).to_string())
out.append("\n  per species:")
out.append(grp[["ref","wgd","median_pident","frac_nonrecip","many_to_one"]].sort_values("wgd").round(3).to_string(index=False))

# method penalty (rbh-besthit) by wgd group
if "rbh" in wang and "besthit" in wang:
    w=wang.copy(); w["wgd"]=w.index.map(lambda r: wgd.get(r,"none")); w["wang_penalty"]=w["rbh"]-w["besthit"]
    sub=w[w.wgd.isin(["post","pre"])]
    out.append("\n=== Wang RBH-besthit penalty by WGD group (post should be >= pre) ===")
    out.append(sub.groupby("wgd").wang_penalty.mean().round(3).to_string())

open(f"{OUTD}/wgd_fungi.txt","w").write("\n".join(out)); print("\n".join(out))

# figure
g=grp.groupby("wgd")[["frac_nonrecip","many_to_one"]].mean()
fig,ax=plt.subplots(1,2,figsize=(9,4))
for a,col in zip(ax,["frac_nonrecip","many_to_one"]):
    gg=g[col]; a.bar(gg.index,gg.values,color=["#d95f02","#1b9e77"]); a.set_title(col); a.grid(alpha=.3,axis="y")
fig.suptitle("Fungi-internal WGD: post-WGD vs pre-WGD paralog ambiguity (focal S. cerevisiae)")
fig.tight_layout(); fig.savefig(f"{OUTD}/fig_wgd_fungi.svg"); fig.savefig(f"{OUTD}/fig_wgd_fungi.pdf")
print(f"\n[S11] wrote -> {OUTD}")
