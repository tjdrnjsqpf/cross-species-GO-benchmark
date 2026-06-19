#!/usr/bin/env python3
"""S15 — unified practical guide table across all 6 tracks: when can you trust borrowed GO?
For each clade x aspect x distance-bin (by median %identity) x method, mean Wang semantic
similarity -> reliability label (reliable>=0.6, caution>=0.45, else unreliable).
Outputs: results/crossclade/guide_table_unified.tsv (+ .md)
"""
import os
import numpy as np, pandas as pd

ROOT="/var2/lsg/Claude_Code/Cross-species-GeneOntology"; OUTD=f"{ROOT}/results/crossclade"
TRACKS={"fish":"zebrafish","mammal":"mouse","plant_rice":"rice",
        "plant_arabidopsis":"arabidopsis","fungi":"yeast","insect":"fruitfly"}
BINS=[0,55,70,100]; BLAB=["far (<55%)","mid (55-70%)","near (>70%)"]

def label(w):
    return "reliable" if w>=0.6 else ("caution" if w>=0.45 else "unreliable")

rows=[]
for t,focal in TRACKS.items():
    sem=pd.read_csv(f"{ROOT}/results/{t}/metrics/semantic_sim.tsv",sep="\t")
    sem["ref"]=sem.label.str.split(".").str[0]; sem["method"]=sem.label.str.split(".").str[1]; sem["evset"]=sem.label.str.split(".").str[2]
    sem=sem[(sem.evset=="all")&(sem.list_id.str.startswith("SYN_"))]
    m=pd.read_csv(f"{ROOT}/results/{t}/metrics/setlevel_metrics.tsv",sep="\t")
    pid=m[m.ic_bin=="all"].groupby("ref").median_pident.first()
    sem["pid"]=sem.ref.map(pid)
    sem["dist"]=pd.cut(sem.pid,bins=BINS,labels=BLAB)
    for (asp,dist,meth),g in sem.groupby(["aspect","dist","method"],observed=True):
        if meth not in ("besthit","rbh"): continue
        w=g.wang_bma.mean()
        rows.append(dict(clade=t,aspect=asp,distance=str(dist),method=meth,
                         wang=round(w,2),n_ref=g.ref.nunique(),reliability=label(w)))
G=pd.DataFrame(rows)
order_asp={"MF":0,"BP":1,"CC":2}; order_d={"near (>70%)":0,"mid (55-70%)":1,"far (<55%)":2}
G["_a"]=G.aspect.map(order_asp); G["_d"]=G.distance.map(order_d)
G=G.sort_values(["clade","_a","_d","method"]).drop(columns=["_a","_d"])
G.to_csv(f"{OUTD}/guide_table_unified.tsv",sep="\t",index=False)

# compact markdown: RBH only (recommended default), aspect x distance grid per clade
rbh=G[G.method=="rbh"].pivot_table(index=["clade","aspect"],columns="distance",
        values="wang",aggfunc="first")[ ["near (>70%)","mid (55-70%)","far (<55%)"] ]
with open(f"{OUTD}/guide_table_unified.md","w") as f:
    f.write("# Practical guide: trust of borrowed-GO enrichment (RBH, Wang similarity)\n\n")
    f.write("reliable >=0.6 | caution 0.45-0.6 | unreliable <0.45\n\n")
    f.write(rbh.round(2).to_markdown())
print(G.to_string(index=False))
print("\n[S15] wrote -> guide_table_unified.tsv + .md")
