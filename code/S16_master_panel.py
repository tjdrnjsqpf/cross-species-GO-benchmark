#!/usr/bin/env python3
"""S16 — one master summary panel (2x2) of the headline findings, for the manuscript.
(a) cross-clade tolerance curves (BP Wang vs %identity)
(b) ID50 + floor by clade (tolerance ranking)
(c) conserved vs taxon-specific (distance-matched, by clade)
(d) orthology method comparison (besthit/rbh/eggnog) Wang per clade — the plant/RBH story
Outputs: results/crossclade/fig_master_panel.{svg,pdf}
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT="/var2/lsg/Claude_Code/Cross-species-GeneOntology"; OUTD=f"{ROOT}/results/crossclade"
TRACKS=["fish","mammal","plant_rice","plant_arabidopsis","fungi","insect"]
clade_color={"fish":"#1b9e77","mammal":"#d95f02","plant_rice":"#7570b3",
             "plant_arabidopsis":"#9e8fd0","fungi":"#e7298a","insect":"#66a61e"}

cc=pd.read_csv(f"{OUTD}/crossclade.tsv",sep="\t")
id50=pd.read_csv(f"{OUTD}/id50.tsv",sep="\t")
cat=pd.read_csv(f"{OUTD}/category.tsv",sep="\t")
cmp=pd.read_csv(f"{OUTD}/eggnog_compare.tsv",sep="\t")

fig,ax=plt.subplots(2,2,figsize=(13,10))

# (a) tolerance curves BP
a=ax[0,0]
ccb=cc[cc.aspect=="BP"]
for t in TRACKS:
    g=ccb[ccb.track==t].groupby("ref").agg(pid=("median_pident","first"),w=("wang_bma","mean")).sort_values("pid",ascending=False)
    if g.empty: continue
    a.plot(g.pid,g.w,marker="o",ms=4,color=clade_color[t],label=t,alpha=.85)
a.invert_xaxis(); a.set_xlabel("median ortholog % identity"); a.set_ylabel("Wang similarity (BP)")
a.set_title("(a) Cross-clade divergence tolerance (RBH, BP)"); a.grid(alpha=.3); a.legend(fontsize=8)

# (b) ID50 + floor (BP)
b=ax[0,1]
ib=id50[id50.aspect=="BP"].set_index("track").reindex(TRACKS)
x=np.arange(len(TRACKS))
b.bar(x-0.2,ib.ID50_rep,0.4,label="ID50_rep (%id at half-collapse)",color="#555")
b.bar(x+0.2,ib.floor*100,0.4,label="floor x100 (far-distance plateau)",color="#bbb")
b.set_xticks(x); b.set_xticklabels(TRACKS,rotation=30,ha="right",fontsize=8)
b.set_ylabel("value"); b.set_title("(b) ID50 & floor — tolerance ranking (BP)"); b.legend(fontsize=8); b.grid(alpha=.3,axis="y")

# (c) conserved vs specific (BP, distance-matched 30-68)
c=ax[1,0]
cm=cat[(cat.aspect=="BP")&(cat.median_pident.between(30,68))]
cons=cm[cm.category=="conserved"].groupby("track").wang_bma.mean().reindex(TRACKS)
spec=cm[cm.category=="specific"].groupby("track").wang_bma.mean().reindex(TRACKS)
c.bar(x-0.2,cons,0.4,label="conserved (metabolism, translation...)",color="#2c7fb8")
c.bar(x+0.2,spec,0.4,label="taxon-specific (development, immunity...)",color="#de2d26")
c.set_xticks(x); c.set_xticklabels(TRACKS,rotation=30,ha="right",fontsize=8)
c.set_ylabel("Wang (BP, 30-68% id)"); c.set_title("(c) Conserved >> taxon-specific (all clades)"); c.legend(fontsize=8); c.grid(alpha=.3,axis="y")

# (d) method comparison Wang per clade
d=ax[1,1]
mp=cmp.pivot_table(index="track",columns="method",values="wang").reindex(TRACKS)
w=0.26
for i,(meth,col) in enumerate([("besthit","#7570b3"),("eggnog","#d95f02"),("rbh","#1b9e77")]):
    if meth in mp: d.bar(x+(i-1)*w,mp[meth],w,label=meth,color=col)
d.set_xticks(x); d.set_xticklabels(TRACKS,rotation=30,ha="right",fontsize=8)
d.set_ylabel("Wang (eval subset)"); d.set_ylim(0.4,0.9)
d.set_title("(d) Orthology method: RBH hurts plants only"); d.legend(fontsize=8); d.grid(alpha=.3,axis="y")

fig.suptitle("Cross-species GO-transfer divergence tolerance — headline results",fontsize=13,y=1.0)
fig.tight_layout()
fig.savefig(f"{OUTD}/fig_master_panel.svg"); fig.savefig(f"{OUTD}/fig_master_panel.pdf")
print(f"[S16] wrote -> {OUTD}/fig_master_panel.svg/.pdf")
