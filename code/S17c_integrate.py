#!/usr/bin/env python3
"""S17c — integrate orthogonal validation: does GO-INDEPENDENT functional conservation
(STRING network; optionally Bgee expression) track the GO-transfer tolerance curve?
Core result: Spearman(GO_transfer_Wang, orthogonal_conservation), per track + pooled.
Figures: S17-fig1 (scatter, the message in one panel), S17-fig2 (both tolerance curves overlaid).
Outputs: results/crossclade/{orthogonal_validation.txt, fig_S17_scatter.*, fig_S17_curves.*}
"""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT=os.environ.get("GOTX_ROOT", os.path.dirname(os.path.dirname(os.path.realpath(__file__)))); OUTD=f"{ROOT}/results/crossclade"
clade_color={"fish":"#1b9e77","mammal":"#d95f02","plant_rice":"#7570b3",
             "plant_arabidopsis":"#9e8fd0","fungi":"#e7298a","insect":"#66a61e"}

# each modality = (own table, value column, label). STRING and Bgee cover DIFFERENT ref subsets
# so they are kept as independent tables (do NOT inner-merge — that loses pairs).
st=pd.read_csv(f"{OUTD}/string_conservation.tsv",sep="\t")
mods=[(st,"string_conservation","STRING neighbour Jaccard (exp+coexp)")]
if "edge_rate" in st.columns:
    mods.insert(0,(st,"edge_rate","STRING edge-conservation rate (exp+coexp)"))
expf=f"{OUTD}/expr_conservation.tsv"
if os.path.exists(expf):
    ex=pd.read_csv(expf,sep="\t")
    mods.append((ex,"expr_conservation","Bgee expression conservation"))

out=[]
for df_m,col,lab in mods:
    d=df_m.dropna(subset=[col,"GO_Wang"])
    out.append(f"\n===== {lab} =====")
    # per-track
    for t in sorted(d.track.unique()):
        s=d[d.track==t]
        if len(s)<4: out.append(f"  {t}: n={len(s)} (too few)"); continue
        rho,p=spearmanr(s[col],s.GO_Wang)
        out.append(f"  {t:18s}: rho={rho:+.3f} p={p:.3f} n={len(s)}")
    rho,p=spearmanr(d[col],d.GO_Wang)
    out.append(f"  POOLED            : rho={rho:+.3f} p={p:.2e} n={len(d)}")
    # also vs %identity (does orthogonal score itself decay with distance?)
    rd,pd_=spearmanr(d[col],d.median_pident)
    out.append(f"  (sanity) {col} vs %identity: rho={rd:+.3f} p={pd_:.2e}")
txt="=== Orthogonal validation: GO-independent functional conservation vs GO-transfer Wang ===\n"+"\n".join(out)
open(f"{OUTD}/orthogonal_validation.txt","w").write(txt); print(txt)

# ---- S17-fig1: scatter (orthogonal x, GO Wang y), colour by track, pooled regression ----
fig,axes=plt.subplots(1,len(mods),figsize=(6.5*len(mods),5),squeeze=False)
for ax,(df_m,col,lab) in zip(axes[0],mods):
    d=df_m.dropna(subset=[col,"GO_Wang"])
    for t in d.track.unique():
        s=d[d.track==t]; ax.scatter(s[col],s.GO_Wang,color=clade_color.get(t,"k"),label=t,s=40,alpha=.8)
    if len(d)>=4:
        rho,p=spearmanr(d[col],d.GO_Wang)
        b,a=np.polyfit(d[col],d.GO_Wang,1); xs=np.linspace(d[col].min(),d[col].max(),50)
        ax.plot(xs,a+b*xs,"k--",alpha=.6)
        ax.set_title(f"{lab}\nSpearman rho={rho:+.2f}, p={p:.1e}, n={len(d)}",fontsize=10)
    ax.set_xlabel(lab); ax.set_ylabel("GO-transfer Wang (BP)"); ax.grid(alpha=.3); ax.legend(fontsize=7)
fig.suptitle("Orthogonal validation: GO-transfer fidelity tracks GO-independent functional conservation")
fig.tight_layout(); fig.savefig(f"{OUTD}/fig_S17_scatter.svg"); fig.savefig(f"{OUTD}/fig_S17_scatter.pdf"); plt.close(fig)

# ---- S17-fig2: both curves vs %identity (GO Wang + STRING conservation), normalised ----
fig,ax=plt.subplots(figsize=(9,5.5))
for t in st.track.unique():
    s=st[st.track==t].sort_values("median_pident",ascending=False)
    if s.empty: continue
    c=clade_color.get(t,"k")
    ax.plot(s.median_pident,s.GO_Wang,"-o",color=c,ms=4,alpha=.85)
    ax.plot(s.median_pident,s.string_conservation,"--s",color=c,ms=3,alpha=.5)
ax.invert_xaxis(); ax.set_xlabel("median ortholog % identity")
ax.set_ylabel("score (solid=GO Wang, dashed=STRING conservation)")
ax.set_title("Both fidelity (GO, solid) and GO-independent conservation (STRING, dashed) decay with distance")
ax.grid(alpha=.3); fig.tight_layout()
fig.savefig(f"{OUTD}/fig_S17_curves.svg"); fig.savefig(f"{OUTD}/fig_S17_curves.pdf"); plt.close(fig)
print(f"\n[S17c] wrote orthogonal_validation.txt + fig_S17_scatter/curves -> {OUTD}")
