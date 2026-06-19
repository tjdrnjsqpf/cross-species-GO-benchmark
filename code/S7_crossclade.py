#!/usr/bin/env python3
"""S7 — cross-clade comparison (the 3-track deliverable).
Overlays fish/mammal/plant tolerance curves normalized on x=%identity, and tests the
plant-WGD hypothesis via a clade x distance interaction regression.
Inputs: each track's results/<track>/metrics/{semantic_sim,setlevel_metrics}.tsv
Outputs: results/crossclade/{crossclade.tsv, fig_crossclade_semantic.png/pdf,
         fig_crossclade_recall.png, regression_clade.txt}
"""
import os
import pandas as pd, numpy as np
import statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.environ.get("GOTX_ROOT", os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
TRACKS = {"fish":"fish", "mammal":"mammal", "plant_rice":"plant", "plant_arabidopsis":"plant",
          "fungi":"fungi", "insect":"insect"}
OUTD = f"{ROOT}/results/crossclade"; os.makedirs(OUTD, exist_ok=True)

rows = []
for track, clade in TRACKS.items():
    base = f"{ROOT}/results/{track}/metrics"
    sem = pd.read_csv(f"{base}/semantic_sim.tsv", sep="\t")
    m = pd.read_csv(f"{base}/setlevel_metrics.tsv", sep="\t")
    pid = m[m.ic_bin=="all"].groupby("ref").median_pident.first()
    rich = m[m.ic_bin=="all"].groupby("ref").ref_richness.first()
    My = m[m.ic_bin=="all"].groupby("ref").My.first()
    sem = sem[sem.list_id.str.startswith("SYN_")].copy()   # controlled synthetic only
    sem["ref"] = sem.label.str.split(".").str[0]
    sem["method"] = sem.label.str.split(".").str[1]
    sem["evset"] = sem.label.str.split(".").str[2]
    sem = sem[(sem.method=="rbh") & (sem.evset=="all")]
    sem["track"] = track; sem["clade"] = clade
    sem["median_pident"] = sem.ref.map(pid)
    sem["My"] = sem.ref.map(My)
    # recall from setlevel
    rec = m[(m.ic_bin=="all")&(m.method=="rbh")&(m.evset=="all")][["list_id","aspect","ref","recall"]]
    sem = sem.merge(rec, on=["list_id","aspect","ref"], how="left")
    rows.append(sem)
df = pd.concat(rows, ignore_index=True)
df.to_csv(f"{OUTD}/crossclade.tsv", sep="\t", index=False)

# ---- aggregate per track/ref/aspect, with bootstrap 95% CI over synthetic lists (FIX 2) ----
RNG = np.random.default_rng(20260616); NBOOT = 1000
def boot_ci(vals):
    vals = np.asarray(vals, float); vals = vals[~np.isnan(vals)]
    if len(vals) < 3: return (np.nan, np.nan)
    means = vals[RNG.integers(0, len(vals), size=(NBOOT, len(vals)))].mean(axis=1)
    return tuple(np.percentile(means, [2.5, 97.5]))

arows=[]
for (track,clade,ref,aspect,pidv),g in df.groupby(["track","clade","ref","aspect","median_pident"]):
    wlo,whi = boot_ci(g.wang_bma.values); rlo,rhi = boot_ci(g.recall.values)
    arows.append(dict(track=track,clade=clade,ref=ref,aspect=aspect,median_pident=pidv,
                      wang=g.wang_bma.mean(), recall=g.recall.mean(), n=len(g),
                      wang_lo=wlo,wang_hi=whi,recall_lo=rlo,recall_hi=rhi))
agg = pd.DataFrame(arows)
agg.to_csv(f"{OUTD}/crossclade_ci.tsv", sep="\t", index=False)

clade_color = {"fish":"#1b9e77","mammal":"#d95f02","plant":"#7570b3","fungi":"#e7298a","insect":"#66a61e"}
track_marker = {"fish":"o","mammal":"s","plant_rice":"^","plant_arabidopsis":"D","fungi":"P","insect":"X"}

def plot(metric, fname, ylab):
    lo_c, hi_c = f"{metric}_lo", f"{metric}_hi"
    fig, axes = plt.subplots(1, 3, figsize=(14,4.6), sharey=True)
    for ax, asp in zip(axes, ["MF","BP","CC"]):
        sub = agg[agg.aspect==asp]
        for track, clade in TRACKS.items():
            t = sub[sub.track==track].sort_values("median_pident", ascending=False)
            if t.empty: continue
            ax.fill_between(t.median_pident, t[lo_c], t[hi_c], color=clade_color[clade], alpha=.12)
            ax.plot(t.median_pident, t[metric], marker=track_marker[track],
                    color=clade_color[clade], label=track, alpha=.85, ms=5)
        ax.set_xlabel("median ortholog % identity"); ax.set_title(asp)
        ax.invert_xaxis(); ax.grid(alpha=.3)
    axes[0].set_ylabel(ylab)
    h,l = axes[1].get_legend_handles_labels()
    axes[2].legend(h,l,fontsize=8,title="track")
    fig.suptitle(f"Cross-clade divergence tolerance ({ylab}, RBH, all-evidence, synthetic)  — x=% identity, shaded=95% bootstrap CI")
    fig.tight_layout()
    fig.savefig(f"{OUTD}/{fname}.svg"); fig.savefig(f"{OUTD}/{fname}.pdf"); plt.close(fig)

plot("wang","fig_crossclade_semantic","Wang semantic similarity")
plot("recall","fig_crossclade_recall","recall")

# ---- WGD hypothesis: clade x distance interaction ----
d = df.dropna(subset=["wang_bma","median_pident"]).copy()
d["pid_z"] = (d.median_pident - d.median_pident.mean())/d.median_pident.std()
out = []
for metric in ["wang_bma","recall"]:
    dd = d.dropna(subset=[metric])
    mod = smf.ols(f"{metric} ~ pid_z * C(clade, Treatment('fish')) + C(aspect)", data=dd).fit()
    out.append(f"\n===== {metric} ~ pid_z * clade + aspect  (n={len(dd)}, R2={mod.rsquared:.3f}) =====")
    out.append(mod.summary().as_text())
    # interaction terms = does the slope vs identity differ by clade? (steeper = larger positive pid_z slope)
    inter = {k:v for k,v in mod.params.items() if "pid_z:" in k}
    out.append("# distance-slope interactions (vs fish baseline): " + ", ".join(f"{k.split('T.')[-1].rstrip(']')}={v:+.3f}" for k,v in inter.items()))
open(f"{OUTD}/regression_clade.txt","w").write("\n".join(out))

# ---- rice <-> arabidopsis concordance (FIX 2): two independent plant focals agree? ----
con=[]
for asp in ["MF","BP","CC"]:
    r=agg[(agg.track=="plant_rice")&(agg.aspect==asp)].sort_values("median_pident")
    a=agg[(agg.track=="plant_arabidopsis")&(agg.aspect==asp)].sort_values("median_pident")
    if len(r)<3 or len(a)<3: continue
    # compare on overlapping %id range via interpolation onto a common grid
    lo=max(r.median_pident.min(),a.median_pident.min()); hi=min(r.median_pident.max(),a.median_pident.max())
    grid=np.linspace(lo,hi,20)
    ri=np.interp(grid,r.median_pident,r.wang); ai=np.interp(grid,a.median_pident,a.wang)
    dw=np.abs(ri-ai)
    rho=np.corrcoef(ri,ai)[0,1]
    con.append(f"  {asp}: overlap {lo:.0f}-{hi:.0f}%id  mean|dWang|={dw.mean():.3f}  max|dWang|={dw.max():.3f}  Pearson r={rho:.3f}")
conc="\n".join(["=== rice <-> arabidopsis concordance (independent plant focals, matched %identity) ==="]+con)
open(f"{OUTD}/plant_concordance.txt","w").write(conc); print(conc); print()

print("=== cross-clade tolerance (Wang, BP) at matched identity ===")
piv = agg[agg.aspect=="BP"].pivot_table(index="track", columns=None, values="wang")
print(agg[agg.aspect=="BP"][["track","ref","median_pident","wang","recall"]].sort_values(["track","median_pident"],ascending=[True,False]).to_string(index=False))
print("\n[S7] wrote figures + regression_clade.txt ->", OUTD)
print("\n".join(out[-1:] if out else []))
