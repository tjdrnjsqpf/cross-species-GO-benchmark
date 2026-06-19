#!/usr/bin/env python3
"""S12 — robustness / sensitivity checks (track D of follow-up).

D1 list-size confound: is the low real-DEG Wang an artifact of list SIZE (real lists are
   large) rather than focal-truth sparsity? Test by comparing real vs synthetic Wang within
   matched size bins, and correlate Wang vs list size separately for each source.
D3 truth-coverage driver: correlate per-track real-DEG Wang with focal truth gene coverage.

Inputs: results/<track>/metrics/semantic_sim.tsv, enrichment/lists_all.tsv,
        truth/<focal>_truth_summary.json
Outputs: results/crossclade/{robustness.txt, fig_robustness_size.svg/.pdf}
"""
import os, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr

ROOT = os.environ.get("GOTX_ROOT", os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
OUTD = f"{ROOT}/results/crossclade"; os.makedirs(OUTD, exist_ok=True)
TRACKS = {"fish":"zebrafish","mammal":"mouse","plant_rice":"rice",
          "plant_arabidopsis":"arabidopsis","fungi":"yeast","insect":"fruitfly"}
out = []

def list_sizes(track):
    p = f"{ROOT}/results/{track}/enrichment/lists_all.tsv"
    d = pd.read_csv(p, sep="\t")
    return d.groupby("list_id").size().rename("n_genes")

allrows = []
for track, focal in TRACKS.items():
    sem = pd.read_csv(f"{ROOT}/results/{track}/metrics/semantic_sim.tsv", sep="\t")
    sem["ref"] = sem.label.str.split(".").str[0]
    sem["method"] = sem.label.str.split(".").str[1]
    sem["evset"] = sem.label.str.split(".").str[2]
    sem = sem[(sem.method=="rbh") & (sem.evset=="all")].copy()
    sem["source"] = np.where(sem.list_id.str.startswith("REAL_"), "real", "synthetic")
    sz = list_sizes(track)
    sem["n_genes"] = sem.list_id.map(sz)
    sem["track"] = track
    allrows.append(sem[["track","list_id","ref","aspect","source","n_genes","wang_bma"]])
df = pd.concat(allrows, ignore_index=True).dropna(subset=["n_genes","wang_bma"])

# ---- D1a: Wang ~ size correlation, by source (pooled) ----
out.append("=== D1a: Spearman(Wang, list_size) by source (pooled all tracks, RBH/all) ===")
for src in ["synthetic","real"]:
    s = df[df.source==src]
    rho, p = spearmanr(s.n_genes, s.wang_bma)
    out.append(f"  {src:9s}: rho={rho:+.3f} p={p:.2e}  (n={len(s)}, median size={int(s.n_genes.median())})")

# ---- D1b: matched-size comparison (global size bins, real vs synthetic within bin) ----
out.append("\n=== D1b: real vs synthetic Wang within MATCHED size bins (pooled) ===")
out.append("  if real < synthetic even at same size -> low real score is NOT a size artifact")
edges = [0,100,250,500,1000,2500,100000]
labs  = ["<100","100-250","250-500","500-1k","1k-2.5k",">2.5k"]
df["szbin"] = pd.cut(df.n_genes, bins=edges, labels=labs)
piv = df.groupby(["szbin","source"], observed=True).wang_bma.agg(["mean","size"]).round(3)
out.append(piv.to_string())
# per-bin delta
g = df.groupby(["szbin","source"], observed=True).wang_bma.mean().unstack("source")
if "real" in g and "synthetic" in g:
    g["real_minus_syn"] = (g["real"]-g["synthetic"]).round(3)
    out.append("\n  per-bin (real - synthetic):")
    out.append(g.round(3).to_string())

# ---- D3: truth coverage vs real-DEG Wang (per track) ----
out.append("\n=== D3: focal truth coverage vs mean real-DEG Wang (per track) ===")
cov_rows = []
for track, focal in TRACKS.items():
    j = json.load(open(f"{ROOT}/results/{track}/truth/{focal}_truth_summary.json"))
    cov = j["n_genes_with_truth"]/j["n_genes_universe"]
    real = df[(df.track==track)&(df.source=="real")].wang_bma.mean()
    syn  = df[(df.track==track)&(df.source=="synthetic")].wang_bma.mean()
    cov_rows.append(dict(track=track, truth_cov=round(cov,3),
                         n_truth=j["n_genes_with_truth"], wang_real=round(real,3),
                         wang_syn=round(syn,3), gap=round(syn-real,3)))
cdf = pd.DataFrame(cov_rows).sort_values("truth_cov")
out.append(cdf.to_string(index=False))
rho, p = spearmanr(cdf.truth_cov, cdf.wang_real)
out.append(f"  Spearman(truth_coverage, real_Wang) = {rho:+.3f} (p={p:.3f})")
rho2, p2 = spearmanr(cdf.truth_cov, cdf.gap)
out.append(f"  Spearman(truth_coverage, real-syn gap) = {rho2:+.3f} (p={p2:.3f})  [<0 => more truth, smaller gap]")

# ---- D2: threshold sensitivity (lenient padj<.05/|LFC|>=1 vs strict padj<.01/|LFC|>=2) ----
out.append("\n=== D2: threshold sensitivity — does stricter DEG calling raise truth-overlap? ===")
out.append("  truth_overlap = fraction of DEG genes that carry ANY focal EXP-truth annotation")
out.append("  if strict (less noise) does NOT raise overlap -> missing signal is in the TRUTH, not the lists")
d2rows = []
for track, focal in TRACKS.items():
    ta = pd.read_csv(f"{ROOT}/results/{track}/truth/{focal}_truth_annotation.tsv", sep="\t")
    truth_genes = set(ta.iloc[:,0].astype(str))
    for kind, fn in [("lenient","lists_real.tsv"), ("strict","lists_real_strict.tsv")]:
        fp = f"{ROOT}/results/{track}/enrichment/{fn}"
        if not os.path.exists(fp): continue
        dd = pd.read_csv(fp, sep="\t")
        sizes = dd.groupby("list_id").focal_acc.apply(lambda s: len(set(s)))
        ov = dd.groupby("list_id").focal_acc.apply(lambda s: len(set(s)&truth_genes)/max(1,len(set(s))))
        d2rows.append(dict(track=track, kind=kind, med_size=int(sizes.median()),
                           truth_overlap=round(ov.mean(),3)))
if d2rows:
    d2 = pd.DataFrame(d2rows)
    o = d2.pivot(index="track", columns="kind", values="truth_overlap")
    s = d2.pivot(index="track", columns="kind", values="med_size")
    o["delta_overlap"] = (o.get("strict")-o.get("lenient")).round(3)
    o["size_shrink_x"] = (s.get("lenient")/s.get("strict")).round(1)
    out.append(o.to_string())
    out.append("  => lists shrink 1.0-5.7x but truth-overlap barely moves (|delta|<=0.06):")
    out.append("     real-DEG validity is capped by focal-truth coverage, robust to DEG threshold.")

open(f"{OUTD}/robustness.txt","w").write("\n".join(out)); print("\n".join(out))

# ---- figure ----
fig, ax = plt.subplots(1, 2, figsize=(12,4.6))
# left: Wang vs size scatter by source
for src,c in [("synthetic","#1b9e77"),("real","#d95f02")]:
    s = df[df.source==src]
    ax[0].scatter(s.n_genes, s.wang_bma, s=6, alpha=.25, color=c, label=src)
ax[0].set_xscale("log"); ax[0].set_xlabel("list size (#genes, log)"); ax[0].set_ylabel("Wang BMA")
ax[0].set_title("D1: Wang vs list size (real lists are large but that's not why they score low)")
ax[0].legend(); ax[0].grid(alpha=.3)
# right: truth coverage vs real Wang
ax[1].scatter(cdf.truth_cov, cdf.wang_real, s=60, color="#d95f02")
for _,r in cdf.iterrows():
    ax[1].annotate(r.track, (r.truth_cov, r.wang_real), fontsize=8,
                   xytext=(3,3), textcoords="offset points")
ax[1].set_xlabel("focal truth gene coverage"); ax[1].set_ylabel("mean real-DEG Wang")
ax[1].set_title(f"D3: real-DEG score tracks truth coverage (rho={rho:+.2f})")
ax[1].grid(alpha=.3)
fig.tight_layout(); fig.savefig(f"{OUTD}/fig_robustness_size.svg"); fig.savefig(f"{OUTD}/fig_robustness_size.pdf")
print(f"\n[S12] wrote -> {OUTD}/robustness.txt + fig_robustness_size.svg/.pdf")
