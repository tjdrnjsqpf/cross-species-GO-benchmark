#!/usr/bin/env python3
"""Revision R1-9 — Fig 1B with 95% bootstrap CIs on ID50 and floor.

Replicates the S10 estimator exactly (constrained logistic / interpolation
fallback, per id50.tsv method_used) and bootstraps the *reported* estimator by
resampling the reference panel (1,000 resamples, same as the pipeline).

Inputs : DATA_DEPOSIT/data/crossclade_ci.tsv  (per-reference mean Wang; fit input)
         DATA_DEPOSIT/data/id50.tsv           (point estimates + method_used)
Outputs: fig1b_ci.tsv, Fig1B_ID50_floor_CI.pdf/.png
"""
import os
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.realpath(__file__))
DATA = os.path.join(HERE, "../../data")
FIGD = os.path.join(HERE, "../../data/revision")
TABD = os.path.join(HERE, "../../data/revision")
NBOOT = 1000
K_UPPER = 5.0
RNG = np.random.default_rng(20260825)

def logistic(x, lo, hi, id50, k):
    return lo + (hi - lo) / (1 + np.exp(-k * (x - id50)))

def interp_id50(x, y):
    o = np.argsort(x); x, y = np.array(x)[o], np.array(y)[o]
    hi = np.mean(np.sort(y)[-2:]); lo = np.min(y); tgt = (hi + lo) / 2
    for i in range(len(x) - 1):
        if (y[i] - tgt) * (y[i + 1] - tgt) <= 0 and y[i + 1] != y[i]:
            return x[i] + (tgt - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i]), lo, hi
    return np.nan, lo, hi

def fit_logistic(x, y):
    lo0, hi0 = float(np.min(y)), float(np.mean(np.sort(y)[-2:]))
    xmn, xmx = float(np.min(x)), float(np.max(x))
    bounds = ([0.0, lo0 - 1e-6, xmn, 1e-3], [max(hi0, lo0) + 1e-6, 1.0, xmx, K_UPPER])
    i0 = interp_id50(x, y)[0]
    p0 = [lo0, hi0, np.clip(i0 if not np.isnan(i0) else np.median(x), xmn, xmx), 0.2]
    try:
        popt, _ = curve_fit(logistic, x, y, p0=p0, bounds=bounds, maxfev=30000)
    except Exception:
        return None
    return popt  # lo, hi, id50, k

cc = pd.read_csv(f"{DATA}/crossclade_ci.tsv", sep="\t")
id50 = pd.read_csv(f"{DATA}/id50.tsv", sep="\t")

rows = []
for _, r in id50.iterrows():
    sub = cc[(cc.track == r.track) & (cc.aspect == r.aspect)]
    x = sub.median_pident.values.astype(float)
    y = sub.wang.values.astype(float)
    n = len(x)

    def estimate(xx, yy):
        if r.method_used == "logistic":
            f = fit_logistic(xx, yy)
            return (f[2], f[0]) if f is not None else (np.nan, np.nan)
        i, lo, _ = interp_id50(xx, yy)
        return (i, lo)

    # sanity: recompute point estimate on full panel
    id50_pt, floor_pt = estimate(x, y)
    boots_id, boots_fl = [], []
    while len(boots_fl) < NBOOT:
        idx = RNG.integers(0, n, n)
        if len(np.unique(x[idx])) < 4:
            continue
        bi, bf = estimate(x[idx], y[idx])
        if not np.isnan(bf):
            boots_fl.append(bf)
            if not np.isnan(bi):
                boots_id.append(bi)
    qi = np.percentile(boots_id, [2.5, 97.5]) if boots_id else (np.nan, np.nan)
    qf = np.percentile(boots_fl, [2.5, 97.5])
    rows.append(dict(track=r.track, aspect=r.aspect, n=n, method_used=r.method_used,
                     ID50_rep=r.ID50_rep, id50_recomputed=round(id50_pt, 1) if not np.isnan(id50_pt) else np.nan,
                     id50_lo=round(qi[0], 1), id50_hi=round(qi[1], 1),
                     floor=r.floor, floor_recomputed=round(floor_pt, 3),
                     floor_lo=round(qf[0], 3), floor_hi=round(qf[1], 3)))

out = pd.DataFrame(rows)
out.to_csv(os.path.join(TABD, "fig1b_ci.tsv"), sep="\t", index=False)
print(out.to_string(index=False))

# consistency check vs pipeline
chk = out.merge(id50[["track", "aspect", "ID50_rep", "floor"]], on=["track", "aspect"], suffixes=("", "_pipe"))
bad = chk[(abs(chk.id50_recomputed - chk.ID50_rep) > 1.5) | (abs(chk.floor_recomputed - chk.floor) > 0.02)]
print(f"\n[check] point-estimate reproduction: {len(chk)-len(bad)}/{len(chk)} within tolerance"
      + ("" if bad.empty else f"\n{bad[['track','aspect','ID50_rep','id50_recomputed','floor','floor_recomputed']]}"))

# ---------------- figure: grouped bars with CI, proof-style ----------------
ORDER = ["mammal", "fish", "insect", "plant_rice", "plant_arabidopsis", "fungi"]
LABEL = {"mammal": "Mammal", "fish": "Fish", "insect": "Insect",
         "plant_rice": "Plant-R", "plant_arabidopsis": "Plant-A", "fungi": "Fungi"}
ASPECTS = ["BP", "MF", "CC"]
ACOL = {"BP": "#4D4D4D", "MF": "#9B72B0", "CC": "#DBA827"}

fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.9))
for ax, (val, lo, hi, ylab, ymax) in zip(axes, [
        ("ID50_rep", "id50_lo", "id50_hi", "ID50 (identity %)", 80),
        ("floor", "floor_lo", "floor_hi", "Floor (Wang)", 0.85)]):
    for j, asp in enumerate(ASPECTS):
        sub = out.set_index(["track", "aspect"])
        xs = np.arange(len(ORDER)) + (j - 1) * 0.26
        vals = [sub.loc[(t, asp), val] for t in ORDER]
        los = [sub.loc[(t, asp), lo] for t in ORDER]
        his = [sub.loc[(t, asp), hi] for t in ORDER]
        err = [np.array(vals) - np.array(los), np.array(his) - np.array(vals)]
        ax.bar(xs, vals, width=0.24, color=ACOL[asp], label=asp, zorder=2)
        ax.errorbar(xs, vals, yerr=err, fmt="none", ecolor="black",
                    elinewidth=0.9, capsize=2.2, zorder=3)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([LABEL[t] for t in ORDER], fontsize=9)
    ax.set_ylabel(ylab, fontsize=10)
    ax.set_ylim(0, ymax)
    ax.grid(axis="y", color="0.9", lw=0.7, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
axes[1].legend(frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(1.01, 1.0))
fig.suptitle("Fig 1B with 95% bootstrap CIs (1,000 resamples of the reference panel)",
             fontsize=10, y=1.0)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIGD, f"Fig1B_ID50_floor_CI.{ext}"), dpi=300, bbox_inches="tight")
print("\nwrote fig1b_ci.tsv + Fig1B_ID50_floor_CI.pdf/.png")
