#!/usr/bin/env python3
"""Revision R1-1 — new Supplementary Fig S1: reference-panel overview.
For each focal species, all references arranged along the divergence axis
(median one-to-one ortholog %identity), with divergence time (My) annotated
and the equidistant richness-contrast group marked.
Input : ../../DATA_DEPOSIT/data/figdata_refmeta.tsv
Output: FigS1_panel_overview.pdf / .png
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.realpath(__file__))
SRC = os.path.join(HERE, "../../data/figdata_refmeta.tsv")
FIGD = os.path.join(HERE, "../../data/revision")

PALETTE = {"mammal": "#D95F02", "fish": "#1B9E77", "insect": "#B8860B",
           "plant_rice": "#E7298A", "plant_arabidopsis": "#3B6CB7", "fungi": "#7570B3"}
ORDER = ["mammal", "fish", "insect", "plant_rice", "plant_arabidopsis", "fungi"]
TITLE = {"mammal": "Mammal — focal: mouse ($\\it{Mus\\ musculus}$)",
         "fish": "Fish — focal: zebrafish ($\\it{Danio\\ rerio}$)",
         "insect": "Insect — focal: fruit fly ($\\it{Drosophila\\ melanogaster}$)",
         "plant_rice": "Plant-R — focal: rice ($\\it{Oryza\\ sativa}$)",
         "plant_arabidopsis": "Plant-A — focal: thale cress ($\\it{Arabidopsis\\ thaliana}$)",
         "fungi": "Fungi — focal: budding yeast ($\\it{Saccharomyces\\ cerevisiae}$)"}

NAMES_TSV = os.path.join(HERE, "../../data/revision/species_names.tsv")
SCI = {}
with open(NAMES_TSV) as _f:
    next(_f)
    for _line in _f:
        _k, _sci, _com = _line.rstrip("\n").split("\t")
        SCI[_k] = _sci

def pretty(ref):
    """Italic scientific name via mathtext (authoritative mapping, species_names.tsv)."""
    sci = SCI[ref]
    return r"$\it{" + sci.replace(" ", r"\ ") + "}$"

rows = list(csv.DictReader(open(SRC), delimiter="\t"))
by_track = {t: sorted([r for r in rows if r["track"] == t],
                      key=lambda r: -float(r["median_pident"])) for t in ORDER}

# optional: per-reference ortholog-identity distribution (server job R1C1)
IQR = {}
for cand in (os.path.join(HERE, "../../data/revision/identity_iqr.tsv"),):
    if os.path.exists(cand):
        for r in csv.DictReader(open(cand), delimiter="\t"):
            IQR[(r["track"], r["ref"])] = (float(r["p5"]), float(r["q25"]),
                                           float(r["q75"]), float(r["p95"]))
        print(f"loaded identity IQR for {len(IQR)} refs from {cand}")
        break
if not IQR:
    print("identity_iqr.tsv not found - drawing medians only (run server job R1C1 for error bars)")

fig, axes = plt.subplots(3, 2, figsize=(13.2, 15.5))
for ax, track in zip(axes.flat, ORDER):
    refs = by_track[track]
    col = PALETTE[track]
    n = len(refs)
    ys = list(range(n, 0, -1))                       # closest at top
    labels = []
    for r, y in zip(refs, ys):
        x = float(r["median_pident"])
        eqd = bool(r["equidistant_group"])
        ax.hlines(y, x, 101, color="0.88", lw=0.8, zorder=1)
        iqr = IQR.get((track, r["ref"]))
        if iqr:
            p5, q25, q75, p95 = iqr
            ax.hlines(y, p5, p95, color=col, lw=0.9, alpha=0.45, zorder=2)
            ax.hlines(y, q25, q75, color=col, lw=3.2, alpha=0.75, zorder=2)
        if eqd:
            ax.scatter([x], [y], s=52, marker="D", facecolor=col,
                       edgecolor="black", linewidth=1.1, zorder=3)
        else:
            ax.scatter([x], [y], s=42, marker="o", facecolor=col,
                       edgecolor="none", zorder=3)
        my = r["My"]
        my_txt = f"{int(float(my))} My" if my not in ("", "NA") else ""
        labels.append(f"{pretty(r['ref'])}  ({my_txt})" if my_txt else pretty(r["ref"]))
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=7.4)
    xlo = 33
    if IQR:
        p5s = [IQR[(track, r["ref"])][0] for r in refs if (track, r["ref"]) in IQR]
        if p5s:
            xlo = min(33, min(p5s) - 2)
    ax.set_xlim(101, xlo)                            # reversed: close -> distant
    ax.set_ylim(0.3, n + 0.9)
    ax.axvline(100, color="0.4", lw=0.8, ls=":")
    ax.text(100, n + 0.72, "Focal", ha="center", va="bottom", fontsize=7, color="0.35")
    ax.set_title(TITLE[track], fontsize=9.3, color="black", loc="left", pad=8)
    ax.set_xlabel("Median ortholog identity (%)", fontsize=8.5)
    ax.tick_params(axis="x", labelsize=8)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)

handles = [Line2D([], [], marker="o", ls="none", mfc="0.45", mec="none", ms=7,
                  label="Reference species (median ortholog identity)"),
           Line2D([], [], marker="D", ls="none", mfc="0.45", mec="black", ms=7,
                  label="Equidistant richness-contrast group")]
if IQR:
    handles += [Line2D([], [], color="0.45", lw=3.2, alpha=0.75,
                       label="Interquartile range of ortholog identities"),
                Line2D([], [], color="0.45", lw=0.9, alpha=0.45,
                       label="5th\u201395th percentile")]
fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=9,
           bbox_to_anchor=(0.5, 0.0))
fig.tight_layout(rect=(0, 0.035, 1, 0.995))
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIGD, f"FigS1_panel_overview.{ext}"), dpi=300)
print("wrote FigS1_panel_overview.pdf/.png ;",
      sum(len(v) for v in by_track.values()), "references,",
      sum(1 for r in rows if r["equidistant_group"]), "equidistant-flagged")
