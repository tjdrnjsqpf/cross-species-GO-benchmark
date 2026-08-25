#!/usr/bin/env python3
"""REVISION R1-C7 — per-category ortholog mapping rate & amino-acid identity.

Reviewer 1, major comment 7: "examine whether ortholog mapping rate and amino
acid identity differ among these process categories (conserved / taxon-specific / other)".

For every track x reference x synthetic list:
  - n_genes        : genes in the list (from lists_synthetic.tsv)
  - n_mapped       : list genes with an RBH ortholog in that reference
  - frac_mapped    : n_mapped / n_genes
  - med_pident_gene: median %identity of the mapped list genes' ortholog pairs
joined with the list's category (ancestry-based, from lists_meta.tsv).

Outputs (results/crossclade/):
  category_mapping.tsv          granular table (one row per track x ref x list)
  category_mapping_summary.txt  per-category means +- SE and simple tests

Run:  python3 R1C7_category_mapping.py            (expects GOTX_ROOT layout)
      GOTX_ROOT=/path/to/project python3 R1C7_category_mapping.py
"""
import os, glob, statistics, math
from collections import defaultdict

ROOT = os.environ.get("GOTX_ROOT", os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
TRACKS = ["mammal", "fish", "insect", "plant_rice", "plant_arabidopsis", "fungi"]
OUTD = f"{ROOT}/results/crossclade"
os.makedirs(OUTD, exist_ok=True)


def load_map(path):
    """focal_acc -> pident (first/best hit, same convention as S3_transfer.load_map)."""
    pid = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            f, p = c[0], float(c[2])
            if f not in pid:
                pid[f] = p
    return pid


rows = []
for track in TRACKS:
    base = f"{ROOT}/results/{track}"
    meta_p = f"{base}/enrichment/lists_meta.tsv"
    lists_p = f"{base}/enrichment/lists_synthetic.tsv"
    if not (os.path.exists(meta_p) and os.path.exists(lists_p)):
        print(f"[R1C7] {track}: enrichment inputs missing, SKIP ({meta_p})")
        continue

    # list_id -> (category, namespace)
    meta = {}
    with open(meta_p) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(header)}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            meta[c[ci["list_id"]]] = (c[ci["category"]], c[ci["namespace"]])

    # list_id -> set(focal_acc)
    members = defaultdict(set)
    with open(lists_p) as fh:
        next(fh)
        for line in fh:
            lid, acc = line.rstrip("\n").split("\t")[:2]
            members[lid].add(acc)

    maps = sorted(glob.glob(f"{base}/mapping/*__*.rbh.tsv"))
    if not maps:
        print(f"[R1C7] {track}: no rbh maps under {base}/mapping, SKIP")
        continue
    print(f"[R1C7] {track}: {len(members)} lists x {len(maps)} refs")

    for mp in maps:
        ref = os.path.basename(mp).split("__")[1].split(".")[0]
        pid = load_map(mp)
        for lid, genes in members.items():
            if lid not in meta:
                continue
            cat, ns = meta[lid]
            hits = [pid[g] for g in genes if g in pid]
            n, m = len(genes), len(hits)
            rows.append((track, ref, lid, ns, cat, n, m,
                         round(m / n, 4) if n else 0.0,
                         round(statistics.median(hits), 2) if hits else ""))

out_p = f"{OUTD}/category_mapping.tsv"
with open(out_p, "w") as fo:
    fo.write("track\tref\tlist_id\tnamespace\tcategory\tn_genes\tn_mapped\tfrac_mapped\tmed_pident_genes\n")
    for r in rows:
        fo.write("\t".join(map(str, r)) + "\n")
print(f"[R1C7] wrote {out_p} ({len(rows)} rows)")

# ---- summary: per-category means +- SE (BP lists, all refs pooled; and per track) ----
def mean_se(vals):
    vals = [v for v in vals if v != ""]
    if not vals:
        return (float("nan"), float("nan"), 0)
    mu = statistics.fmean(vals)
    se = (statistics.stdev(vals) / math.sqrt(len(vals))) if len(vals) > 1 else 0.0
    return (mu, se, len(vals))

lines = ["=== R1-C7: ortholog mapping rate & %identity by category (synthetic BP lists, RBH) ==="]
bp = [r for r in rows if r[3] == "biological_process"]
for scope, sub in [("pooled", bp)] + [(t, [r for r in bp if r[0] == t]) for t in TRACKS]:
    lines.append(f"\n[{scope}]")
    lines.append(f"{'category':10s} {'n':>6s} {'frac_mapped':>16s} {'med_pident':>16s}")
    for cat in ("conserved", "specific", "other"):
        cs = [r for r in sub if r[4] == cat]
        fm, fs, n1 = mean_se([r[7] for r in cs])
        pm, ps, _ = mean_se([float(r[8]) for r in cs if r[8] != ""])
        if n1:
            lines.append(f"{cat:10s} {n1:>6d} {fm:>10.3f}+-{fs:.3f} {pm:>10.2f}+-{ps:.2f}")
# simple two-sample comparison conserved vs specific (pooled BP)
try:
    from scipy import stats as sst
    a = [r[7] for r in bp if r[4] == "conserved"]
    b = [r[7] for r in bp if r[4] == "specific"]
    ai = [float(r[8]) for r in bp if r[4] == "conserved" and r[8] != ""]
    bi = [float(r[8]) for r in bp if r[4] == "specific" and r[8] != ""]
    if min(len(a), len(b), len(ai), len(bi)) >= 8:
        t1 = sst.mannwhitneyu(a, b)
        t2 = sst.mannwhitneyu(ai, bi)
        lines.append(f"\nMann-Whitney conserved vs specific: frac_mapped p={t1.pvalue:.2e}; med_pident p={t2.pvalue:.2e}")
    else:
        lines.append("\n(too few rows for tests)")
except Exception as e:
    lines.append(f"\n(scipy unavailable, skipped tests: {e})")

summ_p = f"{OUTD}/category_mapping_summary.txt"
open(summ_p, "w").write("\n".join(lines) + "\n")
print(f"[R1C7] wrote {summ_p}")
print("\n".join(lines))
