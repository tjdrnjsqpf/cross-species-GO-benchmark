#!/usr/bin/env python3
"""REVISION R1-C1 — per-reference ortholog %identity distribution (IQR).

For the panel-overview figure (new Supp Fig S1) we want error bars on each
reference's identity: the interquartile range (and 5-95 percentile) of the
per-ortholog %identity values behind the reported median.

Input : results/<track>/mapping/<focal>__<ref>.rbh.tsv  (focal_acc, ref_acc, pident)
Output: results/crossclade/identity_iqr.tsv
        (track, ref, n_orthologs, p5, q25, median, q75, p95)

Run:  GOTX_ROOT=$PWD python3 SERVER_rev_scripts/R1C1_identity_iqr.py
"""
import os, glob, statistics

ROOT = os.environ.get("GOTX_ROOT", os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
TRACKS = ["mammal", "fish", "insect", "plant_rice", "plant_arabidopsis", "fungi"]
OUTD = f"{ROOT}/results/crossclade"
os.makedirs(OUTD, exist_ok=True)


def quantile(sorted_vals, q):
    n = len(sorted_vals)
    if n == 1:
        return sorted_vals[0]
    pos = q * (n - 1)
    lo = int(pos)
    hi = min(lo + 1, n - 1)
    return sorted_vals[lo] + (pos - lo) * (sorted_vals[hi] - sorted_vals[lo])


rows = []
for track in TRACKS:
    maps = sorted(glob.glob(f"{ROOT}/results/{track}/mapping/*__*.rbh.tsv"))
    if not maps:
        print(f"[R1C1] {track}: no rbh maps, SKIP")
        continue
    for mp in maps:
        ref = os.path.basename(mp).split("__")[1].split(".")[0]
        pids = []
        seen = set()
        with open(mp) as fh:
            next(fh)
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < 3:
                    continue
                if c[0] in seen:          # first/best hit per focal gene, as in S3
                    continue
                seen.add(c[0])
                pids.append(float(c[2]))
        if not pids:
            continue
        pids.sort()
        rows.append((track, ref, len(pids),
                     round(quantile(pids, 0.05), 1), round(quantile(pids, 0.25), 1),
                     round(statistics.median(pids), 1),
                     round(quantile(pids, 0.75), 1), round(quantile(pids, 0.95), 1)))
    print(f"[R1C1] {track}: {sum(1 for r in rows if r[0]==track)} refs")

out = f"{OUTD}/identity_iqr.tsv"
with open(out, "w") as fo:
    fo.write("track\tref\tn_orthologs\tp5\tq25\tmedian\tq75\tp95\n")
    for r in rows:
        fo.write("\t".join(map(str, r)) + "\n")
print(f"[R1C1] wrote {out} ({len(rows)} rows)")
