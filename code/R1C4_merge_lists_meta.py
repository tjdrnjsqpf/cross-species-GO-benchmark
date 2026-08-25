#!/usr/bin/env python3
"""REVISION R1-C4 — merge the six lists_meta.tsv into one file (adds a track column).

For the new Supplementary Table S4 (gene-list summary) we need the authoritative
per-list metadata (ic_bin tertiles, category, namespace, n_genes).

Run:  GOTX_ROOT=$PWD python3 SERVER_rev_scripts/R1C4_merge_lists_meta.py
Output: results/crossclade/lists_meta_all.tsv  -> return via to_Desktop/revision_R1C4/
"""
import os

ROOT = os.environ.get("GOTX_ROOT", os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
TRACKS = ["mammal", "fish", "insect", "plant_rice", "plant_arabidopsis", "fungi"]
OUTD = f"{ROOT}/results/crossclade"
os.makedirs(OUTD, exist_ok=True)

out = f"{OUTD}/lists_meta_all.tsv"
n = 0
with open(out, "w") as fo:
    header_written = False
    for t in TRACKS:
        p = f"{ROOT}/results/{t}/enrichment/lists_meta.tsv"
        if not os.path.exists(p):
            print(f"[R1C4] {t}: missing {p}, SKIP")
            continue
        with open(p) as fh:
            hdr = fh.readline().rstrip("\n")
            if not header_written:
                fo.write("track\t" + hdr + "\n")
                header_written = True
            for line in fh:
                fo.write(t + "\t" + line)
                n += 1
        print(f"[R1C4] {t}: ok")
print(f"[R1C4] wrote {out} ({n} lists)")
