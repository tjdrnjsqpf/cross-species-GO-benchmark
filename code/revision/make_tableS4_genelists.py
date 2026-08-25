#!/usr/bin/env python3
"""Revision R1-4 — new Supplementary Table S4: gene-list summary.

Part A (per focal species): synthetic-list counts by aspect, size distribution,
IC-tertile composition (authoritative ic_bin from server lists_meta when
available), BP conserved/specific/other counts, large-list and real-DEG counts.
Part B (per Expression Atlas experiment): accession, contrasts, list sizes.

Inputs : DATA_DEPOSIT/data/figdata_setlevel.tsv, category.tsv
         to_Desktop/revision_R1C4/lists_meta_all.tsv   (optional; ic_bin source)
Outputs: tables/tableS4_genelists_summary.tsv, tables/tableS4_realDEG_datasets.tsv
"""
import os, re
import pandas as pd

HERE = os.path.dirname(os.path.realpath(__file__))
ROOT = os.path.join(HERE, "../..")
TABD = os.path.join(HERE, "../../data/revision")
LAB = {"mammal": "Mammal", "fish": "Fish", "insect": "Insect",
       "plant_rice": "Plant-R", "plant_arabidopsis": "Plant-A", "fungi": "Fungi"}
ORDER = list(LAB)

sl = pd.read_csv(f"{ROOT}/data/figdata_setlevel.tsv", sep="\t")
cat = pd.read_csv(f"{ROOT}/data/category.tsv", sep="\t")

# per-list frame (synthetic small/large + real), from setlevel
lists = (sl.groupby(["track", "list_id"])
           .agg(size=("size", "first"), size_class=("size_class", "first"),
                list_source=("list_source", "first")).reset_index())
meta_local = (cat.groupby(["track", "list_id"])
                .agg(category=("category", "first"), namespace=("namespace", "first"))
                .reset_index())
lists = lists.merge(meta_local, on=["track", "list_id"], how="left")

# authoritative ic_bin (server R1C4), with consistency check against local category
icb = None
p = f"{ROOT}/data/revision/lists_meta_all.tsv"
if os.path.exists(p):
    lm = pd.read_csv(p, sep="\t")
    icb = lm[["track", "list_id", "ic_bin", "category", "namespace", "n_genes"]]
    m = lists.merge(icb, on=["track", "list_id"], how="left", suffixes=("", "_srv"))
    chk = m.dropna(subset=["category_srv"])
    mism = ((chk.category != chk.category_srv) | (chk.namespace != chk.namespace_srv)).sum()
    print(f"lists_meta_all loaded: {len(lm)} lists; category/namespace mismatch vs local: {mism}")
    lists = m
else:
    lists["ic_bin"] = None
    print("lists_meta_all.tsv not found - IC-tertile columns marked 'pending' (run server job R1C4)")

NS = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}
lists["aspect"] = lists.namespace.map(NS)

rows = []
for t in ORDER:
    g = lists[lists.track == t]
    syn = g[g.size_class == "small"]
    lg = g[g.size_class == "large"]
    rl = g[g.list_source == "real"]
    r = {"Focal": LAB[t], "Synthetic lists": len(syn)}
    for a in ("BP", "MF", "CC"):
        r[f"Synthetic {a}"] = (syn.aspect == a).sum()
    r["Synthetic size (min-max, median)"] = (
        f"{int(syn['size'].min())}-{int(syn['size'].max())} ({int(syn['size'].median())})")
    for b in ("shallow", "mid", "deep"):
        r[f"IC {b}"] = int((syn.ic_bin == b).sum()) if icb is not None else "pending"
    for c, nm in (("conserved", "BP conserved"), ("specific", "BP taxon-specific"), ("other", "BP other")):
        r[nm] = ((syn.aspect == "BP") & (syn.category == c)).sum()
    r["Large synthetic lists"] = len(lg)
    r["Large size (min-max, median)"] = (
        f"{int(lg['size'].min())}-{int(lg['size'].max())} ({int(lg['size'].median())})")
    r["Real DEG lists"] = len(rl)
    r["Real size (min-max, median)"] = (
        f"{int(rl['size'].min())}-{int(rl['size'].max())} ({int(rl['size'].median())})")
    rows.append(r)
A = pd.DataFrame(rows)
A.to_csv(os.path.join(TABD, "tableS4_genelists_summary.tsv"), sep="\t", index=False)
print(A.to_string(index=False))

# Part B — real-DEG datasets
rl = lists[lists.list_source == "real"].copy()
rl["accession"] = rl.list_id.str.extract(r"REAL_(E-[A-Z]+-\d+)")
B = (rl.groupby(["track", "accession"])
       .agg(n_contrasts=("list_id", "nunique"),
            min_size=("size", "min"), max_size=("size", "max")).reset_index())
B["Focal"] = B.track.map(LAB)
B = B[["Focal", "accession", "n_contrasts", "min_size", "max_size"]]
B.to_csv(os.path.join(TABD, "tableS4_realDEG_datasets.tsv"), sep="\t", index=False)
print(f"\nPart B: {len(B)} Expression Atlas experiments, {B.n_contrasts.sum()} contrasts")
print(B.to_string(index=False))
