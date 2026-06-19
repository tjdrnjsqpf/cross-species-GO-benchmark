#!/usr/bin/env python3
"""S3 — transfer reference GO onto focal genes via the S2 ortholog maps.
For each reference x method x evidence-set:
  - build reference acc->GO (propagated true-path); evidence-set in {all, exp}
  - for each focal gene mapped to a ref acc, borrow that ref acc's GO
Outputs results/transfer/<focal>__<ref>.<method>.<evset>.tsv  (focal_acc, go_id, namespace)
Plus results/transfer/transfer_summary.tsv with mapping + annotation hit-rates.
"""
import os, sys, glob, yaml, argparse
from collections import defaultdict
sys.path.insert(0, os.path.dirname(__file__))
import lib_go

ROOT = "/var2/lsg/Claude_Code/Cross-species-GeneOntology"
import os as _os
OUT = _os.environ.get("GOTX_OUT", ROOT + "/results")


def load_map(path):
    """focal_acc -> ref_acc (first/best). Returns dict + pident list."""
    m = {}; pid = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3: continue
            f, r, p = c[0], c[1], float(c[2])
            if f not in m:           # keep first (best) hit
                m[f] = r; pid[f] = p
    return m, pid


def ref_acc2go(gaf, dag, evset, excl_accs=None, excl_taxon=None):
    keep = lib_go.EXPERIMENTAL if evset == "exp" else None
    cnt = [0]
    pairs = [(g, go, ns) for g, go, ev, ns in lib_go.parse_gaf(
                gaf, evidence_keep=keep, exclude_with_accs=excl_accs,
                exclude_taxon=excl_taxon, count_excluded=cnt)]
    gene2go, go_ns = lib_go.propagate(pairs, dag)
    return gene2go, go_ns, cnt[0]


def load_focal_accs(focal):
    accs = set()
    p = f"{ROOT}/data/proteomes/{focal}.gene2acc.tsv"
    if os.path.exists(p):
        with open(p) as fh:
            next(fh)
            for line in fh:
                accs.add(line.split("\t")[1])
    return accs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=f"{ROOT}/config/track_fish.yaml")
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    focal = cfg["focal"]
    methods = cfg["transfer_methods"]
    refs = [s for s, v in cfg["species"].items() if v["role"] == "reference"]
    outdir = f"{OUT}/transfer"; os.makedirs(outdir, exist_ok=True)

    print("[S3] loading GO DAG ...")
    dag = lib_go.load_godag(f"{ROOT}/data/ontology/go-basic.obo")

    # circularity filter (addendum): drop reference GO whose With/From references focal
    focal_accs = load_focal_accs(focal)
    focal_taxon = f"taxon:{cfg['species'][focal]['taxid']}"
    print(f"[S3] circularity filter: focal accs={len(focal_accs)}, {focal_taxon}")

    summ = open(f"{outdir}/transfer_summary.tsv", "w")
    summ.write("focal\tref\tmethod\tevset\tn_mapped_focal\tn_ref_acc_hit\tannot_hitrate\tn_focal_with_go\tmedian_pident\tn_circular_excluded\n")

    for ref in refs:
        gaf = f"{ROOT}/data/gaf/{ref}.gaf.gz"
        if not os.path.exists(gaf):
            print(f"[S3] {ref}: GAF missing ({gaf}), SKIP"); continue
        a2g = {}; n_excl = {}
        for evset in ("all", "exp"):
            a2g[evset], go_ns, n_excl[evset] = ref_acc2go(gaf, dag, evset, focal_accs, focal_taxon)
            for method in methods:
                mp = f"{OUT}/mapping/{focal}__{ref}.{method}.tsv"
                if not os.path.exists(mp):
                    print(f"[S3] missing map {mp}, skip"); continue
                fmap, fpid = load_map(mp)
                n_mapped = len(fmap)
                ref_hit = set(); focal2go = {}
                pids = []
                for f_acc, r_acc in fmap.items():
                    pids.append(fpid[f_acc])
                    gos = a2g[evset].get(r_acc)
                    if gos:
                        ref_hit.add(r_acc)
                        focal2go[f_acc] = gos
                # write transferred annotation
                outp = f"{outdir}/{focal}__{ref}.{method}.{evset}.tsv"
                with open(outp, "w") as fo:
                    fo.write("focal_acc\tgo_id\tnamespace\n")
                    for f_acc, gos in focal2go.items():
                        for go in gos:
                            fo.write(f"{f_acc}\t{go}\t{go_ns.get(go,'')}\n")
                med_pid = sorted(pids)[len(pids)//2] if pids else 0.0
                hr = len(ref_hit)/n_mapped if n_mapped else 0.0
                summ.write(f"{focal}\t{ref}\t{method}\t{evset}\t{n_mapped}\t{len(ref_hit)}\t{hr:.3f}\t{len(focal2go)}\t{med_pid:.1f}\t{n_excl[evset]}\n")
                print(f"[S3] {ref:12s} {method:7s} {evset:3s}: mapped={n_mapped} annot_hit={hr:.2f} focal_with_go={len(focal2go)} medID={med_pid:.1f} circ_excl={n_excl[evset]}")
    summ.close()
    print(f"[S3] DONE -> {outdir}/transfer_summary.tsv")


if __name__ == "__main__":
    main()
