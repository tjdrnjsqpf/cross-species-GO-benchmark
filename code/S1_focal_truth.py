#!/usr/bin/env python3
"""S1 — build focal TRUTH annotation (EXP-only) with true-path propagation + IC.
Truth gene universe = representative accessions of the focal proteome (S0b output).
Outputs (results/truth/):
  <focal>_truth_annotation.tsv   gene_acc, go_id, namespace   (propagated, EXP-only)
  <focal>_IC.tsv                 go_id, namespace, ic, n_genes
  <focal>_truth_summary.json
"""
import json, os, sys, argparse
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__))
import lib_go

ROOT = "/var2/lsg/Claude_Code/Cross-species-GeneOntology"
import os as _os
OUT = _os.environ.get("GOTX_OUT", ROOT + "/results")


def load_universe(species):
    acc = set()
    with open(f"{ROOT}/data/proteomes/{species}.gene2acc.tsv") as fh:
        next(fh)
        for line in fh:
            acc.add(line.split("\t")[1])
    return acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--focal", default="zebrafish")
    a = ap.parse_args()
    focal = a.focal
    outdir = f"{OUT}/truth"; os.makedirs(outdir, exist_ok=True)

    print(f"[S1] loading GO DAG ...")
    dag = lib_go.load_godag(f"{ROOT}/data/ontology/go-basic.obo")
    universe = load_universe(focal)
    print(f"[S1] focal universe (genes): {len(universe)}")

    gaf = f"{ROOT}/data/gaf/{focal}.gaf.gz"
    raw = []  # (gene, go, aspect) restricted to universe, EXP-only
    ev_counter = Counter()
    for gene, go_id, ev, ns in lib_go.parse_gaf(gaf, evidence_keep=lib_go.EXPERIMENTAL):
        if gene in universe:
            raw.append((gene, go_id, ns)); ev_counter[ev]+=1
    print(f"[S1] EXP annotations in-universe: {len(raw)} ; evidence: {dict(ev_counter)}")

    gene2go, go_ns = lib_go.propagate(raw, dag)
    ic, n_total = lib_go.compute_ic(gene2go)

    # write propagated annotation
    with open(f"{outdir}/{focal}_truth_annotation.tsv","w") as fo:
        fo.write("gene_acc\tgo_id\tnamespace\n")
        for gene, gos in gene2go.items():
            for go_id in gos:
                fo.write(f"{gene}\t{go_id}\t{go_ns.get(go_id,'')}\n")
    with open(f"{outdir}/{focal}_IC.tsv","w") as fo:
        fo.write("go_id\tnamespace\tic\tn_genes\n")
        for go_id,(v,n) in sorted(ic.items(), key=lambda x:-x[1][0]):
            fo.write(f"{go_id}\t{go_ns.get(go_id,'')}\t{v:.4f}\t{n}\n")

    ns_terms = Counter(go_ns[t] for t in go_ns)
    summary = {
        "focal": focal,
        "n_genes_universe": len(universe),
        "n_genes_with_truth": n_total,
        "n_terms_total": len(go_ns),
        "terms_by_namespace": dict(ns_terms),
        "evidence_counts": dict(ev_counter),
    }
    with open(f"{outdir}/{focal}_truth_summary.json","w") as fo:
        json.dump(summary, fo, indent=2)
    print(f"[S1] DONE: {n_total} genes with truth, {len(go_ns)} terms")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
