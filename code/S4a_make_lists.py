#!/usr/bin/env python3
"""S4a — build SYNTHETIC gene lists of interest from focal TRUTH annotation.
Each list = member genes of a chosen GO term (size window), spanning aspects & IC bins.
Tagged conserved(metabolic) vs taxon-specific via GO-name keywords (Primmer liver caveat).
The lists are defined from TRUTH; the TEST is whether transferred annotation reproduces
the truth-enrichment of these lists (independent reference annotation -> not circular).
Outputs:
  results/enrichment/lists_synthetic.tsv   list_id, focal_acc
  results/enrichment/lists_meta.tsv        list_id, go_id, namespace, n_genes, ic, ic_bin, category
  results/enrichment/background.tsv        focal_acc   (measured proteome = focal universe)
"""
import os, sys, argparse
from collections import defaultdict
sys.path.insert(0, os.path.dirname(__file__))
import lib_go

ROOT = "/var2/lsg/Claude_Code/Cross-species-GeneOntology"
import os as _os
OUT = _os.environ.get("GOTX_OUT", ROOT + "/results")
CONSERVED_KW = ["metabolic", "biosynthetic", "catabolic", "glycoly", "respirat",
                "oxidation", "tricarboxylic", "translation", "ribosome", "rrna",
                "trna", "dna replication", "mitochond", "fatty acid", "lipid metabolic"]
SPECIFIC_KW = ["development", "differentiation", "morphogenesis", "fin ", "neural crest",
               "pigment", "swim bladder", "immune", "inflammat", "behavior",
               "regeneration", "axon", "synap", "muscle", "heart", "epitheli"]


def categorize(name):
    n = name.lower()
    if any(k in n for k in CONSERVED_KW): return "conserved"
    if any(k in n for k in SPECIFIC_KW): return "specific"
    return "other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--focal", default="zebrafish")
    ap.add_argument("--min", type=int, default=40)
    ap.add_argument("--max", type=int, default=250)
    ap.add_argument("--per_bin", type=int, default=8, help="terms per (aspect,ic_bin)")
    a = ap.parse_args()

    dag = lib_go.load_godag(f"{ROOT}/data/ontology/go-basic.obo")
    # truth annotation (propagated) -> term -> set(genes); use DIRECT-ish via IC file n_genes
    term_genes = defaultdict(set); term_ns = {}
    with open(f"{OUT}/truth/{a.focal}_truth_annotation.tsv") as fh:
        next(fh)
        for line in fh:
            g, go, ns = line.rstrip("\n").split("\t")
            term_genes[go].add(g); term_ns[go] = ns
    ic = {}
    with open(f"{OUT}/truth/{a.focal}_IC.tsv") as fh:
        next(fh)
        for line in fh:
            go, ns, v, n = line.rstrip("\n").split("\t"); ic[go] = float(v)

    # candidate terms within size window
    cand = [(go, term_ns[go], len(gs)) for go, gs in term_genes.items()
            if a.min <= len(gs) <= a.max and go in ic]
    def spread(tl, k):  # evenly sample k items across an IC-sorted list
        step = max(1, len(tl)//k); return tl[::step][:k]
    # IC bins (tertiles) within each namespace; BP additionally balanced by conserved/specific/other
    chosen = []
    for ns in ("biological_process", "molecular_function", "cellular_component"):
        terms = [(go, n) for go, t_ns, n in cand if t_ns == ns]
        if not terms: continue
        terms_ic = sorted(terms, key=lambda x: ic[x[0]])
        k = len(terms_ic)
        bins = {"shallow": terms_ic[:k//3], "mid": terms_ic[k//3:2*k//3], "deep": terms_ic[2*k//3:]}
        for binname, tl in bins.items():
            if ns == "biological_process":
                by_cat = {"conserved": [], "specific": [], "other": []}
                for go, n in tl:
                    by_cat[lib_go.classify_term(dag, go)].append((go, n))
                for cat, pool in by_cat.items():
                    for go, n in spread(pool, a.per_bin):
                        chosen.append((go, ns, n, ic[go], binname, cat, dag[go].name if go in dag else go))
            else:
                for go, n in spread(tl, a.per_bin):
                    chosen.append((go, ns, n, ic[go], binname,
                                   lib_go.classify_term(dag, go), dag[go].name if go in dag else go))

    os.makedirs(f"{OUT}/enrichment", exist_ok=True)
    fl = open(f"{OUT}/enrichment/lists_synthetic.tsv", "w"); fl.write("list_id\tfocal_acc\n")
    fm = open(f"{OUT}/enrichment/lists_meta.tsv", "w")
    fm.write("list_id\tgo_id\tnamespace\tn_genes\tic\tic_bin\tcategory\tname\n")
    nlist = 0
    for go, ns, n, icv, binname, cat, name in chosen:
        lid = f"SYN_{go.replace(':','')}"
        for g in sorted(term_genes[go]):
            fl.write(f"{lid}\t{g}\n")
        fm.write(f"{lid}\t{go}\t{ns}\t{n}\t{icv:.3f}\t{binname}\t{cat}\t{name}\n")
        nlist += 1
    fl.close(); fm.close()
    # background = focal proteome universe
    with open(f"{OUT}/enrichment/background.tsv", "w") as fb:
        fb.write("focal_acc\n")
        with open(f"{ROOT}/data/proteomes/{a.focal}.gene2acc.tsv") as fh:
            next(fh)
            for line in fh:
                fb.write(line.split("\t")[1] + "\n")
    print(f"[S4a] {nlist} synthetic lists written")
    from collections import Counter
    print("  by namespace:", Counter(c[1] for c in chosen))
    print("  by ic_bin:", Counter(c[4] for c in chosen))
    print("  by category:", Counter(c[5] for c in chosen))


if __name__ == "__main__":
    main()
