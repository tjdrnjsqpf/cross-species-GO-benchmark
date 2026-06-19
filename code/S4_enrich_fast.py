#!/usr/bin/env python3
"""S4 (fast) — vectorized hypergeometric over-representation, equivalent to
clusterProfiler::enricher core: pvalue = hypergeom.sf(k-1, N, M, n), BH over tested terms.
Tested terms = those hit by >=1 list gene AND minGSSize<=M<=maxGSSize (matches enricher).
Validated against clusterProfiler (see scripts/validate_enrich.R).
Output identical schema to S4_enrich.R: label list_id aspect go_id pvalue padj rank
"""
import os, sys, argparse
from collections import defaultdict, Counter
import pandas as pd
import numpy as np
from scipy.stats import hypergeom

ASPECT = {"biological_process":"BP","molecular_function":"MF","cellular_component":"CC"}
MIN_GS, MAX_GS, PCUT = 5, 2000, 0.1

def bh(pvals):
    p = np.asarray(pvals); n = len(p)
    order = np.argsort(p); ranked = p[order]
    adj = ranked * n / (np.arange(n)+1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(adj,0,1)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--lists", required=True)
    ap.add_argument("--background", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    man = pd.read_csv(a.manifest, sep="\t")
    lists = pd.read_csv(a.lists, sep="\t")
    bg = set(pd.read_csv(a.background, sep="\t").iloc[:,0])
    list_genes_full = {lid: (set(g) & bg) for lid, g in lists.groupby("list_id")["focal_acc"]}

    rows = []
    for _, r in man.iterrows():
        label, path = r["label"], r["path"]
        if not os.path.exists(path):
            continue
        ann = pd.read_csv(path, sep="\t")
        ann.columns = ["gene","go","ns"][:ann.shape[1]]
        ann = ann[ann.gene.isin(bg)]
        for nsname, asp in ASPECT.items():
            sub = ann[ann.ns == nsname]
            if sub.empty: continue
            # gene -> terms ; term -> M (size in background)
            gene2terms = defaultdict(list)
            term_M = Counter()
            for g, go in zip(sub.gene.values, sub.go.values):
                gene2terms[g].append(go)
            for g, gos in gene2terms.items():
                for go in set(gos):
                    term_M[go] += 1
            # clusterProfiler convention: background = genes annotated in THIS aspect
            annotated = set(gene2terms.keys())
            N = len(annotated)
            if N == 0: continue
            for lid, genes_full in list_genes_full.items():
                genes = genes_full & annotated      # query genes that are annotated
                n = len(genes)
                if n < 5: continue
                kc = Counter()
                for g in genes:
                    for go in set(gene2terms.get(g, ())):
                        kc[go] += 1
                # tested terms: k>=1 and size filter
                tt = [(go,k) for go,k in kc.items() if MIN_GS <= term_M[go] <= MAX_GS]
                if not tt: continue
                gos = np.array([t[0] for t in tt])
                ks = np.array([t[1] for t in tt])
                Ms = np.array([term_M[go] for go in gos])
                pv = hypergeom.sf(ks-1, N, Ms, n)
                padj = bh(pv)
                rank = pd.Series(pv).rank(method="average").values
                keep = pv < PCUT
                for go,p,pa,rk in zip(gos[keep], pv[keep], padj[keep], rank[keep]):
                    rows.append((label, lid, asp, go, float(p), float(pa), float(rk)))
        print(f"[fast] done {label}", flush=True)
    out = pd.DataFrame(rows, columns=["label","list_id","aspect","go_id","pvalue","padj","rank"])
    out.to_csv(a.out, sep="\t", index=False)
    print(f"[fast] wrote {len(out)} rows -> {a.out}")

if __name__ == "__main__":
    main()
