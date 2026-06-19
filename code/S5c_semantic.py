#!/usr/bin/env python3
"""S5c — set-level SEMANTIC similarity (Wang 2007) between truth-significant and
transferred-significant enriched term sets. Redundancy-robust (trap #5).
Wang method = pure GO-DAG structure (is_a w=0.8, part_of w=0.6); no IC/OrgDb needed.
Set similarity = BMA (best-match average). Sets capped to top-N by p-value for speed.
Output: results/metrics/semantic_sim.tsv  (label,list_id,aspect,n_truth,n_trans,wang_bma)
"""
import os, sys
from collections import deque
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import lib_go

ROOT = "/var2/lsg/Claude_Code/Cross-species-GeneOntology"
import os as _os
OUT = _os.environ.get("GOTX_OUT", ROOT + "/results")
SIG, TOPN = 0.05, 50
W_ISA, W_PARTOF = 0.8, 0.6

def main():
    dag = lib_go.load_godag(f"{ROOT}/data/ontology/go-basic.obo")
    svcache = {}
    def wang(go):
        if go in svcache: return svcache[go]
        if go not in dag:
            svcache[go] = ({go:1.0}, 1.0); return svcache[go]
        sval = {go: 1.0}; dq = deque([go])
        while dq:
            x = dq.popleft(); node = dag[x]
            ps = [(p.id, W_ISA) for p in node.parents]
            rel = getattr(node, "relationship", {})
            ps += [(p.id, W_PARTOF) for p in rel.get("part_of", ())]
            for pid, w in ps:
                cand = w * sval[x]
                if cand > sval.get(pid, 0.0):
                    sval[pid] = cand; dq.append(pid)
        res = (sval, sum(sval.values())); svcache[go] = res; return res
    def sim(a, b):
        sa, SVa = wang(a); sb, SVb = wang(b)
        common = sa.keys() & sb.keys()
        if not common: return 0.0
        return sum(sa[t] + sb[t] for t in common) / (SVa + SVb)
    def bma(s1, s2):
        if not s1 or not s2: return 0.0 if (s1 or s2) else float("nan")
        m = {(a, b): sim(a, b) for a in s1 for b in s2}
        row = sum(max(m[(a, b)] for b in s2) for a in s1)
        col = sum(max(m[(a, b)] for a in s1) for b in s2)
        return (row + col) / (len(s1) + len(s2))

    er = pd.read_csv(f"{OUT}/enrichment/enrich_results.tsv", sep="\t")
    def topset(df):
        d = df[df.padj < SIG].nsmallest(TOPN, "pvalue")
        return list(d.go_id)
    truth = {(l, a): topset(g) for (l, a), g in er[er.label=="truth"].groupby(["list_id","aspect"])}
    rows = []
    for label, grp in er.groupby("label"):
        if label == "truth": continue
        for (lid, asp), tg in grp.groupby(["list_id","aspect"]):
            if (lid, asp) not in truth: continue
            ts = truth[(lid, asp)]; xs = topset(tg)
            if not ts and not xs: continue
            rows.append(dict(label=label, list_id=lid, aspect=asp,
                             n_truth=len(ts), n_trans=len(xs), wang_bma=bma(ts, xs)))
        print(f"[S5c] done {label}", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(f"{OUT}/metrics/semantic_sim.tsv", sep="\t", index=False)
    print(f"[S5c] wrote {len(out)} rows -> results/metrics/semantic_sim.tsv")

if __name__ == "__main__":
    main()
