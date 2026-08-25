#!/usr/bin/env python3
"""S20 — NULL Wang-similarity baselines for interpreting the tolerance FLOOR.

Reviewer point (PI): Wang BMA between two random GO-term sets is NOT 0 — the shared
high-level ancestors in the GO DAG give unrelated sets a non-trivial similarity. So a
far-distance floor of e.g. 0.58 cannot be read as "does not collapse" relative to 0; it
must be read relative to this null. Per track x aspect we compute two baselines:

  null_rand_rand  : Wang BMA of TWO random size-matched term sets   (fully unrelated)
  null_truth_rand : Wang BMA of the REAL truth-enriched set vs a random size-matched set
                    (the direct null for a transferred set at maximal divergence)

Term universe = the pool of terms actually enriched in the TRUTH run for that track x
aspect (the "enrichable" space); set sizes are sampled from the observed truth-enriched
set sizes (top-50, padj<0.05), exactly as the real Wang was scored. Wang(is_a 0.8,
part_of 0.6) + BMA are the same as S5c_semantic.py.

Run (server, after the pipeline has produced results/<track>/enrichment/enrich_results.tsv):
    export GOTX_ROOT=/path/to/repo
    python3 code/S20_null_wang_baseline.py            # all tracks, B=1000
Output: results/baseline/null_wang.tsv  (+ printed table)
"""
import os, sys, argparse, random
from collections import deque
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import lib_go

ROOT = os.environ.get("GOTX_ROOT", os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
TRACKS = ["fish", "mammal", "plant_rice", "plant_arabidopsis", "fungi", "insect"]
DISPLAY = {"fish": "Fish", "mammal": "Mammal", "plant_rice": "Plant-R",
           "plant_arabidopsis": "Plant-A", "fungi": "Fungi", "insect": "Insect"}
SIG, TOPN = 0.05, 50
W_ISA, W_PARTOF = 0.8, 0.6


def build_wang(dag):
    """Return a bma(set1,set2) using the Wang (2007) graph measure, identical to S5c."""
    sv = {}
    def wang(go):
        if go in sv: return sv[go]
        if go not in dag: sv[go] = ({go: 1.0}, 1.0); return sv[go]
        s = {go: 1.0}; dq = deque([go])
        while dq:
            x = dq.popleft(); node = dag[x]
            ps = [(p.id, W_ISA) for p in node.parents]
            rel = getattr(node, "relationship", {}); ps += [(p.id, W_PARTOF) for p in rel.get("part_of", ())]
            for pid, w in ps:
                c = w * s[x]
                if c > s.get(pid, 0.0): s[pid] = c; dq.append(pid)
        sv[go] = (s, sum(s.values())); return sv[go]
    def sim(a, b):
        sa, SVa = wang(a); sb, SVb = wang(b); common = sa.keys() & sb.keys()
        return sum(sa[t] + sb[t] for t in common) / (SVa + SVb) if common else 0.0
    def bma(s1, s2):
        if not s1 or not s2: return float("nan")
        m = {(a, b): sim(a, b) for a in s1 for b in s2}
        row = sum(max(m[(a, b)] for b in s2) for a in s1)
        col = sum(max(m[(a, b)] for a in s1) for b in s2)
        return (row + col) / (len(s1) + len(s2))
    return bma


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=1000, help="random pairs per track x aspect")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--id50", default=f"{ROOT}/results/crossclade/id50.tsv",
                    help="optional floor table to merge (Track/aspect/floor)")
    a = ap.parse_args()

    dag = lib_go.load_godag(f"{ROOT}/data/ontology/go-basic.obo")
    bma = build_wang(dag)

    floors = None
    if os.path.exists(a.id50):
        try:
            floors = pd.read_csv(a.id50, sep="\t")
            floors.columns = [c.lower() for c in floors.columns]
        except Exception:
            floors = None

    rows = []
    for ti, track in enumerate(TRACKS):
        erp = f"{ROOT}/results/{track}/enrichment/enrich_results.tsv"
        if not os.path.exists(erp):
            print(f"[S20] {track}: no enrich_results ({erp}), skip"); continue
        er = pd.read_csv(erp, sep="\t")
        tr = er[er.label == "truth"]
        for ai, asp in enumerate(["BP", "MF", "CC"]):
            ta = tr[tr.aspect == asp]
            if ta.empty: continue
            sets = []
            for lid, g in ta.groupby("list_id"):
                d = g[g.padj < SIG].nsmallest(TOPN, "pvalue")
                if len(d): sets.append(list(d.go_id))
            if not sets: continue
            universe = sorted({t for s in sets for t in s})
            sizes = [len(s) for s in sets]
            rng = random.Random(a.seed + ti * 100 + ai)
            rr, trd = [], []
            for _ in range(a.B):
                n1, n2 = rng.choice(sizes), rng.choice(sizes)
                r1 = rng.sample(universe, min(n1, len(universe)))
                r2 = rng.sample(universe, min(n2, len(universe)))
                rr.append(bma(r1, r2))
                ts = rng.choice(sets)
                rrand = rng.sample(universe, min(len(ts), len(universe)))
                trd.append(bma(ts, rrand))
            rr, trd = np.array(rr, float), np.array(trd, float)
            row = dict(track=track, aspect=asp, universe=len(universe),
                       median_setsize=int(np.median(sizes)),
                       null_rand_rand_mean=round(float(np.nanmean(rr)), 3),
                       null_rand_rand_p95=round(float(np.nanpercentile(rr, 95)), 3),
                       null_truth_rand_mean=round(float(np.nanmean(trd)), 3),
                       null_truth_rand_p95=round(float(np.nanpercentile(trd, 95)), 3))
            if floors is not None and "floor" in floors.columns and "aspect" in floors.columns:
                tcol = "track" if "track" in floors.columns else None
                match = floors[floors.aspect == asp]
                if tcol:
                    match = match[match[tcol].isin([track, DISPLAY[track]])]
                if len(match):
                    fl = float(match.iloc[0]["floor"])
                    row["floor"] = round(fl, 3)
                    row["floor_minus_truthrand"] = round(fl - row["null_truth_rand_mean"], 3)
            rows.append(row)
        print(f"[S20] {track} done", flush=True)

    out = pd.DataFrame(rows)
    od = f"{ROOT}/results/baseline"; os.makedirs(od, exist_ok=True)
    out.to_csv(f"{od}/null_wang.tsv", sep="\t", index=False)
    print(out.to_string(index=False))
    print(f"[S20] wrote {len(out)} rows -> {od}/null_wang.tsv")


if __name__ == "__main__":
    main()
