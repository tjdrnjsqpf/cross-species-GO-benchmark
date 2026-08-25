#!/usr/bin/env python3
"""S21 — RANDOM-ORTHOLOGY baseline. Does CORRECT ortholog transfer beat random?

Reviewer point (PI): to show the transferred signal is real (not just an effect of
annotation density), compare the real ortholog transfer against a density-matched random
control. For each reference we keep exactly the same multiset of reference genes borrowed
(same number of focal genes annotated, same reference-GO multiset) but SHUFFLE which focal
gene receives which reference gene's GO, then re-enrich and score Wang vs the focal truth.
Real >> random confirms that correct orthology carries the signal.

Per track (run with --config, like S3_transfer.py). Self-contained w.r.t. paths: depends
only on lib_go (path-agnostic) + GOTX_ROOT, and shells out to S4_enrich_fast.py (which
takes all paths as arguments). Wang(is_a 0.8, part_of 0.6)+BMA are identical to S5c. BP only.

Run (server):
    export GOTX_ROOT=/path/to/repo            # the root that contains results/, data/ontology/go-basic.obo
    export GOTX_OUT=$GOTX_ROOT/results/<track>   # optional; default derived from config track name
    python3 code/S21_random_orthology_baseline.py --config config/track_fish.yaml --k 5
    # one per track; lower --k or pass --refs a,b,c to shorten runtime
Output: results/<track>/baseline/random_orthology.tsv  (ref, wang_real, wang_random_mean, wang_random_sd, median_pident)
"""
import os, sys, argparse, random, subprocess
from collections import deque
import numpy as np, pandas as pd, yaml
sys.path.insert(0, os.path.dirname(__file__))
import lib_go  # path-agnostic GO utilities (parse_gaf / propagate / load_godag / EXPERIMENTAL)

ROOT = os.environ.get("GOTX_ROOT", os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
HERE = os.path.dirname(os.path.realpath(__file__))
W_ISA, W_PARTOF = 0.8, 0.6
SIG, TOPN = 0.05, 50


# --- helpers inlined from S3_transfer (so paths follow GOTX_ROOT, not a hardcoded ROOT) ---
def load_map(path):
    m = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3: continue
            f, r = c[0], c[1]
            if f not in m: m[f] = r     # keep first (best) hit
    return m


def load_focal_accs(focal):
    accs = set(); p = f"{ROOT}/data/proteomes/{focal}.gene2acc.tsv"
    if os.path.exists(p):
        with open(p) as fh:
            next(fh)
            for line in fh:
                accs.add(line.split("\t")[1])
    return accs


def ref_acc2go(gaf, dag, evset, excl_accs, excl_taxon):
    keep = lib_go.EXPERIMENTAL if evset == "exp" else None
    cnt = [0]
    pairs = [(g, go, ns) for g, go, ev, ns in lib_go.parse_gaf(
        gaf, evidence_keep=keep, exclude_with_accs=excl_accs,
        exclude_taxon=excl_taxon, count_excluded=cnt)]
    return lib_go.propagate(pairs, dag)  # (gene2go, go_ns)


def build_wang(dag):
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
    ap.add_argument("--config", required=True)
    ap.add_argument("--k", type=int, default=5, help="shuffles per reference")
    ap.add_argument("--method", default="rbh")
    ap.add_argument("--evset", default="all")
    ap.add_argument("--refs", default="", help="comma-separated subset of references (default all)")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()

    cfg = yaml.safe_load(open(a.config))
    focal, track = cfg["focal"], cfg["track"]
    OUT = os.environ.get("GOTX_OUT", f"{ROOT}/results/{track}")
    refs = [s for s, v in cfg["species"].items() if v["role"] == "reference"]
    if a.refs:
        want = set(a.refs.split(",")); refs = [r for r in refs if r in want]

    dag = lib_go.load_godag(f"{ROOT}/data/ontology/go-basic.obo")
    bma = build_wang(dag)
    focal_accs = load_focal_accs(focal)
    focal_taxon = f"taxon:{cfg['species'][focal]['taxid']}"
    rng = random.Random(a.seed)

    rdir = f"{OUT}/transfer_random"; os.makedirs(rdir, exist_ok=True)
    manifest = [("truth", f"{OUT}/truth/{focal}_truth_annotation.tsv")]
    label_ref = {}
    for ref in refs:
        gaf = f"{ROOT}/data/gaf/{ref}.gaf.gz"
        mp = f"{OUT}/mapping/{focal}__{ref}.{a.method}.tsv"
        if not (os.path.exists(gaf) and os.path.exists(mp)):
            print(f"[S21] {ref}: missing gaf/map, skip"); continue
        a2g, go_ns = ref_acc2go(gaf, dag, a.evset, focal_accs, focal_taxon)
        fmap = load_map(mp)
        f_list, r_list = list(fmap.keys()), list(fmap.values())
        for k in range(a.k):
            shuf = r_list[:]; rng.shuffle(shuf)
            sm = dict(zip(f_list, shuf))
            lab = f"{ref}.{a.method}.{a.evset}.shuf{k}"
            outp = f"{rdir}/{focal}__{lab}.tsv"
            with open(outp, "w") as fo:
                fo.write("focal_acc\tgo_id\tnamespace\n")
                for f_acc, r_acc in sm.items():
                    gos = a2g.get(r_acc)
                    if gos:
                        for go in gos:
                            fo.write(f"{f_acc}\t{go}\t{go_ns.get(go, '')}\n")
            manifest.append((lab, outp)); label_ref[lab] = ref
        print(f"[S21] {ref}: {a.k} shuffled maps built", flush=True)

    manp = f"{OUT}/enrichment/annot_manifest_random.tsv"
    with open(manp, "w") as fo:
        fo.write("label\tpath\n")
        for lab, p in manifest:
            fo.write(f"{lab}\t{p}\n")
    lists = f"{OUT}/enrichment/lists_all.tsv"
    if not os.path.exists(lists):
        lists = f"{OUT}/enrichment/lists_synthetic.tsv"
    erout = f"{OUT}/enrichment/enrich_random.tsv"
    print(f"[S21] enriching {len(manifest) - 1} shuffled annotations (slow step) ...", flush=True)
    subprocess.run([sys.executable, os.path.join(HERE, "S4_enrich_fast.py"),
                    "--manifest", manp, "--lists", lists,
                    "--background", f"{OUT}/enrichment/background.tsv", "--out", erout],
                   check=True, env={**os.environ, "GOTX_ROOT": ROOT})

    er = pd.read_csv(erout, sep="\t")
    er = er[er.aspect == "BP"]

    def topset(g):
        return list(g[g.padj < SIG].nsmallest(TOPN, "pvalue").go_id)

    truth = {lid: topset(g) for lid, g in er[er.label == "truth"].groupby("list_id")}
    recs = []
    for lab, g in er.groupby("label"):
        if lab == "truth":
            continue
        ws = []
        for lid, gg in g.groupby("list_id"):
            ts = truth.get(lid)
            if not ts:
                continue
            ws.append(bma(ts, topset(gg)))
        if ws:
            recs.append((label_ref[lab], float(np.nanmean(ws))))
    rd = (pd.DataFrame(recs, columns=["ref", "w"]).groupby("ref").w
          .agg(["mean", "std"]).reset_index()
          .rename(columns={"mean": "wang_random_mean", "std": "wang_random_sd"}))

    real = pd.read_csv(f"{OUT}/metrics/semantic_sim.tsv", sep="\t")
    real = real[real.aspect == "BP"].copy()
    real["ref"] = real.label.str.split(".").str[0]
    real["m"] = real.label.str.split(".").str[1]
    real["e"] = real.label.str.split(".").str[2]
    real = (real[(real.m == a.method) & (real.e == a.evset)]
            .groupby("ref").wang_bma.mean().reset_index().rename(columns={"wang_bma": "wang_real"}))

    out = real.merge(rd, on="ref", how="inner")
    for rmp in (f"{ROOT}/results/tables/figdata_refmeta.tsv", f"{ROOT}/data/figdata_refmeta.tsv"):
        if os.path.exists(rmp):
            try:
                out = out.merge(pd.read_csv(rmp, sep="\t")[["ref", "median_pident"]], on="ref", how="left")
            except Exception:
                pass
            break
    if "median_pident" in out:
        out = out.sort_values("median_pident", ascending=False)
    od = f"{OUT}/baseline"; os.makedirs(od, exist_ok=True)
    out.to_csv(f"{od}/random_orthology.tsv", sep="\t", index=False)
    print(out.to_string(index=False))
    print(f"[S21] wrote {len(out)} rows -> {od}/random_orthology.tsv")


if __name__ == "__main__":
    main()
