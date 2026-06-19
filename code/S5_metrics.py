#!/usr/bin/env python3
"""S5 — set-level metrics: transferred enrichment vs truth enrichment, per gene list.
Metrics: precision/recall/F1, Jaccard of significant terms (padj<0.05),
Spearman rho of p-value ranks over commonly-tested terms.
Stratified by aspect x IC-bin (term IC tertiles from focal truth).
Joins explanatory vars: divergence (My), median %identity, reference richness.
Output: results/metrics/setlevel_metrics.tsv
"""
import os, sys, math, yaml
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

ROOT = os.environ.get("GOTX_ROOT", os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import os as _os
OUT = _os.environ.get("GOTX_OUT", ROOT + "/results")
SIG = 0.05

def main():
    cfg = yaml.safe_load(open(_os.environ.get("GOTX_CONFIG", f"{ROOT}/config/track_fish.yaml")))
    sp = cfg["species"]
    focal = cfg["focal"]
    er = pd.read_csv(f"{OUT}/enrichment/enrich_results.tsv", sep="\t")

    # IC per term + tertile bins from truth
    ic = pd.read_csv(f"{OUT}/truth/{focal}_IC.tsv", sep="\t")
    ic_map = dict(zip(ic.go_id, ic.ic))
    q1, q2 = ic.ic.quantile([1/3, 2/3])
    def icbin(go):
        v = ic_map.get(go)
        if v is None: return None
        return "shallow" if v <= q1 else ("mid" if v <= q2 else "deep")
    er["ic_bin"] = er.go_id.map(icbin)

    # reference richness = EXP terms per gene (from transfer_summary focal_with_go / n) -- use ref intrinsic:
    # richness proxy = number of distinct EXP go assigned per mapped gene in rbh.exp transfer
    ts = pd.read_csv(f"{OUT}/transfer/transfer_summary.tsv", sep="\t")
    rich = {}
    for ref in ts.ref.unique():
        sub = ts[(ts.ref==ref)&(ts.method=="rbh")&(ts.evset=="exp")]
        rich[ref] = float(sub.n_focal_with_go.iloc[0]) if len(sub) else 0.0
    pid = {}
    for _,r in ts[(ts.method=="rbh")&(ts.evset=="all")].iterrows():
        pid[r.ref] = r.median_pident

    truth = er[er.label=="truth"]
    truth_idx = {(l,a): g for (l,a),g in truth.groupby(["list_id","aspect"])}

    rows = []
    for label, grp in er.groupby("label"):
        if label == "truth": continue
        ref, method, evset = label.split(".")
        for (lid, asp), tg in grp.groupby(["list_id","aspect"]):
            key = (lid, asp)
            if key not in truth_idx: continue
            trg = truth_idx[key]
            for bin_ in ["all","shallow","mid","deep"]:
                if bin_=="all":
                    tt = trg; xx = tg
                else:
                    tt = trg[trg.ic_bin==bin_]; xx = tg[tg.ic_bin==bin_]
                truth_sig = set(tt[tt.padj<SIG].go_id)
                trans_sig = set(xx[xx.padj<SIG].go_id)
                if not truth_sig and not trans_sig: continue
                inter = truth_sig & trans_sig
                union = truth_sig | trans_sig
                prec = len(inter)/len(trans_sig) if trans_sig else np.nan
                rec  = len(inter)/len(truth_sig) if truth_sig else np.nan
                f1 = (2*prec*rec/(prec+rec)) if (prec and rec and not math.isnan(prec) and not math.isnan(rec) and (prec+rec)>0) else (0.0 if (truth_sig or trans_sig) else np.nan)
                jac = len(inter)/len(union) if union else np.nan
                # spearman over commonly tested terms (pvalue)
                rho = np.nan
                m = pd.merge(tt[["go_id","pvalue"]], xx[["go_id","pvalue"]], on="go_id", suffixes=("_t","_x"))
                if len(m) >= 3:
                    rho = spearmanr(m.pvalue_t, m.pvalue_x).correlation
                rows.append(dict(label=label, ref=ref, method=method, evset=evset,
                    list_id=lid, list_source=("real" if lid.startswith("REAL_") else "synthetic"),
                    aspect=asp, ic_bin=bin_,
                    n_truth_sig=len(truth_sig), n_trans_sig=len(trans_sig), n_overlap=len(inter),
                    precision=prec, recall=rec, f1=f1, jaccard=jac, spearman=rho,
                    My=sp[ref]["timetree_My"], median_pident=pid.get(ref,np.nan),
                    ref_richness=rich.get(ref,np.nan)))
    out = pd.DataFrame(rows)
    os.makedirs(f"{OUT}/metrics", exist_ok=True)
    out.to_csv(f"{OUT}/metrics/setlevel_metrics.tsv", sep="\t", index=False)
    print(f"[S5] wrote {len(out)} rows -> results/metrics/setlevel_metrics.tsv")
    # quick aggregate preview (aspect=all-ic, BP, by ref/method/evset)
    agg = out[out.ic_bin=="all"].groupby(["ref","method","evset","aspect"]).agg(
        f1=("f1","mean"), recall=("recall","mean"), precision=("precision","mean"),
        spearman=("spearman","mean"), n=("f1","size")).reset_index()
    print(agg.to_string(index=False))


if __name__ == "__main__":
    main()
