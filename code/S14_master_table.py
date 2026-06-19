#!/usr/bin/env python3
"""S14 — master per-clade summary table: one row per track consolidating every headline value
(distance range, truth coverage, tolerance, ID50, conserved/specific gap, richness, paralog
signature, real-DEG external validity). Emits TSV + Markdown for the report.
Outputs: results/crossclade/master_table.{tsv,md}
"""
import os, json, glob
import numpy as np, pandas as pd

ROOT=os.environ.get("GOTX_ROOT", os.path.dirname(os.path.dirname(os.path.realpath(__file__)))); OUTD=f"{ROOT}/results/crossclade"
TRACKS={"fish":"zebrafish","mammal":"mouse","plant_rice":"rice",
        "plant_arabidopsis":"arabidopsis","fungi":"yeast","insect":"fruitfly"}

cc  = pd.read_csv(f"{OUTD}/crossclade.tsv", sep="\t")      # per ref/list Wang+recall (rbh,all,SYN)
id50= pd.read_csv(f"{OUTD}/id50.tsv", sep="\t")
cat = pd.read_csv(f"{OUTD}/category.tsv", sep="\t")
nr  = pd.read_csv(f"{OUTD}/nonreciprocity.tsv", sep="\t")

rows=[]
for t,focal in TRACKS.items():
    j=json.load(open(f"{ROOT}/results/{t}/truth/{focal}_truth_summary.json"))
    cov=j["n_genes_with_truth"]/j["n_genes_universe"]
    m=pd.read_csv(f"{ROOT}/results/{t}/metrics/setlevel_metrics.tsv",sep="\t")
    pid=m[m.ic_bin=="all"].groupby("ref").median_pident.first()
    nrefs=pid.size
    # tolerance: BP Wang at nearest ref + mean over refs (rbh/all/SYN)
    ccb=cc[(cc.track==t)&(cc.aspect=="BP")]
    near=ccb.loc[ccb.median_pident.idxmax()] if len(ccb) else None
    wang_near=ccb.groupby("ref").wang_bma.mean().reindex(
        ccb.groupby("ref").median_pident.first().sort_values(ascending=False).index).iloc[0] if len(ccb) else np.nan
    wang_mean=ccb.wang_bma.mean()
    # ID50 (BP) — representative value (logistic if identifiable, else interpolated; FIX 1)
    i=id50[(id50.track==t)&(id50.aspect=="BP")]
    id50bp=i.ID50_rep.iloc[0] if len(i) else np.nan
    id50_method=i.method_used.iloc[0] if len(i) else ""
    ci_lo=i.ID50_CI_low.iloc[0] if len(i) else np.nan
    ci_hi=i.ID50_CI_high.iloc[0] if len(i) else np.nan
    floor=i.floor.iloc[0] if len(i) else np.nan
    plateau=i.plateau.iloc[0] if len(i) else np.nan
    # conserved vs specific gap (BP, distance-matched 30-68%id)
    cm=cat[(cat.track==t)&(cat.aspect=="BP")&(cat.median_pident.between(30,68))]
    cons=cm[cm.category=="conserved"].wang_bma.mean(); spec=cm[cm.category=="specific"].wang_bma.mean()
    # paralog signature: mean best-hit non-reciprocity over refs (computed from mapping files, all 6 tracks)
    fr=[]
    for bh in glob.glob(f"{ROOT}/results/{t}/mapping/{focal}__*.besthit.tsv"):
        rb=bh.replace(".besthit.",".rbh.")
        nb=sum(1 for _ in open(bh))-1
        nrb=(sum(1 for _ in open(rb))-1) if os.path.exists(rb) else 0
        if nb>0: fr.append(1-nrb/nb)
    frac_nr=float(np.mean(fr)) if fr else np.nan
    # real DEG external validity
    sem=pd.read_csv(f"{ROOT}/results/{t}/metrics/semantic_sim.tsv",sep="\t")
    sem["method"]=sem.label.str.split(".").str[1]; sem["evset"]=sem.label.str.split(".").str[2]
    sem=sem[(sem.method=="rbh")&(sem.evset=="all")]
    sem["src"]=np.where(sem.list_id.str.startswith("REAL_"),"real","syn")
    wr=sem[sem.src=="real"].wang_bma.mean(); ws=sem[sem.src=="syn"].wang_bma.mean()
    # thin-truth flag + cross-validation note (FIX 2)
    note=""
    if cov<0.10:
        note=f"THIN TRUTH (cov {cov:.3f}); wide CI"
        if t=="plant_rice":
            note+="; cross-validated by arabidopsis (BP r=0.95, mean dWang=0.07)"
    rows.append(dict(
        clade=t, focal=focal, n_ref=nrefs,
        pid_min=round(pid.min(),1), pid_max=round(pid.max(),1),
        truth_cov=round(cov,3), n_truth=j["n_genes_with_truth"],
        wang_near=round(wang_near,3), wang_mean=round(wang_mean,3),
        ID50_BP=round(id50bp,1), ID50_method=id50_method,
        ID50_CI=f"{ci_lo:.0f}-{ci_hi:.0f}" if not (np.isnan(ci_lo) or np.isnan(ci_hi)) else "NE",
        floor=round(floor,2), plateau=round(plateau,2),
        cons_BP=round(cons,2), spec_BP=round(spec,2), cons_minus_spec=round(cons-spec,2),
        frac_nonrecip=round(frac_nr,2),
        realDEG_wang=round(wr,2), syn_wang=round(ws,2), real_syn_gap=round(ws-wr,2),
        note=note,
    ))
T=pd.DataFrame(rows)
# order clades by ID50 (tolerance)
T=T.sort_values("ID50_BP")
T.to_csv(f"{OUTD}/master_table.tsv",sep="\t",index=False)

# markdown
cols_desc={
 "clade":"track","focal":"focal","n_ref":"#ref","pid_min":"%id min","pid_max":"%id max",
 "truth_cov":"truth cov","n_truth":"#truth genes","wang_near":"Wang(near,BP)","wang_mean":"Wang(mean,BP)",
 "ID50_BP":"ID50(BP)","ID50_method":"ID50 method","ID50_CI":"ID50 95%CI",
 "floor":"floor","plateau":"plateau","cons_BP":"conserved","spec_BP":"specific",
 "cons_minus_spec":"cons−spec","frac_nonrecip":"paralog(nonrecip)",
 "realDEG_wang":"realDEG Wang","syn_wang":"syn Wang","real_syn_gap":"real−syn gap"}
md=T.rename(columns=cols_desc)
with open(f"{OUTD}/master_table.md","w") as f:
    f.write("| "+" | ".join(md.columns)+" |\n")
    f.write("|"+"|".join(["---"]*len(md.columns))+"|\n")
    for _,r in md.iterrows():
        f.write("| "+" | ".join(str(x) for x in r.values)+" |\n")

pd.set_option("display.width",200,"display.max_columns",30)
print(T.to_string(index=False))
print(f"\n[S14] wrote -> {OUTD}/master_table.tsv + master_table.md")
