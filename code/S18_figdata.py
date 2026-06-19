#!/usr/bin/env python3
"""S18 — freeze tidy figure-input tables for local plotting (SERVER_DATA_EXPORT).
Writes deliverable/tables/{figdata_setlevel.tsv, figdata_ic.tsv, figdata_refmeta.tsv}.
UTF-8, tab-sep, header, missing=blank, ASCII names. Metric defs unchanged; columns widened.
"""
import os, glob, gzip
import numpy as np, pandas as pd, yaml

ROOT=os.environ.get("GOTX_ROOT", os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
TAB=f"{ROOT}/deliverable/tables"; os.makedirs(TAB,exist_ok=True)
CLADE={"fish":"fish","mammal":"mammal","plant_rice":"plant","plant_arabidopsis":"plant",
       "fungi":"fungi","insect":"insect"}
TRACKS=list(CLADE)
EXP={"EXP","IDA","IMP","IGI","IPI","IEP","HTP","HDA","HMP","HGI","HEP","TAS","IC"}  # experimental-ish
METHODS={"besthit","rbh","eggnog"}

def list_sizes(track):
    d=pd.read_csv(f"{ROOT}/results/{track}/enrichment/lists_all.tsv",sep="\t")
    return d.groupby("list_id").size()

def load_setlevel(path, track):
    m=pd.read_csv(path,sep="\t")
    m=m[m.method.isin(METHODS)].copy(); m["track"]=track; m["clade"]=CLADE[track]
    return m

def load_sem(path):
    s=pd.read_csv(path,sep="\t")
    parts=s.label.str.split(".",expand=True)
    s["ref"]=parts[0]; s["method"]=parts[1]; s["evset"]=parts[2]
    return s[["ref","method","evset","list_id","aspect","wang_bma"]]

# ---------- A: figdata_setlevel ----------
def list_ic_map(t):
    """list_id -> (list_ic, size_class). SYN_=small (lists_meta), SYNL_=large (lists_large_meta)."""
    ic={}; sc={}
    fm=f"{ROOT}/results/{t}/enrichment/lists_meta.tsv"
    if os.path.exists(fm):
        d=pd.read_csv(fm,sep="\t")
        for lid,v in zip(d.list_id,d.ic): ic[lid]=v; sc[lid]="small"
    fl=f"{ROOT}/results/{t}/enrichment/lists_large_meta.tsv"
    if os.path.exists(fl):
        d=pd.read_csv(fl,sep="\t")
        for lid,v in zip(d.list_id,d.list_ic): ic[lid]=v; sc[lid]="large"
    return ic,sc

def large_sizes(t):
    fl=f"{ROOT}/results/{t}/enrichment/lists_large.tsv"
    if not os.path.exists(fl): return pd.Series(dtype=int)
    d=pd.read_csv(fl,sep="\t"); return d.groupby("list_id").size()

A=[]
for t in TRACKS:
    sz=pd.concat([list_sizes(t), large_sizes(t)])      # small/real + large list sizes
    icmap,scmap=list_ic_map(t)
    # 1) main (besthit/rbh: small synth + real) 2) eggnog_eval (eggnog: small synth)
    #    3) synth_large_eval (besthit/rbh: LARGE synth)
    srcs=[(f"{ROOT}/results/{t}/metrics/setlevel_metrics.tsv",
           f"{ROOT}/results/{t}/metrics/semantic_sim.tsv", {"besthit","rbh"})]
    egg=f"{ROOT}/results/eggnog_eval/{t}/metrics/setlevel_metrics.tsv"
    if os.path.exists(egg):
        srcs.append((egg, f"{ROOT}/results/eggnog_eval/{t}/metrics/semantic_sim.tsv", {"eggnog"}))
    lrg=f"{ROOT}/results/synth_large_eval/{t}/metrics/setlevel_metrics.tsv"
    if os.path.exists(lrg):
        srcs.append((lrg, f"{ROOT}/results/synth_large_eval/{t}/metrics/semantic_sim.tsv", {"besthit","rbh"}))
    for spath,sempath,keepm in srcs:
        m=load_setlevel(spath,t); m=m[(m.ic_bin=="all")&(m.method.isin(keepm))]
        m=m.merge(load_sem(sempath),on=["ref","method","evset","list_id","aspect"],how="left")
        m["size"]=m.list_id.map(sz)
        m["list_ic"]=m.list_id.map(icmap)
        m["size_class"]=np.where(m.list_source=="real","real",m.list_id.map(scmap))
        A.append(m[["track","ref","clade","method","evset","list_source","size_class","list_id","aspect",
                    "n_truth_sig","n_trans_sig","size","list_ic","wang_bma","recall","precision","jaccard",
                    "median_pident","My"]].rename(columns={"n_truth_sig":"n_truth","n_trans_sig":"n_trans"}))
Adf=pd.concat(A,ignore_index=True)
Adf.to_csv(f"{TAB}/figdata_setlevel.tsv",sep="\t",index=False,na_rep="")
print(f"[A] figdata_setlevel.tsv rows={len(Adf)} methods={sorted(Adf.method.unique())} "
      f"size_class={Adf.size_class.value_counts().to_dict()}")

# ---------- B: figdata_ic (IC-bin stratified, aggregated over synthetic lists) ----------
B=[]
for t in TRACKS:
    m=load_setlevel(f"{ROOT}/results/{t}/metrics/setlevel_metrics.tsv",t)
    m=m[(m.ic_bin.isin(["shallow","mid","deep"]))&(m.list_source=="synthetic")]
    g=m.groupby(["track","clade","ref","method","evset","aspect","ic_bin"]).agg(
        recall=("recall","mean"), f1=("f1","mean"), n_terms=("n_truth_sig","mean"),
        median_pident=("median_pident","first")).reset_index()
    g["wang_bma"]=""   # Wang not computed per IC bin (semantic is per-list); left blank by design
    B.append(g)
Bdf=pd.concat(B,ignore_index=True)[["track","ref","clade","method","evset","aspect","ic_bin",
                                    "wang_bma","recall","f1","n_terms","median_pident"]]
Bdf.to_csv(f"{TAB}/figdata_ic.tsv",sep="\t",index=False,na_rep="")
print(f"[B] figdata_ic.tsv rows={len(Bdf)}")

# ---------- C: figdata_refmeta (ref-intrinsic richness from GAFs, ortholog counts, eqd group) ----------
def gaf_path(track, sp, cfg):
    # GAFs are downloaded per-species to data/gaf/<sp>.gaf.gz (shared, ref-intrinsic) -> prefer it.
    v=cfg["species"][sp]
    cands=[f"{ROOT}/data/gaf/{sp}.gaf.gz", f"{ROOT}/data/gaf/{sp}.gaf"]
    if v.get("goa_path"):
        p=v["goa_path"]; cands.append(p if os.path.isabs(p) else f"{ROOT}/{p}")
    for cand in cands:
        if os.path.exists(cand): return cand
    return None

def richness(gaf, exp_only):
    """avg distinct GO terms per annotated gene (direct, NOT propagated). gene = col2 (DB Object ID)."""
    if not gaf or not os.path.exists(gaf): return np.nan
    per={}
    op=gzip.open if gaf.endswith(".gz") else open
    with op(gaf,"rt",errors="replace") as fh:
        for line in fh:
            if line.startswith("!"): continue
            c=line.rstrip("\n").split("\t")
            if len(c)<15: continue
            if c[3].startswith("NOT"): continue
            if exp_only and c[6] not in EXP: continue
            per.setdefault(c[1],set()).add(c[4])
    if not per: return 0.0
    return round(sum(len(v) for v in per.values())/len(per),3)

C=[]
for t in TRACKS:
    cfg=yaml.safe_load(open(f"{ROOT}/config/track_{t}.yaml")); focal=cfg["focal"]
    eqd=set(cfg.get("equidistant_group_names",[]))
    tj=__import__("json").load(open(f"{ROOT}/results/{t}/truth/{focal}_truth_summary.json"))
    cov=round(tj["n_genes_with_truth"]/tj["n_genes_universe"],3)
    m=pd.read_csv(f"{ROOT}/results/{t}/metrics/setlevel_metrics.tsv",sep="\t")
    meta=m[m.ic_bin=="all"].groupby("ref").agg(pid=("median_pident","first"),My=("My","first"))
    for ref in meta.index:
        gaf=gaf_path(t,ref,cfg)
        nb=f"{ROOT}/results/{t}/mapping/{focal}__{ref}.besthit.tsv"
        nr=f"{ROOT}/results/{t}/mapping/{focal}__{ref}.rbh.tsv"
        nob=(sum(1 for _ in open(nb))-1) if os.path.exists(nb) else np.nan
        nor=(sum(1 for _ in open(nr))-1) if os.path.exists(nr) else np.nan
        C.append(dict(track=t,ref=ref,clade=CLADE[t],
                      median_pident=round(meta.loc[ref,"pid"],1),My=meta.loc[ref,"My"],
                      ref_exp_richness=richness(gaf,True),ref_allev_richness=richness(gaf,False),
                      n_orthologs_rbh=nor,n_orthologs_besthit=nob,
                      equidistant_group=(f"{t}_eqd" if ref in eqd else ""),truth_cov=cov))
    print(f"[C] {t}: {len(meta)} refs richness computed")
Cdf=pd.DataFrame(C)
Cdf.to_csv(f"{TAB}/figdata_refmeta.tsv",sep="\t",index=False,na_rep="")
print(f"[C] figdata_refmeta.tsv rows={len(Cdf)}")
print("\n[S18] wrote A/B/C ->", TAB)
