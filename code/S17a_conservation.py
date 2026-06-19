#!/usr/bin/env python3
"""S17a conservation — GO-INDEPENDENT expression conservation (Bgee) per focal-ref pair
(animal tracks fish, mammal). For each RBH ortholog, Spearman of the gene's expression-rank
vector across SHARED Uberon anatomies (cross-species tissue match via Bgee/Uberon). Pair score
= median Spearman over genes with >= MINA shared anatomies. Coverage-robust (value-based, not
set-overlap). Output: results/crossclade/expr_conservation.tsv
"""
import os, gzip
import numpy as np, pandas as pd
from scipy.stats import spearmanr

ROOT=os.environ.get("GOTX_ROOT", os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
BG=f"{ROOT}/data/bgee"; OUTD=f"{ROOT}/results/crossclade"
MINA=5; MINGENES=20
tg=pd.read_csv(f"{BG}/targets.tsv",sep="\t")

def acc2ensg(proteome):
    """UniProt acc -> Ensembl GENE id (acc is the RBH ortholog key). No accset restriction."""
    up=proteome.split("_")[0]
    p=f"{ROOT}/data/proteomes/{proteome}.idmapping.gz"
    if not os.path.exists(p):
        url=f"https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/{up}/{proteome}.idmapping.gz"
        os.system(f"wget -q -O {p} '{url}'")
    m={}
    if not os.path.exists(p) or os.path.getsize(p)<100: return m
    with gzip.open(p,"rt") as fh:
        for line in fh:
            parts=line.rstrip("\n").split("\t")
            if len(parts)!=3: continue
            a,typ,val=parts
            if typ=="Ensembl":
                m.setdefault(a, val.split(".")[0])   # strip Ensembl version suffix (Bgee has none)
    return m

def load_slim_ensg(bgee_name):
    """Ensembl gene id -> {anatomy: rank} (keyed directly by Bgee Gene.ID)."""
    f=f"{BG}/{bgee_name}.slim.tsv"
    if not os.path.exists(f): return None
    d=pd.read_csv(f,sep="\t")
    g={}
    for ge,an,rk in zip(d["Gene.ID"],d["Anatomical.entity.ID"],d["rank"]):
        g.setdefault(ge,{})[an]=rk
    return g

rows=[]
for track in ["fish","mammal"]:
    s=tg[tg.track==track]; fr=s[s.role=="focal"].iloc[0]; focal=fr["name"]
    a2g_f=acc2ensg(fr.proteome); Ef=load_slim_ensg(fr.bgee_name)
    if Ef is None: print(f"[S17a] {track}: focal slim missing"); continue
    print(f"[S17a] {track} focal {focal}: acc->ENSG={len(a2g_f)}, slim ENSG={len(Ef)}")
    for _,r in s[s.role=="reference"].iterrows():
        ref=r["name"]; mp=f"{ROOT}/results/{track}/mapping/{focal}__{ref}.rbh.tsv"
        if not os.path.exists(mp): continue
        a2g_r=acc2ensg(r.proteome); Er=load_slim_ensg(r.bgee_name)
        if Er is None: print(f"[S17a] {track}/{ref}: ref slim missing"); continue
        o=pd.read_csv(mp,sep="\t"); f2r=dict(zip(o.focal_acc.astype(str),o.ref_acc.astype(str)))
        rhos=[]
        for g,gp in f2r.items():
            ge=a2g_f.get(g); gpe=a2g_r.get(gp)
            if ge is None or gpe is None: continue
            ef=Ef.get(ge); er=Er.get(gpe)
            if not ef or not er: continue
            shared=set(ef)&set(er)
            if len(shared)<MINA: continue
            xs=[ef[a] for a in shared]; ys=[er[a] for a in shared]
            rho=spearmanr(xs,ys)[0]
            if not np.isnan(rho): rhos.append(rho)
        if len(rhos)>=MINGENES:
            rows.append(dict(track=track,ref=ref,median_pident=r.median_pident,
                             expr_conservation=float(np.median(rhos)),n_genes=len(rhos)))
            print(f"[S17a] {track}/{ref}: expr_cons={np.median(rhos):.3f} (n={len(rhos)})")
        else:
            print(f"[S17a] {track}/{ref}: too few genes ({len(rhos)})")
ex=pd.DataFrame(rows)
cc=pd.read_csv(f"{OUTD}/crossclade.tsv",sep="\t")
ccb=cc[cc.aspect=="BP"].groupby(["track","ref"]).agg(GO_Wang=("wang_bma","mean"),GO_recall=("recall","mean")).reset_index()
ex=ex.merge(ccb,on=["track","ref"],how="left")
ex.to_csv(f"{OUTD}/expr_conservation.tsv",sep="\t",index=False)
print(f"\n[S17a] wrote {OUTD}/expr_conservation.tsv ({len(ex)} pairs)")
