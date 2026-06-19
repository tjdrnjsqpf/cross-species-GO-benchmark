#!/usr/bin/env python3
"""S17b — GO-INDEPENDENT network conservation (STRING) per focal-reference pair.
Orthogonal validation: does a functional-conservation score that never touches GO track the
GO-transfer tolerance curve? STRING channels restricted to **experiments + coexpression**
(database/textmining excluded -> no annotation leakage).

Per pair: focal gene g (UniProt) with RBH ortholog g' in ref. Project g's STRING neighbours to
ref via orthologs (set P); compare to g''s actual STRING neighbours (set N_r); Jaccard. Pair
score = median Jaccard over genes with >= MINK neighbours on both sides.

Inputs: data/string/{taxid}.protein.{links.detailed,aliases}.v12.0.txt.gz, species_list.tsv,
        results/<track>/mapping/<focal>__<ref>.rbh.tsv, crossclade.tsv (GO Wang/recall)
Output: results/crossclade/string_conservation.tsv
"""
import os, gzip
import numpy as np, pandas as pd, yaml

ROOT=os.environ.get("GOTX_ROOT", os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
SDIR=f"{ROOT}/data/string"; OUTD=f"{ROOT}/results/crossclade"
PRIOR=0.041; THR=0.40; MINK=3   # noisy-OR(exp,coexp) >= THR; need >=MINK neighbours both sides
sl=pd.read_csv(f"{SDIR}/species_list.tsv",sep="\t")

def load_alias(tx):
    """string_id -> UniProt acc (prefer UniProt_AC)."""
    p=f"{SDIR}/{tx}.protein.aliases.v12.0.txt.gz"
    best={}; pri={}
    with gzip.open(p,"rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            parts=line.rstrip("\n").split("\t")
            if len(parts)<3: continue
            sid,alias,source=parts[0],parts[1],parts[2]
            if "UniProt" not in source: continue
            rank=0 if source=="UniProt_AC" else (1 if "UniProt_AC" in source else 2)
            if sid not in pri or rank<pri[sid]:
                best[sid]=alias; pri[sid]=rank
    return best

def load_neighbors(tx, id2acc):
    """UniProt acc -> set(UniProt acc neighbours), using exp+coexp noisy-OR >= THR."""
    p=f"{SDIR}/{tx}.protein.links.detailed.v12.0.txt.gz"
    nb={}
    # columns: protein1 protein2 neighborhood fusion cooccurence coexpression experimental database textmining combined_score
    it=pd.read_csv(p,sep=" ",usecols=["protein1","protein2","coexpression","experimental"],
                   dtype={"protein1":str,"protein2":str,"coexpression":np.int16,"experimental":np.int16},
                   chunksize=2_000_000)
    for ch in it:
        e=(ch.experimental/1000.0-PRIOR).clip(0)/(1-PRIOR)
        c=(ch.coexpression/1000.0-PRIOR).clip(0)/(1-PRIOR)
        comb=1-(1-e)*(1-c)
        ch=ch[comb>=THR]
        for a,b in zip(ch.protein1.map(id2acc),ch.protein2.map(id2acc)):
            if a is None or b is None or a!=a or b!=b: continue
            nb.setdefault(a,set()).add(b); nb.setdefault(b,set()).add(a)
    return nb

def main():
    rows=[]
    # cache per-taxid neighbour graphs
    graph={}
    def get_graph(tx):
        if tx not in graph:
            graph[tx]=load_neighbors(tx, load_alias(tx))
        return graph[tx]
    for track in sl.track.unique():
        s=sl[sl.track==track]; focal_row=s[s.role=="focal"].iloc[0]
        cfg=yaml.safe_load(open(f"{ROOT}/config/track_{track}.yaml")); focal=focal_row["name"]
        ftx=int(focal_row.string_taxid)
        try: Gf=get_graph(ftx)
        except Exception as ex: print(f"[S17b] {track} focal graph fail: {ex}"); continue
        for _,r in s[s.role=="reference"].iterrows():
            ref=r["name"]; rtx=int(r.string_taxid)
            mp=f"{ROOT}/results/{track}/mapping/{focal}__{ref}.rbh.tsv"
            if not os.path.exists(mp): continue
            try: Gr=get_graph(rtx)
            except Exception as ex: print(f"[S17b] {track}/{ref} graph fail: {ex}"); continue
            o=pd.read_csv(mp,sep="\t"); f2r=dict(zip(o.focal_acc,o.ref_acc))
            jac=[]
            # edge-conservation RATE (coverage-robust): of focal high-conf edges whose BOTH
            # endpoints map to ref orthologs, fraction whose ref ortholog-edge is also high-conf.
            edges_mappable=0; edges_conserved=0
            seen=set()
            for g,gp in f2r.items():
                Nf=Gf.get(g); Nr=Gr.get(gp)
                if Nf and Nr:
                    P={f2r[n] for n in Nf if n in f2r}
                    if len(P)>=MINK and len(Nr)>=MINK:
                        inter=len(P&Nr); union=len(P|Nr)
                        if union: jac.append(inter/union)
                # edge rate
                if not Nf: continue
                for n in Nf:
                    if n not in f2r: continue
                    key=(g,n) if g<n else (n,g)
                    if key in seen: continue
                    seen.add(key)
                    edges_mappable+=1
                    np_=f2r[n]
                    if np_ in Gr.get(gp,()):  # ref ortholog-edge present & high-conf
                        edges_conserved+=1
            edge_rate = edges_conserved/edges_mappable if edges_mappable>=30 else np.nan
            if len(jac)>=20 or not np.isnan(edge_rate):
                rows.append(dict(track=track,ref=ref,median_pident=r.median_pident,
                                 string_conservation=float(np.median(jac)) if jac else np.nan,
                                 edge_rate=edge_rate, n_genes=len(jac), n_edges=edges_mappable))
                print(f"[S17b] {track}/{ref}: jaccard={np.median(jac) if jac else float('nan'):.3f} "
                      f"edge_rate={edge_rate:.3f} (n_gene={len(jac)}, n_edge={edges_mappable})")
    cons=pd.DataFrame(rows)
    # join GO-transfer Wang/recall (RBH, all, BP) from crossclade.tsv
    cc=pd.read_csv(f"{OUTD}/crossclade.tsv",sep="\t")
    ccb=cc[cc.aspect=="BP"].groupby(["track","ref"]).agg(
        GO_Wang=("wang_bma","mean"),GO_recall=("recall","mean")).reset_index()
    cons=cons.merge(ccb,on=["track","ref"],how="left")
    cons.to_csv(f"{OUTD}/string_conservation.tsv",sep="\t",index=False)
    print(f"\n[S17b] wrote {OUTD}/string_conservation.tsv ({len(cons)} pairs)")

if __name__=="__main__":
    main()
