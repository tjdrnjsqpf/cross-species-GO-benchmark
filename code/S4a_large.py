#!/usr/bin/env python3
"""S4a_large — LARGE synthetic lists to close the size-confound (SERVER_SYNTH_LARGE).
Existing synthetic lists are <=250 genes; here we add lists in size bins 250-500/500-1k/1k-2.5k/>2.5k
so synthetic overlaps real across ALL sizes. Each list = member genes of a large (=general, low-IC)
GO term; if a (aspect,bin) has no single term, build a composite = union of that aspect's largest
terms until size in bin (composite=TRUE). list_id prefix SYNL_ so it does NOT match the SYN_ filter
used by tolerance/ID50/master_table (kept as list_source=synthetic, excluded from headline curves).
Outputs: results/<track>/enrichment/{lists_large.tsv, lists_large_meta.tsv}
"""
import os, sys, argparse
from collections import defaultdict
import lib_go
ROOT=os.environ.get("GOTX_ROOT", os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
OUT=os.environ.get("GOTX_OUT", ROOT+"/results")
BINS=[(250,500,"250-500"),(500,1000,"500-1k"),(1000,2500,"1k-2.5k"),(2500,10**9,">2.5k")]
PER=3   # single terms per (aspect,bin)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--focal",required=True); a=ap.parse_args()
    term_genes=defaultdict(set); term_ns={}
    with open(f"{OUT}/truth/{a.focal}_truth_annotation.tsv") as fh:
        next(fh)
        for line in fh:
            g,go,ns=line.rstrip("\n").split("\t"); term_genes[go].add(g); term_ns[go]=ns
    ic={}
    with open(f"{OUT}/truth/{a.focal}_IC.tsv") as fh:
        next(fh)
        for line in fh:
            go,ns,v,n=line.rstrip("\n").split("\t"); ic[go]=float(v)
    NS={"biological_process":"BP","molecular_function":"MF","cellular_component":"CC"}
    fl=open(f"{OUT}/enrichment/lists_large.tsv","w"); fl.write("list_id\tfocal_acc\n")
    fm=open(f"{OUT}/enrichment/lists_large_meta.tsv","w")
    fm.write("list_id\tgo_id\tnamespace\taspect\tn_genes\tlist_ic\tsize_class\tcomposite\n")
    nlist=0
    def emit(lid,go,ns,genes,icv,comp):
        nonlocal nlist
        for g in sorted(genes): fl.write(f"{lid}\t{g}\n")
        fm.write(f"{lid}\t{go}\t{ns}\t{NS[ns]}\t{len(genes)}\t{icv:.3f}\tlarge\t{comp}\n"); nlist+=1
    for ns in NS:
        terms=sorted([(go,len(gs)) for go,gs in term_genes.items() if term_ns[go]==ns and go in ic],
                     key=lambda x:x[1])
        for lo,hi,blab in BINS:
            single=[(go,n) for go,n in terms if lo<=n<hi]
            if single:
                step=max(1,len(single)//PER)
                for go,n in single[::step][:PER]:
                    emit(f"SYNL_{go.replace(':','')}",go,ns,term_genes[go],ic[go],"FALSE")
            else:
                # composite: union of largest terms in this aspect until size reaches [lo,hi)
                acc=set(); comps=[]
                for go,n in reversed(terms):           # largest first
                    if len(acc)>=lo: break
                    acc|=term_genes[go]; comps.append(go)
                    if len(acc)>=hi:                    # overshoot -> stop (cap by bin upper)
                        break
                if lo<=len(acc)<( hi if hi<10**9 else 10**12 ) and comps:
                    micv=sum(ic[g] for g in comps)/len(comps)
                    emit(f"SYNL_C_{NS[ns]}_{blab.replace('-','').replace('.','').replace('>','gt')}",
                         "composite:"+("+".join(comps[:5])), ns, acc, micv, "TRUE")
    fl.close(); fm.close()
    print(f"[S4a_large] {a.focal}: {nlist} large synthetic lists -> {OUT}/enrichment/lists_large.tsv")

if __name__=="__main__":
    main()
