#!/usr/bin/env python3
"""S17 resolve — pick, per track, focal + references spanning the %identity range that EXIST in
STRING v12, for orthogonal network-conservation validation. Handles STRING taxid quirks
(yeast S288C 559292 -> 4932). Writes data/string/species_list.tsv (cached probe).
"""
import os, time, urllib.request, yaml
import pandas as pd

ROOT=os.environ.get("GOTX_ROOT", os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
OUT=f"{ROOT}/data/string"; os.makedirs(OUT,exist_ok=True)
BASE="https://stringdb-downloads.org/download/protein.links.detailed.v12.0"
TRACKS=["fish","mammal","plant_rice","plant_arabidopsis","fungi","insect"]
TAXID_FIX={559292:4932}      # STRING uses generic S. cerevisiae taxid
MAX_REFS=10
RNG_PICK=None

def string_taxid(tx): return TAXID_FIX.get(tx,tx)

def exists(tx):
    url=f"{BASE}/{tx}.protein.links.detailed.v12.0.txt.gz"
    req=urllib.request.Request(url,method="HEAD")
    try:
        with urllib.request.urlopen(req,timeout=30) as r: return r.status==200
    except Exception: return False

def pick_spanning(refs_pid, k):
    """refs_pid: list of (ref,pid) sorted; pick k spanning the pid range (incl. extremes)."""
    if len(refs_pid)<=k: return refs_pid
    import numpy as np
    pids=[p for _,p in refs_pid]
    targets=np.linspace(min(pids),max(pids),k)
    chosen=[]; used=set()
    for tgt in targets:
        j=min(range(len(refs_pid)),key=lambda i:(abs(refs_pid[i][1]-tgt), i in used))
        if j not in used: chosen.append(refs_pid[j]); used.add(j)
    return sorted(chosen,key=lambda x:-x[1])

rows=[]
cache={}
for t in TRACKS:
    c=yaml.safe_load(open(f"{ROOT}/config/track_{t}.yaml")); focal=c["focal"]
    m=pd.read_csv(f"{ROOT}/results/{t}/metrics/setlevel_metrics.tsv",sep="\t")
    pid=m[m.ic_bin=="all"].groupby("ref").median_pident.first().to_dict()
    refs=[(r,pid.get(r)) for r,v in c["species"].items() if v.get("role")=="reference" and r in pid]
    refs=sorted([(r,p) for r,p in refs if p is not None],key=lambda x:-x[1])
    sel=pick_spanning(refs, MAX_REFS+4)   # over-pick; some will be 404
    # focal first
    ftx=string_taxid(c["species"][focal]["taxid"])
    if ftx not in cache: cache[ftx]=exists(ftx); time.sleep(0.2)
    rows.append(dict(track=t,name=focal,role="focal",taxid=c["species"][focal]["taxid"],
                     string_taxid=ftx,median_pident=0.0,string_ok=cache[ftx]))
    kept=0
    for r,p in sel:
        if kept>=MAX_REFS: break
        tx=string_taxid(c["species"][r]["taxid"])
        if tx not in cache: cache[tx]=exists(tx); time.sleep(0.2)
        if cache[tx]:
            rows.append(dict(track=t,name=r,role="reference",taxid=c["species"][r]["taxid"],
                             string_taxid=tx,median_pident=round(p,1),string_ok=True))
            kept+=1
    print(f"[S17] {t}: focal {focal} string={cache[ftx]}, refs kept={kept}/{len(sel)} probed")
df=pd.DataFrame(rows)
df.to_csv(f"{OUT}/species_list.tsv",sep="\t",index=False)
print(df[df.role=="reference"].groupby("track").size().to_string())
print(f"\n[S17] wrote {OUT}/species_list.tsv  (focal ok: {df[df.role=='focal'].string_ok.sum()}/6)")
