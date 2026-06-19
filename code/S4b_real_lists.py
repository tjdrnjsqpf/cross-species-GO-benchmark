#!/usr/bin/env python3
"""S4b — REAL DEG lists from EBI Expression Atlas, for ANY focal species (track-parameterized).

Per experiment x contrast: significant genes (adj.p<PADJ & |log2FC|>=LFC), mapped
Atlas GeneName/GeneID -> focal UniProt accession via data/proteomes/<focal>.gene2acc.tsv,
restricted to the track background. Differential experiments are auto-selected from the
Atlas JSON API, prioritising diverse experimentalFactors (Primmer caveat: cover taxon-specific
+ conserved biology). Curated SEED experiments (if any) are tried first.

Usage:
  GOTX_TRACK=mammal GOTX_CONFIG=config/track_mammal.yaml python scripts/S4b_real_lists.py
Outputs: results/<track>/enrichment/{lists_real.tsv, lists_real_meta.tsv}
"""
import os, re, sys, json, urllib.request, urllib.error
from collections import defaultdict
import pandas as pd, yaml

ROOT = "/var2/lsg/Claude_Code/Cross-species-GeneOntology"
TRACK   = os.environ.get("GOTX_TRACK", "fish")
CONFIG  = os.environ.get("GOTX_CONFIG", f"{ROOT}/config/track_{TRACK}.yaml")
OUT     = os.environ.get("GOTX_OUT", f"{ROOT}/results/{TRACK}")
FTP = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments"
APIEXP = "https://www.ebi.ac.uk/gxa/json/experiments"

PADJ = float(os.environ.get("GOTX_DEG_PADJ", 0.05))
LFC  = float(os.environ.get("GOTX_DEG_LFC", 1.0))
MINGENES = 15
TARGET_LISTS = 8        # stop once this many real lists are written
MAX_PER_EXP  = 3        # don't let one experiment dominate (diversity)
MAX_TRY_EXP  = 25       # cap experiments downloaded
# output filename suffix (e.g. "_strict") so sensitivity variants don't clobber the main lists
SUFFIX = os.environ.get("GOTX_DEG_SUFFIX", "")

# Atlas species string per track focal
SPECIES = {
    "fish":"Danio rerio", "mammal":"Mus musculus", "plant_rice":"Oryza sativa",
    "plant_arabidopsis":"Arabidopsis thaliana", "fungi":"Saccharomyces cerevisiae",
    "insect":"Drosophila melanogaster",
}
# curated picks tried first (known diverse biology); auto-fill handles the rest
SEED = {
    "fish": ["E-MTAB-8958","E-MTAB-10068","E-MTAB-9113","E-MTAB-5992"],
}

def load_g2a(focal, bg):
    """lowercased gene/locus -> UniProt acc, restricted to background."""
    p = f"{ROOT}/data/proteomes/{focal}.gene2acc.tsv"
    df = pd.read_csv(p, sep="\t").dropna(subset=["gene", "acc"])
    g2a = {}
    for g, a in zip(df["gene"], df["acc"]):
        a = str(a)
        if a in bg:
            g2a.setdefault(str(g).strip().lower(), a)
    return g2a

def map_gene(g2a, gid, gname):
    for cand in (gname, gid):
        if isinstance(cand, str) and cand.strip():
            a = g2a.get(cand.strip().lower())
            if a:
                return a
    return None

def select_experiments(species):
    """Return ordered list of (accession, factors, desc) differential experiments,
    SEED first, then auto-selected for experimentalFactor diversity."""
    data = json.load(urllib.request.urlopen(APIEXP, timeout=90))
    exps = data["experiments"]
    cands = []
    for e in exps:
        if e.get("species") == species and e.get("experimentType") == "Differential":
            facs = tuple(sorted(e.get("experimentalFactors", [])))
            cands.append((e["experimentAccession"], facs, e.get("experimentDescription", "")))
    seed = [c for c in cands if c[0] in SEED.get(TRACK, [])]
    rest = [c for c in cands if c[0] not in SEED.get(TRACK, [])]
    # greedy diversity: prefer experiments whose factor-set is new
    seen_fac = set(c[1] for c in seed)
    ordered = list(seed)
    for c in rest:                      # first pass: novel factor combos
        if c[1] not in seen_fac:
            ordered.append(c); seen_fac.add(c[1])
    for c in rest:                      # then remainder (for top-up)
        if c not in ordered:
            ordered.append(c)
    return ordered

def fetch_analytics(acc):
    path = f"{ROOT}/data/deg/{acc}-analytics.tsv"
    os.makedirs(f"{ROOT}/data/deg", exist_ok=True)
    if not os.path.exists(path) or os.path.getsize(path) < 1000:
        try:
            urllib.request.urlretrieve(f"{FTP}/{acc}/{acc}-analytics.tsv", path)
        except Exception as e:
            print(f"[S4b] download fail {acc}: {e}"); return None
    try:
        return pd.read_csv(path, sep="\t", low_memory=False)
    except Exception as e:
        print(f"[S4b] parse fail {acc}: {e}"); return None

def main():
    cfg = yaml.safe_load(open(CONFIG)); focal = cfg["focal"]
    species = SPECIES[TRACK]
    bg = set(pd.read_csv(f"{OUT}/enrichment/background.tsv", sep="\t").iloc[:, 0])
    g2a = load_g2a(focal, bg)
    print(f"[S4b] track={TRACK} focal={focal} species='{species}' bg={len(bg)} gene2acc={len(g2a)}")

    cands = select_experiments(species)
    print(f"[S4b] {len(cands)} differential experiments available; targeting {TARGET_LISTS} lists")

    os.makedirs(f"{OUT}/enrichment", exist_ok=True)
    print(f"[S4b] thresholds padj<{PADJ} |log2FC|>={LFC} suffix='{SUFFIX}'")
    fl = open(f"{OUT}/enrichment/lists_real{SUFFIX}.tsv", "w"); fl.write("list_id\tfocal_acc\n")
    fm = open(f"{OUT}/enrichment/lists_real{SUFFIX}_meta.tsv", "w")
    fm.write("list_id\texperiment\tcontrast\tn_deg_raw\tn_mapped\tdescription\n")

    nlist, ntried = 0, 0
    for acc, facs, desc in cands:
        if nlist >= TARGET_LISTS or ntried >= MAX_TRY_EXP:
            break
        ntried += 1
        df = fetch_analytics(acc)
        if df is None or df.shape[1] < 3:
            continue
        cols = df.columns.tolist()
        contrasts = sorted({m.group(1) for c in cols if (m := re.match(r"(.+)\.p-value$", c))})
        gid_col = cols[0]; gname_col = cols[1] if len(cols) > 1 else None
        used_here = 0
        for con in contrasts:
            if nlist >= TARGET_LISTS or used_here >= MAX_PER_EXP:
                break
            pcol, fcol = f"{con}.p-value", f"{con}.log2foldchange"
            if pcol not in df or fcol not in df:
                continue
            sig = df[(pd.to_numeric(df[pcol], errors="coerce") < PADJ) &
                     (pd.to_numeric(df[fcol], errors="coerce").abs() >= LFC)]
            accs = set()
            for _, r in sig.iterrows():
                a = map_gene(g2a, r[gid_col], r[gname_col] if gname_col else None)
                if a:
                    accs.add(a)
            if len(accs) >= MINGENES:
                lid = f"REAL_{acc}_{con}"
                for a in sorted(accs):
                    fl.write(f"{lid}\t{a}\n")
                d = (desc or "+".join(facs))[:120].replace("\t", " ")
                fm.write(f"{lid}\t{acc}\t{con}\t{len(sig)}\t{len(accs)}\t{d}\n")
                nlist += 1; used_here += 1
                print(f"[S4b] {lid}: raw_deg={len(sig)} mapped={len(accs)}  ({'+'.join(facs)[:50]})")
    fl.close(); fm.close()
    print(f"[S4b] {nlist} real DEG lists written -> {OUT}/enrichment/lists_real{SUFFIX}.tsv")
    if nlist == 0:
        print("[S4b] WARNING: no lists passed filter (check gene2acc matching)")

if __name__ == "__main__":
    main()
