#!/usr/bin/env python3
"""S13b — build OG-membership-strict ortholog mapping (3rd transfer method 'eggnog') and set up a
scoped evaluation directory per track.

eggnog mapping := besthit pairs whose focal & reference proteins are CO-MEMBERS of an eggNOG
orthologous group (shared OG id at any sub-root taxonomic level). This isolates orthology
stringency: besthit(loose) vs rbh(reciprocity-strict) vs eggnog(OG-curated-strict), all sharing
the besthit %identity so the divergence axis is identical.

Sets up results/eggnog_eval/<track>/ with symlinked truth + synthetic lists + background, copied
besthit/rbh subset maps, new eggnog maps, and a temp config (methods=besthit,rbh,eggnog) so the
existing S3->S4->S5/S5c can run on the subset only.
Outputs also: results/crossclade/eggnog_retention.txt (eggnog vs rbh retention of besthit pairs)
"""
import os, glob, yaml, shutil
import pandas as pd

ROOT = "/var2/lsg/Claude_Code/Cross-species-GeneOntology"
ANN  = f"{ROOT}/results/eggnog_eval/annotations"
OUTX = f"{ROOT}/results/crossclade"; os.makedirs(OUTX, exist_ok=True)

# representative reference subset per track (span %identity; must exist in besthit + be emapper'd)
SUBSET = {
  "fish": ["carp","medaka","human","celegans"],
  "mammal": ["rat","dog","xenopus","gar"],
  "plant_rice": ["o_nivara","sorghum","maize","chlamydomonas"],
  "plant_arabidopsis": ["brassica_rapa","tomato","selaginella","moss"],
  "fungi": ["s_uvarum","z_rouxii","k_lactis","pombe"],
  "insect": ["d_yakuba","anopheles","tribolium","apis"],
}

def load_og_sets(name):
    """protein acc -> set of sub-root OG ids (excludes @1|root)."""
    p = f"{ANN}/{name}.emapper.annotations"
    if not os.path.exists(p):
        return None
    og = {}
    with open(p) as fh:
        header = None
        for line in fh:
            if line.startswith("#query"):
                header = line.rstrip("\n").lstrip("#").split("\t"); continue
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if header:
                d = dict(zip(header, f)); q = d.get("query", f[0]); ogs = d.get("eggNOG_OGs","")
            else:
                q = f[0]; ogs = f[4] if len(f) > 4 else ""
            ids = set()
            for tok in ogs.split(","):
                tok = tok.strip()
                if not tok or "@1|root" in tok:
                    continue
                ids.add(tok.split("@")[0])
            if ids:
                og[q] = ids
    return og

def main():
    rep = ["track\tref\tn_besthit\tn_rbh\tn_eggnog\trbh_frac\teggnog_frac"]
    for track, refs in SUBSET.items():
        cfg = yaml.safe_load(open(f"{ROOT}/config/track_{track}.yaml"))
        focal = cfg["focal"]
        focal_og = load_og_sets(focal)
        if focal_og is None:
            print(f"[S13b] {track}: focal {focal} not emapper'd yet, skip"); continue
        base = f"{ROOT}/results/{track}"
        evald = f"{ROOT}/results/eggnog_eval/{track}"
        os.makedirs(f"{evald}/mapping", exist_ok=True)
        os.makedirs(f"{evald}/enrichment", exist_ok=True)
        # symlink shared inputs (truth, synthetic lists as lists_all, background)
        def link(src, dst):
            if os.path.islink(dst) or os.path.exists(dst): os.remove(dst)
            os.symlink(src, dst)
        if os.path.exists(f"{base}/truth"):
            if os.path.islink(f"{evald}/truth"): os.remove(f"{evald}/truth")
            if not os.path.exists(f"{evald}/truth"): os.symlink(f"{base}/truth", f"{evald}/truth")
        link(f"{base}/enrichment/lists_synthetic.tsv", f"{evald}/enrichment/lists_all.tsv")
        link(f"{base}/enrichment/background.tsv", f"{evald}/enrichment/background.tsv")

        kept_refs = []
        for ref in refs:
            bh = f"{base}/mapping/{focal}__{ref}.besthit.tsv"
            rb = f"{base}/mapping/{focal}__{ref}.rbh.tsv"
            ref_og = load_og_sets(ref)
            if not os.path.exists(bh) or ref_og is None:
                print(f"[S13b] {track}/{ref}: missing besthit or OG, skip"); continue
            b = pd.read_csv(bh, sep="\t")
            # eggnog filter: focal & ref share an OG
            def shares(r):
                fo = focal_og.get(r.focal_acc); ro = ref_og.get(r.ref_acc)
                return bool(fo and ro and (fo & ro))
            mask = b.apply(shares, axis=1)
            egg = b[mask]
            egg.to_csv(f"{evald}/mapping/{focal}__{ref}.eggnog.tsv", sep="\t", index=False)
            # copy besthit + rbh subset maps into eval
            shutil.copy(bh, f"{evald}/mapping/{focal}__{ref}.besthit.tsv")
            nrb = 0
            if os.path.exists(rb):
                shutil.copy(rb, f"{evald}/mapping/{focal}__{ref}.rbh.tsv")
                nrb = len(pd.read_csv(rb, sep="\t"))
            rep.append(f"{track}\t{ref}\t{len(b)}\t{nrb}\t{len(egg)}\t{nrb/max(1,len(b)):.3f}\t{len(egg)/max(1,len(b)):.3f}")
            kept_refs.append(ref)
            print(f"[S13b] {track}/{ref}: besthit={len(b)} rbh={nrb} eggnog={len(egg)} "
                  f"(eggnog/besthit={len(egg)/max(1,len(b)):.2f})")
        # temp config: subset species + methods=besthit,rbh,eggnog
        sub_species = {focal: cfg["species"][focal]}
        for ref in kept_refs:
            sub_species[ref] = cfg["species"][ref]
        newcfg = dict(cfg); newcfg["species"] = sub_species
        newcfg["transfer_methods"] = ["besthit","rbh","eggnog"]
        with open(f"{ROOT}/config/track_{track}_eggnog.yaml","w") as fo:
            yaml.safe_dump(newcfg, fo, sort_keys=False)
    open(f"{OUTX}/eggnog_retention.txt","w").write("\n".join(rep))
    print("\n".join(rep))
    print(f"\n[S13b] wrote eval dirs + temp configs + {OUTX}/eggnog_retention.txt")

if __name__ == "__main__":
    main()
