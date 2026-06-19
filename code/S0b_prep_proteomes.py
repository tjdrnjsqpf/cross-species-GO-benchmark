#!/usr/bin/env python3
"""S0b — collapse each UniProt proteome to one protein per gene.
Gene = UniProt 'GN=' token in the FASTA header; keep the longest sequence per gene.
Entries lacking GN are kept individually (keyed by accession).
Outputs: <species>.prot.fasta (header = accession) + <species>.gene2acc.tsv (gene\tacc\tlen).
"""
import gzip, os, re, sys, glob

PDIR = "/var2/lsg/Claude_Code/Cross-species-GeneOntology/data/proteomes"
GN_RE = re.compile(r"\bGN=(\S+)")
ACC_RE = re.compile(r"^(?:sp|tr)\|([^|]+)\|")


def read_fasta(path):
    op = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    acc=None; gene=None; seq=[]
    with op as fh:
        for line in fh:
            if line.startswith(">"):
                if acc:
                    yield acc, gene, "".join(seq)
                h = line[1:].strip()
                m = ACC_RE.match(h); acc = m.group(1) if m else h.split()[0]
                mg = GN_RE.search(h); gene = mg.group(1) if mg else None
                seq=[]
            else:
                seq.append(line.strip())
    if acc:
        yield acc, gene, "".join(seq)


def main(species):
    src = os.path.join(PDIR, f"{species}.fasta.gz")
    best = {}   # key -> (acc, gene, seq)
    n_in = 0
    for acc, gene, seq in read_fasta(src):
        n_in += 1
        key = f"GN:{gene}" if gene else f"AC:{acc}"
        if key not in best or len(seq) > len(best[key][2]):
            best[key] = (acc, gene if gene else acc, seq)
    out_fa = os.path.join(PDIR, f"{species}.prot.fasta")
    out_map = os.path.join(PDIR, f"{species}.gene2acc.tsv")
    with open(out_fa,"w") as fo, open(out_map,"w") as fm:
        fm.write("gene\tacc\tlen\n")
        for acc, gene, seq in best.values():
            fo.write(f">{acc}\n{seq}\n")
            fm.write(f"{gene}\t{acc}\t{len(seq)}\n")
    print(f"{species}: {n_in} entries -> {len(best)} genes  ({out_fa})")


if __name__ == "__main__":
    targets = sys.argv[1:] or [os.path.basename(p)[:-9] for p in glob.glob(os.path.join(PDIR,"*.fasta.gz"))]
    for s in targets:
        main(s)
