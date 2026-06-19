"""Shared GO/GAF utilities for the cross-species GO-transfer benchmark.
ID space = UniProt accession (GAF DB_Object_ID for GOA/QuickGO files).
"""
import gzip, math, sys
from collections import defaultdict

EXPERIMENTAL = {"EXP","IDA","IPI","IMP","IGI","IEP","HTP","HDA","HMP","HGI","HEP"}
ASPECT = {"P":"BP","F":"MF","C":"CC"}  # GAF aspect -> namespace label


def open_maybe_gz(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_gaf(path, evidence_keep=None, drop_not=True,
              exclude_with_accs=None, exclude_taxon=None, count_excluded=None):
    """Yield (gene_acc, go_id, evidence, aspect_label) from a GAF 2.x file.
    gene_acc = column 2 (DB_Object_ID) = UniProt accession for GOA/QuickGO.
    evidence_keep: set of evidence codes to keep (None = all).
    Circularity filter (addendum): drop rows whose With/From (col 8) references the FOCAL
    species (its accessions or taxon) -- those are GO propagated FROM focal by orthology,
    so transferring them back would be circular and fake a high-%identity plateau.
      exclude_with_accs: set of focal UniProt accessions
      exclude_taxon: e.g. 'taxon:7955'
      count_excluded: optional 1-element list; [0] is incremented per dropped circular row
    """
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("!") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            qualifier, go_id, evidence, aspect = c[3], c[4], c[6], c[8]
            if drop_not and "NOT" in qualifier.split("|"):
                continue
            if evidence_keep is not None and evidence not in evidence_keep:
                continue
            if aspect not in ASPECT:
                continue
            if exclude_with_accs or exclude_taxon:
                wf = c[7]
                circ = False
                if exclude_taxon and exclude_taxon in wf:
                    circ = True
                elif exclude_with_accs and wf:
                    for tok in wf.replace(",", "|").split("|"):
                        acc = tok.split(":")[-1]
                        if acc in exclude_with_accs:
                            circ = True; break
                if circ:
                    if count_excluded is not None:
                        count_excluded[0] += 1
                    continue
            yield c[1], go_id, evidence, ASPECT[aspect]


def load_godag(obo_path):
    """Load GO DAG with is_a + part_of ancestor lookups via goatools."""
    from goatools.obo_parser import GODag
    dag = GODag(obo_path, optional_attrs={"relationship"}, prt=None)
    return dag


def ancestors(dag, go_id):
    """All ancestors over is_a + part_of (true-path), excluding self."""
    if go_id not in dag:
        return set()
    seen, stack = set(), [go_id]
    while stack:
        t = stack.pop()
        node = dag[t]
        parents = set(p.id for p in node.parents)  # is_a
        rel = getattr(node, "relationship", {})
        parents |= set(p.id for p in rel.get("part_of", set()))
        for p in parents:
            if p not in seen:
                seen.add(p); stack.append(p)
    seen.discard(go_id)
    return seen


def propagate(annot_pairs, dag):
    """annot_pairs: iterable of (gene, go_id, aspect). Returns dict gene->set(go_id)
    with true-path propagation; also returns go->namespace map for kept terms."""
    gene2go = defaultdict(set)
    go_ns = {}
    cache = {}
    for gene, go_id, asp in annot_pairs:
        if go_id not in dag:
            continue
        if go_id not in cache:
            cache[go_id] = ancestors(dag, go_id)
        terms = {go_id} | cache[go_id]
        gene2go[gene].update(terms)
        for t in terms:
            go_ns.setdefault(t, dag[t].namespace if t in dag else asp)
    return gene2go, go_ns


# --- conserved vs taxon-specific classification (Primmer caveat, design principle #3) ---
# Ancestry-based: a term inherits the category of its high-level ancestors.
SPECIFIC_ROOTS = {  # lineage-variable / multicellular / taxon-specific processes
    "GO:0032502",  # developmental process
    "GO:0048856",  # anatomical structure development
    "GO:0002376",  # immune system process
    "GO:0007610",  # behavior
    "GO:0032501",  # multicellular organismal process
    "GO:0022414",  # reproductive process
    "GO:0007154",  # cell communication
    "GO:0023052",  # signaling
}
CONSERVED_ROOTS = {  # deeply conserved housekeeping / primary metabolism
    "GO:0008152",  # metabolic process
    "GO:0006412",  # translation
    "GO:0006259",  # DNA metabolic process
    "GO:0007049",  # cell cycle
    "GO:0006091",  # generation of precursor metabolites and energy
    "GO:0042254",  # ribosome biogenesis
    "GO:0016070",  # RNA metabolic process
}

def classify_term(dag, go_id):
    """Return 'specific' | 'conserved' | 'other' by ancestry.
    Specific roots take precedence (Primmer's risky taxon-specific processes)."""
    if go_id not in dag:
        return "other"
    anc = ancestors(dag, go_id) | {go_id}
    if anc & SPECIFIC_ROOTS:
        return "specific"
    if anc & CONSERVED_ROOTS:
        return "conserved"
    return "other"


def compute_ic(gene2go):
    """Information content per term from a gene2go map (already propagated).
    IC = -log2(n_genes_with_term / n_total_genes)."""
    n_total = len(gene2go)
    cnt = defaultdict(int)
    for go_set in gene2go.values():
        for t in go_set:
            cnt[t] += 1
    ic = {}
    for t, n in cnt.items():
        ic[t] = (-math.log2(n / n_total), n)
    return ic, n_total
