#!/bin/bash
# S2 — ortholog mapping focal -> each reference via DIAMOND.
#   method 1 (besthit): focal->ref top hit
#   method 2 (rbh):     reciprocal best hit (focal<->ref)
# Params logged per Primmer Box 6.
set -eu
ROOT=/var2/lsg/Claude_Code/Cross-species-GeneOntology
OUT=${GOTX_OUT:-$ROOT/results}
PDIR=$ROOT/data/proteomes
MDIR=$OUT/mapping
DBDIR=$MDIR/db
mkdir -p "$DBDIR"
DIAMOND=/var2/lsg/miniforge3/envs/gotx/bin/diamond
THREADS=${THREADS:-32}
EVALUE=1e-5; QCOV=50; SCOV=50; SENS=--more-sensitive
FMT="6 qseqid sseqid pident length evalue bitscore"

FOCAL=${1:-zebrafish}
shift || true
REFS="$@"
[ -z "$REFS" ] && REFS="medaka chicken mouse human fruitfly yeast"

makedb () { local s=$1
  [ -f "$DBDIR/$s.dmnd" ] || { echo "[db ] $s"; $DIAMOND makedb --in "$PDIR/$s.prot.fasta" -d "$DBDIR/$s" --quiet; }
}
blast () { local q=$1 db=$2 out=$3
  [ -f "$out" ] || $DIAMOND blastp -q "$PDIR/$q.prot.fasta" -d "$DBDIR/$db" -o "$out" \
      -p $THREADS -e $EVALUE --query-cover $QCOV --subject-cover $SCOV $SENS \
      --max-target-seqs 1 --outfmt $FMT --quiet
}

echo "[S2] focal=$FOCAL refs=$REFS  ($(date -u))"
makedb "$FOCAL"
for r in $REFS; do makedb "$r"; done

LOG=$MDIR/S2_params.log
echo "diamond=$($DIAMOND version 2>/dev/null) evalue=$EVALUE qcov=$QCOV scov=$SCOV sens=more-sensitive date=$(date -u +%F)" > "$LOG"

for r in $REFS; do
  echo "[map] $FOCAL -> $r"
  fwd=$MDIR/${FOCAL}__${r}.fwd.m8
  rev=$MDIR/${FOCAL}__${r}.rev.m8
  blast "$FOCAL" "$r" "$fwd"     # focal query vs ref db
  blast "$r" "$FOCAL" "$rev"     # ref query vs focal db (for RBH)

  # besthit: focal_acc ref_acc pident evalue bitscore  (top hit already, max-target-seqs 1)
  awk -F'\t' 'BEGIN{OFS="\t"; print "focal_acc","ref_acc","pident","evalue","bitscore"}
       {print $1,$2,$3,$5,$6}' "$fwd" > "$MDIR/${FOCAL}__${r}.besthit.tsv"

  # rbh: keep pairs where focal->ref best == ref->focal best (reciprocal)
  awk -F'\t' 'NR==FNR{rev[$1]=$2; next} ($2 in rev) && (rev[$2]==$1){
        print $1"\t"$2"\t"$3"\t"$5"\t"$6}' "$rev" "$fwd" \
      | (echo -e "focal_acc\tref_acc\tpident\tevalue\tbitscore"; cat) > "$MDIR/${FOCAL}__${r}.rbh.tsv"

  nb=$(($(wc -l < "$MDIR/${FOCAL}__${r}.besthit.tsv")-1))
  nr=$(($(wc -l < "$MDIR/${FOCAL}__${r}.rbh.tsv")-1))
  echo "     besthit=$nb  rbh=$nr"
done
echo "[S2] DONE ($(date -u))"
