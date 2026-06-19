#!/bin/bash
# S0 — data acquisition for Fish track. UniProt-centric. Logs sizes + access date (Primmer Box 6).
set -u
ROOT="${GOTX_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
PDIR=$ROOT/data/proteomes; GDIR=$ROOT/data/gaf; ODIR=$ROOT/data/ontology
PBASE=https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota
GBASE=https://ftp.ebi.ac.uk/pub/databases/GO/goa
QGO=https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch
STAMP=$(date -u +%Y-%m-%dT%H:%M:%SZ)
MAN=$ROOT/data/MANIFEST.tsv
echo -e "file\tsource_url\tbytes\taccess_utc" > $MAN

dl () { # url outpath
  local url=$1 out=$2
  if [ -s "$out" ]; then echo "[skip] $out exists"; else
    echo "[get ] $url"
    curl -fSL --retry 3 --retry-delay 5 -o "$out" "$url" || { echo "[FAIL] $url"; return 1; }
  fi
  local sz=$(stat -c%s "$out" 2>/dev/null || echo 0)
  echo -e "$(basename $out)\t$url\t$sz\t$STAMP" >> $MAN
  echo "[ ok ] $out  ($sz bytes)"
}

# proteome ID map: name -> UP_taxid
declare -A P=( [zebrafish]=UP000000437_7955 [medaka]=UP000001038_8090 \
  [chicken]=UP000000539_9031 [mouse]=UP000000589_10090 [human]=UP000005640_9606 \
  [fruitfly]=UP000000803_7227 [yeast]=UP000002311_559292 )

echo "===== ontology ====="
dl "https://current.geneontology.org/ontology/go-basic.obo" "$ODIR/go-basic.obo"

echo "===== proteomes (canonical fasta + gene2acc) ====="
for s in "${!P[@]}"; do
  up=${P[$s]}; updir=${up%%_*}
  dl "$PBASE/$updir/${up}.fasta.gz"   "$PDIR/${s}.fasta.gz"
  dl "$PBASE/$updir/${up}.gene2acc.gz" "$PDIR/${s}.gene2acc.gz"
done
# focal idmapping (for DEG-list ID -> UniProt later)
dl "$PBASE/UP000000437/UP000000437_7955.idmapping.gz" "$PDIR/zebrafish.idmapping.gz"

echo "===== GAFs (GOA species dirs) ====="
declare -A G=( [zebrafish]=ZEBRAFISH/goa_zebrafish.gaf.gz [chicken]=CHICKEN/goa_chicken.gaf.gz \
  [mouse]=MOUSE/goa_mouse.gaf.gz [human]=HUMAN/goa_human.gaf.gz \
  [fruitfly]=FLY/goa_fly.gaf.gz [yeast]=YEAST/goa_yeast.gaf.gz )
for s in "${!G[@]}"; do dl "$GBASE/${G[$s]}" "$GDIR/${s}.gaf.gz"; done

echo "===== medaka via QuickGO (no dedicated GOA dir) ====="
if [ ! -s "$GDIR/medaka.gaf.gz" ]; then
  echo "[get ] QuickGO taxonId=8090"
  curl -fSL --retry 3 -H "Accept:text/gaf" \
    "$QGO?taxonId=8090&taxonUsage=exact&downloadLimit=500000" -o "$GDIR/medaka.gaf" \
    && gzip -f "$GDIR/medaka.gaf" \
    && { sz=$(stat -c%s "$GDIR/medaka.gaf.gz"); echo -e "medaka.gaf.gz\t$QGO?taxonId=8090\t$sz\t$STAMP" >> $MAN; echo "[ ok ] medaka.gaf.gz ($sz)"; } \
    || echo "[FAIL] medaka QuickGO"
else echo "[skip] medaka.gaf.gz exists"; fi

echo "===== DONE $(date -u) ====="
echo "manifest: $MAN"
