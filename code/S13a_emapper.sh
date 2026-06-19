#!/bin/bash
# S13a — run eggNOG-mapper (diamond mode) on focal + a representative reference subset per track,
# to obtain per-protein eggNOG orthologous-group (OG) membership. OGs are then used (S13b) to
# build an OG-membership-strict ortholog mapping (besthit ∩ shared-OG) = 3rd transfer method.
# Heavy/long: run in background. Skips proteomes already done. one-per-gene .prot.fasta inputs.
set -u
ROOT=/var2/lsg/Claude_Code/Cross-species-GeneOntology
EMAPPER=/var2/lsg/miniforge3/envs/eggnog/bin/emapper.py
DB=$ROOT/data/eggnog_db
CPU=${CPU:-16}
PAR=${PAR:-4}          # concurrent emapper jobs (PAR*CPU threads total)
OUTD=$ROOT/results/eggnog_eval/annotations
mkdir -p "$OUTD"
LOG=$ROOT/logs/s13a_emapper.log

# 6 focal + 4 representative refs/track spanning the %identity range (deduped across tracks)
PROTEOMES="
zebrafish mouse rice arabidopsis yeast fruitfly
carp medaka human celegans
rat dog xenopus gar
o_nivara sorghum maize chlamydomonas
brassica_rapa tomato selaginella moss
s_uvarum z_rouxii k_lactis pombe
d_yakuba anopheles tribolium apis
"

run_one(){
  local p=$1
  local fa=$ROOT/data/proteomes/$p.prot.fasta
  local ann=$OUTD/$p.emapper.annotations
  if [ ! -f "$fa" ]; then echo "  [skip] no fasta $p" >> "$LOG"; return; fi
  if [ -f "$ann" ] && [ $(wc -l < "$ann") -gt 5 ]; then echo "  [done] $p" >> "$LOG"; return; fi
  echo "[$(date +%H:%M:%S)] emapper START $p ($(grep -c '>' "$fa") prot)" >> "$LOG"
  $EMAPPER -i "$fa" -o "$p" --output_dir "$OUTD" --temp_dir "$OUTD" \
    -m diamond --data_dir "$DB" --cpu "$CPU" \
    --no_file_comments --override --dmnd_iterate no >> "$LOG" 2>&1 \
    && echo "[$(date +%H:%M:%S)]   [ok] $p" >> "$LOG" || echo "[$(date +%H:%M:%S)]   [FAIL] $p" >> "$LOG"
}
export -f run_one; export EMAPPER DB CPU OUTD LOG ROOT

echo "[$(date +%H:%M:%S)] S13a START cpu=$CPU par=$PAR" | tee -a "$LOG"
echo "$PROTEOMES" | tr ' ' '\n' | grep -v '^$' | xargs -P "$PAR" -I{} bash -c 'run_one "$@"' _ {}
echo "[$(date +%H:%M:%S)] S13a ALL DONE" | tee -a "$LOG"
