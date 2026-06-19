#!/bin/bash
# S13 orchestrator — waits for eggNOG DB, smoke-tests emapper on yeast, then runs the full
# eggNOG method-comparison: S13a (emapper all) -> S13b (maps+configs) -> S13c (eval) -> S13d (compare).
set -u
ROOT="${GOTX_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
GOTX=python3
EMAPPER=emapper.py
DB=$ROOT/data/eggnog_db
cd $ROOT
S=logs/s13_orchestrator.log
say(){ echo "[$(date +%H:%M:%S)] $*" >> $S; }
: > $S
say "waiting for eggNOG DB download"
while ! grep -q "EGGNOG MANUAL DOWNLOAD DONE" logs/eggnog_download.log 2>/dev/null; do sleep 30; done
say "DB ready: $(du -sh $DB | cut -f1)"
[ -f "$DB/eggnog.db" ] && [ -f "$DB/eggnog_proteins.dmnd" ] || { say "FATAL: db files missing"; exit 1; }

# smoke test on yeast (skip if already produced)
mkdir -p results/eggnog_eval/annotations
if [ ! -s results/eggnog_eval/annotations/yeast.emapper.annotations ]; then
  say "smoke test: emapper on yeast"
  $EMAPPER -i data/proteomes/yeast.prot.fasta -o yeast --output_dir results/eggnog_eval/annotations \
    --temp_dir results/eggnog_eval/annotations -m diamond --data_dir "$DB" --cpu 16 \
    --no_file_comments --override --dmnd_iterate no >> logs/s13a_emapper.log 2>&1
  [ -s results/eggnog_eval/annotations/yeast.emapper.annotations ] || { say "FATAL: emapper smoke test no output"; exit 1; }
fi
say "smoke ok ($(wc -l < results/eggnog_eval/annotations/yeast.emapper.annotations) lines)"

say "S13a: emapper on all proteomes (parallel, long)"
CPU=16 PAR=4 bash scripts/S13a_emapper.sh >> $S 2>&1
say "S13b: build eggnog maps + temp configs"
$GOTX scripts/S13b_eggnog_map.py >> $S 2>&1
say "S13c: transfer+enrich+metrics for besthit/rbh/eggnog (subset)"
bash scripts/S13c_run_eggnog_eval.sh fish mammal plant_rice plant_arabidopsis fungi insect >> $S 2>&1
say "S13d: compare methods"
$GOTX scripts/S13d_eggnog_compare.py >> $S 2>&1
say "S13 ALL DONE"
