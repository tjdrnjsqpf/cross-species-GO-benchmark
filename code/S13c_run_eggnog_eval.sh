#!/bin/bash
# S13c — run transfer->enrich->metrics for the eggNOG eval (subset refs, methods besthit/rbh/eggnog)
# Reuses existing truth + synthetic lists + background (symlinked by S13b). Does NOT re-run S0-S2.
set -eu
ROOT="${GOTX_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
GOTX=python3
cd $ROOT
for track in "$@"; do
  CFG=$ROOT/config/track_${track}_eggnog.yaml
  [ -f "$CFG" ] || { echo "[S13c] no temp config for $track, run S13b first"; continue; }
  FOCAL=$($GOTX -c "import yaml;print(yaml.safe_load(open('$CFG'))['focal'])")
  export GOTX_OUT=$ROOT/results/eggnog_eval/$track GOTX_CONFIG=$CFG GOTX_TRACK=$track GOTX_FOCAL=$FOCAL
  mkdir -p "$GOTX_OUT"/{transfer,metrics}
  LOG=$ROOT/logs/eggnog_eval_${track}.log
  echo "=== S13c eval $track (focal $FOCAL) ===" | tee "$LOG"
  $GOTX scripts/S3_transfer.py --config "$CFG" >> "$LOG" 2>&1
  { echo -e "label\tpath"; echo -e "truth\t$GOTX_OUT/truth/${FOCAL}_truth_annotation.tsv";
    for f in "$GOTX_OUT"/transfer/${FOCAL}__*.*.*.tsv; do b=$(basename "$f" .tsv); echo -e "${b#${FOCAL}__}\t$f"; done
  } > "$GOTX_OUT"/enrichment/annot_manifest.tsv
  $GOTX scripts/S4_enrich_fast.py --manifest "$GOTX_OUT"/enrichment/annot_manifest.tsv \
    --lists "$GOTX_OUT"/enrichment/lists_all.tsv --background "$GOTX_OUT"/enrichment/background.tsv \
    --out "$GOTX_OUT"/enrichment/enrich_results.tsv >> "$LOG" 2>&1
  $GOTX scripts/S5_metrics.py  >> "$LOG" 2>&1
  $GOTX scripts/S5c_semantic.py >> "$LOG" 2>&1
  echo "=== DONE eval $track ===" | tee -a "$LOG"
done
