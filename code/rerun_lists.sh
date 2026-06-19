#!/bin/bash
# Re-run only the list->enrichment->metrics->figures stages (S4a..S6) for a track.
# Use when synthetic-list definition changed but S0-S3 (truth/mapping/transfer) are unchanged.
set -eu
ROOT="${GOTX_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
CFG=$1
GOTX=python3
RS=Rscript
read TRACK FOCAL EQD < <($GOTX - "$CFG" <<'PY'
import sys, yaml
c=yaml.safe_load(open(sys.argv[1])); print(c["track"], c["focal"], ",".join(c.get("equidistant_group_names",[])))
PY
)
export GOTX_OUT=$ROOT/results/$TRACK GOTX_CONFIG=$CFG GOTX_TRACK=$TRACK GOTX_FOCAL=$FOCAL GOTX_EQD=$EQD
LOG=$ROOT/logs/rerun_${TRACK}.log
echo "=== rerun S4a..S6 $TRACK (focal $FOCAL) ===" | tee "$LOG"
$GOTX scripts/S4a_make_lists.py --focal $FOCAL >> "$LOG" 2>&1
cp "$GOTX_OUT"/enrichment/lists_synthetic.tsv "$GOTX_OUT"/enrichment/lists_all.tsv
# append real DEG lists if present (fish)
[ -f "$GOTX_OUT"/enrichment/lists_real.tsv ] && tail -n +2 "$GOTX_OUT"/enrichment/lists_real.tsv >> "$GOTX_OUT"/enrichment/lists_all.tsv || true
$GOTX scripts/S4_enrich_fast.py --manifest "$GOTX_OUT"/enrichment/annot_manifest.tsv \
  --lists "$GOTX_OUT"/enrichment/lists_all.tsv --background "$GOTX_OUT"/enrichment/background.tsv \
  --out "$GOTX_OUT"/enrichment/enrich_results.tsv >> "$LOG" 2>&1
$GOTX scripts/S5_metrics.py >> "$LOG" 2>&1
$GOTX scripts/S5b_regression.py >> "$LOG" 2>&1
$GOTX scripts/S5c_semantic.py >> "$LOG" 2>&1
$RS scripts/S6_figures.R >> "$LOG" 2>&1
echo ">> S4a category balance:" | tee -a "$LOG"; grep -E "category|namespace|synthetic lists" "$LOG" | tail -4
echo "=== DONE $TRACK ===" | tee -a "$LOG"
