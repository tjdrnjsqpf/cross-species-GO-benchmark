#!/bin/bash
# run_track.sh <config.yaml> — full pipeline for one track into results/<track>/.
# Data (proteomes/gaf) shared in data/. Reuses existing per-species files.
set -eu
ROOT=/var2/lsg/Claude_Code/Cross-species-GeneOntology
CFG=$1
GOTX=/var2/lsg/miniforge3/envs/gotx/bin/python
RS=/var2/lsg/miniforge3/envs/gotx-r/bin/Rscript

# read track/focal/refs/eqd from config via python
read TRACK FOCAL EQD REFS < <($GOTX - "$CFG" <<'PY'
import sys, yaml
c=yaml.safe_load(open(sys.argv[1]))
refs=[s for s,v in c["species"].items() if v["role"]=="reference"]
eqd=c.get("equidistant_group_names", [])
print(c["track"], c["focal"], ",".join(eqd), " ".join(refs))
PY
)
export GOTX_OUT=$ROOT/results/$TRACK
export GOTX_CONFIG=$CFG
export GOTX_TRACK=$TRACK GOTX_FOCAL=$FOCAL GOTX_EQD=$EQD
mkdir -p "$GOTX_OUT"/{truth,mapping,transfer,enrichment,metrics,figures,tables}
LOG=$ROOT/logs/run_${TRACK}.log
echo "=== TRACK=$TRACK FOCAL=$FOCAL EQD=$EQD ===" | tee "$LOG"
echo "REFS=$REFS" | tee -a "$LOG"

# NOTE: steps write full output to $LOG via >> (no pipe), so set -e catches real failures
# (a piped `| grep` would mask a non-zero exit and silently truncate results).
echo ">> S0 download" | tee -a "$LOG";        $GOTX scripts/S0_download.py "$CFG" >> "$LOG" 2>&1
echo ">> S0b prep proteomes" | tee -a "$LOG"; $GOTX scripts/S0b_prep_proteomes.py $FOCAL $REFS >> "$LOG" 2>&1
echo ">> S1 focal truth" | tee -a "$LOG";     $GOTX scripts/S1_focal_truth.py --focal $FOCAL >> "$LOG" 2>&1
echo ">> S2 mapping" | tee -a "$LOG";         THREADS=${THREADS:-48} bash scripts/S2_map.sh $FOCAL $REFS >> "$LOG" 2>&1
echo ">> S3 transfer" | tee -a "$LOG";        $GOTX scripts/S3_transfer.py --config "$CFG" >> "$LOG" 2>&1
echo ">> S4a synthetic lists" | tee -a "$LOG";$GOTX scripts/S4a_make_lists.py --focal $FOCAL >> "$LOG" 2>&1
cp "$GOTX_OUT"/enrichment/lists_synthetic.tsv "$GOTX_OUT"/enrichment/lists_all.tsv
# annotation manifest (truth + all transfers)
{ echo -e "label\tpath"; echo -e "truth\t$GOTX_OUT/truth/${FOCAL}_truth_annotation.tsv";
  for f in "$GOTX_OUT"/transfer/${FOCAL}__*.*.*.tsv; do b=$(basename "$f" .tsv); echo -e "${b#${FOCAL}__}\t$f"; done
} > "$GOTX_OUT"/enrichment/annot_manifest.tsv
echo ">> S4 enrichment" | tee -a "$LOG"
$GOTX scripts/S4_enrich_fast.py --manifest "$GOTX_OUT"/enrichment/annot_manifest.tsv \
  --lists "$GOTX_OUT"/enrichment/lists_all.tsv --background "$GOTX_OUT"/enrichment/background.tsv \
  --out "$GOTX_OUT"/enrichment/enrich_results.tsv >> "$LOG" 2>&1
# sanity: enrich_results must contain every manifest label (else a step silently truncated)
NMAN=$(($(wc -l < "$GOTX_OUT"/enrichment/annot_manifest.tsv)-1))
NRES=$(cut -f1 "$GOTX_OUT"/enrichment/enrich_results.tsv | tail -n +2 | sort -u | wc -l)
echo "   enrich labels: $NRES / manifest $NMAN" | tee -a "$LOG"
echo ">> S5 metrics" | tee -a "$LOG";    $GOTX scripts/S5_metrics.py >> "$LOG" 2>&1
echo ">> S5b regression" | tee -a "$LOG";$GOTX scripts/S5b_regression.py >> "$LOG" 2>&1
echo ">> S5c semantic" | tee -a "$LOG";  $GOTX scripts/S5c_semantic.py >> "$LOG" 2>&1
echo ">> S6 figures" | tee -a "$LOG";     $RS scripts/S6_figures.R >> "$LOG" 2>&1
echo "=== DONE $TRACK -> $GOTX_OUT ===" | tee -a "$LOG"
