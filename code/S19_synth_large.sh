#!/bin/bash
# S19 — generate large synthetic lists and compute their set-level + Wang metrics in an isolated
# eval dir (does NOT touch existing per-track results). Reuses S4_enrich_fast/S5_metrics/S5c_semantic.
set -eu
ROOT="${GOTX_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"; cd $ROOT
G=python3
declare -A FOCAL=( [fish]=zebrafish [mammal]=mouse [plant_rice]=rice [plant_arabidopsis]=arabidopsis [fungi]=yeast [insect]=fruitfly )
for t in "$@"; do
  f=${FOCAL[$t]}; M=$ROOT/results/$t; E=$ROOT/results/synth_large_eval/$t
  L=$ROOT/logs/s19_$t.log; echo "=== S19 $t (focal $f) ===" | tee "$L"
  GOTX_OUT=$M $G scripts/S4a_large.py --focal $f >> "$L" 2>&1
  mkdir -p "$E"/enrichment "$E"/metrics
  for d in truth transfer; do [ -e "$E/$d" ] || ln -s "$M/$d" "$E/$d"; done
  $G scripts/S4_enrich_fast.py --manifest "$M"/enrichment/annot_manifest.tsv \
     --lists "$M"/enrichment/lists_large.tsv --background "$M"/enrichment/background.tsv \
     --out "$E"/enrichment/enrich_results.tsv >> "$L" 2>&1
  GOTX_OUT=$E GOTX_CONFIG=$ROOT/config/track_$t.yaml GOTX_FOCAL=$f $G scripts/S5_metrics.py  >> "$L" 2>&1
  GOTX_OUT=$E GOTX_CONFIG=$ROOT/config/track_$t.yaml GOTX_FOCAL=$f $G scripts/S5c_semantic.py >> "$L" 2>&1
  nl=$(tail -n +2 "$M"/enrichment/lists_large.tsv | cut -f1 | sort -u | wc -l)
  echo "  [$t] large lists=$nl  setlevel rows=$(($(wc -l < "$E"/metrics/setlevel_metrics.tsv)-1))" | tee -a "$L"
done
echo "S19 ALL DONE"
