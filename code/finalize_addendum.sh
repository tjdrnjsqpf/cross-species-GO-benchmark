#!/bin/bash
# Final batch after all 6 tracks complete:
#  - fish: re-run list->S6 to include real DEG (now copied in) + fig2/fig4 fixes
#  - other tracks: regenerate S6 only (fig2 dodge fix; fig4 auto-skips w/o real)
#  - cross-clade: S7 (6 tracks), S9 (category), S10 (ID50), S11 (fungi WGD)
set -u
ROOT="${GOTX_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
GOTX=python3
RS=Rscript
cd $ROOT
declare -A FOCAL=( [fish]=zebrafish [mammal]=mouse [plant_rice]=rice [plant_arabidopsis]=arabidopsis [fungi]=yeast [insect]=fruitfly )
declare -A EQD=( [fish]=chicken,mouse,human [mammal]=human,cow,dog,pig [plant_rice]=arabidopsis,tomato,soybean,poplar [plant_arabidopsis]=rice,maize [fungi]=pombe,aspergillus,neurospora [insect]=apis,tribolium,anopheles )

echo "### fish: re-run lists (include real DEG) + figs"
bash scripts/rerun_lists.sh config/track_fish.yaml

echo "### regenerate S6 for other tracks (fig2 fix)"
for t in mammal plant_rice plant_arabidopsis fungi insect; do
  export GOTX_OUT=$ROOT/results/$t GOTX_CONFIG=$ROOT/config/track_$t.yaml \
         GOTX_TRACK=$t GOTX_FOCAL=${FOCAL[$t]} GOTX_EQD=${EQD[$t]}
  $RS scripts/S6_figures.R >> logs/rerun_$t.log 2>&1 && echo "  [$t] S6 ok" || echo "  [$t] S6 FAIL"
done

echo "### cross-clade analyses (6 tracks)"
$GOTX scripts/S7_crossclade.py   2>&1 | tail -2
$GOTX scripts/S9_category.py     2>&1 | grep -vE "load obo|nodes imported" | tail -3
$GOTX scripts/S10_id50.py        2>&1 | grep -vE "load obo|nodes imported" | tail -3
$GOTX scripts/S11_wgd_fungi.py   2>&1 | tail -3
echo "### FINALIZE DONE"
