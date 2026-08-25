# Baseline analyses for reviewer feedback (S20, S21)

Two null/baseline analyses that contextualise the Wang-similarity results. Run on the
**server** (they need `go-basic.obo` + the finished per-track pipeline outputs), after
`run_track.sh` has produced `results/<track>/`.

Prereqs (already present on the analysis server):
- `data/ontology/go-basic.obo`
- `results/<track>/enrichment/enrich_results.tsv`, `.../metrics/semantic_sim.tsv`,
  `.../mapping/`, `.../transfer/`, `.../truth/<focal>_truth_annotation.tsv`
- conda env with goatools, pandas, numpy, scipy, pyyaml (same as the main pipeline)

```bash
export GOTX_ROOT=/path/to/Cross-species-GeneOntology   # repo root
```

## S20 — null Wang baselines (floor interpretation)  [PI point 2]
Wang BMA of random GO-term sets is not 0; this quantifies the null so the floor can be
read against it (not against 0).
```bash
python3 code/S20_null_wang_baseline.py            # all 6 tracks, B=1000 random pairs
```
Output `results/baseline/null_wang.tsv`: per track x aspect —
`null_rand_rand_*` (two random sets), `null_truth_rand_*` (truth vs random; the direct
null for a far-distance transferred set), and, if `results/crossclade/id50.tsv` exists,
`floor` and `floor_minus_truthrand`. **Interpretation:** a floor is "does not collapse"
only if it sits clearly above `null_truth_rand_mean`.

## S21 — random-orthology baseline (signal is real)  [PI point 3]
Shuffles which focal gene receives which reference gene's GO (density-matched control),
re-enriches, and scores Wang vs truth. Real >> random = correct orthology carries signal.
```bash
# per track config; --k shuffles per reference (start small), --refs to subset
python3 code/S21_random_orthology_baseline.py --config config/track_fish.yaml   --k 5
python3 code/S21_random_orthology_baseline.py --config config/track_mammal.yaml --k 5
# ... one per track (fish, mammal, insect, plant_rice, plant_arabidopsis, fungi)
```
Output `results/<track>/baseline/random_orthology.tsv`: per reference —
`wang_real`, `wang_random_mean`, `wang_random_sd`, `median_pident`.

Return `results/baseline/null_wang.tsv` and every `results/<track>/baseline/random_orthology.tsv`
to the Desktop for figure/text integration.
```
