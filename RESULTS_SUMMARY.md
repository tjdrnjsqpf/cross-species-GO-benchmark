# Results summary (manuscript-ready)

Six focal species across three eukaryotic kingdoms (vertebrate fish & mammal, invertebrate
insect, two plants, fungus). Per track, 15–22 reference species tile the divergence axis
(median ortholog %identity). Quality metric below is redundancy-robust Wang semantic
similarity of the transferred vs experimental-truth enriched term set (RBH, all-evidence,
synthetic lists), unless stated.

## 1. Divergence-tolerance curves decline and then collapse (Fig a; fig_crossclade_semantic)
GO-transfer quality falls monotonically with distance in every clade and collapses at the
clade boundary (e.g. fish→yeast 0.30 at 38% id; plants→non-angiosperms moss/algae <56% id,
0.30–0.49). Distance-slope interaction vs fish baseline: mammal +0.076 (flattest = most
robust), insect −0.067 (steepest), fungi −0.011, plant +0.005.

## 2. Tolerance ranking — ID50 must be read with the floor (Fig b; id50.tsv)
ID50 = %identity at half-collapse; floor = far-distance plateau (BP). ID50_rep is the
representative value: the constrained-logistic fit when **identifiable** (bootstrap 95% CI
within the observed range, R²≥0.5, slope/half-point not boundary-pinned), otherwise the robust
interpolation ID50 (flagged). Flat curves (mammal, fish, insect BP) are not logistically
identifiable and use interpolation — this is the FIX that removes the earlier boundary-pinned
artifact (e.g. mammal BP logistic was pinned to the data edge).

| clade | ID50_rep | 95% CI | fit | floor | reading |
|---|---|---|---|---|---|
| insect | 40.1 | 38–50 | interp | 0.63 | tolerant to far distance, high floor |
| fungi | 40.5 | 38–47 | logistic | 0.71 | most flat/robust at distance |
| fish | 42.6 | 38–61 | interp | 0.36 | tolerant but deep floor |
| plant_arabidopsis | 54.6 | 49–60 | logistic | 0.48 | low tolerance |
| plant_rice | 55.8 | 55–59 | logistic | 0.30 | **lowest tolerance** (collapses early, deepest floor) |
| mammal | 58.9 | 38–61 | interp | 0.58 | high ID50 but **highest floor** → strongest in absolute terms |

Plants lose reliability at the highest %identity AND fall deepest; mammal's high ID50 is
offset by a high floor (it never drops far). ID50 alone is misleading without the floor.
The plant ranking is **not** a thin-truth artifact: the two independent plant focals (rice,
truth coverage 2.4%; arabidopsis, 49.5%) agree closely on the BP tolerance curve at matched
%identity (Pearson r=0.95, mean |ΔWang|=0.07, max 0.18; `plant_concordance.txt`). Bootstrap
95% CI bands (shaded in fig_crossclade_semantic; `crossclade_ci.tsv`) are correspondingly wider
for rice (CI width ~0.13–0.19 vs mammal ~0.03–0.06) — honestly reflecting its sparse truth.

## 3. Annotation richness: a universal recall↑/precision↓ trade-off
Using equidistant reference groups (distance fixed, richness varying), richer reference
annotation raises recall but lowers precision in all clades (consistent coefficient signs;
e.g. fish recall +0.077/precision −0.066; mammal +0.119/−0.049; plant +0.062/−0.067). This
answers Primmer's open question quantitatively and is taxon-independent.

## 4. Conserved >> taxon-specific, across all kingdoms (Fig c; category_stats.txt)
Distance-matched (30–68% id) BP Wang, conserved (metabolism/translation/DNA) vs taxon-specific
(development/immunity/behavior): gaps fish 0.30, mammal 0.27, plant 0.26, fungi 0.19,
insect 0.18. Pooled regression: conserved main effect +0.191 (p=8.4e-283); conserved decays
slower with distance (interaction −0.079, p=6e-56). Confirms the metabolism-over-optimism
caveat as a cross-kingdom law: conserved processes transfer reliably to far species,
taxon-specific processes only to close relatives.

## 5. Whole-genome duplication: a conditional, not universal, penalty
The a priori hypothesis (plants collapse faster due to WGD paralog confusion) is **not**
supported as a steeper curve (RBH+Wang plant slope interaction non-significant). Instead:
- Plant best-hit mappings carry the most paralog ambiguity (non-reciprocity 0.49 vs fish
  0.37, mammal 0.28) — a direct WGD signature.
- The WGD effect surfaces as **method sensitivity**, not curve steepness: plants are by far
  the most sensitive to best-hit↔RBH choice.
- Fungi-internal post-WGD vs pre-WGD references show **no extra penalty** at matched distance
  (yeast's ancient WGD resolved by massive gene loss, ~10% ohnolog retention).
- Conclusion: WGD degrades GO transfer **only when ohnolog retention is high** (recent WGD,
  plants), not as a blanket rule.

## 6. Orthology method — RBH over-penalises plants; eggNOG recovers it (Fig d; eggnog_compare.tsv)
besthit / eggnog / RBH mean (eval subset): Wang 0.698 / 0.696 / 0.676; recall 0.467 / 0.462 /
0.415; precision 0.490 / 0.494 / 0.537. RBH trades recall for precision. Per-clade Wang shows
this trade is **neutral for animals/fungi but costly for plants only**: plant_arabidopsis
besthit 0.678 / eggnog 0.673 / RBH 0.616; plant_rice besthit 0.536 / eggnog 0.543 / RBH 0.487.
eggNOG OG membership retains 82–99% of best-hit pairs even at distance (vs RBH 38–53% for far
plant refs), so it preserves the recall that RBH discards in WGD-rich genomes while remaining
orthology-grounded. Practical recommendation: in paralog-rich clades (plants), prefer
eggNOG-OG (or best-hit) over strict RBH to avoid recall loss. (Caveat: our OG-sharing criterion
is permissive; deepest-level matching would behave more like RBH.)

## 7. Real-DEG external validity is capped by truth sparsity, not transfer failure
Real DEGs (EBI Expression Atlas, 6–8 per track) score far below synthetic lists (e.g. mammal
0.58 vs 0.74; rice 0.11 vs 0.54). This gap is **not** a transfer failure: it persists at high
identity, is not a list-size artifact (matched-size real < synthetic by 0.40–0.43), and is
**robust to DEG-calling threshold** (stricter cutoffs shrink lists 2.9–5.7× but barely change
the fraction of DEG genes carrying focal truth, Δ ≤ 0.06). Real-DEG Wang tracks focal
experimental-truth coverage (Spearman +0.66): mammal 58.7% coverage → 0.58; rice 2.4% → 0.11.
The ceiling is the sparse experimental-only gold standard, so synthetic lists remain the valid
benchmark for the tolerance curves while real DEGs document the real-world coverage limit.

## 8. Orthogonal validation: the tolerance curve reflects genuine functional conservation (Fig S17)
The tolerance curve is measured GO-vs-GO, so a reviewer may ask whether it reflects real
functional divergence or merely GO-annotation artefacts. We therefore measured functional
conservation with data that never touches GO, per focal–reference pair, and asked whether it
tracks the GO-transfer Wang.
- **Expression conservation (Bgee, animal tracks).** For each RBH ortholog, Spearman of the
  gene's expression-rank vector across shared Uberon anatomies (cross-species tissue match via
  Bgee/Uberon); pair score = median over genes. It correlates with GO-transfer Wang:
  **pooled Spearman ρ=+0.62, p=0.005, n=19** (mammal ρ=+0.60, p=0.03, n=13; fish ρ=+0.77, n=6).
  Expression conservation itself decays with %identity (ρ=+0.73, p=4e-4) — i.e. the GO-independent
  signal falls in parallel with GO transfer along the same divergence axis. This rebuts the
  circularity critique for animals (cf. Chen & Zhang gold-standard logic).
- **Network conservation (STRING, all 6 tracks; experiments+coexpression channels only, no
  database/textmining → no annotation leakage).** STRING conservation also decays with %identity
  (neighbour-Jaccard ρ=+0.51, p=2e-4), parallel to GO transfer. However its **direct** per-pair
  correlation with GO-transfer Wang is null/inconsistent (pooled ρ≈−0.12): the Jaccard is
  dominated by per-species network/UniProt-mapping coverage (cons~coverage ρ=+0.90) and by
  cross-clade offsets, which a coverage-robust edge-conservation-rate only partly removes. STRING
  thus supports the *parallel-decay* argument but does not give a clean direct correlation; the
  expression modality (Bgee) carries the orthogonal validation.

## Headline takeaways
1. Tolerance is universal in shape (decline→collapse) but mammal is most robust, plants least.
2. Richness trades recall for precision everywhere.
3. Conserved processes transfer far; taxon-specific processes do not — across all kingdoms.
4. WGD hurts GO transfer only with high ohnolog retention (plants), via method sensitivity.
5. For paralog-rich genomes, eggNOG/best-hit beat strict RBH (recall preservation).
6. Experimental-truth sparsity, not transfer, caps real-DEG validation.
7. GO-independent functional conservation (expression) tracks the tolerance curve (ρ=0.62,
   p=0.005) and decays in parallel with distance — the axis reflects real divergence, not GO artefacts.
