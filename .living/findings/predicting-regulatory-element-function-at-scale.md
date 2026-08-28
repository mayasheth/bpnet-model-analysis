---
topic: predicting-regulatory-element-function-at-scale
description: Building and benchmarking models that predict regulatory element function, enhancer–gene links, and variant effects genome-wide (the question behind consortium-scale efforts such as IGVF).
created: 2026-07-21
last_updated: 2026-08-24
status: active
---

# Predicting regulatory element function at scale

_Seed topic (Engreitz Lab). No findings recorded yet — crystallize-findings will
append `F-NNN` entries here as analyses produce them. (Broad question behind
IGVF-style consortium work — the slug avoids naming the consortium.)_

### Open Questions
- How well do predictive models of element function generalize across cell types?
- What experimental data most improves genome-wide prediction accuracy?

---

## F-001: Sequence adds ~2.2x more predictive information over accessibility for p300 than for H3K27ac
**Status:** preliminary
**Claim:** In K562, on DNase candidate elements, models predicting p300 gain +0.301 Pearson from adding sequence to an accessibility-only model (0.305 -> 0.606 on the top signal quintile), whereas models predicting H3K27ac gain only +0.125 (0.543 -> 0.668). H3K27ac is the more predictable target overall (0.668 vs 0.606 multimodal), but a larger share of its predictability is carried by chromatin accessibility alone.
**Implications:** Target choice matters for sequence-interpretability work independently of how well a model scores. A high-performing H3K27ac model is substantially reading the accessibility track, so attributions from it will reflect less sequence-specific information than the equivalent p300 model, despite the better headline correlation. For motif-syntax questions, p300 remains the better substrate. Also a caution against ranking targets by overall correlation: the accessibility-only control is what makes the comparison interpretable.
**Tags:** h3k27ac, p300, sequence-vs-accessibility, multimodal, interpretability, k562, target-choice

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-24 | job 40783892 | ENCSR000AKP (H3K27ac), ENCSR000EGE (p300), K562 ATAC | 2026_0824_H3K27ac_model | Top-quintile Pearson, 5 folds pooled: p300 ATAC-only 0.305, p300 multimodal 0.606; H3K27ac ATAC-only 0.543, H3K27ac multimodal 0.668, H3K27ac sequence-only 0.357 | supports |
| 2026-08-26 | job 40917535 | ENCSR000AKP 5-prime target, K562 ATAC | 2026_0824_H3K27ac_model | Re-derived on the corrected 5-prime target, 5 folds, top quintile: H3K27ac ATAC-only 0.5508 -> multimodal 0.6861, so sequence adds +0.135 (was +0.125 on the superseded fragment target). p300 remains +0.301, so the ratio is ~2.2x rather than ~2.4x. Conclusion unchanged. | refines |

### Open Questions
- The p300 models were trained on ENCSR000EGE peak summits but evaluated on element centers, which likely understates p300 and so probably widens rather than closes the gap. Does a matched element-centered p300 model change the magnitude?
- Does the wider (+/-1000 bp) counting window raise the *sequence* contribution for H3K27ac, or only the accessibility contribution? This decides whether the wider window helps the syntax question at all.
- Does the pattern hold in a second cell type (GM12878 data already exists for p300), and for other coactivator vs histone-mark pairs?
- Would a sequence-only model trained on the residual of an accessibility-only model isolate the sequence contribution more cleanly than comparing two independently trained models?

---

## F-002: Overall correlation ranks coactivator/histone targets backwards; residual-beyond-accessibility reverses it
**Status:** preliminary
**Claim:** Evaluated identically on K562 DNase candidate elements, an H3K27ac sequence+ATAC model reaches higher OVERALL correlation than the equivalent p300 model (0.809 vs 0.791) but much lower RESIDUAL correlation beyond an ATAC-only baseline (0.475 vs 0.654), and less than half the incremental variance explained (R^2 +0.099 vs +0.265). The ranking of the two targets reverses depending on which metric is used. Separately, an independently-trained sequence-only H3K27ac model captures almost none of the ATAC residual (residual r = 0.100), i.e. it is largely redundant with accessibility rather than complementary to it.
**Implications:** Overall correlation is not a valid basis for choosing a prediction target when accessibility is available as an input, because it is dominated by the dead-vs-active contrast that accessibility already resolves. Any claim of the form "target X is more predictable than target Y" needs an accessibility-only control and a residual metric. For sequence-interpretability work specifically, p300 carries substantially more information that accessibility cannot supply, so it remains the better substrate despite scoring lower on the conventional metric. The sequence-only redundancy result further implies that a jointly-trained or independently-trained sequence model should not be expected to discover the accessibility complement on its own — it has to be trained on the residual explicitly.
**Tags:** h3k27ac, p300, residual, evaluation, sequence-vs-accessibility, metric-choice, interpretability, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-24 | job 40783892 | ENCSR000AKP, ENCSR000EGE, K562 ATAC | 2026_0824_H3K27ac_model | Top-quintile Pearson: p300 ATAC-only 0.305 -> multimodal 0.606 (+0.301); H3K27ac ATAC-only 0.543 -> multimodal 0.668 (+0.125) | supports |
| 2026-08-25 | jobs 40799753, 40803766 | ENCSR000AKP, ENCSR000EGE, K562 ATAC | 2026_0824_H3K27ac_model | Residual r beyond ATAC-only model, 5 folds pooled: p300 multimodal 0.654 (incr R^2 +0.265); H3K27ac multimodal 0.475 (+0.099); H3K27ac sequence-only 0.100 (+/-0.388 worse than ATAC alone). Overall r was 0.791 vs 0.809 in the opposite direction. | refines |

### Open Questions
- p300's residual r is high even where ATAC is nearly right (0.260 in the lowest |residual| quintile vs 0.066 for H3K27ac), so its sequence contribution is broadly useful rather than concentrated at the extremes. Is that a property of coactivator binding, or of p300 having been trained on peak summits?
- Does explicitly training a sequence model on `observed - atac_pred` recover a complement that the independently-trained sequence model missed, for either target?
- Does the H3K27ac residual become more predictable with a wider receptive field, given acetylation domains extend kilobases?
- Does the sequence-only redundancy hold across cell types, or is it a K562-specific artifact of ATAC and H3K27ac both tracking the same K562 activity axis?

---

## F-003: Cross-cell-type H3K27ac transfer is strongly asymmetric; sequence is weak everywhere but its apparent collapse is largely target-cell-type difficulty
**Status:** supported
**Claim:** K562-trained models evaluated on GM12878 candidate elements (top signal quintile) retain very different fractions of their in-cell-type performance: sequence-only falls 0.357 -> 0.153 (43% retained), ATAC-only 0.543 -> 0.473 (87% retained), sequence+ATAC 0.668 -> 0.522 (78% retained). The multimodal advantage over accessibility alone shrinks from +0.125 in K562 to +0.049 in GM12878.
**Implications:** Taken with F-002 (an independently-trained sequence model captures almost none of the ATAC residual, r = 0.100), the sequence contribution to H3K27ac prediction looks largely cell-type-specific rather than a generalizable sequence grammar. Accessibility transfers well because it is a direct measurement of the state rather than an inference from sequence. This bears directly on the stated goal of predicting H3K27ac "in a generalizable way": the current sequence-only model does not meet it, and roughly 60% of the multimodal model's edge over accessibility does not survive a cell-type change. It also means in-cell-type performance is a poor guide to generalization here, so any architecture work should be scored on transfer, not only on held-out chromosomes.
**Tags:** h3k27ac, transferability, gm12878, k562, sequence-vs-accessibility, generalization

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-25 | job 40874283 | ENCSR000AKP (K562), ENCFF645BAL+ENCFF865OOP (GM12878), K562+GM12878 ATAC | 2026_0824_H3K27ac_model | Top-quintile Pearson on GM12878: sequence 0.1533, ATAC-only 0.4727, multimodal 0.5219, vs K562 in-cell-type 0.357 / 0.543 / 0.668 | supports |
| 2026-08-26 | job 40876333 | GM12878 H3K27ac replicates | 2026_0824_H3K27ac_model | GM12878 top-quintile ceiling 0.8321 vs K562 0.7601, so GM12878 is the EASIER target. Normalized as fraction of own ceiling: sequence 47% (K562) -> 18% (GM12878); ATAC-only 71% -> 57%; multimodal 88% -> 63%. Transfer failure is not explained by target difficulty. | refines |
| 2026-08-27 | jobs 41010962, 41011219 | K562 + GM12878 ATAC and H3K27ac | 2026_0824_H3K27ac_model | Model-free ATAC-H3K27ac coupling is lower in GM12878 (top-quintile r 0.409 vs 0.510), which fully accounts for the ATAC-only model's transfer drop (0.542 -> 0.467; model/coupling ratio 1.06 -> 1.14, i.e. no degradation). Sequence-only falls 0.360 -> 0.146 with no such explanation. Refines the claim: accessibility transfers, only sequence collapses. | refines |
| 2026-08-27 | jobs 41025926-41025940, 41107130 | K562 + GM12878 H3K27ac 5-prime, both directions | 2026_0824_H3K27ac_model | RECIPROCAL: GM-trained models score HIGHER on K562 than in GM12878 for every modality, ceiling-normalised (sequence 34% vs 23%, ATAC 58% vs 52%, multimodal 65% vs 61%). Sequence-only retention is 38% for K562->GM12878 but 147% for GM12878->K562. The forward-only reading of a sequence transfer failure does not survive; GM12878 is the harder target. Sequence remains weak in absolute terms in all four settings (15-41% of ceiling), and the sequence margin shrinks on transfer in both directions (+0.138->+0.044 and +0.086->+0.061). | contradicts |

### Open Questions
- ~~The GM12878 ceiling has not been computed~~ **RESOLVED 2026-08-26**: GM12878's ceiling is HIGHER (0.8321 vs 0.7601 top quintile), so it is the easier target and the transfer failure is understated by the raw numbers, not overstated. Remaining nit: the transfer used the fragment-extended GM12878 target while the ceiling is on the 5' one; equivalence is assumed from K562 (0.760 vs 0.761) rather than measured in GM12878.
- GM12878 ATAC is a different experiment at different depth, which could independently depress ATAC-only and multimodal transfer.
- Does the reverse direction (train GM12878, test K562) show the same asymmetry, as it did for p300?
- Would residual training — fitting `observed - atac_pred` explicitly — produce a sequence component that transfers better, given it would be forced onto signal accessibility cannot supply?
- p300 sequence models transfer at 0.277 (from the earlier transferability work). Is the H3K27ac sequence collapse worse than p300's, evaluated identically?
