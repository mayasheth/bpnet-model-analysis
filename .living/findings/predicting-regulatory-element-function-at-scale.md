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

## F-001: Sequence adds ~2.4x more predictive information over accessibility for p300 than for H3K27ac
**Status:** preliminary
**Claim:** In K562, on DNase candidate elements, models predicting p300 gain +0.301 Pearson from adding sequence to an accessibility-only model (0.305 -> 0.606 on the top signal quintile), whereas models predicting H3K27ac gain only +0.125 (0.543 -> 0.668). H3K27ac is the more predictable target overall (0.668 vs 0.606 multimodal), but a larger share of its predictability is carried by chromatin accessibility alone.
**Implications:** Target choice matters for sequence-interpretability work independently of how well a model scores. A high-performing H3K27ac model is substantially reading the accessibility track, so attributions from it will reflect less sequence-specific information than the equivalent p300 model, despite the better headline correlation. For motif-syntax questions, p300 remains the better substrate. Also a caution against ranking targets by overall correlation: the accessibility-only control is what makes the comparison interpretable.
**Tags:** h3k27ac, p300, sequence-vs-accessibility, multimodal, interpretability, k562, target-choice

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-24 | job 40783892 | ENCSR000AKP (H3K27ac), ENCSR000EGE (p300), K562 ATAC | 2026_0824_H3K27ac_model | Top-quintile Pearson, 5 folds pooled: p300 ATAC-only 0.305, p300 multimodal 0.606; H3K27ac ATAC-only 0.543, H3K27ac multimodal 0.668, H3K27ac sequence-only 0.357 | supports |

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
