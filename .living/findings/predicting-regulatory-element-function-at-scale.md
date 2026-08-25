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
