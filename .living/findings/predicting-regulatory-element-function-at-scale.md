---
topic: predicting-regulatory-element-function-at-scale
description: Building and benchmarking models that predict regulatory element function, enhancer–gene links, and variant effects genome-wide (the question behind consortium-scale efforts such as IGVF).
created: 2026-07-21
last_updated: 2026-09-05
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
**Status:** supported
**Claim:** Evaluated identically on K562 DNase candidate elements, an H3K27ac sequence+ATAC model reaches higher OVERALL correlation than the equivalent p300 model (0.809 vs 0.791) but much lower RESIDUAL correlation beyond an ATAC-only baseline (0.475 vs 0.654), and less than half the incremental variance explained (R^2 +0.099 vs +0.265). The ranking of the two targets reverses depending on which metric is used. Separately, an independently-trained sequence-only H3K27ac model captures almost none of the ATAC residual (residual r = 0.100; 0.149 [0.076, 0.222] re-measured per-fold on the 5-prime target). **This is a property of the training objective, not of sequence.** A sequence model trained explicitly on `observed - atac_pred` reaches residual r = 0.459 [0.421, 0.496], and the jointly-trained multimodal model reaches 0.551 [0.510, 0.592]. Sequence therefore carries substantial information that accessibility cannot supply; a model trained to predict total H3K27ac simply spends its capacity on the accessibility-correlated component instead. **Replicated independently in GM12878** (0.133 -> 0.372 for sequence-only; multimodal 0.449): the residual objective GAINS +0.310 (K562) and +0.239 (GM12878) for a sequence-blind input but COSTS -0.0369 and -0.0350 when accessibility is already an input, negative in all ten folds across both cell types. An ATAC-input model trained on the residual -- a negative control asked to predict an ATAC model's own errors from the same ATAC input -- scores -0.003 (K562) and +0.010 (GM12878), both CIs including zero, so the metric is not manufacturing the effect.
**Implications:** Overall correlation is not a valid basis for choosing a prediction target when accessibility is available as an input, because it is dominated by the dead-vs-active contrast that accessibility already resolves. Any claim of the form "target X is more predictable than target Y" needs an accessibility-only control and a residual metric. For sequence-interpretability work specifically, p300 carries substantially more information that accessibility cannot supply, so it remains the better substrate despite scoring lower on the conventional metric. The sequence-only redundancy result was initially read as implying that neither an independently-trained nor a jointly-trained sequence model would discover the accessibility complement on its own. **Half of that is wrong.** It holds for the independently-trained sequence model, but the jointly-trained multimodal model captures MORE of the residual (0.551) than explicit residual training does (0.459) — joint training is the better way to reach the complement, not a failure mode. What genuinely does not work is training sequence against the total signal and expecting the complement to fall out. **The effect is a property of the training objective, not of a cell type**: the multimodal cost replicates to within 0.002 across two cell types that differ in ATAC-H3K27ac coupling (0.510 vs 0.409), inter-replicate ceiling (0.760 vs 0.832), accessibility library and element derivation. Practically: use joint multimodal training for prediction, and reserve residual training for attribution work, where forcing the model onto accessibility-independent signal is the goal rather than a cost. It should not need re-testing per cell type.
**Tags:** h3k27ac, p300, residual, evaluation, sequence-vs-accessibility, metric-choice, interpretability, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-24 | job 40783892 | ENCSR000AKP, ENCSR000EGE, K562 ATAC | 2026_0824_H3K27ac_model | Top-quintile Pearson: p300 ATAC-only 0.305 -> multimodal 0.606 (+0.301); H3K27ac ATAC-only 0.543 -> multimodal 0.668 (+0.125) | supports |
| 2026-08-25 | jobs 40799753, 40803766 | ENCSR000AKP, ENCSR000EGE, K562 ATAC | 2026_0824_H3K27ac_model | Residual r beyond ATAC-only model, 5 folds pooled: p300 multimodal 0.654 (incr R^2 +0.265); H3K27ac multimodal 0.475 (+0.099); H3K27ac sequence-only 0.100 (+/-0.388 worse than ATAC alone). Overall r was 0.791 vs 0.809 in the opposite direction. | refines |
| 2026-08-29 | jobs 41232315, 41233130 | ENCSR000AKP 5-prime, K562 ATAC | 2026_0824_H3K27ac_model | Residual models scored at last (trained 2026-08-27, never evaluated). Per-fold mean [95% CI] on the 5-prime target: residual r sequence-only 0.149 [0.076, 0.222], EXPLICIT RESIDUAL 0.459 [0.421, 0.496], multimodal 0.551 [0.510, 0.592]. Incremental R^2 over ATAC-only: -0.418, +0.071 [0.058, 0.085], +0.102 [0.087, 0.117]; on the top signal quintile -0.156, +0.115 [0.075, 0.154], +0.169 [0.131, 0.207]. Artifact ruled out: r(residual output, atac_pred) = -0.090 [-0.104, -0.076] and partialling atac_pred out of both sides leaves 0.459 unchanged. Sequence carries a real accessibility-independent component; the independently-trained model just does not surface it. Required fixing 2.4, which double-subtracted atac_pred from offset-trained models. | contradicts |
| 2026-09-01 | jobs 41427509-41427523, 41514616 | ENCSR000AKC 5-prime, GM12878 ATAC | 2026_0824_H3K27ac_model | REPLICATION in GM12878. Residual r: sequence-on-total 0.133 [0.090, 0.176], sequence-on-residual 0.372 [0.334, 0.410], multimodal-on-total 0.449 [0.414, 0.485], multimodal-on-residual 0.414 [0.378, 0.451], ATAC control 0.010 [-0.019, 0.040]. Paired: +0.239 [0.220, 0.258] p<1e-4 for sequence, -0.0350 [-0.0440, -0.0261] p=4e-4 for multimodal. Matches K562 (+0.310, -0.0369) in sign and, for the multimodal cost, to within 0.002. | supports |

### Open Questions
- p300's residual r is high even where ATAC is nearly right (0.260 in the lowest |residual| quintile vs 0.066 for H3K27ac), so its sequence contribution is broadly useful rather than concentrated at the extremes. Is that a property of coactivator binding, or of p300 having been trained on peak summits?
- ~~Does explicitly training a sequence model on `observed - atac_pred` recover a complement that the independently-trained sequence model missed?~~ **RESOLVED 2026-08-29 for H3K27ac: yes, decisively.** Residual r 0.459 [0.421, 0.496] vs 0.149 [0.076, 0.222] for the independently-trained sequence model; non-overlapping CIs. Ruled out the metric's built-in artifact: because `true_resid = obs - atac_pred`, anything anti-correlated with `atac_pred` scores positive for free, but the residual model's output is only -0.090 [-0.104, -0.076] correlated with `atac_pred`, and partialling `atac_pred` out of both sides leaves the score unchanged (0.459 -> 0.459). Confirmed against the observed signal too: adding its output to `atac_pred` lifts held-out r from 0.815 to 0.856 [0.840, 0.873], incremental R^2 +0.071 [0.058, 0.085] (+0.115 [0.075, 0.154] on the top signal quintile). Still open for p300.
- ~~Given the multimodal model beats explicit residual training, is residual training useful at all here?~~ **RESOLVED 2026-09-01**: its value is interpretability, not prediction. The multimodal model beats it in both cell types (+0.093 K562, +0.077 GM12878), and the residual objective actively costs a multimodal model ~0.035. Use joint training to predict, residual training to attribute.
- Does the residual model TRANSFER better than the multimodal one? The residual design never transfers the ATAC->H3K27ac mapping -- accessibility is re-measured locally and only the sequence complement moves -- so it may retain more across cell types. Running as of 2026-09-01 in two variants: an upper bound using the target's own ATAC model as offset, and the deployable case where that offset model is also transferred (in the real application there is ATAC but no H3K27ac in the target).
- Does the H3K27ac residual become more predictable with a wider receptive field, given acetylation domains extend kilobases?
- Does the sequence-only redundancy hold across cell types, or is it a K562-specific artifact of ATAC and H3K27ac both tracking the same K562 activity axis?

---

## F-003: Cross-cell-type H3K27ac transfer is strongly asymmetric; sequence is weak everywhere but its apparent collapse is largely target-cell-type difficulty
**Status:** supported
**Claim:** K562-trained models evaluated on GM12878 candidate elements (top signal quintile) retain very different fractions of their in-cell-type performance: sequence-only falls 0.357 -> 0.153 (43% retained), ATAC-only 0.543 -> 0.473 (87% retained), sequence+ATAC 0.668 -> 0.522 (78% retained). The multimodal advantage over accessibility alone shrinks from +0.125 in K562 to +0.049 in GM12878.
**Implications:** ~~Taken with F-002 (an independently-trained sequence model captures almost none of the ATAC residual, r = 0.100), the sequence contribution looks largely cell-type-specific rather than a generalizable sequence grammar.~~ **That cross-reference is stale — F-002 has since been corrected.** The r = 0.100 figure reflects the training objective, not the information in sequence: trained on the residual directly, the same sequence input reaches 0.459 (K562) and 0.372 (GM12878). And on transfer (2026-09-01), every sequence-containing model beats a transferred ATAC-only model by 0.04-0.06 on the top quintile in both directions, so the sequence contribution does survive a cell-type change. What remains true is that sequence is weak in absolute terms everywhere and that sequence-ONLY transfer is useless (top-quintile 0.147 and 0.317, against 0.477 and 0.541 for transferred ATAC-only): accessibility must be measured in the target cell type. Accessibility transfers well because it is a direct measurement of the state rather than an inference from sequence. This bears directly on the stated goal of predicting H3K27ac "in a generalizable way": a sequence-only model does not meet it, and part of the multimodal model's edge over accessibility does not survive a cell-type change — though the edge that does survive is substantial (0.04-0.06 top-quintile over transferred ATAC-only, both directions). It also means in-cell-type performance is a poor guide to generalization here, so any architecture work should be scored on transfer, not only on held-out chromosomes.
**Tags:** h3k27ac, transferability, gm12878, k562, sequence-vs-accessibility, generalization

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-25 | job 40874283 | ENCSR000AKP (K562), ENCFF645BAL+ENCFF865OOP (GM12878), K562+GM12878 ATAC | 2026_0824_H3K27ac_model | Top-quintile Pearson on GM12878: sequence 0.1533, ATAC-only 0.4727, multimodal 0.5219, vs K562 in-cell-type 0.357 / 0.543 / 0.668 | supports |
| 2026-08-26 | job 40876333 | GM12878 H3K27ac replicates | 2026_0824_H3K27ac_model | GM12878 top-quintile ceiling 0.8321 vs K562 0.7601, so GM12878 is the EASIER target. Normalized as fraction of own ceiling: sequence 47% (K562) -> 18% (GM12878); ATAC-only 71% -> 57%; multimodal 88% -> 63%. Transfer failure is not explained by target difficulty. | refines |
| 2026-08-27 | jobs 41010962, 41011219 | K562 + GM12878 ATAC and H3K27ac | 2026_0824_H3K27ac_model | Model-free ATAC-H3K27ac coupling is lower in GM12878 (top-quintile r 0.409 vs 0.510), which fully accounts for the ATAC-only model's transfer drop (0.542 -> 0.467; model/coupling ratio 1.06 -> 1.14, i.e. no degradation). Sequence-only falls 0.360 -> 0.146 with no such explanation. Refines the claim: accessibility transfers, only sequence collapses. | refines |
| 2026-08-27 | jobs 41025926-41025940, 41107130 | K562 + GM12878 H3K27ac 5-prime, both directions | 2026_0824_H3K27ac_model | RECIPROCAL: GM-trained models score HIGHER on K562 than in GM12878 for every modality, ceiling-normalised (sequence 34% vs 23%, ATAC 58% vs 52%, multimodal 65% vs 61%). Sequence-only retention is 38% for K562->GM12878 but 147% for GM12878->K562. The forward-only reading of a sequence transfer failure does not survive; GM12878 is the harder target. Sequence remains weak in absolute terms in all four settings (15-41% of ceiling), and the sequence margin shrinks on transfer in both directions (+0.138->+0.044 and +0.086->+0.061). | contradicts |
| 2026-09-01 | jobs 41559350, 41559736, 41562791 | K562 + GM12878 5-prime, both directions | 2026_0824_H3K27ac_model | DEPLOYMENT transfer, top signal quintile, everything transferred and only target ATAC used: K562->GM12878 ATAC-only 0.477, sequence-only 0.147, multimodal 0.520, residual-multimodal 0.512; GM12878->K562 0.541 / 0.317 / 0.602 / 0.600. Sequence adds 0.04-0.06 over transferred ATAC-only in both directions. Residual and multimodal are indistinguishable on the top quintile (paired p = 0.17-0.76); the all-elements difference is significant and sign-flipping by variant but lives in the dead-vs-active contrast. Transferring the ATAC baseline rather than fitting it locally costs only ~0.017. | refines |

### Open Questions
- ~~The GM12878 ceiling has not been computed~~ **RESOLVED 2026-08-26**: GM12878's ceiling is HIGHER (0.8321 vs 0.7601 top quintile), so it is the easier target and the transfer failure is understated by the raw numbers, not overstated. Remaining nit: the transfer used the fragment-extended GM12878 target while the ceiling is on the 5' one; equivalence is assumed from K562 (0.760 vs 0.761) rather than measured in GM12878.
- GM12878 ATAC is a different experiment at different depth, which could independently depress ATAC-only and multimodal transfer.
- Does the reverse direction (train GM12878, test K562) show the same asymmetry, as it did for p300?
- ~~Would residual training produce a sequence component that transfers better?~~ **RESOLVED 2026-09-01: no, not usefully.** On the top signal quintile the residual and multimodal models transfer indistinguishably in both directions (paired p = 0.17-0.76). The residual design does win when the accessibility baseline is fitted in the TARGET cell type rather than transferred (+0.005 to +0.011 on all elements, p < 0.01 both directions), but that requires target H3K27ac, which the deployment scenario lacks. Practical consequence: transfer the multimodal model. Secondary consequence worth acting on: a shallow target H3K27ac library, enough only to fit an ATAC-only model, would make the residual architecture preferable.
- p300 sequence models transfer at 0.277 (from the earlier transferability work). Is the H3K27ac sequence collapse worse than p300's, evaluated identically?

---

## F-004: Predicted H3K27ac does not improve ABC's CRISPR-benchmark performance, and the failure is compressed dynamic range rather than inaccurate prediction
**Status:** established
**Claim:** In K562, substituting model-predicted H3K27ac for observed H3K27ac in ABC's `activity_base` gives a CRISPR-benchmark AUPRC of 0.482 [0.437, 0.530] at best, against a floor of 0.457 (ATAC-only activity) and a ceiling of 0.519 (ATAC x observed H3K27ac) — no resolvable improvement over the floor. The two deployment-realistic arms, GM12878-trained models applied to K562, score 0.455 and 0.452, at or below the floor. Predicted and observed H3K27ac nonetheless agree at Spearman 0.794 over all 153,545 candidate regions, so the predictions are not inaccurate; their dynamic range is compressed exactly where the benchmark is decided. On the regions carrying a CRISPR-regulated pair the predicted activity term's p99/p50 ratio is 4.03 against 6.00 observed, its top-decile mean/median ratio 4.23 against 6.46, and agreement with observed falls to Spearman 0.663.
**Implications:** Correlation with observed H3K27ac is a poor proxy for downstream utility, and the two can be improved independently: every architecture change that raised top-quintile Pearson left this benchmark unmoved. The actionable target is the *spread* of the predicted activity term on strongly acetylated elements, which no architecture change tried so far addresses. Note also that the ceiling is only 0.062 above the floor — observed H3K27ac itself buys ABC very little in K562 — so this experiment had limited room from the start, and a negative result here bounds the value of the whole predicted-activity idea rather than only of this model.
**Update 2026-09-06:** this finding is specific to the H3K27ac target. F-008 shows a
predicted *p300* track clears the same floor by +0.055 [+0.034, +0.074] on a paired
bootstrap. Note also that F-004's conclusion rests on unpaired per-predictor CIs, which
the p300 run demonstrates cannot resolve differences of this size; the H3K27ac deltas
should be re-tested with the paired bootstrap in `scripts/4.17`.

**Tags:** h3k27ac, abc, crispr-benchmark, dynamic-range, downstream-evaluation, k562, negative-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-04 | job 42158462 | EPCrisprBenchmark_ensemble_data_GRCh38, scE2G intGENCODEv43 universes | 2026_0824_H3K27ac_model | AUPRC: ceiling 0.519, best predicted 0.482, floor 0.457, GM12878-trained deployment arms 0.455/0.452, sequence-only arms 0.393/0.275, distance-to-TSS 0.435 | supports |
| 2026-09-04 | job 42173879 | 153,545 ABC candidate regions | 2026_0824_H3K27ac_model | Spearman(obs, pred) 0.794 all regions / 0.663 on regions of regulated pairs; activity p99/p50 4.03 predicted against 6.00 observed; 12 of 429 regulated pairs caught by observed and missed by predicted | supports |

---

## F-005: In-cell-type architecture gains for H3K27ac do not transfer, and the better in-cell-type models have the larger transfer drops
**Status:** established
**Claim:** Across a 2x2 of receptive field (~1.1 vs ~4.2 kb) and accessibility representation (flat vs five fragment-size channels), no architecture beats the simplest by a resolvable margin when transferred between K562 and GM12878. Paired within-fold differences in top-quintile Pearson against narrow+flat transferred: wide+flat +0.0115 (*p*=0.14) and +0.0053 (*p*=0.53); fragment channels **-0.0012** (*p*=0.86) and **-0.0029** (*p*=0.75); wide+fragments +0.0159 (*p*=0.066). Fragment channels gain +0.016 (*p*=0.002) in-cell type in K562 and are exactly null transferred in both directions. The transfer drop is larger for the richer architectures: local minus transferred is +0.055/+0.085 for narrow+flat against +0.059/+0.107 for wide+flat and +0.061/+0.103 for narrow+fragments (all *p*<0.006). Transferring the best model buys only +0.013 to +0.036 over using the target cell type's own ATAC-only model (0.518 and 0.584), against locally trained models reaching 0.602 and 0.733.
**Implications:** In-cell-type ranking is not the deployment ranking, so any architecture change intended for cross-cell-type use must be scored transferred before adoption. Fragment-length structure in particular appears to be a property of a specific ATAC library as much as of chromatin, which no paired in-cell-type test can detect. On the mechanistic metric the shortfall is larger still: transferred `residual_pearson` is 0.212-0.262 where a local model reaches 0.406-0.481, so roughly half of what a model extracts beyond accessibility does not survive the move.
**Tags:** h3k27ac, transfer, deployment, architecture, fragment-channels, receptive-field, gm12878, k562, negative-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-05 | job 42179731 | K562 ENCSR000AKP + GM12878 ENCSR000AKC H3K27ac, 5 folds each direction | 2026_0824_H3K27ac_model | Top-quintile Pearson transferred: 0.531/0.542/0.530/0.547 (K562->GM) and 0.614/0.620/0.611 (GM->K562) against target ATAC-only floors of 0.518 and 0.584 | supports |

---

## F-006: p300 is not predictable at the elements where an H3K27ac model fails, so it is not usable as an auxiliary target
**Status:** established
**Claim:** At the elements a K562 H3K27ac multimodal model most under-predicts, observed EP300 ChIP signal is elevated 4.89x over its genome-wide median but its own **input control** is elevated 2.10x, so real enrichment is about 2.3x rather than the 6.6x that peak overlap implies. Existing p300-target models predict 1.83x (multimodal) and 1.55x (ATAC-only) there — below the input control. The same p300 models also over-predict the accessible-but-unmarked tail (1.92x and 1.79x where observed p300 is 1.36x and its control 0.50x), i.e. they fail in the same direction on the same elements.
**Implications:** Neither a multi-task model with H3K27ac and p300 heads nor stacking predicted p300 as an input feature can recover these elements, because the auxiliary target is not learnable exactly where it would need to be. More generally, a peak-overlap enrichment is not sufficient evidence that an auxiliary target carries usable signal: without the ChIP input control the premise here looked roughly 3x stronger than it is. This does not affect F-001 — p300 remains the better substrate for sequence *attribution*, which is a claim about where motif work should be done, not about what to predict.
**Tags:** p300, h3k27ac, multi-task, auxiliary-target, input-control, negative-result, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-05 | job 42182188 | ENCSR000EGE p300 + control, 153,545 ABC candidate regions | 2026_0824_H3K27ac_model | Fold elevation at H3K27ac under-predicted 1%: observed H3K27ac 31.33, predicted 5.81, observed p300 4.89, p300 input control 2.10, predicted p300 1.83 (multimodal) / 1.55 (ATAC-only) | refutes |

---

## F-007: The elements an H3K27ac model over-predicts are accessible-but-unacetylated, and CTCF binding is a passenger rather than the cause
**Status:** established
**Claim:** In K562, over the 153,545 ABC candidate regions, the elements a sequence + ATAC H3K27ac model most over-predicts are accessible (ATAC 9.9 against 0.9 RPM typical), GC-rich (0.593 against 0.490), CpG-rich (o/e 0.551 against 0.277) and promoter-enriched (31.6% against 13.1%), and carry no more H3K27ac than a random element (0.94x typical RPKM) with IgG *below* background (0.33 against 0.48). Pooled CTCF signal is enriched 2.6x there, but this is a minority effect and is not causal: the tail's per-element CTCF median is 0.72 RPKM against 0.55 typical, only 15.7% of the 1% tail (26.5% of the 5% tail) exceeds the 90th percentile of typical-stratum CTCF, and splitting the tail at that threshold leaves the CTCF-low 73.5% over-predicted *more* than the CTCF-high 26.5% (median residualised error 0.821 against 0.759; predicted fold-elevation 3.83 against 1.63). Within accessibility deciles, Spearman rho(error, GC) is 0.24-0.37 in every decile against rho(error, CTCF) of 0.02-0.23 (median 0.084). The opposite tail is unambiguous: EP300 3.4x, H3K4me1 3.3x, H3K27ac 14x, H3K27me3 depleted to 0.31x and no element overlapping an H3K27me3 peak.
**Implications:** The failure mode is "accessible but unacetylated" as a general phenotype, of which CTCF sites are one instance covering about a quarter. Interventions aimed specifically at CTCF -- motif channels for CTCF, insulator annotations -- would address a minority of the problem. GC content being the stronger correlate points instead at the training-element composition: the negative pool is a uniform genome sample at mean GC 0.389 against 0.466-0.593 for candidate elements, so a sequence branch can satisfy much of its training objective with a GC detector, which is consistent with both this correlation and its measured lack of dynamic range within candidate elements.
**Tags:** h3k27ac, error-analysis, ctcf, gc-content, attributable-fraction, training-composition, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-05 | jobs 42179002, 42179552, 42181781 | K562 CTCF/EP300/H3K4me1/H3K27me3/H3K27ac/IgG BAMs and peak calls | 2026_0824_H3K27ac_model | Stratum-level RPKM: over-predicted 5% CTCF 4.69 against 1.79 typical, IgG 0.39 against 0.48; under-predicted 1% EP300 2.90, H3K4me1 4.35, H3K27ac 17.52, H3K27me3 0.31 | supports |
| 2026-09-05 | job 42221086 | per-element CTCF/IgG counts over 153,545 ABC regions | 2026_0824_H3K27ac_model | CTCF-high = 26.5% of the over-predicted 5%; median error 0.759 CTCF-high against 0.821 CTCF-low; rho(err,CTCF) median 0.084 within ATAC deciles against rho(err,GC) 0.24-0.37 | refutes |

---

## F-008: A sequence + ATAC model predicting p300 improves ABC's CRISPR-benchmark performance, and observed p300 is a better activity term than observed H3K27ac
**Status:** established
**Claim:** In K562, substituting model-predicted p300 for the H3K27ac term in ABC's `activity_base` beats the ATAC-only floor by **+0.055 AUPRC [+0.034, +0.074]** on a paired bootstrap over the shared 10,342-pair benchmark set (466 regulated), with the sign preserved in 100% of 2,000 resamples, and is statistically indistinguishable from *measured* H3K27ac (-0.007 [-0.032, +0.017]). Predicted p300 used as the entire activity term also clears the floor (+0.045 [+0.012, +0.077], 99.6%). Separately and independently of any model, **observed p300 beats observed H3K27ac** as the activity term by **+0.038 [+0.019, +0.058]** (99.9%), widening the headroom over the floor from 0.062 to 0.099. The gain requires sequence: the multimodal p300 model beats the ATAC-only p300 model by +0.051 [+0.032, +0.069] (100%), while the ATAC-only p300 model sits at the floor (+0.004 [-0.012, +0.021]).
**Implications:** This is the first predicted-activity arm in the project to clear the floor with a resolvable margin, and it reverses the headline negative from F-004, which was specific to the H3K27ac target. Two separable causes contribute and should not be conflated: p300 is a better activity concept (more headroom), and the model captures more of the headroom it is given (55% against a non-resolvable fraction for H3K27ac). The sequence-branch result is the mirror image of H3K27ac, where sequence-containing arms scored 0.393 and 0.275 against a 0.457 floor — same architecture, same accessibility input, same benchmark, so **target choice rather than architecture was the binding constraint**. Methodologically: the unpaired per-predictor CIs the pipeline emits (~±0.05) resolve none of these differences and every arm's interval overlaps every other's; only the paired bootstrap over the shared pair set resolves them. The H3K27ac conclusion in F-004 was drawn from unpaired intervals and should be re-tested paired.
**Caveats:** K562-trained and K562-evaluated, so this is not yet a transferability result; no GM12878 p300 model exists. The mechanism is inferred rather than measured — whether predicted p300 succeeds by being less compressed than predicted H3K27ac has not been tested with the `4.9` dynamic-range diagnostic. Does not contradict F-006, which concerns p300's learnability at the specific elements where H3K27ac prediction fails and correctly killed the multi-head and stacking designs; p300 as the *primary* target is a different question.
**Tags:** p300, abc, crispr-benchmark, target-choice, paired-bootstrap, k562, positive-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-06 | ABC 42221146, benchmark 42233184 | EPCrisprBenchmark_ensemble_data_GRCh38, scE2G intGENCODEv43 universes | 2026_0824_H3K27ac_model | AUPRC: observed p300 0.557, observed H3K27ac 0.519, predicted p300 (seq+ATAC) 0.512, predicted p300 alone 0.502, predicted p300 (ATAC-only) 0.461, floor 0.457, distToTSS 0.435 | supports |
| 2026-09-06 | job 42233184 + `4.17` paired bootstrap, 2,000 resamples | shared 10,342 element-gene pairs | 2026_0824_H3K27ac_model | Paired deltas: obs p300 - obs H3K27ac +0.038 [+0.019, +0.058]; pred p300 - floor +0.055 [+0.034, +0.074]; pred p300 - obs H3K27ac -0.007 [-0.032, +0.017]; multimodal - ATAC-only p300 +0.051 [+0.032, +0.069] | supports |
| 2026-09-06 | cross-run anchor check | same two anchor arms in both pipeline runs | 2026_0824_H3K27ac_model | Floor 0.4573 vs 0.4572 and observed H3K27ac 0.5192 vs 0.5192 across independent runs, licensing the 0.512-against-0.482 comparison | supports |

**Update 2026-09-06 (same day): the gain does not transfer, and this finding is in-cell-type only.**
A GM12878-trained multimodal p300 model applied to K562 scores **+0.009 AUPRC [-0.004, +0.022]**
over the ATAC-only floor on the same paired bootstrap and the same 10,342-pair set — an interval
spanning zero, i.e. indistinguishable from using the target cell type's own ATAC. The transfer
drop is **-0.046 [-0.061, -0.030]** (sign kept 100%), removing 83% of the in-cell-type gain, and
the transferred arm sits **-0.053 [-0.074, -0.032]** below measured H3K27ac. Predicted-alone
transferred is +0.005 [-0.019, +0.028].

Deployment therefore fails for **both** targets — H3K27ac deployment arms 0.455 and 0.452, p300
0.466, against a 0.457 floor — so the barrier is transfer rather than target choice, which is the
same conclusion F-005 reached for architecture changes by an independent route. The claims that
survive are: observed p300 beats observed H3K27ac as an activity term (+0.038, a fact about the
assay), and a K562-trained p300 model matches measured H3K27ac *in K562*.

**Missing control.** The drop conflates transfer with the possibility that the GM12878 p300 model
is weaker in absolute terms: different experiment (ENCSR000DZG against ENCSR000EGE), different
project, never scored in its own cell type. The local-versus-transferred pairing used throughout
the H3K27ac transfer matrix is required before "transfer destroys the gain" can be claimed rather
than "the gain does not transfer".

| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-06 | ABC 42269593, benchmark 42275735, bootstrap `4.17` | 9 arms on one shared 10,342-pair set | 2026_0824_H3K27ac_model | transferred - floor +0.009 [-0.004, +0.022]; local - transferred +0.046 [+0.030, +0.061]; transferred - observed H3K27ac -0.053 [-0.074, -0.032] | refutes |
