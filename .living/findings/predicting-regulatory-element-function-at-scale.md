---
topic: predicting-regulatory-element-function-at-scale
description: Building and benchmarking models that predict regulatory element function, enhancer-gene links, and variant effects genome-wide (the question behind consortium-scale efforts such as IGVF).
created: 2026-07-21
last_updated: 2026-09-05
status: active
---

# Predicting regulatory element function at scale

_Seed topic (Engreitz Lab). No findings recorded yet, crystallize-findings will
append `F-NNN` entries here as analyses produce them. (Broad question behind
IGVF-style consortium work, the slug avoids naming the consortium.)_

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
**Implications:** Overall correlation is not a valid basis for choosing a prediction target when accessibility is available as an input, because it is dominated by the dead-vs-active contrast that accessibility already resolves. Any claim of the form "target X is more predictable than target Y" needs an accessibility-only control and a residual metric. For sequence-interpretability work specifically, p300 carries substantially more information that accessibility cannot supply, so it remains the better substrate despite scoring lower on the conventional metric. The sequence-only redundancy result was initially read as implying that neither an independently-trained nor a jointly-trained sequence model would discover the accessibility complement on its own. **Half of that is wrong.** It holds for the independently-trained sequence model, but the jointly-trained multimodal model captures MORE of the residual (0.551) than explicit residual training does (0.459), joint training is the better way to reach the complement, not a failure mode. What genuinely does not work is training sequence against the total signal and expecting the complement to fall out. **The effect is a property of the training objective, not of a cell type**: the multimodal cost replicates to within 0.002 across two cell types that differ in ATAC-H3K27ac coupling (0.510 vs 0.409), inter-replicate ceiling (0.760 vs 0.832), accessibility library and element derivation. Practically: use joint multimodal training for prediction, and reserve residual training for attribution work, where forcing the model onto accessibility-independent signal is the goal rather than a cost. It should not need re-testing per cell type.
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
**Implications:** ~~Taken with F-002 (an independently-trained sequence model captures almost none of the ATAC residual, r = 0.100), the sequence contribution looks largely cell-type-specific rather than a generalizable sequence grammar.~~ **That cross-reference is stale, F-002 has since been corrected.** The r = 0.100 figure reflects the training objective, not the information in sequence: trained on the residual directly, the same sequence input reaches 0.459 (K562) and 0.372 (GM12878). And on transfer (2026-09-01), every sequence-containing model beats a transferred ATAC-only model by 0.04-0.06 on the top quintile in both directions, so the sequence contribution does survive a cell-type change. What remains true is that sequence is weak in absolute terms everywhere and that sequence-ONLY transfer is useless (top-quintile 0.147 and 0.317, against 0.477 and 0.541 for transferred ATAC-only): accessibility must be measured in the target cell type. Accessibility transfers well because it is a direct measurement of the state rather than an inference from sequence. This bears directly on the stated goal of predicting H3K27ac "in a generalizable way": a sequence-only model does not meet it, and part of the multimodal model's edge over accessibility does not survive a cell-type change, though the edge that does survive is substantial (0.04-0.06 top-quintile over transferred ATAC-only, both directions). It also means in-cell-type performance is a poor guide to generalization here, so any architecture work should be scored on transfer, not only on held-out chromosomes.
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

## F-004: Predicted H3K27ac from an ATAC-input model does not improve ABC's CRISPR-benchmark performance, and the failure is compressed dynamic range rather than inaccurate prediction
**Status:** established for ATAC-input models; NOT general, see F-018 (2026-09-15); transfer arms confirmed at the floor by paired test (2026-09-17)
**Claim:** In K562, substituting model-predicted H3K27ac for observed H3K27ac in ABC's `activity_base` gives a CRISPR-benchmark AUPRC of 0.482 [0.437, 0.530] at best, against a floor of 0.457 (ATAC-only activity) and a ceiling of 0.519 (ATAC x observed H3K27ac), no resolvable improvement over the floor. The two deployment-realistic arms, GM12878-trained models applied to K562, score 0.455 and 0.452, at or below the floor. Predicted and observed H3K27ac nonetheless agree at Spearman 0.794 over all 153,545 candidate regions, so the predictions are not inaccurate; their dynamic range is compressed exactly where the benchmark is decided. On the regions carrying a CRISPR-regulated pair the predicted activity term's p99/p50 ratio is 4.03 against 6.00 observed, its top-decile mean/median ratio 4.23 against 6.46, and agreement with observed falls to Spearman 0.663.
**Implications:** Correlation with observed H3K27ac is a poor proxy for downstream utility, and the two can be improved independently: every architecture change that raised top-quintile Pearson left this benchmark unmoved. The actionable target is the *spread* of the predicted activity term on strongly acetylated elements, which no architecture change tried so far addresses. Note also that the ceiling is only 0.062 above the floor, observed H3K27ac itself buys ABC very little in K562, so this experiment had limited room from the start, and a negative result here bounds the value of the whole predicted-activity idea rather than only of this model.
**Update 2026-09-06:** this finding is specific to the H3K27ac target. F-008 shows a
predicted *p300* track clears the same floor by +0.055 [+0.034, +0.074] on a paired
bootstrap. Note also that F-004's conclusion rests on unpaired per-predictor CIs, which
the p300 run demonstrates cannot resolve differences of this size; the H3K27ac deltas
should be re-tested with the paired bootstrap in `scripts/4.17`.

**Update 2026-09-15: the paired re-test was run and this finding SURVIVES for ATAC-input
models, but it does NOT generalise to the model's accessibility input, see F-018.** Paired
on 10,342 shared element-gene pairs, the ATAC-input multimodal arm beats the floor by only
+0.0121 [-0.0053, +0.0285] (sign kept 91%) in the geomean form and +0.0248 [-0.0034,
+0.0518] (95.3%) as activity alone, both spanning zero. So the negative result was not an
artefact of unpaired intervals. What does clear the floor is the same architecture retrained
with DNase as its accessibility input: +0.0403 [+0.0230, +0.0563] with ATAC in the activity
slot and +0.0839 [+0.0613, +0.1048] with DNase in it (F-018). F-004's scope is therefore
"predicted H3K27ac from an ATAC-input model", not "predicted H3K27ac".


**Update 2026-09-17: the two transfer arms were re-tested with the paired bootstrap, and they
are genuinely AT the floor, not merely unresolvable.** The original 0.455 and 0.452 were judged
on unpaired per-predictor CIs, which F-008 showed cannot resolve effects of this size. Paired on
the same 10,342 element-gene pairs, GM12878-trained applied to K562 gives **-0.0028 [-0.0132,
+0.0078] (sign kept 73.2%)** for the ATAC-input model and **-0.0048 [-0.0216, +0.0110] (72.7%)**
for the multimodal one, against the ATAC-only floor. Pairing shrank the intervals from about
+/-0.05 to +/-0.01 and they still straddle zero, so the arms sit on the floor rather than above
or resolvably below it.
**What pairing DOES resolve is the transfer penalty itself**, which the unpaired intervals could
not: in-cell minus transferred is **+0.0169 [+0.0041, +0.0297] (99.6%)** for the multimodal model
and +0.0053 [-0.0013, +0.0119] (94.5%) for the ATAC-input one. So transfer costs the multimodal
model a resolvable amount downstream, and that is the first paired downstream confirmation of
the transfer failure F-005, F-016 and F-017 found upstream. `results/f004_transfer_paired.tsv`.
Note the AUPRC values on the shared pair set differ from `performance_summary.txt` by up to
0.011 (floor 0.4680 here against 0.4573 there); `4.17` recomputes average precision on the
common set by design, so use it for deltas and the pipeline's own table for absolute values.

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
**Claim:** At the elements a K562 H3K27ac multimodal model most under-predicts, observed EP300 ChIP signal is elevated 4.89x over its genome-wide median but its own **input control** is elevated 2.10x, so real enrichment is about 2.3x rather than the 6.6x that peak overlap implies. Existing p300-target models predict 1.83x (multimodal) and 1.55x (ATAC-only) there, below the input control. The same p300 models also over-predict the accessible-but-unmarked tail (1.92x and 1.79x where observed p300 is 1.36x and its control 0.50x), i.e. they fail in the same direction on the same elements.
**Implications:** Neither a multi-task model with H3K27ac and p300 heads nor stacking predicted p300 as an input feature can recover these elements, because the auxiliary target is not learnable exactly where it would need to be. More generally, a peak-overlap enrichment is not sufficient evidence that an auxiliary target carries usable signal: without the ChIP input control the premise here looked roughly 3x stronger than it is. This does not affect F-001, p300 remains the better substrate for sequence *attribution*, which is a claim about where motif work should be done, not about what to predict.
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
**Implications:** This is the first predicted-activity arm in the project to clear the floor with a resolvable margin, and it reverses the headline negative from F-004, which was specific to the H3K27ac target. Two separable causes contribute and should not be conflated: p300 is a better activity concept (more headroom), and the model captures more of the headroom it is given (55% against a non-resolvable fraction for H3K27ac). The sequence-branch result is the mirror image of H3K27ac, where sequence-containing arms scored 0.393 and 0.275 against a 0.457 floor, same architecture, same accessibility input, same benchmark, so **target choice rather than architecture was the binding constraint**. Methodologically: the unpaired per-predictor CIs the pipeline emits (~+/-0.05) resolve none of these differences and every arm's interval overlaps every other's; only the paired bootstrap over the shared pair set resolves them. The H3K27ac conclusion in F-004 was drawn from unpaired intervals and should be re-tested paired.
**Caveats:** K562-trained and K562-evaluated, so this is not yet a transferability result; no GM12878 p300 model exists. The mechanism is inferred rather than measured, whether predicted p300 succeeds by being less compressed than predicted H3K27ac has not been tested with the `4.9` dynamic-range diagnostic. Does not contradict F-006, which concerns p300's learnability at the specific elements where H3K27ac prediction fails and correctly killed the multi-head and stacking designs; p300 as the *primary* target is a different question.
**Tags:** p300, abc, crispr-benchmark, target-choice, paired-bootstrap, k562, positive-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-06 | ABC 42221146, benchmark 42233184 | EPCrisprBenchmark_ensemble_data_GRCh38, scE2G intGENCODEv43 universes | 2026_0824_H3K27ac_model | AUPRC: observed p300 0.557, observed H3K27ac 0.519, predicted p300 (seq+ATAC) 0.512, predicted p300 alone 0.502, predicted p300 (ATAC-only) 0.461, floor 0.457, distToTSS 0.435 | supports |
| 2026-09-06 | job 42233184 + `4.17` paired bootstrap, 2,000 resamples | shared 10,342 element-gene pairs | 2026_0824_H3K27ac_model | Paired deltas: obs p300 - obs H3K27ac +0.038 [+0.019, +0.058]; pred p300 - floor +0.055 [+0.034, +0.074]; pred p300 - obs H3K27ac -0.007 [-0.032, +0.017]; multimodal - ATAC-only p300 +0.051 [+0.032, +0.069] | supports |
| 2026-09-06 | cross-run anchor check | same two anchor arms in both pipeline runs | 2026_0824_H3K27ac_model | Floor 0.4573 vs 0.4572 and observed H3K27ac 0.5192 vs 0.5192 across independent runs, licensing the 0.512-against-0.482 comparison | supports |

**Update 2026-09-06 (same day): the gain does not transfer, and this finding is in-cell-type only.**
A GM12878-trained multimodal p300 model applied to K562 scores **+0.009 AUPRC [-0.004, +0.022]**
over the ATAC-only floor on the same paired bootstrap and the same 10,342-pair set, an interval
spanning zero, i.e. indistinguishable from using the target cell type's own ATAC. The transfer
drop is **-0.046 [-0.061, -0.030]** (sign kept 100%), removing 83% of the in-cell-type gain, and
the transferred arm sits **-0.053 [-0.074, -0.032]** below measured H3K27ac. Predicted-alone
transferred is +0.005 [-0.019, +0.028].

Deployment therefore fails for **both** targets, H3K27ac deployment arms 0.455 and 0.452, p300
0.466, against a 0.457 floor, so the barrier is transfer rather than target choice, which is the
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

---

## F-009: p300 transfer is strongly asymmetric, a K562-trained model retains most of its advantage in GM12878, a GM12878-trained model retains none in K562, and both fit their own cell type equally well
**Status:** established
**Claim:** Scoring both multimodal p300 models in both cell types on each cell type's own EP300, top-quintile Pearson: K562-trained is 0.619 locally and **0.507** transferred to GM12878; GM12878-trained is 0.608 locally and **0.312** transferred to K562; the ATAC-only p300 floors are 0.306 (K562) and 0.299 (GM12878). Paired within fold: both local models beat their own floor by the same margin (+0.313 [+0.286, +0.340] and +0.309 [+0.263, +0.355]), so neither model is weaker at home. But K562->GM12878 retains **+0.207 [+0.158, +0.257]** over the target's floor (83% of the local advantage, accessibility residual *r* +0.139 [+0.067, +0.210]), while GM12878->K562 retains **+0.006 [-0.038, +0.050]**, *p*=0.72, indistinguishable from the target's own accessibility model, with a **negative** residual (-0.088 [-0.140, -0.035]), i.e. worse than the accessibility baseline at what the baseline already does.
**Implications:** The training cell type is a first-class design variable: identical architecture and target give +0.207 or +0.006 over the floor depending only on where the model was trained. This retires two earlier readings, that the p300 benchmark gain simply does not transfer (F-008 update), and that the GM12878 p300 model might be weaker in absolute terms. Critically, **the CRISPR benchmark can only test the failing direction**: CRISPR data exists only for K562, so the downstream metric is structurally unable to evaluate GM12878-as-target, which is the direction that works. Any deployment claim for predicted p300 therefore rests on correlation metrics until a benchmark in a second cell type exists, a limitation of the evaluation, not of the model. The H3K27ac explanation for its own transfer asymmetry (F-003: GM12878 is the harder cell type) cannot apply here, since both cell types are equally predictable locally; the remaining candidate is that ENCSR000EGE supports a more portable sequence model than ENCSR000DZG despite equal local fit.
**Caveats:** Two cell types only, so "K562-trained models are portable" is not separable from "K562 transfers to GM12878". The metric is top-quintile Pearson on p300, not benchmark AUPRC, so the magnitudes are not comparable to F-008's numbers.
**Tags:** p300, transfer, asymmetry, deployment, training-cell-type, k562, gm12878, benchmark-limitation

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-06 | job 42280519 (`2.27`) | K562 EP300 ENCSR000EGE + GM12878 EP300 ENCSR000DZG, candidate elements in both, 5 folds | 2026_0824_H3K27ac_model | local 0.619 / 0.608 against floors 0.306 / 0.299; transferred 0.507 (K562->GM) and 0.312 (GM->K562); drops -0.102 and -0.307 | supports |

**Update 2026-09-08: training-signal volume does not explain the asymmetry.** Retraining the
K562 p300 model on GM12878's exact budget (30.0M mapped reads from 51.1M, 21,068 peaks from
28,532) leaves the transfer intact: +0.202 [+0.147, +0.257] over GM12878's ATAC-only floor
against the full-depth model's +0.207 [+0.158, +0.257], a paired difference of -0.005
[-0.021, +0.010] (*p*=0.40). In-cell it costs -0.020 [-0.040, +0.001] (*p*=0.056), so the
depth-matched model is barely worse at home and is not simply a weaker model. A K562 model on
GM12878's budget transfers at +0.202; the GM12878 model on that same budget transfers at
+0.006.

Peak-set geometry is also excluded: both training sets are summit-centred fixed-width (316 bp
K562, 350 bp GM12878), the trainer centres on `start + summit`, and GM12878's summits are the
tighter of the two (|summit - midpoint| p90 0 bp against 6 bp).

The asymmetry is therefore a property of K562 and GM12878 as training cell types that survives
equalising volume. Two cell types cannot separate "K562-trained models are portable in general"
from "K562 happens to transfer to GM12878", so a third cell type is required rather than
merely desirable.

| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-08 | build 42400234, training 42402641-49, scoring 42440125 | K562 EP300 subsampled to GM12878's read and peak budget | 2026_0824_H3K27ac_model | depth-matched transferred +0.202 vs full-depth +0.207, paired -0.005 (*p*=0.40); in-cell 0.599 vs 0.619 | refutes |

---

## F-010: Swapping DNase for ATAC as the model's accessibility input beats ATAC in both cell types in-cell; its transfer claim was an artefact of an ATAC-only floor and is withdrawn
**Status:** established in-cell; the transfer claim is SUPERSEDED by F-016 (2026-09-14)
**Claim:** Retraining the multimodal H3K27ac model with a DNase 5'-insertion track in place of the ATAC one gains **+0.037 [+0.019, +0.056]** top-quintile Pearson in-cell K562 (*p*=0.005) and **+0.086 [+0.076, +0.095]** in-cell GM12878 (*p*<1e-4). Transferred K562 to GM12878 it gains **+0.043 [+0.003, +0.084]** over the ATAC-input model (*p*=0.041) and clears GM12878's own ATAC-only floor by **+0.056 [+0.014, +0.099]** (*p*=0.020), where the ATAC-input model manages +0.013 (*p*=0.25). Absolute top-quintile: in-cell 0.736 against 0.699 (K562) and 0.671 against 0.585 (GM12878); transferred 0.574 against 0.531, floor 0.518. On the accessibility residual the transferred DNase model scores 0.236 [0.211, 0.262] against the transferred ATAC model's 0.041 [0.003, 0.079] and the floor's 0.026.
**Implications:** This is the first architecture or input change in the project to survive transfer with a resolvable margin, against a background where the receptive field, fragment channels, the sequence gate, the asymmetric loss and GC-matched negatives all failed to travel (F-005) and where predicted-H3K27ac ABC arms never cleared the benchmark floor (F-004). It also supplies the premise the ATAC-to-DNase converter needed: DNase is the better representation to be in, the advantage survives transfer, and DNase has a reproducible base-resolution profile to target (0.848 and 0.686 at 1 bp against H3K27ac's 0.21 and 0.18). The residual numbers locate the effect: transferred, the ATAC-input model adds essentially nothing beyond accessibility while the DNase-input model adds real signal.
**REVISED 2026-09-14, see F-016.** The transfer claim in the title and in the sentence above -- "clears GM12878's own ATAC-only floor by +0.056 [+0.014, +0.099]" -- does not survive an assay-matched floor. That floor was an ATAC-only model; a DNase-input model can only be deployed where DNase exists, so the comparator must also be a DNase-only model. GM12878's DNase-only floor is 0.615 against its ATAC-only floor's 0.518, and the transferred K562 DNase model scores 0.570, i.e. **-0.045 [-0.080, -0.010] (*p*=0.023) BELOW** the assay-matched floor. The in-cell gains (+0.037 K562, +0.086 GM12878) are unaffected, as is the premise it supplied to the converter, since both are within-cell-type comparisons against the same-cell ATAC arm. What is withdrawn is "the first H3K27ac model to survive transfer".
**Caveats:** The comparison is depth-confounded, but against the result rather than for it: DNase is the shallower input in both cell types (301.1M usable reads against ATAC's 545.7M in K562; 53.8M against 571.4M in GM12878) and GM12878's DNase is the weaker library by its own count ceiling (0.445 against K562's 0.925). The advantage is nonetheless larger where the handicap is larger (+0.086 at 10.6x shallower against +0.037 at 1.8x), so depth-matching would be expected to widen the gap. Not yet tested downstream: no ABC or CRISPR-benchmark arm has been run with a DNase-input model.
**Tags:** dnase, atac, accessibility-input, transfer, deployment, h3k27ac, k562, gm12878, positive-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-08 | tracks 42491539, training 42499405-19, scoring 42514074 | K562 DNase ENCSR000EOT+ENCSR000EKS, GM12878 DNase ENCSR000EMT, H3K27ac targets in both | 2026_0824_H3K27ac_model | in-cell +0.037 (K562) and +0.086 (GM12878); transferred +0.043 over ATAC and +0.056 over the target's own ATAC-only floor | supports |

---

## F-011: A sequence + ATAC model predicts the DNase base-resolution profile at 95% of the inter-replicate ceiling in K562, and the advantage over the raw ATAC track survives transfer
**Status:** established (correlation metrics only; no downstream test)
**Claim:** Trained in K562 on ATAC-derived candidate elements to predict pooled stranded DNase 5' insertions, the converter reaches a raw within-element shape correlation of **0.792** at 1 bp on the top signal quintile, against an inter-replicate ceiling of 0.834, **95.1% of ceiling**, where the untouched observed ATAC track sits at 0.292 (35%). Paired within fold, converter minus observed ATAC is **+0.500 [+0.492, +0.508]** (*p*<1e-4). Applied unchanged to GM12878 the converter reaches **0.490** against that cell type's lower ceiling of 0.684 (71.6%), where observed ATAC is 0.206 (30%), a paired gain of **+0.284 [+0.279, +0.289]** (*p*<1e-4). An ATAC-only arm, same target, same elements, sequence branch removed, reaches 0.563 in K562 and 0.318 transferred, so sequence contributes **+0.229** in-cell and **+0.172** on transfer rather than the model merely rescaling its own input. On counts the sequence margin over the ATAC-only arm is +0.064 [+0.047, +0.082] top-quintile Pearson (*p*=0.0005).
**Implications:** This supplies the second premise the converter needed. F-010 established that DNase is the better representation to be in and that the advantage survives transfer; this establishes that a model with no DNase at inference can produce a base-resolution DNase track that is far more DNase-like than the ATAC track it would replace, and that the gain travels. The sequence margin over the ATAC-only arm retains 91% of its in-cell size on transfer (27.6 vs 25.1 percentage points of ceiling), against a project background where in-cell gains have repeatedly been exactly null transferred (F-005, fragment channels). Count loss weight is not a live variable here: clw=1 and clw=10 differ by 0.001 in shape at every bin size, so the project-standard clw=10 is kept.
**Caveats:** **Correlation only, no downstream test exists.** Whether feeding the painted track to the H3K27ac model improves it is untested, and so is the DNase input itself (F-010 carries the same gap). The absolute transfer number is much weaker than in-cell: 71.6% of ceiling against 95.1%, and 0.490 raw. The reported ceiling is computed from two replicates (142M reads) while the target is the 3-BAM pool (301M), so the true ceiling is higher and every "% of ceiling" here is optimistic. The paired CIs are very tight because fold-to-fold variation in a mean over ~10k elements is small; the real uncertainty is the third cell type, not fold sampling, and F-009's warning applies unchanged, two cell types cannot separate "K562-trained converters are portable" from "K562 transfers to GM12878". Elements are ATAC-derived by design so the converter never trains on regions chosen by the signal it predicts; this departs from every other K562 model in the project, justified by element derivation making no difference on the H3K27ac task (*p*=0.83).
**Tags:** converter, atac-to-dnase, profile, base-resolution, transfer, deployment, k562, gm12878, positive-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-10 | targets 42654625 (`0.34`), training 42713703-13 + 42722426-56 (`1.24`), scoring 42728292 (`2.34`), baseline 42728298 (`0.35`) | K562 ATAC 5' input, pooled stranded DNase target (ENCSR000EOT+ENCSR000EKS), ATAC-derived K562 elements, 5 folds; applied to GM12878 ATAC + ENCSR000EMT DNase | 2026_0824_H3K27ac_model | 1 bp top-quintile shape vs observed DNase: K562 ATAC 0.292 -> converter 0.792 (ceiling 0.834); GM12878 ATAC 0.206 -> converter 0.490 (ceiling 0.684); paired +0.500 and +0.284, both *p*<1e-4 | supports |

**Methodological note.** The metric here is 0.25's raw `shape_corr`, not 2.15's
`profile_pearson`. The latter is Gaussian-smoothed (`kernel_sigma=7, kernel_width=81`) and is
not on the same scale as the 1 bp inter-replicate ceiling: on the same fold-0 model it reads
0.774 where the raw 1 bp value is 0.799 and the ceiling is 0.831. `2.31` therefore computes
the model score and the ceiling on the same held-out windows with the same estimator, and
`0.35` supplies the comparator the ceiling framing omits, fraction of the DNase ceiling says
how close the converter gets to DNase, but the deployment question is whether it beats the
ATAC track already in the input slot, which is a different and lower bar that had to be
measured rather than assumed.


---

## F-012: Real DNase beats real ATAC as the ABC activity term by +0.071, and neither ATAC-to-DNase converter recovers it
**Status:** established
**Claim:** Scored as the whole ABC activity term on the K562 CRISPR benchmark, 10,342
element-gene pairs with 466 regulated and every arm on the identical pair set, real DNase
reaches AUPRC **0.5280** against real ATAC's **0.4572**, a paired-bootstrap delta of
**+0.0709 [+0.0476, +0.0938]** with the sign kept in 100% of resamples. DNase predicted from
ATAC by the K562-trained converter reaches 0.4780, **+0.0208 [-0.0049, +0.0441]** over real
ATAC, which does not resolve; the GM12878-trained converter applied to K562 reaches 0.4707,
**+0.0136 [-0.0118, +0.0364]**, which also does not resolve. Real DNase still beats the better
converter by **+0.0501 [+0.0298, +0.0714]** (100% sign kept), so the shortfall is resolvable
even though the gain is not. The two converters are indistinguishable from each other
(+0.0072 [-0.0123, +0.0273]).
**Implications:** Swapping the accessibility assay is a larger downstream lever than anything
the modelling has produced: +0.071 against +0.055 for predicted p300 (F-008) and against no
predicted-H3K27ac arm ever clearing the floor (F-004). But the converter captures at most 29%
of that gain and not resolvably, so F-011's profile result does not transfer into downstream
value here. **This benchmark cannot test what the converter was built for.** ABC sums the
activity track over each region and discards the profile, so it scores the converter's COUNTS,
and its counts are not good enough to substitute for the assay. The test that exercises the
profile is the H3K27ac model, whose accessibility branch reads base resolution.
**Caveats:** All four arms are 5-prime bigwigs so that no arm differs by counting path (D-21),
which means these numbers are not the July run's: real ATAC reproduces it (0.4572 against
0.4573) but real DNase does not (0.5280 against 0.5763), so the July DHS number must not be
quoted beside these. Training cell type makes no difference here, but the GM12878 converter
also trained on a differently derived element set (D-19), so "no transfer penalty" and "no
element-derivation penalty" are not separated. The converters' painted tracks are the
`4.1`-style flat-within-region kind, which is correct for ABC and useless as a base-resolution
input.
**Tags:** dnase, atac, converter, abc, crispr, activity-term, deployment, negative-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-11 | ABC 42798274, CRISPR 42801835, paired `4.17` | K562 ATAC 5-prime, K562 DNase 5-prime (ENCSR000EOT+ENCSR000EKS), converted DNase from `1.24` and `1.25` | 2026_0824_H3K27ac_model | real DNase 0.5280, converters 0.4780 and 0.4707, real ATAC 0.4572; paired +0.0709, +0.0208, +0.0136 | supports for DNase, refutes for the converters |


---

## F-013: Most of DNase's input advantage is sub-250 bp structure, and the converted track is worse than the ATAC track it would replace
**Status:** established
**Claim:** Four accessibility inputs, H3K27ac target, identical element set, 5 folds, paired
within fold against the ATAC-input model. Top-quintile Pearson: real DNase **0.728**
(+0.0379 [+0.0213, +0.0546], *p*=0.0032, replicating F-010's +0.037), real DNase smoothed at
250 bp **0.704** (+0.0141 [-0.0057, +0.0338], *p*=0.12, NOT resolvable), ATAC **0.690**, and
the converted DNase track **0.680** (**-0.0099 [-0.0175, -0.0023]**, *p*=0.023, i.e.
resolvably WORSE than ATAC). Destroying structure finer than 250 bp while preserving
magnitude costs real DNase **+0.0239 [+0.0027, +0.0451]** (*p*=0.035) of its +0.0379
advantage, and what remains no longer separates from ATAC. Real DNase beats the converted
track by +0.0478 [+0.0384, +0.0573] (*p*=0.0001). On the accessibility residual the converted
track scores 0.133 [0.113, 0.153] against real DNase's 0.456 and smoothed DNase's 0.413, and
its incremental R2 over the ATAC baseline is negative, -0.020 [-0.026, -0.013].
**Implications:** The shape-versus-magnitude question is answered in favour of shape: the
majority of DNase's advantage as an input lives in base-resolution structure, and magnitude
alone does not resolvably beat ATAC. That is the branch of the control (fig18) that keeps a
converter conceptually alive. But this converter does not deliver it. Reproducing the DNase
profile at 95% of its inter-replicate ceiling (F-011) while performing WORSE than raw ATAC
downstream means the painted track carries a defect that outweighs its shape fidelity, and
the named candidate is magnitude: the 1.54x inflation of the lowest observed-signal quintile
that survived thresholding (D-22). Shape fidelity is necessary and demonstrably not
sufficient. Separately, `profile_pearson` is 0.063-0.064 across ALL four arms, so the H3K27ac
profile head learns nothing regardless of what accessibility it is given, which is the
premise the DNase-auxiliary-target proposal rests on.
**Caveats:** The smoothed arm's interval is wide, [-0.0057, +0.0338] against real DNase's
+0.0379, so "most of the advantage is shape" is the point estimate and the interval does not
exclude magnitude carrying much of it. 250 bp is one filter width, chosen because observed
ATAC and DNase already agree about shape there (r = 0.92 against 0.29 at 1 bp); a width sweep
would localise the scale that matters and has not been run. The converted arm confounds shape
fidelity with the painted track's magnitude distortion and its 34% genome-wide total, so its
negative result does not isolate either. K562 only; no transfer arm for the converted or
smoothed inputs.
**Tags:** dnase, atac, converter, accessibility-input, shape, magnitude, h3k27ac, control, negative-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-12 | control track `0.36`, painting 42923575 + `4.22`, training 43080124-33 (`1.26`), scoring 43118504 (`2.36`) | K562 H3K27ac 5' target under four accessibility inputs: ATAC, real DNase, DNase smoothed 250 bp, converted DNase | 2026_0824_H3K27ac_model | top-quintile Pearson 0.690 / 0.728 / 0.704 / 0.680; smoothing costs +0.0239 of DNase's +0.0379; converted arm -0.0099 against ATAC | supports shape, refutes this converter |


---

## F-014: The converter fails as an input in all four transfer directions, and whether DNase wins on shape or magnitude depends on the cell type
**Status:** established
**Claim:** Four accessibility inputs x two training cell types, each scored in both cell
types, 5 folds, paired within fold against the ATAC-input arm of the same direction.
Top-quintile Pearson:

| input | K562->K562 | K562->GM12878 | GM->GM | GM->K562 |
|---|---|---|---|---|
| ATAC | 0.690 | 0.527 | 0.579 | **0.607** |
| real DNase | **0.728** (+0.038) | **0.570** (+0.043) | **0.663** (+0.084) | 0.474 (-0.133) |
| DNase smoothed 250 bp | 0.704 (+0.014, ns) | 0.503 (-0.024, ns) | 0.662 (+0.083) | 0.309 (-0.297) |
| converted, power gamma=1.2 | 0.679 (-0.011) | 0.469 (-0.058) | 0.558 (-0.021) | 0.484 (-0.123) |

**The converted input is resolvably worse than plain ATAC in every direction**: -0.0105
[-0.0202, -0.0007] in-cell K562, -0.0579 [-0.1005, -0.0153] transferred to GM12878, -0.0206
[-0.0313, -0.0100] in-cell GM12878. **Whether DNase's advantage is shape or magnitude is
cell-type dependent**: smoothing at 250 bp costs real DNase +0.0239 [+0.0027, +0.0451]
(*p*=0.035) in K562 and +0.0677 [+0.0533, +0.0822] (*p*=0.0002) transferred to GM12878, but
**+0.0011 [-0.0069, +0.0091] (*p*=0.72) in-cell GM12878**, i.e. nothing at all. **Transfer is
asymmetric in the same direction as p300** (F-009): K562-trained models keep the DNase
advantage in GM12878, while every DNase-family input transferred GM12878->K562 lands far below
ATAC, and the plain ATAC-input model is the best transferred model in that direction.
**Implications:** The converter is finished as a model input. Two independent attempts at the
magnitude defect diagnosed in F-013, quantile mapping (0.37) and a power transform fitted
out-of-cell (0.39, 0.40), both left the downstream number unmoved (-0.0099 -> -0.0105), so the
magnitude distortion was not what was costing it. Reproducing the DNase profile at 95% of its
inter-replicate ceiling (F-011) is not sufficient to substitute for the assay, and no further
post-hoc transform is worth trying. F-013's headline, that DNase wins on sub-250 bp structure,
survives only where DNase is deep enough to HAVE reliable shape: GM12878's DNase is 53.8M reads
against K562's 301.1M with a profile ceiling of 0.686 against 0.848, and there the entire
advantage is magnitude. F-013 was K562-only and flagged that as a caveat; the caveat turned out
to be load-bearing.
**Caveats:** The GM12878->K562 collapse is large enough (-0.13 to -0.30) to be worth a
sanity check that it is not an artefact of the GM12878 models' own weakness rather than of
direction; the in-cell GM12878 arms are strong (0.663), which argues against that, but a third
cell type is the only clean test. Two cell types still cannot separate "K562 is a good training
cell type" from "K562->GM12878 is a good pair", and this is now the second finding pointing at
that gap (F-009 was the first). The converted arm used gamma fitted in GM12878 and applied to
both, which is deployment-legal but leaves K562's dynamic range at 24x against real DNase's 39x.
**Tags:** dnase, atac, converter, accessibility-input, transfer, asymmetry, shape, magnitude, negative-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-13 | transform 43300449/51 (`0.39`, `0.40`), training 43302465-77 (`1.26`, `1.27`), scoring 43337247 (`2.37`) | K562 and GM12878 H3K27ac under ATAC / real DNase / DNase smoothed 250 bp / converted-power inputs, both training cell types, both targets | 2026_0824_H3K27ac_model | converted input below ATAC in all four directions; smoothing costs +0.024 in K562 and +0.001 in GM12878; GM12878->K562 collapses for every DNase input | refutes the converter, qualifies F-013 |

---

## F-015: A DNase profile head is a learnable auxiliary task and improves H3K27ac counts, but only on overall correlation; the top-quintile gain does not survive an epoch-matched control
**Status:** established
**Claim:** In K562, swapping the profile head's target from H3K27ac to pooled stranded DNase while the counts head keeps H3K27ac turns a dead auxiliary task into a live one: validation profile Pearson rises from 0.060-0.068 to 0.530-0.561 in all five folds, and the head abandons H3K27ac shape entirely (`profile_pearson` 0.004 on held-out elements against the baseline's 0.064). Against an **epoch-matched** baseline the counts head gains **+0.0026 [+0.0006, +0.0046]** (*p*=0.022) in overall Pearson, with `incremental_r2` 0.005 [0.002, 0.007] and `incremental_r2_topq` 0.015 [0.001, 0.030] both clearing zero where the control's do not (0.000 and 0.005, both spanning zero). **The top-quintile Pearson gain is NOT resolvable against that control: +0.0076 [-0.0078, +0.0230], *p*=0.24.** It measures +0.0109 (*p*=0.043) against the early-stopped baseline, but roughly a third of that is the epoch difference. Needs no DNase at inference: DNase is a training-time target only.
**Implications:** The premise behind the profile head's redesign is confirmed. H3K27ac's 1 bp inter-replicate ceiling is 0.21 and `profile_pearson` sits at 0.063-0.064 across all four accessibility inputs (F-013), so that head was consuming capacity for a near-zero gradient, exactly as D-3's consequences note predicted; DNase at a 0.848 ceiling gives it a real task, and the trunk does get better at counts. But the effect is small and, on this project's headline metric, below resolution: top-quintile Pearson 0.701 against the control's 0.693 and the baseline's 0.690. It is an order of magnitude below swapping the accessibility INPUT to real DNase (+0.038, F-010). Treat this as a cheap architectural improvement worth keeping, not as progress on the gap the DNase input opens.
**The epoch confound was real and is now excluded.** Early stopping watches the total validation loss, which for the multi-task arm includes the reweighted DNase profile term, so the arms stopped on different objectives and the multi-task arm trained longer in every fold (53-99 epochs against 32-55). Given the same 100-epoch budget the baseline changes by **+0.0000 [-0.0008, +0.0009]** (*p*=0.92) overall and +0.0033 (*p*=0.40) top-quintile, and its best checkpoint still lands at epoch 33-47, so it had already converged where early stopping put it. Longer training was a consequence of the different loss, not the cause of the gain. This is the general lesson: when an intervention changes the loss, it changes the stopping rule too, and an epoch-matched arm is the only way to attribute the result.
**Caveats:** `--profile-loss-weight` was set to 0.0561, the measured H3K27ac/DNase per-window depth ratio (mean total 51.1 against 910.2, 17.8x), so the profile term starts at the baseline's magnitude and `count_loss_weight` keeps its meaning. MNLL scales with the target's read depth, so an unweighted swap would have divided the effective count weight by about 18 and made a loss on counts uninterpretable. 0.0561 is the smallest defensible weight, and a larger one has not been tried, so the ceiling on this effect is unknown. K562 only, chosen because F-014 showed GM12878's DNase has no reliable sub-250 bp structure to learn. `profile_*` columns compare every arm against H3K27ac, so the multi-task arm's drop there is by construction and carries no information.
**Tags:** h3k27ac, dnase, multi-task, auxiliary-target, profile-head, k562, architecture, early-stopping, methodology

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-14 | training 43416027-31 (`1.28`), scoring 43428399 (`2.39`) | K562 H3K27ac counts + K562 pooled stranded DNase profile, ATAC input, DNase-derived elements, 5 folds | 2026_0824_H3K27ac_model | validation profile *r* 0.530-0.561 against 0.060-0.068; top-quintile 0.701 against the early-stopped baseline's 0.690, paired +0.0109 (*p*=0.043) | supports |
| 2026-09-14 | control 43440391-96 (`1.29`), 3-arm scoring 43451521 (`2.39`) | same baseline given a 100-epoch budget, checkpoint selection unchanged | 2026_0824_H3K27ac_model | control moves +0.0000 (*p*=0.92) overall, so the epoch gap explains nothing; but multi-task minus control is +0.0026 (*p*=0.022) overall and only +0.0076 (*p*=0.24) top-quintile | qualifies |

---

## F-016: No H3K27ac model transfers across cell types once the floor uses the same assay as the model, in any of six directions over three cell types
**Status:** established
**Claim:** Three cell types with real DNase (K562, GM12878, THP-1), each scored against its OWN locally trained DNase-only model as the floor, 5 folds, paired within fold. Top-quintile Pearson:

| target | local floor (DNase only) | local sequence+DNase | transferred in |
|---|---|---|---|
| K562 | 0.607 | **0.728** (+0.121) | GM12878 0.474 (**-0.133**), THP-1 0.428 (**-0.178**) |
| GM12878 | 0.615 | **0.663** (+0.048) | K562 0.570 (**-0.045**), THP-1 0.599 (-0.016, *p*=0.11) |
| THP-1 | 0.660 | **0.717** (+0.057) | GM12878 0.590 (**-0.070**), K562 0.272 (**-0.388**) |

**All six transferred arms land at or below the target's own DNase-only floor**, five of them resolvably (*p* from 0.0003 to 0.023) and the sixth indistinguishable from it. Meanwhile a locally trained sequence+DNase model beats that floor in all three cell types (+0.121, +0.048, +0.057, all *p*<0.006), so sequence does add real information wherever it can be trained locally.
**Implications:** The deployment claim this project has been pursuing does not hold for H3K27ac. If a new cell type has DNase, a DNase-only model trained on that cell type's own H3K27ac beats every sequence+DNase model trained elsewhere. F-010's transfer result came from comparing a DNase-input model against an ATAC-only floor, and an ATAC-only floor is much the weaker comparator (GM12878 0.518 against 0.615); the general rule is that **a floor must use the same assay as the model it bounds**, or the comparison credits the assay swap to the architecture. Note what this does NOT say: neither floor is a true deployment baseline, since both need target-cell H3K27ac to train. The honest deployment comparator is the raw accessibility track used directly as the activity proxy, which is what ABC does and what F-012 measured (+0.0709 for real DNase over real ATAC). That remains the only DNase result in this project with a downstream effect.
**It also answers the question F-009 and F-014 were both stuck on, in the negative.** Neither "K562 is a good training cell type" nor "K562->GM12878 is a good pair" survives. K562 transfers acceptably into GM12878 (-0.045) and collapses into THP-1 (-0.388), so it is not a good universal trainer; and K562 is simultaneously the hardest target (-0.133 and -0.178 coming in). The best-transferring pair is GM12878<->THP-1 (-0.016 and -0.070), which is not the pair lineage would predict, K562 and THP-1 both being myeloid leukaemia lines against GM12878's lymphoblastoid.
**Caveats:** **The mechanism was partly depth, see F-017.** Per-window DNase input depth differs sharply (mean total over 1 kb windows: K562 894.6, THP-1 281.4, GM12878 126.3). This finding originally argued depth was unlikely to be the cause because the mismatch ratios do not order with transfer quality (K562->GM12878 is a 7.1x mismatch at -0.045 while K562->THP-1 is 3.2x at -0.388). **That reasoning was wrong**: depth-matching the libraries (F-017) removed the largest collapse entirely, taking GM12878->K562 from -0.133 to +0.003 (*p*=0.68). It did not change the conclusion of this finding, since five of six directions stay resolvably below the floor and the sixth only reaches parity, but the numbers quoted above are the unmatched ones and F-017's supersede them as the measure of what transfer can do. `k562_to_thp1` has a very wide interval ([0.120, 0.425]) so it is unstable across folds rather than uniformly bad, and its collapse should not be quoted as a point estimate. THP-1 has one DNase replicate, so it has no DNase shape ceiling; its H3K27ac is deeper than K562's (208.2 against 51.1 per window) while its DNase is shallower, so it is not simply the hardest cell type. Element sets are each cell type's own DNase-derived set, so transferred arms are scored on regions chosen by the target's accessibility, which is the deployment-correct choice but means the three columns are not the same regions.
**Tags:** h3k27ac, dnase, transfer, deployment, floor-definition, thp1, k562, gm12878, negative-result, methodology

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-14 | training 43442172-84 (`1.30`), floors 43445448-533 (`1.22`/`1.23` MODE atac), scoring 43458705 (`2.41`) | K562 / GM12878 / THP-1 DNase inputs and H3K27ac targets, each cell type's own DNase-derived elements, 5 folds | 2026_0824_H3K27ac_model | all six transferred arms at or below the target's own DNase-only floor; local sequence+DNase beats that floor in all three | refutes the transfer claim in F-010 |

---

## F-017: Depth-matching the DNase libraries removes the largest transfer collapse and explains the asymmetry, but transfer still loses to a locally trained model
**Status:** established
**Claim:** Thinning K562 and THP-1 DNase by binomial subsampling to GM12878's 50.3M main-chromosome reads (`0.41`, achieved within 0.01%) brings per-window means from 894.6 / 281.4 / 126.3 to 153.1 / 172.1 / 126.3, a 7.1x spread down to 1.4x. Element sets are not recalled, so the regions scored are identical to F-016's. Top-quintile Pearson, transferred, with the paired delta against the target's own DNase-only floor:

| direction | test/train depth | F-016 | depth-matched |
|---|---|---|---|
| GM12878->K562 | 7.1x deeper | 0.474 (-0.133) | 0.573 (**+0.003, *p*=0.68**) |
| THP-1->K562 | 3.2x deeper | 0.428 (-0.178) | 0.530 (-0.041) |
| GM12878->THP-1 | 2.2x deeper | 0.590 (-0.070) | 0.613 (-0.042) |
| K562->GM12878 | 7.1x shallower | 0.570 (-0.045) | 0.576 (-0.039) |
| THP-1->GM12878 | 2.2x shallower | 0.599 (-0.016, ns) | 0.561 (-0.054) |
| K562->THP-1 | 3.2x shallower | 0.272 (-0.388) | 0.358 (-0.297) |

**The asymmetry that F-009 and F-014 both ended on was library depth.** GM12878->K562 was the severe collapse; depth-matched it reaches parity with K562's own floor. A model trained on a shallow library and then handed a 7x deeper one is out of distribution on its own input. The mild direction, K562->GM12878, was NOT depth and is unchanged (-0.045 to -0.039). **Transfer nonetheless still fails**: five of six directions remain resolvably below the target's floor and the sixth only reaches parity, while a locally trained sequence+DNase model beats that floor in all three cell types (+0.136, +0.048, +0.064, all *p*<0.003).
**Implications:** F-016's deployment conclusion stands and its mechanism is now half-explained. Accessibility library depth is a model input property that must be matched before any cross-cell-type comparison, which generalises the fragment-channel lesson (F-005) from fragment-length structure to read depth. Any future transfer claim in this project should depth-match first, or a collapse will be misattributed to the cell types. Note also that thinning COSTS the deep cell type real performance in its own right: K562's floor falls 0.607 to 0.570 and its local model 0.728 to 0.707, so depth-matching is the right control for a transfer question and the wrong choice for a production model, which should use all the reads it has.
**Caveats:** **The direction of the mismatch looks to matter more than its size, with one exception, and I would not lean on the rule.** The three directions tested on DEEPER input than they trained on all recovered (+0.099, +0.102, +0.023 absolute) while the shallower ones did not (+0.006, -0.038, +0.086). K562->THP-1 breaks it, and that arm also carries the widest interval in the panel ([0.257, 0.459]), so it is unstable across folds rather than a clean counterexample. F-016's original caveat argued depth was unlikely because the ratios do not order with transfer quality; the ratios still do not order with it, and the effect was real anyway, so ratio-ordering was the wrong test of the hypothesis. Thinning is one random draw at seed 0; a second seed would bound the sampling noise and has not been run. GM12878 is unthinned by construction, so its arms are shared between the two panels and its in-cell numbers are identical, which is what makes this a one-variable change but also means the panel has no cell type where the thinning was validated against an independent deeper replicate.
**Tags:** dnase, transfer, deployment, read-depth, library-properties, thinning, k562, gm12878, thp1, methodology

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-15 | thinning 43483926 (`0.41`/`0.42`), training 43485132-186 (`1.26`/`1.30` with ACC_BW), scoring 43537443 (`2.41` DEPTH_MATCHED=1) | K562 and THP-1 DNase thinned to GM12878's 50.3M reads, same H3K27ac targets and element sets as F-016 | 2026_0824_H3K27ac_model | GM12878->K562 recovers from -0.133 to +0.003; five of six directions still below the target floor | qualifies F-016 |

---

## F-018: Predicted H3K27ac from a DNase-input model clears the ABC benchmark floor and adds on top of real DNase, the first predicted-activity arm in this project to do either
**Status:** established
**Claim:** In K562, the multimodal H3K27ac model retrained with DNase in place of ATAC as its accessibility input (F-010) was predicted genome-wide with per-fold leakage assembly and RC averaging over the 153,545 ABC candidate regions (0 dropped), then substituted into ABC's activity term. Paired bootstrap over 10,342 element-gene pairs scored by every arm, 2,000 resamples, 466 regulated (4.51%):

| contrast | delta AUPRC | 95% CI | sign kept |
|---|---|---|---|
| DNase-input predictor - ATAC-input predictor, geomean with ATAC | **+0.0282** | [+0.0164, +0.0400] | 100% |
| DNase-input predictor - ATAC-input predictor, prediction alone | **+0.0287** | [+0.0107, +0.0467] | 100% |
| DNase x predicted H3K27ac - real DNase alone | **+0.0238** | [+0.0051, +0.0429] | 99.2% |
| DNase x predicted H3K27ac - ATAC-only floor | **+0.0839** | [+0.0613, +0.1048] | 100% |
| ATAC x predicted H3K27ac - ATAC-only floor | **+0.0403** | [+0.0230, +0.0563] | 100% |

Absolute AUPRC on the shared pair set: DNase x predicted 0.5519, observed H3K27ac 0.5296, real DNase alone 0.5280, predicted alone 0.5215, ATAC x predicted 0.5083, ATAC-input arms 0.4928 and 0.4801, floor 0.4680.
**Implications:** Two things follow that the project did not have before. **First, predicted H3K27ac is downstream-useful after all**, but only from a model that reads DNase: the identical architecture reading ATAC cannot separate from the floor even under the paired test (+0.0121, CI spanning zero, F-004's update). The input assay, not the architecture and not the statistics, is what moved the benchmark. **Second, the prediction adds to the assay rather than merely recovering it.** F-012 established real DNase alone as the best activity term in the project (+0.0709 over real ATAC); putting predicted H3K27ac on top of real DNase beats real DNase alone by +0.0238 [+0.0051, +0.0429]. So the model contributes information the accessibility track does not carry, which is the first direct evidence for that claim anywhere in this work, every earlier version having been a correlation argument.
**This is NOT a deployment result, and the scarcity runs the wrong way.** DNase is available for FEWER cell types than ATAC, which is the entire reason the application target is ATAC-only and the reason an ATAC-to-DNase converter was attempted and closed (F-014). So a DNase-input model is LESS deployable than the ATAC-input one it beats. F-018 says what information improves H3K27ac prediction enough to matter downstream; it does not supply a route to a new cell type that has ATAC alone. The only deployment-legal candidate still untested is the multi-task model (F-015), which uses DNase as a training TARGET and reads ATAC alone at inference; its top-quintile gain was +0.0026, but F-004 established that this benchmark and top-quintile Pearson disagree, so it has not been ruled out.
**Caveats:** **The comparison against the observed-H3K27ac arm is NOT attributable and must not be quoted as beating the ceiling.** The DNase x predicted arm exceeds it by +0.0222 [+0.0041, +0.0401], but that arm counts reads through count_bam/count_tagalign while every bigwig arm goes through count_bigwig, so a difference of this size sits inside the counting-path confound. The contrasts in the table above are each matched on counting path by construction (`4.6` records which pairs are legitimate). **This is the "you have DNase" scenario, not the ATAC-only deployment target**: the model needs DNase at inference, so it does not address the application constraint that motivated the converter. **It is also in-cell K562 only.** F-016 and F-017 showed no H3K27ac model transfers across cell types once the floor is assay-matched, so a transferred version of this arm should be expected to fail and has not been run. F-004's named mechanism, compressed dynamic range of the predicted activity term on regulated-pair regions (p99/p50 4.03 against 6.00 observed), has NOT been re-measured for the DNase-input prediction, so whether the gain comes from better spread or from something else is untested.
**Tags:** h3k27ac, dnase, abc, crispr-benchmark, downstream-evaluation, activity-term, k562, positive-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-15 | prediction 43542121 (`4.4`), ABC 43544032 (`4.5`), benchmark 43549686 (`4.7`), paired bootstrap (`4.17`) | EPCrisprBenchmark_ensemble_data_GRCh38; 153,545 ABC candidate regions, region set 7d5995ce shared by all arms | 2026_0824_H3K27ac_model | DNase x predicted H3K27ac 0.5519 against a 0.4680 floor, +0.0839 paired; beats real DNase alone by +0.0238 | refutes F-004's generality |

---

## F-019: A sequence+ATAC model's predicted H3K27ac clears the ABC benchmark floor, DNase supervision is not what does it, and the longer-training explanation does not survive a run-variance check
**Status:** floor-clearing established; **the epoch attribution is WITHDRAWN, 2026-09-17, see F-020**
**Claim:** In K562, two models that need ONLY ATAC at inference were put through ABC and the CRISPR benchmark. Both are mode `multimodal`, i.e. sequence AND ATAC in, H3K27ac out; "ATAC-only" here means no DNase is required at prediction time, NOT the accessibility-only architecture that `atac` mode denotes in the trainer. They are: the multi-task model whose profile head trains on DNase (F-015, `1.28`) and the epoch-matched plain baseline that trains to the same 100-epoch budget with no DNase anywhere (`1.29`). Paired bootstrap, 2,000 resamples of the shared pair set:

| contrast | delta AUPRC | 95% CI | sign kept |
|---|---|---|---|
| multi-task - epoch-matched baseline, geomean | +0.0047 | [-0.0035, +0.0122] | 87.5% |
| multi-task - epoch-matched baseline, prediction alone | **+0.0129** | [+0.0022, +0.0240] | 99.3% |
| **epoch-matched baseline - early-stopped baseline** | **+0.0086** | [+0.0005, +0.0175] | 98.2% |
| epoch-matched baseline - ATAC-only floor | **+0.0208** | [+0.0039, +0.0373] | 99.2% |
| multi-task - ATAC-only floor | **+0.0255** | [+0.0068, +0.0439] | 99.7% |

**Implications:** **An ATAC-only-at-inference model can clear the benchmark floor.** That is new, and it is the only result in this project that does so under the actual deployment constraint, F-018's winner requiring DNase at prediction time. But the attribution is not what the multi-task hypothesis predicted: against the early-stopped baseline the multi-task arm gains +0.0133, and +0.0086 of that is simply the longer training the changed loss caused. DNase supervision itself is resolvable only in the prediction-alone form (+0.0129) and not in the geomean form (+0.0047). One speculative reading of that split, untested: the geomean multiplies the prediction by real ATAC, which compresses differences between predictions, so the prediction-alone arm is the purer test.
**The incidental result is arguably the more useful one: training longer is free, and nobody was doing it.** Every predicted-activity arm behind F-004 was early-stopped at patience 10, typically stopping near epoch 35. Running the identical model to 100 epochs moves the benchmark by +0.0086 [+0.0005, +0.0175] and takes it over the floor. **The same change is worth exactly nothing upstream**: F-015 measured the epoch-matched baseline at +0.0000 [-0.0008, +0.0009] (*p*=0.92) overall Pearson, with its best checkpoint still landing at epoch 33-47. So a model that is converged by its own validation loss, and indistinguishable on the project's headline metric, is still measurably better downstream. This is a third instance of the benchmark and top-quintile Pearson disagreeing (F-004 and F-018 being the others) and the first where the benchmark rewards something the upstream metric scores at zero.
**Caveats:** **The epoch comparison conflates the budget with run-to-run variance, and this is its main weakness.** `1.29` is a SEPARATE training run, not a continuation of the baseline, and cuDNN convolution is not bit-reproducible across nodes, which is why this project's regression gates use a tolerance rather than byte-identity (decision 2026-09-03). With one run per condition, +0.0086 [+0.0005, +0.0175] cannot be cleanly attributed to the extra epochs. The one piece of evidence that does isolate the budget: within the 100-epoch run the best validation checkpoint landed at epoch 40, 38 and 47 in folds 1, 2 and 4, later than where the early-stopped run ended in those folds, so the budget demonstrably changes which checkpoint is selected. Suggestive, not sufficient. **This was tested on 2026-09-17 and the attribution did not survive; see the update below and F-020.** Also note the epoch effect is modest in absolute terms. Checkpoint selection is unchanged in the control, so this is not overfitting-by-another-name: `1.29` keeps best-validation-loss selection and only prevents the run ending early. All arms here are geomean(real ATAC tagAligns, predicted bigwig) or a single predicted bigwig, so they are matched on counting path; the observed-H3K27ac arm is not and is reference only. K562 in-cell only. Whether the multi-task arm's prediction-alone advantage survives in another cell type is untested and, given F-016 and F-017, should not be assumed.

**Update 2026-09-17: the epoch attribution is WITHDRAWN. `1.31.epoch_budget_vs_run_variance.py`.**
The proposed within-run checkpoint test is impossible, because training saves only two
checkpoints per fold, best-validation and final. What the logs do support is comparing the two
runs at the SAME epoch, where the budget cannot yet have acted, so every difference there is
run variance. Run variance at epochs 20/25/30 across all five folds is **sd 0.1753 validation
MNLL and sd 0.00180 validation count Pearson**, against a scored best-checkpoint difference
between the two runs of **-0.0636 MNLL and +0.00003 count Pearson**. The claimed effect is 2.8x
smaller than the noise on MNLL and 60x smaller on count Pearson; the worst single matched-epoch
draw, fold 1 at epoch 25, is -0.4905 MNLL, roughly 8x the whole effect. The budget was also not
binding in every fold: baseline fold 1 ran to epoch 55 with its save at 45, LATER than the
100-epoch run's save at 40, so the picture of a baseline uniformly cut short is wrong.
Separately, and with zero run variance because both checkpoints come from one run, the
100-epoch run's final epoch is WORSE than its own best-validation checkpoint on all five folds
(MNLL +0.018 to +0.096, count Pearson -0.00046 to -0.00154), so "more gradient steps are
better" is not the mechanism either. **The +0.0086 downstream number itself is not refuted**,
and the benchmark has disagreed with upstream metrics three times, but it can no longer be
called the epoch budget. Resolving it needs 2-3 seeds per condition scored downstream, which
makes it no longer cheap. **Do not adopt 100 epochs as a default.**

**Tags:** h3k27ac, abc, crispr-benchmark, multi-task, early-stopping, deployment, atac-only, k562, methodology

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-16 | predictions 43803158 and 43832450 (`4.4`), ABC 43808834 and 43837720 (`4.5`), benchmark 43854329 (`4.7`), paired bootstrap (`4.17`) | EPCrisprBenchmark_ensemble_data_GRCh38; 153,545 ABC candidate regions, region set 7d5995ce | 2026_0824_H3K27ac_model | multi-task 0.4935 and epoch-matched baseline 0.4888 against a 0.4680 floor; DNase supervision itself +0.0047 (ns) in geomean, +0.0129 alone | qualifies F-015, extends F-004 |

---

## F-020: Run-to-run variance between two trainings of the same configuration is large enough to swallow several effects this project has treated as real, and it is measurable from the training logs alone
**Status:** established
**Claim:** Two runs of `multimodal5p_accs5p_hw500_clw10`, the `1.11` baseline and the 100-epoch
`1.29`, differ only in the epoch budget. At epochs 20, 25 and 30 neither run has stopped, so the
budget cannot have affected either one and any difference at a matched epoch is pure run variance
from initialisation, batch shuffling and non-deterministic cuDNN convolution. Across all five
folds:

| quantity, at matched epochs | mean | sd | max abs |
|---|---|---|---|
| validation MNLL | -0.0389 | **0.1753** | 0.4905 |
| validation count Pearson | -0.00034 | **0.00180** | 0.00413 |

For scale, the effects this project has recently reported and acted on include +0.0026 overall
Pearson for the multi-task profile head (F-015) and +0.0000 for the epoch budget (F-019). A
second, independent observation from the same logs: within the 100-epoch run the final epoch is
worse than that run's own best-validation checkpoint on all five folds, MNLL +0.018 to +0.096
and count Pearson -0.00046 to -0.00154, which is a zero-variance measurement because both
checkpoints come from one weights lineage.
**Implications:** **One run per condition cannot resolve a count-Pearson effect of order 0.002,
and this project has been running one per condition.** The immediate casualty is F-019's epoch
attribution, withdrawn. The general rule is that any single-run architecture or hyperparameter
comparison whose effect is smaller than the matched-epoch sd needs seeds before it is adopted,
and `1.31` costs nothing to run because it reads logs that every training already writes. It
should be run against any new adopt-or-not arm before the arm is adopted. Note also what this
does NOT say: it does not refute a downstream effect, because the benchmark and the upstream
metrics have disagreed in both directions three times (F-004, F-018, F-019). It says the
upstream metrics cannot be used to attribute one.
**Caveats:** These are the metrics logged on the validation split during training, not the
held-out test-fold overall and top-quintile Pearson that the evaluation scripts report, so the
sd is the same order as the project's reported effects rather than directly comparable to them.
Two runs give a 15-point sd estimate across folds and epochs, which is enough to establish the
order of magnitude and not much more. Checkpoint selection is best-validation-loss in both runs,
where the loss is MNLL plus the weighted count MSE, so the saved checkpoint is not always the
MNLL minimum; `1.31` reads the LAST saved epoch, which is what the checkpoint file holds.
**Tags:** methodology, reproducibility, run-variance, early-stopping, statistical-power, negative-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-17 | `1.31.epoch_budget_vs_run_variance.py` over the `1.11` and `1.29` training logs, 5 folds each | K562 H3K27ac 5' targets, validation split of each chromosome-holdout fold | 2026_0824_H3K27ac_model | matched-epoch run variance sd 0.1753 MNLL / 0.00180 count Pearson, against a best-checkpoint budget difference of -0.0636 / +0.00003; final epoch worse than best-val on 5 of 5 folds | withdraws F-019's epoch attribution |

---

## F-021: The best ATAC-input model on the transferred CRISPR benchmark predicts p300, and its transfer failure shows up as an accessibility-dependent error that ABC's qnorm and the CRISPR element subset largely manufacture
**Status:** the ranking is established and survives the F-025 retrain; **the accessibility-slope MECHANISM is withdrawn, 2026-09-17, see the update below**. The transferred arm here was trained on a BPNet model's predicted plus strand (F-025); retrained on the observed target it scores **+0.0135 over the floor rather than the +0.0091 quoted below**, with a transfer penalty of +0.0414 rather than +0.0458. Quote the corrected numbers.
**Claim:** Of every arm this project has benchmarked whose model reads ATAC (or ATAC+sequence)
and is applied to a cell type it was not trained on, **the best is a GM12878-trained multimodal
p300 model applied to K562**: AUPRC 0.4771 against a 0.4680 ATAC-only floor, +0.0091 [-0.0042,
+0.0222], sign kept 91.4%. It is the only transferred ATAC-input arm above the floor at all; the
two H3K27ac-target transfer arms are at it (-0.0028 and -0.0048, F-004 update 2026-09-17). The
comparison is architecturally clean: same PyTorch `multimodal_bpnet` code, same negatives and
genome, both arms predicted from the SAME K562 ATAC track, geomean with the same ATAC tagAligns.

**Transfer is what costs it, and that cost is fully resolvable**: in-cell minus transferred is
+0.0458 [+0.0296, +0.0610], sign kept 100%, which is 83% of the in-cell arm's entire +0.0549
margin over the floor. In-cell, predicted p300 reaches 0.5229 against an observed-H3K27ac anchor
of 0.5296.

**The error has a specific shape.** Both arms share the identical ATAC track, so stratifying by
ATAC decile is legitimate and the ratio of the two predicted activity tracks isolates what
transfer changed:

| ATAC decile | n positives / 861 | transferred / in-cell | transferred / observed p300 | in-cell / observed p300 | mean percentile shift, negatives |
|---|---|---|---|---|---|
| 1 | 5 | 0.64 | 0.45 | 0.70 | -0.028 |
| 5 | 23 | 0.88 | 0.59 | 0.67 | -0.009 |
| 9 | 100 | 1.10 | 0.83 | 0.76 | +0.018 |
| 10 | 75 | **2.16** | **1.12** | **0.52** | **+0.041** |

Transfer multiplies predicted activity by 0.64 at the least accessible decile and by 2.16 at the
most accessible one, a 3.4x swing that is monotone in accessibility. **The in-cell model
suppresses the top accessibility decile and the transferred model does not**: at decile 10 the
in-cell model predicts 0.52x observed p300 while the transferred model predicts 1.12x, so in
level terms the transferred model is the more accurate of the two exactly where it loses.
**The damage is promotion of negatives, not demotion of positives.** Entering the top 10% of
pairs, transfer admits 121 negatives against in-cell's 102, while 25 positives leave and only 5
arrive. Mean percentile shift on positives is -0.0020 against +0.0001 on negatives overall, but
+0.041 on negatives in the top ATAC decile.
**Why suppressing the top decile is right for the benchmark: the positive rate is not monotone
in accessibility.** Regulated pairs peak at decile 9 (100 of 861) and fall at decile 10 (75 of
861). A model whose predictions keep climbing with accessibility therefore ranks decile-10
negatives above decile-9 positives. Since ABC's qnorm is rank-based, level accuracy buys nothing
and this ordering error is the whole cost.

**The counts metrics see this, but only if measured in the right place.**

| p300, evaluated in K562 | in-cell | transferred | loss |
|---|---|---|---|
| overall Pearson, genome-wide held-out folds | 0.794 | 0.679 | 0.115 |
| top-quintile Pearson, genome-wide | 0.619 | 0.312 | 0.307 |
| Pearson on CRISPR-tested elements | 0.725 | 0.489 | **0.236** |
| Pearson on regulated-pair elements only | 0.747 | 0.487 | **0.260** |

The genome-wide overall Pearson understates the transfer damage on the elements the benchmark is
decided on by roughly a factor of two, and the transferred model's incremental R2 over
accessibility is NEGATIVE (-0.170 overall, -0.284 top-quintile) while its top-quintile Pearson,
0.312, is indistinguishable from the accessibility-only model's 0.306.
**Implications:** Three targeted strategies follow from the shape of the error rather than from
guesswork, in increasing cost. First, the error is monotone in accessibility, so a monotone
recalibration of the predicted track against target-cell-type ATAC deciles could remove most of
it, and the correction can be fitted in the SOURCE cell type where both assays exist, which keeps
it deployment-legal. Second, the transferred model normalises K562 ATAC with GM12878's stored
statistics (acc_mean 4.552, acc_std 1.360 against K562's 4.155 and 1.313), a known and free
confound to remove before concluding anything about the slope; F-017 established that input-side
matching, in that case depth, removed a -0.133 transfer collapse. Third, the residual objective was the obvious
candidate for this channel, but **checking before testing it downstream showed transfer had
already been measured upstream and is null**: `deploy_gm_to_k562` gives residual multimodal
0.823 overall / 0.600 top-quintile against plain multimodal's 0.828 / 0.602, floor 0.802 /
0.541. It also does not remove the accessibility channel so much as relocate it, since
prediction is residual plus an accessibility-only offset model that carries its own slope
error on transfer. Worth a downstream test only because the benchmark has dissociated from
Pearson three times, and worth stating as expected-null.
**Caveats:** The +0.0091 over the floor is not resolvable, so the headline is "best of the
transferred arms", not "works". The p300 and H3K27ac families were scored in separate comparison
runs, so cross-target AUPRCs are comparable only through the shared July anchors and are not
paired against each other; within-run pairings above are. The H3K27ac transfer arms show the same
pattern far more weakly (0.75 to 1.08 across deciles, against p300's 0.64 to 2.16) and, notably,
their counts accuracy on CRISPR elements does NOT degrade on transfer (0.672 to 0.691 overall,
0.515 to 0.621 on regulated-pair elements) even though both arms sit at the floor, so counts
accuracy and benchmark utility are dissociated there in the opposite direction. K562 is the only
cell type with CRISPR data, so "transfer" here always means into K562 and one source cell type.
Accessibility deciles are computed on the CRISPR-tested elements, not genome-wide. The GM12878
p300 model is from `2026_0606_GM12878_transferability/GM12878_multimodal_BPNet/models/atac`,
outside this analysis directory, and its `training_target.json` was never written; the target was
confirmed from that project's `config/input_data_gm12878_multimodal.json` (GM12878 EP300 peaks
ENCFF926AKK, signal ENCFF960OFK/ENCFF941MGK) and from the predicted track correlating 0.913 with
the K562 p300 prediction against 0.678 with the GM12878 H3K27ac prediction.

**Update 2026-09-17, same day: the accessibility-slope mechanism is WITHDRAWN. The models
barely have the defect; ABC's qnorm and the CRISPR subset manufacture it.** The decile table
above was computed on `normalized_h3k27ac_enh` in the ABC putative predictions, restricted to
CRISPR-tested element-gene pairs. Measuring the SAME transferred-over-in-cell ratio at two
earlier points in the pipeline (`4.28`, `4.29`) decomposes it:

| where the ratio is measured | decile 1 -> 10 | spread |
|---|---|---|
| raw predicted bigwigs, all 153,459 candidate regions | 0.475 -> 0.516 | **x1.20** |
| post-qnorm activity, all 153,349 elements | 0.806 -> 1.235 | **x1.85** |
| post-qnorm activity, CRISPR-tested pairs (the table above) | 0.64 -> 2.16 | **x3.40** |

On a log scale the model contributes about 15% of the swing, ABC's rank normalization about
another 35%, and restriction to the CRISPR-tested subset the remaining 50%. **The two models'
raw predictions differ by a nearly constant factor of about 0.48 with only 1.20x accessibility
dependence**, so "the transferred model reads accessibility too steeply" and "the in-cell model
suppresses the top decile" are not supported. Post-qnorm values are also coarsely quantized at
low accessibility (medians land on multiples of about 0.1024) and deciles 2, 4, 5 and 6 give a
ratio of exactly 1.0000 because both arms map to the same reference value, so the two arms are
indistinguishable there by construction.
**What survives.** The AUPRC numbers, the +0.0458 [+0.0296, +0.0610] transfer penalty, the
promotion of accessible negatives in the top decile, and the counts-accuracy table are all
untouched, because none of them depend on the mechanism claim. What changes is the explanation:
the transferred arm loses because its activity distribution has a different SHAPE from the
in-cell arm's, and rank normalization against a shared reference converts a difference in shape
into an accessibility-dependent reordering, which the benchmark then evaluates on the
accessible, gene-proximal elements where the remapping bites hardest.
**Consequence for the proposed fixes.** The monotone recalibration in the implications below
**cannot work as described**: a correction fitted on raw predictions in the source cell type
cannot repair a defect the raw predictions barely have. Any recalibration has to target the
post-qnorm activity, which means it is a statement about ABC's normalization rather than about
the model, and ABC re-derives that normalization per arm. The accessibility-input
renormalisation was still worth doing for its own reason, the measured 0.73 sd against 0.38 sd
centring asymmetry, and it moves the raw ratio from 0.475-0.516 to 0.669-0.648 without
flattening the spread (x1.32 against x1.20).
**Method lesson, recorded because it generalises.** A ratio measured on a benchmark's own input,
over the benchmark's own element subset, is not a property of the model that produced it. Two
normalizations and one subset selection sat between the model and the number, and each was
worth a factor. Measure at the model output first, then add one pipeline stage at a time.

**Tags:** p300, h3k27ac, abc, crispr-benchmark, transfer, accessibility, calibration, deployment, k562, methodology

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-17 | `4.25.transfer_error_anatomy.py`, `4.17` paired bootstrap over the 2026_0906_p300_all and 2026_0904_predicted_activity comparisons | EPCrisprBenchmark_ensemble_data_GRCh38, 10,342 element-gene pairs, 466 regulated; ABC putative predictions for 5 arms | 2026_0824_H3K27ac_model | transferred p300 multimodal 0.4771 vs 0.4680 floor (+0.0091 ns); transfer penalty +0.0458 [+0.0296, +0.0610]; transferred/in-cell activity ratio 0.64 to 2.16 across ATAC deciles; CRISPR-element Pearson 0.725 to 0.489 | extends F-005, F-008, F-016 |

---

## F-022: ABC's rank qnorm makes the predicted activity track's SCALE unusable as a lever, so any fix that only rescales a prediction is null by construction
**Status:** established, structurally and empirically
**Claim:** `run_qnorm(qnorm_method="rank")` in ABC's `neighborhoods.py` maps each element's
WITHIN-ARM rank through a linear interpolation onto a reference distribution, separately for
promoters and nonpromoters. The post-qnorm activity is therefore a function of the predicted
track's per-element ORDERING and of nothing else. Any transformation of the track that preserves
per-element order is discarded exactly.
Tested on the 2026-09-17 accessibility renormalisation, which changed the p300 tracks' scale by
x1.5 (in-cell) and x1.8 (transferred) while moving the per-element ranking by only 1.6% and 3.7%
(Spearman 0.9969 and 0.9853 against the originals):

| p300 arm | post-qnorm activity median, before -> after | activity Spearman | ABC.Score Spearman |
|---|---|---|---|
| in-cell | 0.5372 -> 0.5374 | 0.9947 | 0.9918 |
| transferred | 0.5369 -> 0.5369 | 0.9806 | 0.9831 |

**The scale change vanished to four decimal places.** What survives into ABC is only the small
rank perturbation, and the transferred arm's ranking moved about twice as far as the in-cell
arm's, consistent with its larger input-centring error (0.73 sd against 0.38 sd).
**Implications:** **This retires a whole class of proposed interventions.** Anything that
rescales, recentres, gamma-corrects or otherwise monotonically remaps a predicted activity track
is null in ABC by construction, and does not need an experiment to rule out. That covers the
accessibility renormalisation tested here, the monotone recalibration proposed and dropped in
F-021, and retrospectively explains F-013, where quantile mapping fixed a converter track's
marginal distribution and moved nothing downstream: quantile mapping is monotone, so ABC could
not see it. It also explains why F-004's dynamic-range diagnosis, correct as a description, was
never actionable in the form it was stated; the TODOLIST already recorded the principle as "ABC
qnorm removes scale by construction, so rank is the only channel available", but it was not being
applied to new proposals.
**What this leaves as real levers**, for a predicted activity track: the per-element RANKING, and
the SHAPE of the marginal distribution only insofar as two arms' differing shapes are remapped
differently against a shared reference. The geomean arms have one additional channel, because
`geomean(real ATAC, predicted)` combines values before qnorm sees the product, so there the
predicted track's magnitude does enter; the prediction-alone arms have no such channel.
**Caveats:** Monotone in the per-element COUNTED value, which is what qnorm consumes. A monotone
transformation of the per-base bigwig does not give a monotone transformation of the per-element
sum, so a per-base change can still reorder elements; the renormalisation here was a per-base
change and did move 1.9% of the transferred arm's ranking. `separate_promoters=True` is the
default, so the mapping is monotone within the promoter and nonpromoter classes but differs
between them, and a change that moves elements across that boundary is not covered. Measured on
the p300 arms only, though the argument is about ABC's code rather than about any model.
**Tags:** abc, qnorm, methodology, crispr-benchmark, dynamic-range, negative-result, calibration

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-17 | `4.26`/`4.4` renormalised predictions, ABC 44056013, `4.30.qnorm_invariance.py` (44066859) | 153,447 and 153,349 ABC candidate elements; 10.1M element-gene pairs per arm | 2026_0824_H3K27ac_model | x1.5 and x1.8 track rescaling produced post-qnorm activity medians identical to 4 dp (0.5372->0.5374, 0.5369->0.5369); only a 1.9% rank shift survived | explains F-013, retires F-021's recalibration proposal |

---

## F-023: Both cheap fixes proposed for the transferred ATAC-input arm are null downstream, and renormalising the accessibility input is mildly harmful in the geomean form
**Status:** established
**Claim:** F-021 proposed two low-cost interventions for the transferred p300 arm. Both were run
through ABC and the CRISPR benchmark and scored against the arms they were meant to improve in a
single 14-arm comparison, so every delta below is paired on the same 10,342 element-gene pairs
(466 regulated).

**Renormalising the accessibility input on the prediction regions** (`4.26` gives 3.7232 /
1.1428 for K562 ATAC over the ABC candidate regions, against the models' stored 4.155 / 1.313
in-cell and 4.552 / 1.360 transferred; both arms re-predicted):

| contrast | delta AUPRC | 95% CI | sign kept |
|---|---|---|---|
| transferred, geomean form | -0.0004 | [-0.0064, +0.0048] | 52.2% |
| **in-cell, geomean form** | **-0.0058** | **[-0.0085, -0.0029]** | **100.0%** |
| transferred, prediction alone | -0.0015 | [-0.0121, +0.0089] | 60.7% |
| in-cell, prediction alone | -0.0027 | [-0.0122, +0.0067] | 72.2% |

**The residual objective** (`residual5p_multimodal_hw500_clw10` and its GM12878 twin, each with
its own accessibility-only offset model, prediction = residual + offset):

| contrast | delta AUPRC | 95% CI | sign kept |
|---|---|---|---|
| transferred residual - transferred plain | -0.0011 | [-0.0114, +0.0096] | 59.3% |
| in-cell residual - in-cell plain | -0.0006 | [-0.0093, +0.0086] | 54.9% |
| transferred residual - ATAC-only floor | -0.0059 | [-0.0215, +0.0095] | 78.1% |

**Implications:** Neither intervention helps, and the accessibility renormalisation should NOT be
adopted: the only resolvable effect it has is to make the in-cell p300 geomean arm 0.0058 worse,
with 100% sign retention. **The one asymmetry consistent with F-022's stated exception is that
the resolvable harm lands in a geomean arm and not in a prediction-alone arm.** `geomean(real
ATAC, predicted)` combines values before qnorm sees the product, so the predicted track's
magnitude does enter there, and rescaling it by 1.5x shifts the balance between the two factors.
The transferred geomean arm is nonetheless null despite a larger 1.8x rescaling, so magnitude
alone does not predict the sign or size, and this should be treated as an observation rather
than a mechanism.
**F-021's headline is unchanged by either.** The transferred p300 arm sits at +0.0087 [-0.0048,
+0.0219] over the floor after renormalisation against +0.0091 [-0.0042, +0.0222] before, still
the best transferred ATAC-input arm in the project and still not resolvably above the floor.
**Both results were predicted before they were run, from different evidence**: the
renormalisation from F-022, because ABC's rank qnorm discards scale and only 1.6-3.7% of
per-element ranks moved; the residual objective from `deploy_gm_to_k562`, which had already
measured it upstream on transfer at 0.600 against 0.602 top-quintile Pearson. Two correct
predictions of a null is weak evidence that the reasoning is sound, and cheap: the expensive
half of this was the code to make a residual model predictable at all.
**Caveats:** K562 in-cell and one transfer direction, as always with this benchmark. Every arm
here is geomean(real ATAC tagAligns, predicted bigwig) or a single predicted bigwig, so they are
matched on counting path; the two July anchors are not and are reference only. The residual arms
use the SOURCE cell type's offset model, which is the only deployment-legal choice, so "residual
objective on transfer" here means the whole residual+offset object transferred, not a
target-fitted offset. The prediction-alone residual arms were built and scored but the decile
anatomy (`4.25`) was not re-run on any of the new arms, because F-022 implies the post-qnorm
activity is nearly unchanged (Spearman 0.98-0.99) and the AUPRC deltas confirm it; if that
anatomy is ever wanted it is one scoring pass.
**Tags:** p300, h3k27ac, abc, crispr-benchmark, transfer, calibration, residual-objective, negative-result, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-17 | predictions 44049689/44049692/44050508/44050541 (`4.4`), ABC 44056013 and 44056017 (`4.5`), benchmark 44078638 (`4.7`), paired bootstrap `4.17` | EPCrisprBenchmark_ensemble_data_GRCh38; 14 arms on 10,342 shared element-gene pairs, 466 regulated | 2026_0824_H3K27ac_model | accnorm -0.0004 transferred / -0.0058 in-cell geomean; residual -0.0011 transferred / -0.0006 in-cell; transferred p300 still +0.0087 over floor | closes F-021's two proposed fixes |

---

## F-024: Observed p300 is a resolvably better ABC activity term than observed H3K27ac, which raises the project's ceiling by 60% and leaves two equal, resolvable gaps
**Status:** established
**Claim:** Paired on the same 10,342 element-gene pairs, with ATAC in the other half of the
geomean throughout:

| arm | AUPRC | vs ATAC-only floor |
|---|---|---|
| ATAC only (floor) | 0.4680 | - |
| predicted p300, GM12878 -> K562 | 0.4771 | +0.0091 [-0.0042, +0.0222] |
| predicted H3K27ac, K562 in-cell | 0.4801 | +0.0121 [-0.0053, +0.0285] |
| predicted p300, K562 in-cell | 0.5229 | +0.0549 |
| observed H3K27ac | 0.5296 | +0.0616 |
| **observed p300** | **0.5673** | **+0.0993 [+0.0786, +0.1198]** |

**Observed p300 beats observed H3K27ac by +0.0377 [+0.0187, +0.0583], sign kept 99.9%.** The
activity term's ceiling is therefore +0.0993 over the floor, not the +0.0616 that the
H3K27ac anchor implies, and the project has been measuring itself against a ceiling 60% too low.
**Implications:** **Two gaps remain, they are almost exactly equal, and both are resolvable**,
which is new; nearly every contrast in this project until now has been unresolvable against this
benchmark.

| gap | size | resolvable? |
|---|---|---|
| in-cell predicted p300 against observed p300 | **0.0444** [-0.0633, -0.0264] | yes, 100% |
| transfer, in-cell against GM12878 -> K562 | **0.0458** [+0.0296, +0.0610] | yes, 100% (F-021) |

In-cell prediction captures 55% of the available headroom, and the transferred arm 9%. So the
in-cell problem is NOT nearly solved, which is what the H3K27ac ceiling made it look like.
**The strategic consequence is that p300 should be the default target.** Predicted p300 in-cell
clears the floor by +0.0549 while predicted H3K27ac does not clear it resolvably at all
(+0.0121, CI spans zero, F-004). The project is organised around H3K27ac by inheritance rather
than by evidence, and the better activity target has been sitting in the p300 arms since F-008.
**It also sets a minimum effect size worth chasing.** With 466 regulated pairs the paired CIs
run about +/-0.010 to +/-0.020, so an intervention expected to move less than roughly 0.015
cannot be evaluated here whatever its merit. Both gaps above are comfortably above that; the
epoch budget (+0.0086) and the multi-task auxiliary loss (+0.0047) never were.
**Caveats:** K562 only, as always. The observed-p300 arm is geomean(real ATAC tagAligns,
observed p300 BAMs) and so shares the counting path of the predicted arms, unlike the July
anchors; this contrast is therefore legitimate, and it is the reason the number is quotable
where a comparison against the July floor would not be. A higher ceiling does not mean a
predicted track can reach it: observed p300 has information no ATAC-input model can recover.
Whether p300's advantage holds in another cell type is untested, and no GM12878 CRISPR benchmark
exists to test it.
**Tags:** p300, h3k27ac, abc, crispr-benchmark, activity-term, ceiling, strategy, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-17 | `4.17` over the 2026_0906_p300_all comparison | EPCrisprBenchmark_ensemble_data_GRCh38; 10,342 shared element-gene pairs, 466 regulated | 2026_0824_H3K27ac_model | observed p300 0.5673 against observed H3K27ac 0.5296 (+0.0377) and floor 0.4680 (+0.0993); in-cell predicted p300 short of its ceiling by 0.0444 | reframes F-004, F-008, F-021 |

---

## F-025: The GM12878 p300 model was trained on a BPNet model's PREDICTED plus strand paired with the OBSERVED minus strand, so every result that uses it is compromised
**Status:** established; the provenance error is real, but the retrain (2026-09-19) shows its downstream cost is only +0.0044 and **the "conflates transfer with target corruption" claim below is WITHDRAWN**
**Claim:** `2026_0606_GM12878_transferability/GM12878_multimodal_BPNet/models/atac`, the only
GM12878 p300 model in the project, was trained with
`signal_plus_bw = ENCFF960OFK_plus.bw` and `signal_minus_bw = ENCFF941MGK_minus.bw`. Those two
files are not the two strands of one experiment. Both belong to **annotation ENCSR038OGP, a
BPNet-model annotation** from the Kundaje lab, not to the EP300 experiment ENCSR000DZG, and
their ENCODE `output_type` fields are:

| accession | output_type | used as |
|---|---|---|
| ENCFF960OFK | **predicted** signal profile (plus strand) | plus strand of the target |
| ENCFF941MGK | observed signal profile (minus strand) | minus strand of the target |

The annotation contains the matching observed plus strand, **ENCFF557UDP**, one row away in
the same file list. So the model's target was a BPNet model's OUTPUT on one strand and real
data on the other.

**Confirmed numerically.** Rebuilding the GM12878 EP300 5' plus track from the actual
experiment's BAMs (ENCSR000DZG: ENCFF515HYM, ENCFF215GSQ) with `0.45` reproduces the correct
observed track exactly, and neither matches what was used:

| comparison, chr8 EP300 peaks, 1 bp | Pearson |
|---|---|
| `0.45` rebuild vs ENCFF557UDP (observed plus) | **1.0000** |
| `0.45` rebuild vs ENCFF960OFK (predicted plus) | 0.2734 |
| ENCFF557UDP vs ENCFF960OFK | 0.2734 |

Means agree at 0.0607 and maxima at 7.0 for the two observed tracks, while the predicted track
peaks at 2.7, which is the smoothing a profile model produces.
**Implications:** **F-021's headline arm is this model.** "The best ATAC-input model on the
transferred CRISPR benchmark predicts p300", AUPRC 0.4771 against a 0.4680 floor, was produced
by a model whose training target was half model-output. The AUPRC is still a real measurement
of that predictor, but it cannot be described as a GM12878 p300 model, and **the +0.0458
transfer penalty conflates transfer with target corruption** and should not be quoted as a
transfer effect. F-023's transferred p300 arms, including both accnorm arms, inherit the same
model and the same caveat; its conclusion that the two fixes are null is unaffected, since
those were within-model comparisons.
**What is NOT affected.** Everything K562-trained, which is F-004, F-008, F-018, F-019, F-020
and F-022, because the K562 p300 target comes from the experiment's own BAMs. The 2026-09-08
depth-subsampling test trained on K562 and tested into GM12878, so its MODEL is clean; whether
its GM12878 evaluation target was this same file has not been checked and should be.
**How it happened, and the general lesson.** ENCODE annotations of type `BPNet-model` publish
observed and predicted profiles side by side with near-identical names, differing only by the
words "observed" and "predicted" in `output_type`, which does not appear in the filename. The
project's own convention of naming local copies `<accession>_plus.bw` discards exactly the
field that distinguishes them. **Record `output_type` and `dataset` alongside any accession,
and never take a signal track from an `/annotations/` dataset when an `/experiments/` one
exists.**
**Caveats:** The model is a real predictor and its downstream numbers are real; what is wrong
is the label and the interpretation, not the arithmetic. Only the plus strand is predicted, so
roughly half the target is genuine, which may be why the model works at all. Retraining
GM12878 p300 on the corrected target is now cheap, since `0.45` has already built the correct
tracks for the multi-cell-type panel.

**Update 2026-09-19: the model was retrained on the corrected target and benchmarked. The
provenance bug is real but its downstream cost is small, and F-021's headline SURVIVES.**
Retrained with every hyperparameter, the peaks, negatives, accessibility track and folds
copied from the original run, so the target was the only difference. Predicted with the same
accessibility track and the same default normalisation as the original arm. Paired on 10,342
element-gene pairs:

| contrast | delta AUPRC | 95% CI | sign kept |
|---|---|---|---|
| corrected - corrupted | +0.0044 | [-0.0042, +0.0123] | 84.3% |
| **corrected - ATAC-only floor** | **+0.0135** | [-0.0003, +0.0270] | **97.3%** |
| corrupted - floor (the F-021 number) | +0.0091 | [-0.0042, +0.0222] | 91.4% |
| in-cell - corrected (transfer penalty) | **+0.0414** | [+0.0256, +0.0558] | 100% |

**A CLAIM IN THIS FINDING IS WITHDRAWN.** It said the +0.0458 transfer penalty "conflates
transfer with target corruption and should not be quoted as a transfer effect". Fixing the
target moved that penalty only to +0.0414, so corruption accounted for roughly 0.004 of it
and the penalty is a transfer effect after all. The caution was reasonable before the test
and wrong after it; quote +0.0414 and treat the transferred p300 arm's transfer loss as real.
**What stands.** The provenance error itself, the identity problem (the old arm cannot be
called a GM12878 p300 model), and the lesson about `output_type` on annotation datasets.
F-021's ranking also stands: the transferred p300 arm is still the best transferred
ATAC-input arm in the project, now at +0.0135 over the floor rather than +0.0091, and still
not resolvable because the interval grazes zero at -0.0003.
**The direction of the correction is itself another upstream/downstream dissociation.** The
corrected model is WORSE on its own validation, count Pearson 0.75-0.78 against the original's
0.80-0.82, because the original's predicted plus strand was smoothed and therefore easier to
fit. It is nonetheless BETTER downstream, by +0.0044. **No upstream metric could have caught
this bug, and the one that exists pointed the wrong way.** That is the fifth time the
benchmark and the upstream metrics have disagreed in this project (F-004, F-018, F-019,
F-022's retrospective reading of F-013, and now this).
**Magnitude check, for calibration.** The corrected and corrupted predictions correlate at
Spearman 0.939 over the 153,545 candidate regions, a mean per-element rank shift of 7.2%.
Per F-022 only that reordering can reach ABC, since qnorm discards the 1.18x scale change.
For reference the accessibility renormalisation shifted 3.7% of ranks and moved AUPRC by
-0.0004, so +0.0044 from twice the reordering is the expected order of magnitude.

**Tags:** p300, gm12878, data-provenance, transfer, methodology, negative-result, encode

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-18 | `0.45` rebuild (44164220), `0.46` construction control (44181538), ENCODE metadata for ENCSR038OGP | GM12878 EP300 ENCSR000DZG BAMs; annotation ENCSR038OGP bigwigs; 1,052 chr8 peaks | 2026_0824_H3K27ac_model | rebuild matches observed plus at r=1.0000 and the used track at r=0.2734; used track is output_type "predicted signal profile (plus strand)" | invalidates the identity of F-021's transferred arm |

---

## F-026: Training p300 across four cell types beats single-source transfer and is the first model in this project to clear an assay-matched floor on transfer
**Status:** established upstream; the advantage over single-source DOES NOT REPLICATE downstream, see F-027 (2026-09-22)
**Claim:** A multimodal p300 model trained on A549, GM12878, HepG2 and MCF-7 pooled, holding
K562 out entirely, was scored on K562 p300 over K562's own EP300 peaks with K562 ATAC as the
accessibility input. Baseline is K562's own accessibility-only p300 model, the assay-matched
floor F-016 requires. Five chromosome folds, paired within fold.

| arm | overall Pearson | top-quintile |
|---|---|---|
| K562 in-cell (upper reference) | 0.691 [0.671, 0.710] | 0.482 [0.430, 0.534] |
| **multi-cell, K562 held out** | **0.574 [0.538, 0.609]** | **0.317 [0.222, 0.412]** |
| GM12878 single-source, target-normalised | 0.530 [0.499, 0.561] | 0.266 [0.178, 0.355] |
| GM12878 single-source, own statistics | 0.518 [0.489, 0.547] | 0.253 [0.162, 0.345] |
| K562 accessibility-only FLOOR | 0.478 [0.450, 0.507] | 0.256 [0.165, 0.346] |

| paired contrast | overall | top-quintile |
|---|---|---|
| **multi-cell - single-source, same normalisation** | **+0.0436** (*p*=0.0092) | **+0.0508** (*p*=0.0004) |
| multi-cell - single-source, own statistics | +0.0557 (*p*=0.0042) | +0.0636 (*p*=0.0003) |
| **multi-cell - assay-matched floor** | **+0.0954** (*p*=0.0002) | **+0.0614** (*p*=0.036) |
| single-source - assay-matched floor | +0.0396 (*p*=0.016) | **-0.0022** (*p*=0.91) |
| in-cell - multi-cell | +0.1169 (*p*=0.0004) | +0.1650 (*p*=0.006) |

**Implications:** **F-016 found no model beating an assay-matched floor in six transfer
directions, and single-source transfer fails here too on top-quintile (-0.0022, *p*=0.91).
Pooling four source cell types clears that floor (+0.0614, *p*=0.036).** This is the first
positive transfer result in the project. In headroom terms, against the in-cell model as the
attainable ceiling, single-source recovers 4% of the top-quintile gap and multi-cell 27%; on
overall Pearson, 19% against 45%.
The comparison that carries the claim is the policy-matched one. A multi-cell model has no
single training statistic, because accessibility is standardised per cell type so library
depth cannot enter the input, so it MUST be normalised on the target. Scoring the
single-source model under the same policy costs it nothing to arrange and removes
normalisation as an explanation: the gap narrows from +0.0557 to +0.0436 overall, so about a
fifth of the raw difference was policy and the rest is the training set.
**Caveats:** **Pooling more cell types also pools more data, and this design does not
separate them.** The multi-cell model saw 34,669 training peaks against GM12878's 14,602, so
"four cell types" is confounded with "2.4x the peaks". The 2026-09-08 subsampling test is
weak evidence against volume mattering for p300 transfer, since depth-matching K562 to
GM12878's budget left transfer unchanged (+0.202 against +0.207), but that varied read depth
rather than peak count and was single-source. **The clean control is a single-source model
trained on a peak count matched to the pool**, which no single cell type here can supply;
the nearest feasible version is pooling two cell types and then four, and seeing whether the
gain tracks cell-type count or peak count.
Upstream only. The multi-cell model has not been through ABC or the CRISPR benchmark, and
this project has had the benchmark and Pearson disagree five times, so the downstream result
should not be assumed from this. One held-out cell type, K562, because it is the only one
with an assay-matched accessibility-only p300 floor; the other four leave-one-out models can
be scored against each other but not against a floor without 15 more training runs. GM12878
enters the pool with the F-025-corrected target. Single run per condition, and F-020 puts
run-to-run variance at sd 0.0018 count Pearson, well below these effects.
**Tags:** p300, transfer, multi-cell-type, deployment, k562, positive-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-19 | training 44301326 (25 tasks), eval 44326695 (`2.15` via `2.43`) | K562 EP300 peaks (ENCSR000EGE inliers); panel A549/GM12878/HepG2/MCF-7, 34,669 pooled training peaks | 2026_0824_H3K27ac_model | multi-cell 0.317 top-quintile against single-source 0.266 and floor 0.256; +0.0508 paired (*p*=0.0004) and +0.0614 over the floor (*p*=0.036) | first transfer result to clear an assay-matched floor; extends F-016 |

## F-027: The multi-cell-type advantage does not survive the CRISPR benchmark, and both transferred arms sit at the floor
**Status:** established
**Claim:** The F-026 multi-cell p300 model (A549+GM12878+HepG2+MCF-7, K562 held out) was
painted onto K562's ABC candidate regions, run through ABC, and benchmarked against the K562
CRISPR data alongside a policy-matched single-source comparator. Upstream the multi-cell
model beat that comparator by +0.0508 top-quintile Pearson (*p*=0.0004). Downstream the same
contrast is **-0.0020 [-0.0128, +0.0079]**, with the sign retained in 62.5% of resamples.
Paired bootstrap, 2,000 resamples, 10,342 element-gene pairs scored by every arm, 466
regulated (4.51%).

| arm | AUPRC (shared pair set) |
|---|---|
| ATAC x observed p300 (ceiling) | 0.5673 |
| ATAC x predicted p300, K562 in-cell | 0.5229 |
| GM12878 single-source, target-derived norm | 0.4829 |
| GM12878 single-source, own statistics | 0.4815 |
| **multi-cell, K562 held out** | **0.4810** |
| ATAC only (floor) | 0.4680 |

| paired contrast | delta AUPRC | 95% CI | sign kept |
|---|---|---|---|
| **multi-cell - single-source, same normalisation** | **-0.0020** | [-0.0128, +0.0079] | 62.5% |
| multi-cell - floor | +0.0130 | [-0.0041, +0.0287] | 94.2% |
| single-source (target norm) - floor | +0.0149 | [+0.0012, +0.0278] | 98.1% |
| target-derived norm - own statistics | +0.0014 | [-0.0042, +0.0074] | 67.5% |
| K562 in-cell - floor | +0.0549 | [+0.0343, +0.0741] | 100% |
| observed p300 - floor | +0.0993 | [+0.0786, +0.1198] | 100% |

**Implications:** **This is the sixth time the CRISPR benchmark has disagreed with a
correlation metric, and the first time it has erased a result that looked like the project's
headline.** F-026's upstream claim stands as measured; what fails is the inference from it.
The two transferred arms are indistinguishable from each other and barely separated from the
floor, at roughly a quarter of the in-cell model's margin and an eighth of the observed-p300
ceiling. **Pooling source cell types is not, on this evidence, a route to deployable
transfer**, and the 2-vs-4 cell ladder proposed in F-026 to break the data-volume confound is
no longer worth running for its own sake: there is no downstream effect left to attribute.
The normalisation policy is separately worth nothing (+0.0014, sign kept 67.5%). That does
NOT make the policy-matched comparator wasted work -- it is the measurement that rules
normalisation out as an explanation for the upstream gap, and without it the downstream null
would be ambiguous between "multi-cell does not help" and "the two arms were normalised
differently".
**Caveats:** One benchmark, one held-out cell type, and AUPRC on 466 regulated pairs is a
coarse instrument: the floor-to-ceiling span is only 0.0993, so a 0.002 difference is well
inside the noise of what this benchmark can resolve. The arms share one pair set and one
resample per draw, so the deltas are paired, but the CIs still overlap heavily for every
transferred arm. **This does not overturn F-024** (observed p300 beats observed H3K27ac as
the activity term); that contrast is between OBSERVED tracks and is unaffected.
**Tags:** p300, transfer, multi-cell-type, crispr-benchmark, k562, negative-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-22 | predictions 44736068/44736070, ABC 44748720, benchmark 44763411, bootstrap `4.17` | K562 ABC candidate regions (153,545) and K562 CRISPR pairs; multi-cell panel A549/GM12878/HepG2/MCF-7 | 2026_0824_H3K27ac_model | multi-cell 0.4810 AUPRC against single-source 0.4829 and floor 0.4680; -0.0020 [-0.0128, +0.0079] paired | contradicts the downstream reading of F-026; extends the Pearson/benchmark disagreement to six |

## F-028: Most of the gap between predicted and observed p300 in ABC comes from how the CRISPR-tested elements are ranked
**Status:** established
**Claim:** Oracle substitution. For each model, a hybrid activity track was built (`4.31`)
that gives the 3,406 ABC candidate regions overlapping a CRISPR-tested element their
OBSERVED p300 within-class rank, and keeps the model's predicted values at the other 150,139
regions. Ranks rather than values, because ABC's rank qnorm keeps only ranks (F-022). The
hybrid was run through ABC and the benchmark alongside the all-predicted and all-observed
arms, on one shared pair set (10,342 pairs, 466 regulated), paired bootstrap 2,000 resamples.

| arm | AUPRC |
|---|---|
| C: observed p300 everywhere | 0.5673 |
| B': in-cell model, observed ranks at tested elements | 0.5504 |
| B: multi-cell model, observed ranks at tested elements | 0.5502 |
| A': in-cell model, predicted everywhere | 0.5229 |
| A: multi-cell model, predicted everywhere | 0.4810 |
| ATAC only (floor) | 0.4680 |

| contrast | delta AUPRC | 95% CI | share of gap |
|---|---|---|---|
| multi-cell: B - A (tested elements) | +0.0692 | [+0.0511, +0.0868] | 80% |
| multi-cell: C - B (all other elements) | +0.0171 | [+0.0047, +0.0312] | 20% |
| in-cell: B' - A' (tested elements) | +0.0276 | [+0.0126, +0.0430] | 62% |
| in-cell: C - B' (all other elements) | +0.0169 | [+0.0052, +0.0289] | 38% |

**Implications:** **The activity term fails mainly on the ordering of the elements the
benchmark scores, not on the surrounding elements that enter ABC's denominator and qnorm
reference.** A better activity model can in principle close most of the gap, so work aimed
at element ranking (the ranking loss, better inputs) is aimed at the right quantity. Both
oracle arms land at the same 0.550 although the models start 0.042 apart: **the whole
difference between the in-cell and the transferred model is in how they rank the tested
elements**; on the other 97% of candidates they cost ABC the same ~0.017. That ~0.017 is a
floor on what a model that gets only the tested elements right can reach.
**Caveats:** An oracle, not a model: it uses observed p300 at the tested elements and is
never deployable. The tested elements are, by construction, where AUPRC is computed, so B - A
is the upper bound on what fixing their ranks can buy, not an estimate of what a real model
will. Two models; the decomposition could differ for an activity track with a very different
error structure.
**Tags:** p300, abc, crispr-benchmark, oracle, diagnostic, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-23 | oracle tracks `4.31`, ABC 44913144, benchmark 44930707, bootstrap `4.17` | K562 ABC candidate regions and K562 CRISPR pairs | 2026_0824_H3K27ac_model | tested-element ranks carry 80% (transfer) and 62% (in-cell) of the predicted-to-observed gap | locates the F-027 loss; supports ranking-targeted work |

## F-029: In ENCODE-rE2G, observed p300 adds +0.12 AUPRC over ATAC alone and predicted p300 adds nothing, because the predictions carry little beyond ATAC
**Status:** established for the rE2G result; the explanation is supported, not established
**Claim:** ENCODE-rE2G (ATAC, powerlaw contact) retrained on the K562 CRISPR data only, with
every model on the ALL-features table (10 features ATAC-only, 14 with an activity track) and
the activity slot filled by each p300 track. All five datasets share the July ABC candidate
region set (`4.32`), so the comparison is paired. Leave-one-chromosome-out predictions,
paired bootstrap over all 10,356 K562 pairs (471 regulated), scored exactly as rE2G's
`compare_all_models.py` scores them (`4.34`; reproduces rE2G's table to within 0.0014).

| model | AUPRC | delta vs ATAC-only | 95% CI | sign kept |
|---|---|---|---|---|
| ATAC + observed p300 | 0.6660 | **+0.1225** | [+0.0884, +0.1517] | 100% |
| ATAC + observed H3K27ac | 0.6266 | +0.0831 | [+0.0537, +0.1093] | 100% |
| ATAC + predicted p300, multi-cell | 0.5516 | +0.0082 | [-0.0106, +0.0230] | 78% |
| ATAC + predicted p300, K562 in-cell | 0.5484 | +0.0050 | [-0.0289, +0.0365] | 61% |
| ATAC only | 0.5435 | | | |

Observed p300 beats observed H3K27ac by +0.0394 [+0.0107, +0.0679], **replicating F-024 in a
second framework.**

**Why the predicted tracks add nothing (`4.35`).** Partial rank correlation of each activity
track with observed p300 after removing ATAC, over all 153,545 candidate regions, within
ABC's classes:

| track | rho with obs p300 | rho with ATAC | partial rho given ATAC (non-promoter) | (promoter) |
|---|---|---|---|---|
| ATAC itself | 0.628 | 1 | 0 | 0 |
| observed H3K27ac (reference) | 0.698 | 0.590 | 0.522 | 0.600 |
| predicted p300, K562 in-cell | 0.675 | 0.675 | 0.437 | 0.368 |
| predicted p300, multi-cell | 0.623 | 0.770 | 0.281 | 0.318 |
| predicted p300, GM12878 target norm | 0.573 | 0.793 | 0.157 | 0.280 |
| predicted p300, GM12878 own stats | 0.562 | 0.798 | 0.130 | 0.230 |

(non-promoter rows shown for the first two columns; 132,186 non-promoter and 21,359
promoter regions.)

**Implications:** **The transferred predictions correlate with observed p300 no better than
raw ATAC does (0.56-0.62 against 0.63) and correlate with ATAC at 0.77-0.80: they are mostly
ATAC re-encoded.** What they carry beyond ATAC is partial rho 0.13-0.28, against 0.44 for the
in-cell model and 0.52 for observed H3K27ac. A supervised model that already has ATAC
features has little to gain from them, which is the rE2G result. The partial correlation is
a leak-free upstream metric that orders the arms the way the benchmark does at the extremes
(in-cell far above transfer); within the transferred arms it separates multi-cell from
single-source where the benchmark cannot resolve a difference (F-027).
**It does not explain everything.** The in-cell model's partial rho (0.437 non-promoter) is
not far below observed H3K27ac's (0.522), yet in rE2G one adds +0.005 and the other +0.083.
So agreement with observed p300 beyond ATAC, at the element, is not the whole of what rE2G
extracts from a real activity track. The largest gap between the two is at PROMOTERS (0.368
against 0.600), and rE2G uses promoter activity of the target gene as a feature; that is the
leading hypothesis, testable by ablating rE2G's promoter-activity features.
**Caveats:** One cell type, one CRISPR dataset. The partial correlation controls for ATAC at
the element only, not for rE2G's other features (neighbourhood sums, contact, distance), so
some of the in-cell model's beyond-ATAC signal may be redundant with those. Comparing a
predicted-p300 track with observed H3K27ac against observed p300 mixes assays; H3K27ac is a
reference for "what a real activity track looks like", not a like-for-like comparator.
**Tags:** p300, re2g, supervised-e2g, crispr-benchmark, k562, negative-result, metric

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-23 | rE2G configs `4.32`, driver 44914391 (`4.33`), bootstrap `4.34`, partial correlation `4.35` | K562 CRISPR ensemble (10,356 pairs), K562 ABC candidate regions | ENCODE_rE2G results/2026_0923_p300_features; 2026_0824_H3K27ac_model | observed p300 +0.1225, predicted p300 +0.0050 to +0.0082 over ATAC-only rE2G; transferred predictions partial rho 0.13-0.28 given ATAC | replicates F-024; predicted activity not useful in a supervised framework as built |

## F-030: A within-cell-type ranking loss gives the transferred p300 model its first resolvable benchmark gain, peaking at weight 10
**Status:** supported (one training run per weight)
**Claim:** The multi-cell p300 model (K562 held out) retrained with a RankNet pairwise term added
to the counts MSE (`_pairwise_rank_loss`, pairs within a cell type, `1.36`), five folds per
weight, everything else identical to the MSE-only comparator. Paired bootstrap on one pair set
per framework (ABC: 10,342 pairs, 466 regulated; ENCODE-rE2G: 10,356 pairs, 471 regulated).

| rank weight | ABC delta vs MSE-only | rE2G delta vs MSE-only | best epoch | val count Pearson |
|---|---|---|---|---|
| 0 (MSE only) | | | 23.6 | 0.7129 |
| 2 | +0.0040 [-0.0013, +0.0095] | **+0.0143** [+0.0048, +0.0250] | 25.4 | 0.7149 |
| 10 | **+0.0087** [+0.0006, +0.0179] | **+0.0172** [+0.0072, +0.0293] | 23.0 | 0.7225 |
| 20 | +0.0034 [-0.0028, +0.0099] | -0.0060 [-0.0178, +0.0053] | 25.6 | 0.7275 |
| 50 | +0.0020 [-0.0068, +0.0117] | +0.0032 [-0.0076, +0.0153] | 22.0 | 0.7193 |
| 100 | -0.0003 [-0.0094, +0.0098] | +0.0028 [-0.0084, +0.0157] | 14.4 | 0.7198 |

Weight 10 against the ATAC-only floor: +0.0217 [+0.0049, +0.0374] in ABC and +0.0255
[+0.0066, +0.0413] in rE2G; the MSE-only model clears neither (+0.0130 [-0.0041, +0.0287];
+0.0082 [-0.0106, +0.0230]). Weights 20-100 are below weight 10 in both frameworks (rE2G
-0.0140 to -0.0232, all resolvable; ABC -0.0053 to -0.0090, resolvable only for weight 100).
**Implications:** **The first training change for the transferred model that moves the benchmark
in both frameworks, and the first transferred p300 model to clear the floor in both.** It is
small (about +0.01 AUPRC, a quarter of the 0.042 in-cell gap) and does not touch the TF-identity
errors of F-032: the weight-10 model's GATA and CEBP residual coefficients match the MSE-only
model's. Weight 20 had the best validation Pearson and is among the worst downstream: the
upstream metric picked the wrong model again. Partial rho given ATAC (4.35) PREDICTED
NO GAIN (0.271 and 0.265 against 0.281 over all non-promoter candidates, with rho with ATAC rising 0.770 ->
0.783 and 0.800); the gain sits at the tested elements (weight 10's tested-element partial rho
0.266 -> 0.295), which the all-candidate screen cannot see. Report 2 was updated to say a flat
screen score does not rule out a gain.
**Caveats:** One training run per weight. The non-monotone rE2G curve (weight 20 below 50 and
100) is the size of the run-to-run variation F-020 predicts, so the weight-2/10 gains may be
partly the luck of those runs; a seed replicate of weight 10 would settle it. Weight 100 stops
~9 epochs earlier than the rest, a training-length confound (F-019) for that row only. The
predicted-p300-alone arms do not agree (weight 2 -0.0071, weight 10 +0.0103, both n.s.).
**Tags:** p300, ranking-loss, transfer, multicell, abc, re2g, crispr-benchmark, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-24 | training 44916433 (w 2, 10); ABC 45081577, CRISPR 45101720, rE2G 45081983 | K562 CRISPR (ABC and rE2G pair sets) | 2026_0824_H3K27ac_model; ENCODE_rE2G 2026_0923_p300_features | w10 +0.0087 ABC, +0.0172 rE2G over MSE-only | first resolvable transfer gain |
| 2026-09-26 | sweep training 45134216 / 45233340 (w 20, 50, 100); ABC 45253419, CRISPR 45395938, rE2G 45261496 | same | same | gain peaks at w10, falls back above it | bounds the effect at ~+0.01 |

## F-031: Power-law neighbour information raises partial correlation with p300 but adds nothing to ABC, even fitted on K562's own p300
**Status:** established for ABC; not tested in rE2G
**Claim:** For each ABC candidate region, the power-law-weighted (ABC's gamma 1.024, scale 5.96,
5 kb floor) mean and log-sum of the predictions and ATAC of every other region within 1 Mb (median
149 neighbours), fitted to K562 observed p300 with 5-fold chromosome CV (`4.36`), linearly and by
gradient boosting, against the same fit without neighbour terms. Two graphs: within ABC class, and
across classes with terms split by the neighbour's class. Reordered tracks (`4.37`) run through ABC
and the benchmark.

| non-promoter partial rho given ATAC | boosting, no neighbours | + same-class neighbours | + cross-class |
|---|---|---|---|
| K562 in-cell model | 0.437 | 0.500 | 0.516 |
| multi-cell model | 0.291 | 0.377 | 0.397 |
| ATAC only | -0.014 | 0.119 | 0.204 |

| ABC delta | same-class | cross-class |
|---|---|---|
| multi-cell: with vs without neighbours | -0.0016 [-0.0104, +0.0070] | -0.0010 [-0.0089, +0.0077] |
| in-cell: with vs without neighbours | -0.0030 [-0.0091, +0.0029] | +0.0039 [-0.0029, +0.0115] |
| ATAC + neighbour ATAC vs ATAC-only floor | | -0.0041 [-0.0147, +0.0067] |

**Implications:** **Neighbour information is real for predicting p300 and useless to ABC.** The
likely reason is ABC's own formula: an element's score is its activity x contact over the sum of
activity x contact of every candidate near the gene, so boosting elements in active neighbourhoods
raises the denominator by the same route. A graph network feeding ABC's activity term is therefore
not worth building. Within-locus ranking (`local_rho`) improved slightly in every fit, so this is
not smoothing. **The partial-correlation screen failed a second time, in the opposite direction to
F-030**: +0.09 partial rho, zero AUPRC. Refitting to K562 labels without neighbours also bought
nothing (+0.0031 multi-cell, -0.0083 in-cell).
**Caveats:** An upper bound by design (fitted on K562 labels); a deployable version would be
weaker. rE2G, whose logistic regression could weight neighbour features separately from activity,
was not tested; it already has neighbourhood features.
**Tags:** p300, graph, neighbours, abc, negative-result, metric, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-25 | `4.36` same-class and cross-class; ABC 45233335 / 45141432; CRISPR 45254263 / 45395936 | K562 candidate regions, K562 CRISPR | 2026_0824_H3K27ac_model | +0.06 to +0.22 partial rho, no AUPRC change | closes the neighbour line for ABC |

## F-032: The transferred p300 model mis-ranks elements by which TFs are active in K562, and GATA elements carry about half the transfer gap
**Status:** established
**Claim:** Per candidate region, within ABC class, residual r = rank(observed p300) - rank(predicted)
in percentile points, regressed on each JASPAR 2024 motif (FIMO p < 1e-4, centre +/- 250 bp) and on
ENCODE K562 peak overlaps, partialling out ATAC rank (+ square), GC, CpG o/e and log distance to TSS
(`4.41`, 153,545 regions). Transfer-specific = r_multicell - r_incell.

| non-promoter | regions | in-cell | multi-cell | transfer-specific |
|---|---|---|---|---|
| GATA2 | 12,532 | -2.65 | +6.87 | **+9.52** |
| GATA1::TAL1 | 17,154 | -3.18 | +5.02 | **+8.20** |
| AP-1 (FOS::JUN) | 15,266 | +1.25 | +2.96 | +1.71 |
| STAT5a::STAT5b | 19,124 | -0.53 | +1.24 | +1.77 |
| CEBPB | 2,630 | +0.39 | -4.90 | **-5.29** |
| ESR2 | 18,752 | +0.25 | -4.47 | **-4.71** |
| EP300 peak | 19,660 | +12.81 | +21.14 | +8.33 |
| H3K4me1 peak | 52,470 | +3.63 | +8.60 | +4.97 |

The top seven transfer-specific motifs are all GATA-family (GATA2-6, GATA1::TAL1, TRPS1);
promoters show the same at about half the size (GATA5 +4.72). The ranking-loss model (F-030) has
the same coefficients as the MSE-only one (GATA2 +6.79 against +6.87).

GATA oracle (`4.42`): the multi-cell track with the 43,513 GATA-motif regions (28.3% of
candidates) given another track's within-class ranks, against the same number of ATAC-decile- and
class-matched non-GATA regions given the same ranks, ABC and benchmark, paired against the
multi-cell model:

| ranks from | GATA regions | matched non-GATA | GATA-specific |
|---|---|---|---|
| in-cell model | +0.0186 [+0.0108, +0.0266] | +0.0083 [+0.0003, +0.0156] | +0.0103 [-0.0001, +0.0213] |
| observed p300 | +0.0414 [+0.0295, +0.0538] | +0.0197 [+0.0097, +0.0296] | **+0.0217** [+0.0060, +0.0385] |

(in-cell model over multi-cell: +0.0419 [+0.0278, +0.0556].)
**Implications:** **The missing signal is TF identity.** The transferred model has learned that
the training cells' TFs recruit p300 (CEBP from HepG2 and A549, ESR from MCF-7) and applies that
in K562, where they are weak, while never having seen GATA1/TAL1 recruit p300, since no panel cell
is erythroid. The in-cell model shows neither pattern, so the motifs are not hard to predict; the
training cells taught the wrong mapping. Ranking GATA elements as observed p300 does recovers about
twice what the same number of matched elements does, and the GATA-specific part (+0.022) is about
half the transfer gap; ranking them as the in-cell model does gives +0.010, a quarter. **This
explains why no intervention so far moved the transferred model by more than ~0.01: none gave it
information about the target cell's TFs.**
**Caveats:** FIMO hits are sequence matches, not binding; eight near-identical GATA matrices define
the GATA set. An oracle: not deployable. One cell type; the analogous test for another held-out
cell type needs a CRISPR benchmark there.
**Tags:** p300, transfer, tf-identity, gata, motifs, error-analysis, oracle, abc, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-25 | `4.41` job 45293000 | K562 candidate regions, JASPAR 2024, ENCODE K562 peaks | 2026_0824_H3K27ac_model | GATA +9.5, CEBP/ESR -4.7 to -5.3 transfer-specific | names the missing factor |
| 2026-09-27 | `4.42`; ABC 45396268 / 45440340; CRISPR 45437389 / 45540585 | K562 CRISPR | 2026_0824_H3K27ac_model | GATA-specific +0.0217 (observed ranks), +0.0103 (in-cell ranks) | GATA is about half the transfer gap |

## F-033: ATAC and sequence recover the TF activities that matter for transfer to K562, but not the whole motif ranking
**Status:** supported
**Claim:** In each of the five p300-panel cell types, over the union of their ABC candidate regions
(`4.43`), the TRUTH is each motif's coefficient for p300 beyond ATAC (from that cell's own p300),
against four estimates from ATAC and sequence only: within-cell ATAC association, cell-specific
association (minus the other cells' mean), a discriminative score (open here, not elsewhere), and a
chromVAR-style deviation.

| estimate vs truth (beyond ATAC), Spearman across ~780 motifs | K562 | GM12878 | HepG2 | A549 | MCF-7 |
|---|---|---|---|---|---|
| within-cell ATAC association | 0.573 | 0.513 | 0.422 | 0.508 | 0.120 |
| chromVAR deviation | 0.063 | 0.213 | 0.245 | 0.423 | -0.020 |

For the K562-minus-panel contrast, the best estimate reaches Spearman 0.32 with the true contrast
(7 of the top 20 shared). The true contrast itself correlates 0.63 with F-032's transfer-specific
residuals. Every estimator places the GATA family high in K562 and low in the four panel cells
(chromVAR z +29.4 in K562, -4.7 to -14.7 elsewhere), and CEBP high in HepG2 (+10.2) and low in K562
(-7.2). CTCF is the top motif by within-cell ATAC association in K562 and near the bottom of the
truth (rank 741 of 782); the discriminative and chromVAR estimators place it at 311-664, the
cell-specific one only at 15-52.
**Implications:** **The TF activities that drive the transfer errors can be read from the target
cell's ATAC, so telling a model about them is deployable.** Raw accessibility strength is the wrong
estimator (it puts CTCF first). The specificity-based estimators remove CTCF but measure "different
from the reference panel", so a new cell type similar to a training cell would have its important
TFs cancelled; the proposed alternative is to learn the mapping from each motif's ATAC-only
features to its p300 importance across the panel, evaluated leave-one-cell-out. Test 1 (a model
with a TF-activity input channel) is running with oracle weights first.
**Caveats:** Five cell types. Truth and estimates share the region set and the FIMO calls, so
FIMO's limits apply to both. The per-cell top-20 overlap varies from 0-2 (HepG2, MCF-7) to 15-18.
**Tags:** p300, tf-activity, chromvar, atac, motifs, transfer, deployment

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-26 | `4.43` job 45437518; panel ABC regions 2026_0924_p300_panel_obs | K562, GM12878, HepG2, A549, MCF-7 ATAC + p300 bigwigs | 2026_0824_H3K27ac_model | GATA and CEBP recovered; whole ranking -0.09 to 0.57 | TF-activity input is deployable in principle |

## F-034: Learning motif p300 importance from ATAC-only features ranks motifs better overall but loses the cell-specific part; plain ATAC association minus CTCF keeps it best
**Status:** supported
**Claim:** Target: F-033's truth (each motif's coefficient for p300 beyond ATAC, standardised within
cell). Features per (cell, motif), from that cell's ATAC and the sequence only, no motif identity:
ATAC association, top-decile association, promoter fraction, prevalence, CTCF co-occurrence,
footprint depth and flank height from base-level ATAC at genome-wide FIMO hits, motif GC.
Gradient boosting, leave one cell out (`4.47`). Two corrections (`4.48`): STACK = least-squares
mix of the learned score and ATAC association, weights fitted on inner leave-one-out predictions;
DEPARTURE = mean truth over the training cells plus a model trained to predict each cell's
departure from the other training cells' mean. "Departure" columns score (estimate - prior)
against (truth - prior), prior = mean truth over the other four cells: the cell-specific part.

| estimate | mean Spearman vs truth | mean departure Spearman | K562 departure | K562-minus-panel vs true contrast | GATA in K562 (truth +3.88) |
|---|---|---|---|---|---|
| raw ATAC association | 0.427 | **0.413** | **0.767** | **+0.445** | **+3.22** |
| learned (`4.47`) | 0.658 | 0.300 | 0.361 | +0.096 | +1.15 |
| stack: learned + association | 0.654 | 0.371 | 0.518 | +0.153 | +2.05 |
| departure: prior + model | **0.693** | 0.114 | 0.178 | -0.154 | -0.09 |

Per held-out cell, learned vs raw association: K562 0.731 vs 0.573, GM12878 0.734 vs 0.513, HepG2
0.709 vs 0.422, A549 0.637 vs 0.508, MCF-7 0.481 vs 0.120; top-20 overlap is lower for learned in
K562, GM12878 and A549. Learned pushes CTCF to ranks 737-779 without being told to. Permutation
importance: association 0.29, CTCF co-occurrence 0.19, top-decile 0.18, promoter fraction 0.16,
footprint depth 0.13. Stack weights: learned 0.49-0.61, association 0.14-0.26.
**Implications:** **The learned estimators gain Spearman by learning which motifs matter for p300
beyond ATAC in ANY cell** (AP-1 tops the learned K562 list), which is the part the transferred model
already has. The cell-specific part, the part F-032 says is missing, is carried best by the target
cell's own ATAC association, whose only large K562 error is CTCF (rank 1; truth rank 741). With
CTCF removed its K562 top 10 is the GATA family, TRPS1 and GATA1::TAL1. That is the deployable
weighting taken forward (`4.45 --weighting assoc_noctcf`). No estimator works in MCF-7.
**Caveats:** Five cell types, so every leave-one-out model trains on four. The departure model's
failure may be that four cells are too few to learn departures from, not that the features lack
them. FIMO calls shared by truth and estimates.
**Tags:** p300, tf-activity, motifs, atac, estimator, deployment, gradient-boosting

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-28 | `4.47` job 45613770; `4.48` job 45709699 | five p300-panel cell types, ATAC + p300 + FIMO | 2026_0824_H3K27ac_model | learned 0.66 mean Spearman but K562 contrast +0.10; raw association +0.45 | use ATAC association minus CTCF |

## F-035: Telling the transferred model which TFs matter (oracle weights) matches the ranking loss in rE2G, barely moves ABC, and shrinks the GATA error by a third
**Status:** supported (one training run)
**Claim:** The multi-cell p300 model (K562 held out) retrained with one extra input channel, a
TF-activity track: at each base, the sum over genome-wide FIMO hits of that motif's truth weight in
the cell (F-033, from the cell's own p300), rescaled to genome-wide sd 1 (`4.44`-`4.46`, `1.37`).
Predicted in K562 with K562's own oracle track. Everything else identical to the MSE-only
comparator. Paired bootstrap as in F-030.

| motifact_oracle minus | ABC | rE2G |
|---|---|---|
| MSE-only multi-cell | +0.0036 [-0.0044, +0.0126] | **+0.0171** [+0.0061, +0.0288] |
| ranking loss w=10 | -0.0052 [-0.0110, +0.0004] | -0.0002 [-0.0100, +0.0092] |
| ATAC-only floor | **+0.0165** [+0.0002, +0.0316] | **+0.0253** [+0.0060, +0.0413] |
| in-cell model | **-0.0384** [-0.0536, -0.0231] | +0.0203 [-0.0085, +0.0488] |

Validation count Pearson 0.7343 against 0.7129; best epoch 29.4 against 23.6. Oracle tops per
cell are the expected regulators (K562 GATA, GM12878 SPIB/ETS/RUNX, HepG2 CEBP, A549 AP-1, MCF-7
ESR). F-032's residual screen on this model (`4.41 --extra`), non-promoter, percentile points at
matched ATAC:

| motif | in-cell | multi-cell | rank w=10 | + oracle channel |
|---|---|---|---|---|
| GATA family (mean of 8) | -2.65 | +5.65 | +5.63 | **+3.15** |
| GATA1::TAL1 | -3.18 | +5.02 | +4.97 | +3.18 |
| CEBPB | +0.39 | -4.90 | -5.09 | -2.52 |
| ESR2 | +0.25 | -4.47 | -3.96 | -1.05 |

The channel model's transfer-specific pattern correlates 0.87 with the MSE-only model's across
motifs; mean |coefficient| rises slightly (0.77 -> 0.82).
**Implications:** **Even with the true answer the model uses TF information only partly.** It
removes about a third of the GATA under-ranking and half the CEBP/ESR over-ranking, and adds
small errors elsewhere; the ABC gain is a sixth of what directly re-ranking the GATA elements
gave (F-032, +0.022). In rE2G it equals the ranking loss, and the two were combined
(rank10_motifact_oracle) to test whether they add. The upstream metric improved more than any
previous change (+0.021 Pearson) and ABC barely moved: the upstream metric misled a third time.
**The bottleneck is now how the model uses the channel, not whether the information exists.**
Likely reasons: only four training cells, so the channel can act as a cell fingerprint rather than
local evidence; one summed track mixes positive and negative motifs; the weights include
cell-invariant importance (AP-1) the model already has from sequence.
**Caveats:** One training run; the channel model trained ~6 epochs longer, and training length
alone moved ABC by +0.0086 in F-019, so the ABC delta may be partly length. Oracle weights: not
deployable.
**Tags:** p300, tf-activity, input-channel, transfer, oracle, abc, re2g, crispr-benchmark, k562, gata

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-28 | `1.37` oracle; predict 45613275; ABC 45613278; CRISPR 45709851; bootstrap 45731270; rE2G 45616410, bootstrap 45709881; `4.41` 45732088 | K562 CRISPR (ABC and rE2G pair sets), K562 candidate regions | 2026_0824_H3K27ac_model; ENCODE_rE2G 2026_0923_p300_features | rE2G +0.0171, ABC +0.0036 over MSE-only; GATA residual +5.65 -> +3.15 | information helps; model under-uses it |

## F-036: The deployable TF-activity channel (ATAC association, CTCF removed) is the best transferred p300 model so far in both frameworks; adding the ranking loss to the oracle channel does not add
**Status:** REVISED by F-038: seed-averaged, the channel gain is about +0.006 ABC and the assoc_noctcf vs oracle difference is noise
**Claim:** Two follow-ups to F-035, same design and pair sets. `assoc_noctcf`: the channel built
from each cell's own within-cell ATAC association (F-033), CTCF-family motifs set to 0, in
training and at K562 prediction (`4.45 --weighting assoc_noctcf`, `1.37 WEIGHTING=assoc_noctcf`).
`rank10_motifact_oracle`: the oracle channel plus the ranking loss w=10 (`1.36 TAG=motifact_oracle`).

| minus | ABC | rE2G |
|---|---|---|
| assoc_noctcf - MSE-only multi-cell | **+0.0100** [+0.0031, +0.0177] | **+0.0221** [+0.0120, +0.0331] |
| assoc_noctcf - oracle channel | +0.0065 [-0.0002, +0.0131] | +0.0047 [-0.0047, +0.0138] |
| assoc_noctcf - ATAC-only floor | **+0.0230** [+0.0063, +0.0393] | **+0.0300** [+0.0091, +0.0475] |
| rank10 + oracle - MSE-only | +0.0028 [-0.0069, +0.0130] | **+0.0122** [+0.0005, +0.0238] |
| rank10 + oracle - rank10 | -0.0059 [-0.0135, +0.0014] | -0.0054 [-0.0146, +0.0032] |
| rank10 + oracle - oracle channel | -0.0008 [-0.0080, +0.0060] | -0.0052 [-0.0143, +0.0040] |

AUPRC: ABC assoc_noctcf 0.4910 (rank10 0.4897, oracle 0.4845, MSE-only 0.4810, in-cell 0.5229);
rE2G 0.5735 (rank10 0.5689, oracle 0.5688, MSE-only 0.5514). Best epoch / val count Pearson:
MSE-only 23.6 / 0.7129, oracle 29.4 / 0.7344, assoc_noctcf 31.2 / 0.7285, rank10 + oracle
29.2 / 0.7390.
**Implications:** **A TF-activity channel built from the target cell's own ATAC, with no p300
and no reference panel, gives the largest transfer gain yet** (+0.010 ABC, +0.022 rE2G), at
least as large as the oracle channel's. One reading: its training weights come from the same
estimator used at prediction, so there is no train/test mismatch, whereas the oracle weights
are noisy regression coefficients in training and exact in K562. The ranking loss and the
channel do not add (in both frameworks the combination is at or below either alone), so they
likely fix the same errors. Still about three quarters of the in-cell gap remains in ABC.
**Caveats:** One run per model, and F-020's run-to-run noise is about the size of the
assoc_noctcf vs oracle difference; a seed replicate is needed before claiming it beats the
oracle. The channel models trained 6-8 epochs longer than MSE-only (F-019 training-length
confound).
**Tags:** p300, tf-activity, input-channel, deployment, transfer, abc, re2g, crispr-benchmark, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-28 | tracks 45718829; training 45723634 / 45718831; predict 45732090 / 45732092; ABC 45732096; CRISPR 45759984; rE2G 45739352; bootstraps 45793992 / 45793995 | K562 CRISPR (ABC and rE2G pair sets) | 2026_0824_H3K27ac_model; ENCODE_rE2G 2026_0923_p300_features | assoc_noctcf +0.0100 ABC, +0.0221 rE2G over MSE-only | deployable channel works; rank + channel do not add |

## F-037: Teaching the model to use the oracle channel (additive head, cross-cell paired loss) shrinks the GATA error further but does not move the benchmark beyond the channel alone
**Status:** supported (one training run per arm)
**Claim:** Three arms on the oracle channel (`1.38`), otherwise identical to it. HEAD: the channel
bypasses the trunk and enters the counts through 4 weights (positive and negative motif weight,
each also times softplus(ATAC)), zero-initialised (`--motif-head`). PAIR: the channel in the trunk
plus a cross-cell loss on 20,000 pooled training-peak windows extracted in all four training cell
types: predicted minus observed log-count difference between two cell types, centred within each
cell-type pair, weight 10 (`--pair-loss-weight`). HEADPAIR: both.

| arm | ABC vs oracle channel | ABC vs MSE-only | rE2G vs oracle channel | rE2G vs MSE-only |
|---|---|---|---|---|
| head | +0.0019 [-0.0074, +0.0105] | **+0.0055** [+0.0001, +0.0108] | +0.0003 [-0.0117, +0.0125] | **+0.0176** [+0.0077, +0.0282] |
| pair | +0.0032 [-0.0093, +0.0152] | +0.0068 [-0.0048, +0.0182] | +0.0037 [-0.0102, +0.0172] | **+0.0210** [+0.0082, +0.0337] |
| headpair | +0.0069 [-0.0044, +0.0178] | **+0.0105** [+0.0009, +0.0196] | +0.0053 [-0.0070, +0.0180] | **+0.0226** [+0.0110, +0.0358] |

Against the deployable `assoc_noctcf` channel (F-036), every arm is within +/-0.005 in both
frameworks (headpair +0.0004 ABC, +0.0005 rE2G). AUPRC: ABC headpair 0.4914, assoc_noctcf 0.4910,
pair 0.4877, head 0.4864, oracle 0.4845, MSE-only 0.4810, in-cell 0.5229.

F-032's residual screen (`4.41 --extra`), non-promoter, percentile points at matched ATAC:

| | in-cell | MSE-only | oracle | assoc_noctcf | head | pair | headpair |
|---|---|---|---|---|---|---|---|
| GATA family | -2.64 | +5.69 | +3.15 | +4.77 | +5.39 | **+2.14** | +4.88 |
| CEBP family | -0.09 | -3.61 | -2.15 | -2.14 | -2.67 | -2.74 | -2.83 |
| ESR1/2 | +0.28 | -3.89 | -0.10 | -1.99 | -3.08 | -0.69 | -3.19 |
| mean abs coefficient | 0.60 | 0.77 | 0.82 | 0.80 | 0.69 | **0.60** | 0.67 |

Learned head weights (final): positive motifs +0.09 to +0.30, negative motifs -0.49 to -1.24,
negative x ATAC +0.25 to +0.91. Validation count Pearson: head 0.712, pair 0.726, headpair 0.711
(oracle channel 0.734, MSE-only 0.713).
**Implications:** **The paired loss does what it was designed to do, and the benchmark does not
follow.** It gives the smallest GATA error of any transferred model and an overall residual size
equal to the in-cell model's, yet its ABC gain over the channel alone is +0.003 (n.s.). The head
alone barely changes the GATA error, so a fixed additive use of TF activity is not enough: the
channel has to interact with sequence and ATAC in the trunk. Every TF-channel variant, oracle or
deployable, now sits in one band (ABC 0.485-0.491, rE2G 0.569-0.574) about +0.005-0.010 ABC and
+0.017-0.023 rE2G above MSE-only, with differences among them inside one-run noise. The GATA
residual and the benchmark have come apart: residuals are over all 153k candidates, the benchmark
over 3,406 tested elements, the same gap that sank the partial-correlation screen (F-030, F-031).
**Caveats:** One run per arm; the ~0.005 spread among channel models is F-020-sized. Oracle
weights only; the paired loss on the deployable channel is untested.
**Tags:** p300, tf-activity, input-channel, paired-loss, motif-head, transfer, abc, re2g, k562, gata

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-29 | `1.38` 45793888; predict 45826172/176/242; ABC 45826245; CRISPR 45842255; rE2G 45829844; bootstraps 45856314/5; `4.41` 45826253 | K562 CRISPR, K562 candidate regions | 2026_0824_H3K27ac_model; ENCODE_rE2G 2026_0923_p300_features | pair: GATA +5.69 -> +2.14; ABC/rE2G within noise of the channel alone | channel models have plateaued at one run's resolution |

## F-038: Seed replicates put run-to-run noise at sd ~0.002 AUPRC in ABC and ~0.005 in rE2G; the TF-activity channel's real gain is about +0.006 ABC and +0.012-0.014 rE2G, and the frozen two-stage model matches it
**Status:** established (three runs per model)
**Claim:** Two extra trainings (`1.37 SEED=1,2`, `--seed`: torch, numpy and batch order; same data)
of the MSE-only multi-cell model, the oracle channel and the assoc_noctcf channel, plus the frozen
two-stage models (`1.39`, DNA-free 5,409-parameter correction on the original MSE-only model).
One ABC pass and one rE2G pass; paired bootstrap with seed-averaged groups (`4.17`/`4.34 --group`:
mean AUPRC of the three runs inside every draw).

| | ABC runs (orig, s1, s2) | ABC mean, sd | rE2G runs | rE2G mean, sd |
|---|---|---|---|---|
| MSE-only | 0.4810, 0.4832, 0.4842 | 0.4828, 0.0017 | 0.5516, 0.5577, 0.5643 | 0.5579, 0.0063 |
| oracle channel | 0.4845, 0.4853, 0.4885 | 0.4861, 0.0021 | 0.5687, 0.5747, 0.5722 | 0.5719, 0.0030 |
| assoc_noctcf channel | 0.4910, 0.4880, 0.4872 | 0.4887, 0.0020 | 0.5735, 0.5706, 0.5642 | 0.5695, 0.0047 |

| seed-averaged | ABC | rE2G |
|---|---|---|
| oracle - MSE-only | +0.0033 [-0.0029, +0.0112] | **+0.0140** [+0.0056, +0.0232] |
| assoc_noctcf - MSE-only | **+0.0059** [+0.0001, +0.0132] | **+0.0116** [+0.0038, +0.0202] |
| assoc_noctcf - oracle | +0.0026 [-0.0012, +0.0063] | -0.0024 [-0.0073, +0.0027] |
| MSE-only - ATAC floor | +0.0148 [-0.0021, +0.0303] | +0.0144 [-0.0036, +0.0300] |
| assoc_noctcf - ATAC floor | **+0.0207** [+0.0043, +0.0362] | **+0.0260** [+0.0068, +0.0414] |
| assoc_noctcf - in-cell | **-0.0341** [-0.0484, -0.0196] | +0.0210 [-0.0078, +0.0497] |
| two-stage oracle - its frozen base | **+0.0042** [+0.0003, +0.0086] | **+0.0146** [+0.0062, +0.0238] |
| two-stage assoc_noctcf - its frozen base | +0.0027 [-0.0017, +0.0069] | **+0.0138** [+0.0068, +0.0222] |

Same-configuration replicate differences reach +0.0032 (ABC) and +0.0126 (rE2G, resolvable on its
own: seed 2 minus original MSE-only).
**Implications:** **Single-run differences under about 0.005 ABC or 0.013 rE2G are not evidence.**
The TF-activity channel's gain is real in rE2G in both weightings and marginal in ABC (+0.006,
about a seventh of the 0.041 in-cell gap); F-036's "+0.010, beats the oracle" was the original
assoc_noctcf run being the luckiest of three. The deployable and oracle weightings are
indistinguishable. The frozen two-stage model, whose channel correction cannot see DNA, gains as
much as the channel trained end to end, so the DNA-shortcut hypothesis does not explain the
plateau. Every channel variant (F-035 to F-037) now sits in one band, and earlier single-run rank
results (F-030, +0.0087 ABC) should be read with the same noise in mind. Deciding between
candidates of this size needs at least three runs each.
**Caveats:** Three runs per model; sd from three values is itself imprecise. The two-stage models
were compared with their one frozen base, not seed-averaged.
**Tags:** p300, seeds, run-to-run-variance, tf-activity, two-stage, abc, re2g, crispr-benchmark, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-30 | training 45897299/321/357/363/369/503 (+45957358/360 after a GPU fault), 45903507/509; ABC 45897566; CRISPR 45998015; rE2G 45980539; bootstraps 46073940/42 | K562 CRISPR (ABC and rE2G pair sets) | 2026_0824_H3K27ac_model; ENCODE_rE2G 2026_0923_p300_features | channel +0.006 ABC, +0.012-0.014 rE2G seed-averaged; run sd 0.002 / 0.005 | resolves F-036; channel gain is small and real |

## F-039: ATAC fragment features predict H3K27ac beyond ATAC three times as well as p300, with the opposite sign, mostly because H3K27ac-high elements sit in broad accessible domains; the weights carry between K562 and GM12878
**Status:** supported (two cell types, one library each)
**Claim:** 4.51's per-element fragment features (non-promoter candidates, >= 10 fragments;
121,219 in K562, 124,613 in GM12878) scored against three targets after removing ATAC rank,
rank^2, GC, CpG o/e and log TSS distance (`4.55`): p300 over the candidate region, H3K27ac over
the region, and H3K27ac within +/- 1 kb of its centre (5' ChIP-seq bigwigs from the model inputs).

| partial Spearman (K562 / GM12878) | p300 | H3K27ac region | H3K27ac +/- 1 kb |
|---|---|---|---|
| di_flank (300-450 bp, 150-600 bp out) | +0.03 / -0.06 | **+0.47 / +0.42** | +0.43 / +0.38 |
| mono_flank (150-250 bp, 100-450 bp out) | +0.04 / -0.03 | **+0.42 / +0.33** | +0.31 / +0.20 |
| flank_ratio | -0.12 / -0.11 | +0.31 / +0.34 | +0.19 / +0.20 |
| short_frac | +0.13 / +0.11 | -0.27 / -0.34 | -0.16 / -0.21 |
| nfr (< 120 bp at centre) | +0.19 / +0.13 | -0.16 / -0.27 | -0.06 / -0.15 |
| positioning | -0.01 / -0.01 | -0.22 / -0.21 | -0.22 / -0.20 |
| joint score, weights from the other cell | 0.224 / 0.162 | **0.521 / 0.476** | 0.454 / 0.397 |
| joint score, in-cell weights | 0.248 / 0.179 | 0.533 / 0.484 | 0.468 / 0.407 |

Sign agreement between the cells: 7 of 9 features for p300, 8 of 9 for H3K27ac. The ATAC-only
Spearman is lower for H3K27ac (0.43 / 0.48 region) than for p300 (0.63 / 0.70). Residual p300
and residual H3K27ac (+/- 1 kb) correlate +0.31 in K562 and +0.13 in GM12878.
Average fragment profiles (`4.56`, fig35; 3,000 regions per group, top vs bottom residual
quintile, candidate ATAC matched): the p300-high group has a taller, sharper central peak of
short fragments (K562 249 vs 171 short fragments per region). The H3K27ac-high group has NO
central peak: short and mono-nucleosome fragments are spread over +/- 1 kb, with short-fragment
maxima near +/- 450 bp, while the H3K27ac-low group is one narrow isolated peak at the centre.
Same in GM12878.
**Implications:** **The fragment features mostly measure how accessibility is spread around
the element, not nucleosome packing as such.** At the same ATAC over the candidate, H3K27ac is
high where the candidate sits inside a broad accessible domain (neighbouring open sites and
the nucleosomes between them) and low at a narrow isolated site; p300 is high where the
central site itself is sharper and more open. The features are a large, cell-portable
predictor of H3K27ac and a small one for p300. The joint score's weights, fitted in one cell, lose only 0.01 in the other,
so at the element level, ranked within each cell, the signal is not library-specific between
these two libraries. That contrasts with F-005, where per-base fragment-size input channels
gained in-cell for H3K27ac and were null transferred; the difference may be the within-cell
rank normalisation. For the p300 activity model the fragment channels (1.40) have less to add.
For an H3K27ac target they could add a lot, if given as within-cell normalised features.
**Caveats:** Partly mechanical: H3K27ac is a broad mark counted over the candidate, so a
candidate inside a broad domain collects signal from its neighbours; ATAC within +/- 1 kb may
explain much of it without fragment lengths (4.57 `tx_breadth` arm tests this). Only K562 and GM12878 have H3K27ac here, so it
cannot be checked in MCF-7, where 4.51's p300 signs flipped. One ATAC library set per cell.
Element-level association, not a benchmark result.
**Tags:** h3k27ac, p300, nucleosome, fragment-length, atac, transfer, k562, gm12878

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-30 | 4.51 features (nucsig array); 4.55 job 46100582 | K562, GM12878 ATAC PE BAMs; H3K27ac and p300 5' bigwigs; union ABC candidates | 2026_0824_H3K27ac_model | joint partial rho H3K27ac 0.48-0.52 vs p300 0.16-0.22, weights transfer | nucleosome signal real for H3K27ac, small for p300 |

## F-040: Nucleosome-length ATAC fragments are not longer next to H3K27ac-high or p300-high elements, even at the nucleosome adjacent to the peak; flanks of weak inactive peaks are ~2 bp longer
**Status:** established (K562, three libraries, two window sizes)
**Claim:** Hypothesis (Maya): flanking nucleosomes of acetylated elements give one-nucleosome
ATAC fragments ~10-15 bp longer. K562 ATAC MACS2 peaks (q < 0.01, > 1 kb from a TSS, no other
confident peak within 2 kb; 91,528) grouped by ATAC, H3K27ac (+/- 1 kb) or p300, and CTCF
(motif within 100 bp + ChIP ENCFF216XRV); fragments counted only if they lie entirely in a flank
(MACS2 bound + 50 bp out to 1 or 2 kb), from the three ENCSR868FGK libraries (`4.59`).

| mean length of 140-260 bp fragments, 95% CI | whole 1 kb flank | nearest flank (50-150 bp out), 2 kb run |
|---|---|---|
| weak ATAC, H3K27ac low | 205.7 [205.5, 205.9] | 207.6 [207.1, 208.1] |
| strong ATAC, H3K27ac low | 203.9 [203.6, 204.2] | 205.5 [204.5, 206.5] |
| strong ATAC, H3K27ac high | 203.7 [203.6, 203.9] | 204.2 [203.7, 204.7] |
| any ATAC, H3K27ac high | 204.0 [203.8, 204.1] | 204.7 [204.3, 205.1] |
| CTCF-bound, H3K27ac low | 203.7 [203.5, 203.9] | 203.2 [202.6, 203.6] |

The mono-nucleosome mode is 206-209 bp in every group; p300 groups give the same pattern. The
two-nucleosome peak (392-409 bp) is too broad to resolve the repeat length (intervals ~10 bp).
Flank GC differs by 0.01-0.03 between groups. The three libraries superimpose in every group
(fig37). In the standard view (all fragments within +/- 1 kb, peak included) the groups differ
strongly (CTCF: tallest sub-nucleosomal peak; weak peaks: most nucleosomal); with the peak
excluded they nearly converge.
**Implications:** **No support for longer nucleosomal fragments at active elements.** Matched on
ATAC, H3K27ac-high flanks are 0.2 bp (whole flank) to 1.3 bp (nearest nucleosome) SHORTER than
H3K27ac-low flanks, against a hypothesised +10-15 bp that intervals of +/- 0.5 bp would have
resolved. The window was not too wide: the effect is absent at the adjacent nucleosome too. The
category differences in standard fragment-length plots come from the peak itself (how open it
is and what occupies it). Two small real effects: flanks of weak inactive peaks are ~2 bp longer
at every distance, and CTCF-adjacent nucleosomes give the shortest fragments (203 bp), rising
with distance, consistent with tightly phased nucleosomes next to CTCF.
**Caveats:** K562 only. The outermost distance bin of each window is biased short (a fragment
must end inside the window). ATAC measures Tn5 access to linkers, not nucleosome footprint
directly; MNase-seq would test footprint length.
**Tags:** h3k27ac, p300, nucleosome, fragment-length, atac, ctcf, k562, negative-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-10-01 | 4.59 jobs 46141966/67, 46143202/03, 46144822/23 | K562 ATAC PE BAMs (ENCFF077FBI, ENCFF128WZG, ENCFF534DCE); MACS2 peaks; H3K27ac, p300, CTCF ChIP | 2026_0824_H3K27ac_model | H3K27ac-high flanks 0.2-1.3 bp shorter, not 10-15 bp longer | hypothesis not supported |

## F-041: TF lists from RNA do not beat the ATAC-derived TF list; RNA + ATAC ties it and RNA alone is worse
**Status:** established (three runs per model)
**Claim:** TF-activity channel weighted by RNA (4.54 -> 4.45; ENCODE Gingeras poly(A) RNA-seq, family
maximum of log TPM) alone (`rnafam`) or z-summed with ATAC association (`rnaatac`), three training
runs each (`1.37 SEED`), seed-averaged paired bootstrap against the F-038 groups.

| seed-averaged | ABC | rE2G |
|---|---|---|
| rnaatac - MSE-only | +0.0058 [-0.0001, +0.0134] | **+0.0091** [+0.0020, +0.0178] |
| rnaatac - assoc_noctcf (ATAC list) | -0.0001 [-0.0037, +0.0034] | -0.0025 [-0.0093, +0.0045] |
| rnafam - MSE-only | +0.0014 [-0.0038, +0.0089] | +0.0027 [-0.0040, +0.0108] |
| rnafam - assoc_noctcf | **-0.0045** [-0.0085, -0.0005] | **-0.0089** [-0.0136, -0.0039] |

Means: ABC rnafam 0.4842, rnaatac 0.4886, assoc_noctcf 0.4887, MSE-only 0.4828; rE2G rnafam
0.5606, rnaatac 0.5670, assoc_noctcf 0.5694.
**Implications:** **In bulk ENCODE data, expression adds nothing to the ATAC-derived TF list,
and expression alone is worse.** The motif-level proxy ranking (RNA 0.28, RNA + ATAC 0.41, ATAC
0.44) predicted the order correctly this time. The DeepSCENIC-style per-TF factorised model
(option 2) and the MotifCompendium remapping were conditional on RNA helping; deprioritise both.
**Caveats:** Motif -> TF mapping by JASPAR name plus co-occurrence families, not MotifCompendium;
bulk RNA from different ENCODE experiments than the ATAC.
**Tags:** p300, tf-activity, rna, deepscenic, abc, re2g, crispr-benchmark, k562, negative-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-10-01 | training 46085027-39; ABC 46085041; CRISPR 46162345; rE2G 46157330; bootstraps 46200986, 46201552 | K562 CRISPR (ABC and rE2G pair sets) | 2026_0824_H3K27ac_model | rnaatac = assoc_noctcf; rnafam -0.005/-0.009 | RNA does not help |

## F-042: Fragment-size ATAC input channels do not help the transferred p300 model (three runs)
**Status:** established (three runs per model)
**Claim:** Multi-cell p300 model with 5' ATAC insertions split into four fragment-size channels
(sub <= 139, mono 140-329, di 330-620, poly >= 621 bp) plus total (`frag5p`), against the same
model with total 5' insertions only (`flat5p`), three runs each (`1.40`).

| seed-averaged | ABC | rE2G |
|---|---|---|
| frag5p - flat5p | -0.0046 [-0.0127, +0.0033] | -0.0098 [-0.0197, +0.0007] |
| flat5p - MSE-only (coverage input) | +0.0041 [-0.0028, +0.0104] | +0.0027 [-0.0085, +0.0127] |

Means: ABC flat5p 0.4869 (sd 0.0042), frag5p 0.4823; rE2G flat5p 0.5606 (sd 0.0151), frag5p 0.5508.
**Implications:** **Fragment-size channels are null to slightly negative on transfer**, as F-005
found for H3K27ac. Together with F-039/F-040, fragment-length structure does not carry
deployable p300 signal. The 5' insertion input is at least as good as coverage. Close this line.
**Caveats:** flat5p's rE2G runs vary widely (0.547-0.577).
**Tags:** p300, fragment-length, input-channel, transfer, abc, re2g, crispr-benchmark, k562, negative-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-10-01 | channels 46088402/07/09; training 46088410-36; ABC 46150705; CRISPR 46171084; rE2G 46157330; bootstraps 46200988, 46201552 | K562 CRISPR | 2026_0824_H3K27ac_model | frag5p - flat5p -0.005 ABC, -0.010 rE2G | null |

## F-043: In ABC the ATAC term costs observed p300 0.025 and observed H3K27ac 0.017, and p300 x H3K27ac adds nothing over p300; in rE2G the combination leads (+0.012, unresolved); the deployable combination and the fragment correction of transferred H3K27ac do not help
**Status:** supported (single runs)
**Claim:** All tracks painted over the July K562 candidates (4.57, 4.58). ABC (`k27frag` run,
10,342 pairs) and rE2G (10,356 pairs).

| observed | ABC with ATAC | ABC without ATAC | rE2G |
|---|---|---|---|
| H3K27ac | 0.5133 | 0.5305 | 0.6196 |
| p300 | 0.5657 | **0.5909** | 0.6602 |
| p300 x H3K27ac (sqrt of qnormed product) | 0.5519 | **0.5911** | **0.6718** |

ABC: p300 x H3K27ac with ATAC - p300 -0.0138 [-0.0230, -0.0042]; without ATAC - with ATAC
+0.0392 [+0.0235, +0.0552]. No-ATAC controls: p300 alone - with ATAC +0.0252 [+0.0074, +0.0425];
H3K27ac alone - with ATAC +0.0172 [+0.0022, +0.0315]; p300 x H3K27ac without ATAC - p300 without
ATAC +0.0001 [-0.0182, +0.0183]. rE2G: p300 x H3K27ac - p300 +0.0116 [-0.0031, +0.0262];
- H3K27ac +0.0522 [+0.0312, +0.0738]. Painting costs 0.016 ABC / 0.007 rE2G for H3K27ac
against the read-count arms, so compare within the painted set.

Predicted, deployable: predicted p300 (MSE-only multi-cell) x corrected transferred H3K27ac -
predicted p300 alone -0.0104 [-0.0208, +0.0008] ABC, -0.0071 [-0.0224, +0.0085] rE2G; with the
assoc_noctcf p300 -0.0151 [-0.0247, -0.0053] ABC.

Transferred H3K27ac correction (fitted in GM12878): ABC recalibration - repaint -0.0002,
+ fragment and +/- 1 kb ATAC features -0.0055 to -0.0072 (worse); rE2G recalibration + 0.0207
[+0.0077, +0.0367] over repaint, fragments then -0.0074; the best corrected track (0.5411) only
reaches rE2G's ATAC-only floor (0.5435), and the uncorrected one is below it (-0.0230).
**Implications:** **In ABC the whole no-ATAC gain is from removing ATAC, not from combining
marks**: observed p300 alone without ATAC (0.591) equals the combination, and the ATAC term in
ABC's geometric mean costs both observed marks (+0.025 p300, +0.017 H3K27ac when dropped). In
rE2G, which learns its own weights, the combination is the best observed activity term (+0.012
over p300, not resolved). Dropping ATAC is a free change for the predicted arms too and is being
tested (F-044). The deployable
version does not work: predicted H3K27ac transferred from GM12878 carries no information beyond
ATAC (it is at or below the ATAC floor), so multiplying it into predicted p300 only adds noise.
The element-level gains of the correction (Spearman 0.732 -> 0.787, F-039) did not reach the
benchmark. A two-target model is worth building only if its H3K27ac head transfers far better
than the GM12878 H3K27ac model does.
**Caveats:** Single runs; the observed arms are ceilings, not deployable. p300 and H3K27ac weighted
equally in the product; other weightings untested.
**Tags:** p300, h3k27ac, activity-combination, abc, re2g, crispr-benchmark, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-10-01 | 4.57 46105628; 4.58 46130383; ABC 46130384; CRISPR 46201297; rE2G 46157330; bootstraps 46200994, 46201552 | K562 CRISPR (ABC and rE2G pair sets) | 2026_0824_H3K27ac_model; ENCODE_rE2G 2026_0923_p300_features | obs p300 x H3K27ac no-ATAC ABC 0.591; rE2G 0.672; predicted combos below p300 alone | observed complementary; deployable not |

## F-044: Predicted p300 needs ABC's ATAC term (dropping it costs 0.021, the opposite of observed p300); observed H3K4me1 is no substitute for H3K27ac in ABC or rE2G
**Status:** supported (predicted arms 3 seeds; observed arms single runs)
**Claim:** ABC (`noatac` run, 10,342 pairs, 466 regulated), seed-averaged over 3 runs for the
predicted arms; rE2G on the shared July candidate set (10,356 pairs, 471 regulated).

| activity term (ABC) | with ATAC | without ATAC | without - with |
|---|---|---|---|
| predicted p300, MSE multi-cell (mean of 3) | 0.4828 | 0.4613 | -0.0215 [-0.0383, -0.0047] |
| predicted p300, `assoc_noctcf` (mean of 3) | 0.4887 | 0.4750 | -0.0137 [-0.0311, +0.0024] |
| observed p300 | 0.5673 | 0.5895 | +0.0222 [+0.0040, +0.0393] |
| observed H3K27ac | 0.5167 | 0.5383 | +0.0216 [+0.0057, +0.0372] |
| observed H3K4me1 | 0.4572 | 0.4053 | -0.0519 [-0.0799, -0.0229] |

ATAC only 0.4680. Predicted p300 without ATAC - ATAC only -0.0067 [-0.0336, +0.0202]. Observed
H3K4me1 with ATAC - ATAC only -0.0108 [-0.0382, +0.0168]; - observed H3K27ac -0.0595
[-0.0873, -0.0317].

rE2G, H3K4me1 in the H3K27ac slot (all 14 ATAC + H3K27ac features computed from H3K4me1):
0.5681; - ATAC only (0.5434) +0.0247 [-0.0053, +0.0498]; - H3K27ac (0.6266) -0.0585
[-0.0953, -0.0239]; - observed p300 (0.6660) -0.0979 [-0.1307, -0.0640].

rE2G, H3K4me1 as two EXTRA features beside the 10 ATAC-only features (`ATAC_H3K4me1sep_powerlaw`,
4.61; enhancer and promoter H3K4me1, nothing else changed): 0.5638; - ATAC only +0.0204
[+0.0061, +0.0310] (sign kept 99.6%); - H3K4me1 swap -0.0042 [-0.0248, +0.0173]; - H3K27ac
-0.0628 [-0.0919, -0.0331]; - observed p300 -0.1021 [-0.1307, -0.0730]. Full-model
coefficients: H3K4me1 at E +0.28, at P -0.24 (ATAC at E +0.30, at P -1.39).
**Implications:** Removing ATAC from ABC's activity term helps only an activity track that
carries information ATAC does not. Observed p300 and H3K27ac gain about 0.022 each (F-043
replicated in a fresh run); predicted p300 loses 0.014 to 0.021, and without ATAC it falls to
the ATAC-only floor. So the predicted track's ranking is mostly ATAC re-encoded (as 4.35
measured), and ABC's geometric mean with real ATAC is what keeps it above the floor. Keep the
ATAC term for every deployable arm; no no-ATAC re-scoring.
Observed K562 H3K4me1 is a poor activity term: in ABC it adds nothing over ATAC and alone is
0.063 below it; in rE2G, where it is weighted freely, it adds a resolved +0.020 over ATAC as two
separate features, the same as swapping it in for H3K27ac, and stays 0.063 below H3K27ac. It
cannot replace H3K27ac, and the most an H3K4me1 feature is worth on this benchmark is ~0.02. This sets a high bar for training H3K4me1 models (the
observed mark is the ceiling for any prediction of it).
**Caveats:** Observed arms single runs. rE2G has no model with H3K27ac AND H3K4me1, so whether
H3K4me1 adds anything once H3K27ac is present is untested (irrelevant to deployment, where
neither is measured). rE2G coefficients are on its own feature scaling; compare signs, not sizes.
**Tags:** p300, h3k4me1, h3k27ac, activity-term, abc, re2g, crispr-benchmark, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-10-01 | ABC 46228377; CRISPR 46228383; bootstrap 46228386; rE2G 46228392; bootstrap 46228394 (`log/noatac_launch.sh`); separate features 4.61 46240698, rE2G 46246447, bootstrap 46246449 (`log/k4sep_launch.sh`) | K562 CRISPR (ABC and rE2G pair sets) | 2026_0824_H3K27ac_model; ENCODE_rE2G 2026_0923_p300_features | predicted p300 no-ATAC -0.021 ABC; H3K4me1 ABC 0.457 vs ATAC 0.468; rE2G swap 0.568, separate 0.564 vs ATAC 0.543, H3K27ac 0.627 | ATAC term needed for predicted; H3K4me1 weak |

## F-045: Most of H3K4me1 is not ATAC, the remainder is reproducible, shared across cells and tracks p300 more than H3K27ac does; other cells' H3K4me1 predicts a held-out cell's p300 nearly as well as other cells' p300
**Status:** supported (one analysis, five cells)
**Claim:** `4.60` over 4.43's union candidate set, each cell's own non-promoter regions
(122,645 to 126,762 per cell). All partial Spearman, controls = ATAC rank + rank^2, GC, CpG o/e,
log distance to TSS. H3K4me1 from `0.47` (stranded 5', +/- 250 bp and +/- 1 kb).

| cell | rho with ATAC | R2 on controls | replicate rho of residual | partial rho with p300 | H3K27ac partial rho with p300 |
|---|---|---|---|---|---|
| K562 | 0.447 | 0.333 | 0.700 | 0.542 | 0.376 |
| GM12878 | 0.483 | 0.266 | 0.809 | 0.272 | 0.156 |
| HepG2 | 0.482 | 0.380 | 0.920 | 0.244 | - |
| A549 | 0.456 | 0.415 | 0.868 | 0.352 | - |
| MCF-7 | 0.453 | 0.300 | 0.856 | 0.424 | - |

(+/- 250 bp window; +/- 1 kb is within 0.05 everywhere, with higher replicate rho, 0.80-0.97.)

Across cells, over regions in both: H3K4me1 residuals correlate 0.43-0.65 (500 bp; 0.47-0.68 at
1 kb), against 0.17-0.35 for p300 residuals and 0.37-0.56 for raw ATAC rank. Leave-one-cell-out:
the mean H3K4me1 residual of the other four cells has partial rho with the held-out cell's p300
residual 0.160-0.271 (500 bp), against 0.243-0.370 for the other four cells' p300; equal in A549
(0.271 vs 0.264) and MCF-7 (0.249 vs 0.248).
**Implications:** H3K4me1 passes the pre-training gate: 58-73% of its rank variance is beyond
ATAC and composition, that part is reproducible (0.70-0.97 between replicates), and it tracks p300
beyond ATAC more strongly than H3K27ac does in both cells that have H3K27ac. It is also the most
cell-shared of the signals measured, which is what a model trained on four cells needs to
transfer. But the shared part is, by construction, not what makes regulation cell-specific, and
the benchmark already says observed K562 H3K4me1 adds nothing over ATAC in ABC and cannot
replace H3K27ac in rE2G (F-044). So a well-transferring H3K4me1 model would predict something
with a low benchmark ceiling. The separate-feature rE2G arm (F-044) then
showed observed H3K4me1 adds +0.020 [+0.006, +0.031] over ATAC: real but small, and the ceiling
for any predicted H3K4me1.
**Caveats:** No input-control normalisation, and composition controls cover GC and CpG only:
part of the cross-cell agreement could be shared ChIP background (mappability, copy number in
cancer lines) rather than biology. Replicate agreement cannot separate the two. Promoters
excluded throughout. Partial correlations over all candidates, not over the CRISPR-tested
elements.
**Tags:** h3k4me1, p300, h3k27ac, beyond-atac, cross-cell, transfer, panel

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-10-01 | 0.47 46238895; 4.60 46239127 | ENCODE H3K4me1, ATAC, p300 in K562, GM12878, HepG2, A549, MCF-7 | 2026_0824_H3K27ac_model | H3K4me1 residual rep rho 0.70-0.97, p300 partial 0.24-0.54, cross-cell 0.43-0.65, LOO 0.16-0.27 | signal beyond ATAC, shared across cells |

## F-046: Adding three ENCODE tissues (stomach, transverse colon, ovary) to multi-cell p300 training does not help K562 (+0.0015, three runs)
**Status:** supported (3 runs per arm)
**Claim:** 1.41: the hold-out-K562 MSE-only multi-cell model, identical to 1.37 WEIGHTING=none
except for the panel, now seven training cell types (four lines + three tissues; 46,828 training
peaks). Tissue EP300 is 76 bp PE, built read 1 only (0.45); tissue ATAC 91-106M filtered reads,
EP300 enrichment at own peaks 11.7-12.7x on chr8, level with A549 / HepG2 (12.2-12.4x). ABC
(`tissue` run, 10,342 pairs, 466 regulated), activity = geomean with observed ATAC (F-044).

| arm (mean of 3 runs) | AUPRC | runs |
|---|---|---|
| multi-cell, MSE only, 4 lines (`mse`) | 0.4828 | 0.4810, 0.4832, 0.4842 |
| + 3 tissues (`tissue`) | 0.4843 | 0.4856, 0.4793, 0.4881 |
| assoc_noctcf TF channel (`ctcf`) | 0.4887 | 0.4910, 0.4880, 0.4872 |
| ATAC only | 0.4680 | |

tissue - mse +0.0015 [-0.0077, +0.0104] (sign kept 64%); tissue - ctcf -0.0044 [-0.0162,
+0.0061]; tissue - ATAC only +0.0163 [-0.0034, +0.0348]; mse - ATAC only +0.0148 [-0.0021, +0.0303].
**Implications:** More and more TF-diverse training cell types do not move K562, including the one
tissue with GATA activity (stomach, 4.49/4.52). The run-to-run spread with tissues (sd 0.0045) is
larger than the effect. With F-034 to F-038 and F-041, panel composition joins the TF-list and
channel variants that plateau within run noise near 0.483-0.489, about 0.015-0.02 above ATAC
alone. Stop adding cell types to this panel.
**Caveats:** One biological replicate per tissue and fewer peaks (3.4-5.3k training peaks each);
tissues are mixed cell populations, so their GATA signal may come from a minor cell type. Only
the MSE-only model was trained on the tissue panel; the TF-channel model was not. ABC only, no
rE2G.
**Tags:** p300, multi-cell, panel, tissues, transfer, abc, crispr-benchmark, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-10-02 | 0.45 46270676; 1.41 46270680-82; preds 46286321/22/41; ABC 46286343; CRISPR 46286345; bootstrap 46286348 (`log/tissue_launch.sh`, `log/tissue_score_launch.sh`) | ENCODE EP300 + ATAC stomach, transverse colon, ovary; K562 CRISPR | 2026_0824_H3K27ac_model | tissue 0.4843 vs mse 0.4828, +0.0015 [-0.0077, +0.0104] | null |

## F-047: Simple bias-corrected footprints from Xu K562 multiome pseudobulk detect CTCF binding but not GATA1, and their motif ranking is dominated by GC-rich motifs
**Status:** supported (one dataset, one scoring method)
**Claim:** `2026_1002_sce2g_tf_activity/scripts/1.1` (TOBIAS/HINT-style: core vs 30 bp flanks,
canonical 8-mer Tn5 bias from pooled insertions in 200 bp bins with >= 20 insertions) over all
4.44 FIMO hits with >= 10 insertions (172,082,584 hits; Xu pseudobulk, chr1 alone 30.6M fragments).
`1.2` validation, AUROC of the footprint score for ChIP-bound vs unbound hits WITHIN deciles of
the hit's insertion count:

| motif | ChIP | hits | bound | AUROC, insertions | AUROC, FS within decile | top decile |
|---|---|---|---|---|---|---|
| CTCF (MA0139.2) | CTCF ENCFF519CXF | 313,773 | 40,537 | 0.892 | 0.694 | 0.754 |
| GATA1 (MA0035.5) | GATA1 (K562_GATA1_BPNet peaks) | 67,741 | 270 | 0.982 | 0.419 | 0.413 |
| GATA1::TAL1 (MA0140.3) | GATA1 | 184,480 | 2,157 | 0.986 | 0.497 | 0.495 |

Mean FS per motif (hits with >= 50 insertions, 782 motifs): top 20 are GC/CpG-rich (TCFL5,
ZBTB33, KLF15, TFDP1, ETS, E2F8, NRF1); CTCF ranks 166, GATA1 765, GATA1::TAL1 533, and the
liver factor HNF1A (not expressed in K562) ranks 94.
**Implications:** The score is a real footprint for a long-residence factor (CTCF) but carries
no binding information for GATA, the factor family behind about half the transfer gap (F-032),
and the motif ranking shows residual sequence bias. In this form it cannot supply the missing
element-level GATA evidence; it was not gated against p300 or put into scE2G.
**Caveats:** Crude bias model (8-mer, estimated on accessible DNA rather than naked DNA, so
accessibility-linked composition leaks in). GATA1 ChIP set is small (4,603 peaks). A learned
Tn5 bias model and multiscale footprints (PRINT / scPrinter) are untested; GATA footprints are
weak in ATAC generally, so a better model may still fail for GATA.
**Tags:** footprinting, atac, multiome, tn5-bias, gata, ctcf, k562, sce2g

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-10-02 | 1.1 bias 46346748, score 46346750; 1.2 46349267 | Xu K562 multiome fragments; K562 CTCF and GATA1 ChIP | 2026_1002_sce2g_tf_activity | CTCF AUROC 0.69 within accessibility deciles; GATA1 0.42-0.50; GC-rich motifs top the ranking | footprints fail for GATA |

## F-048: TOBIAS footprints on Xu K562 pseudobulk detect GATA1 binding beyond accessibility (AUROC 0.67), unlike the homemade score, but its per-motif ranking keeps the GC-rich bias
**Status:** supported (one dataset, one method)
**Claim:** TOBIAS 0.17.5 (`2026_1002_sce2g_tf_activity/scripts/2.2`: ATACorrect `--read_shift 0 0`
on the 2.1 insertion BAM -> ScoreBigwig -> BINDetect, 879 JASPAR 2024 motifs) over scE2G's
155,134 candidate regions. `2.3` validation, same test as F-047: AUROC of the TOBIAS footprint
score for ChIP-bound vs unbound sites WITHIN deciles of accessibility (uncorrected insertions in
the site centre +/- 100 bp), 95% interval from 200 within-decile site bootstraps:

| motif | ChIP | sites | bound | AUROC, insertions | AUROC, score within decile | F-047 homemade |
|---|---|---|---|---|---|---|
| CTCF (MA0139.2) | CTCF ENCFF519CXF | 83,492 | 28,726 | 0.777 | 0.752 [0.748, 0.756] | 0.694 |
| GATA1 (MA0035.5) | GATA1 (K562_GATA1_BPNet peaks) | 14,903 | 263 | 0.953 | 0.674 [0.636, 0.702] | 0.419 |
| GATA1::TAL1 (MA0140.3) | GATA1 | 42,262 | 2,165 | 0.954 | 0.674 [0.660, 0.687] | 0.497 |

BINDetect's bound call covers 82% / 91% of ChIP-bound GATA1 / GATA1::TAL1 sites vs 9% / 14% of
unbound ones. Mean footprint score per motif (879 motifs): top 20 are again GC/CpG-rich (KLF15,
ZBED4, ETS, NRF1, ZBTB33, SP); CTCF ranks 77, GATA1::TAL1 274, GATA1 534, below HNF4A (516).
**Implications:** A published bias model recovers site-level GATA evidence that the homemade
score lacked, so feature B is viable at the site level: aggregate per element (e.g. max or sum
of GATA/TAL1/NFE2 site scores) and take it to the beyond-ATAC p300 gate. Do not use
per-motif mean scores (BINDetect's ranking) as TF activity: they are composition-biased.
**Caveats:** Site sets differ from F-047 (TOBIAS's own motif scan inside candidate regions vs
4.44 genome-wide FIMO hits with >= 10 insertions), so the comparison is not hit-matched. GATA1
ChIP set is small (263 bound sites in the MA0035.5 test). Accessibility deciles use a 200 bp
window; residual within-decile accessibility could still contribute. scPrinter untested.
**Tags:** footprinting, atac, multiome, tn5-bias, tobias, gata, ctcf, k562, sce2g

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-10-02 | 2.1 insertion BAM 46354090/46356012; 2.2 TOBIAS 46373951; 2.3 46380940 | Xu K562 multiome fragments; scE2G candidate regions; K562 CTCF and GATA1 ChIP | 2026_1002_sce2g_tf_activity | GATA1 0.674 [0.636, 0.702], GATA1::TAL1 0.674 [0.660, 0.687], CTCF 0.752 within accessibility deciles; GC-rich motifs top per-motif ranking | footprints work for GATA at site level |

## F-049: Per-element TOBIAS footprints of erythroid TFs track observed K562 p300 beyond ATAC and beyond motif presence (partial rho 0.09-0.10); feature B passes its gate
**Status:** supported (one cell, one dataset; gate only, not yet in scE2G)
**Claim:** `2026_1002_sce2g_tf_activity/scripts/2.4` over scE2G's 135,965 non-promoter candidate
elements (Xu run), features from the F-048 TOBIAS sites. Partial Spearman with observed K562 p300
(ENCSR000EGE, +/- 250 bp), 95% intervals from 200 element bootstraps. Controls Z = Xu ATAC
(+ square), GC, CpG o/e, log distance to TSS; the sensitivity run adds ENCODE K562 ATAC (+ square)
to soak up Xu pseudobulk noise.

TEST 2, elements containing the motif, footprint max beyond Z + motif count + motif strength:

| motif set | elements | Z = Xu ATAC | Z + ENCODE ATAC |
|---|---|---|---|
| erythroid (GATA1, GATA1::TAL1, NFE2, MAF::NFE2, TAL1::TCF3) | 85,720 | 0.086 [0.080, 0.092] | 0.104 [0.098, 0.110] |
| GATA1::TAL1 | 32,914 | 0.096 [0.083, 0.104] | 0.116 [0.104, 0.126] |
| NFE2 | 28,477 | 0.067 [0.055, 0.079] | 0.092 [0.080, 0.105] |
| GATA1 | 12,743 | 0.007 [-0.010, 0.025] | 0.032 [0.015, 0.050] |
| HNF4A (control, not expressed) | 33,124 | 0.031 [0.020, 0.042] | 0.019 [0.009, 0.031] |
| CTCF (control) | 52,898 | -0.004 [-0.014, 0.005] | -0.044 [-0.054, -0.035] |

TEST 1, all elements (no site = lowest rank), Z + ENCODE ATAC: erythroid fp_sum 0.120
[0.113, 0.125], erythroid motif count alone 0.102 [0.098, 0.107], TF-agnostic max footprint
0.020 [0.014, 0.025]. With Xu ATAC alone in Z, ENCODE ATAC itself adds 0.248 [0.242, 0.254]:
the Xu pseudobulk ATAC is noisy, and a second ATAC measurement is worth more than any footprint.
**Implications:** Feature B passes: the erythroid footprint carries p300 information that neither
ATAC nor motif presence gives, about 0.085 above the non-expressed-TF control, roughly a quarter
of observed H3K27ac's 0.376 (F-045, different region set). Take erythroid fp_sum / fp_max into
scE2G as external features against the baseline. CTCF footprints point the other way (insulator
sites carry less p300), so keep them separate rather than in a pooled score.
**Caveats:** Gate is correlation with p300, not CRISPR: the scE2G paired bootstrap decides. One
cell; the erythroid motif set is chosen with K562 knowledge, so a transfer version needs
expression-weighted motif sets. The HNF4A control is not zero, so part of each motif's effect is
generic footprint shape. Single TOBIAS run, default parameters.
**Tags:** footprinting, tobias, p300, beyond-atac, gate, gata, nfe2, k562, sce2g

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-10-03 | 2.4 46419803 (Z = Xu ATAC), 46420329 (`--z-encode-atac`) | TOBIAS sites (F-048), scE2G candidates, K562 p300 ENCSR000EGE, ENCODE K562 ATAC | 2026_1002_sce2g_tf_activity | erythroid footprint beyond motif 0.086 [0.080, 0.092] / 0.104 [0.098, 0.110]; HNF4A control 0.031 / 0.019 | supports feature B |

## F-050: TOBIAS erythroid footprints as scE2G features do not improve K562 CRISPR AUPRC (-0.004 [-0.009, +0.002]), despite passing the p300 gate
**Status:** supported (one cell, one footprint method, one model family)
**Claim:** scE2G `multiome_powerlaw_n6` (v3 features) retrained on the Xu K562 run with the
F-049 element features added as external features (`scripts/4.0`, `4.1`; results
`$SCRATCH/sce2g_runs/2026_1003_k562_xu_featB_tobias`). Paired bootstrap (`scripts/4.2`, 2,000
resamples of the 10,342 CRISPR pairs incl. the 1,075 with no candidate element at score 0; 466
positives) of held-out-chromosome predictions:

| model | AUPRC [95% CI] | vs baseline |
|---|---|---|
| baseline (retrained; predictions identical to the first baseline run) | 0.632 [0.587, 0.677] | |
| + erythroid motif count | 0.630 [0.585, 0.675] | -0.0022 [-0.0037, -0.0007] |
| + motif count + erythroid TOBIAS footprint sum | 0.628 [0.582, 0.674] | -0.0036 [-0.0090, +0.0019] |

Footprint beyond motif count: -0.0014 [-0.0067, +0.0040]. The model does use the footprint
(coefficient +0.30, with the motif count turning to -0.22), and ARC.E2G.Score's coefficient
falls from 1.73 to 1.59, so the footprint mostly re-expresses what ARC-E2G already carries.
scE2G's own table reports 0.629 / 0.626 / 0.625 (a bootstrap mean, so ~0.003 lower).
**Implications:** Passing the p300 gate (partial rho ~0.1, F-049) was not enough to move CRISPR
AUPRC: an element feature that small, layered on ARC-E2G (which already contains ABC activity
and Kendall), is below what 466 positives can detect, if it helps at all. Feature B in this form
is closed for TOBIAS; the PRINT version gets the same test before closing the footprint line.
Reopen only with a feature whose beyond-ATAC p300 partial is several times larger, or as an
input to the activity term itself (ABC activity) rather than an additive scE2G feature.
**Caveats:** One model family (multiome v3, logistic, 6 features); the erythroid motif set is
chosen with K562 knowledge; features enter linearly and unscaled; the CRISPR set's positives are
concentrated in a few loci, so element-level gains elsewhere would not show.
**Tags:** footprinting, tobias, sce2g, crispr, external-features, null, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-10-03 | setup 46441027; driver 46443820; bootstrap 46450290 | Xu K562 multiome, K562 CRISPR ensemble benchmark | 2026_1002_sce2g_tf_activity | fp model - baseline -0.0036 [-0.0090, +0.0019]; fp - motif -0.0014 [-0.0067, +0.0040] | null for feature B (TOBIAS) |

## F-051: PRINT (scPrinter) footprints beat TOBIAS on CTCF ChIP but not GATA, pass the p300 gate, and are also null as scE2G features (-0.002 [-0.006, +0.003]); feature B is closed
**Status:** supported (one cell, pretrained PRINT TFBS model, same sites as TOBIAS)
**Claim:** PRINT TF binding scores (`2026_1002_sce2g_tf_activity/scripts/3.1`; pretrained
`TFBS_model_py.pt`, hg38 v2 Tn5 bias, scE2G's 7,821 filtered Xu cells as one pseudobulk) at the
304,473 BINDetect sites of the seven 2.4 motifs, each scored at its 10 bp centre on its strand
with motif score fixed at 1 (footprint only). Same tests as F-048 to F-050, via `--score print`:

| test | PRINT | TOBIAS |
|---|---|---|
| ChIP AUROC within accessibility deciles, CTCF | 0.851 [0.848, 0.854] | 0.752 [0.748, 0.756] |
| GATA1 | 0.605 [0.573, 0.644] | 0.674 [0.636, 0.702] |
| GATA1::TAL1 | 0.564 [0.549, 0.576] | 0.674 [0.660, 0.687] |
| p300 partial, erythroid fp beyond motif, Z = Xu ATAC | 0.111 [0.104, 0.117] | 0.086 [0.080, 0.092] |
| same, Z + ENCODE ATAC | 0.126 [0.119, 0.132] | 0.104 [0.098, 0.110] |
| HNF4A control (Xu / + ENCODE) | 0.052 / 0.053 | 0.031 / 0.019 |
| scE2G CRISPR AUPRC, fp model - baseline | -0.0020 [-0.0061, +0.0026] | -0.0036 [-0.0090, +0.0019] |
| fp model - motif-count model | +0.0002 [-0.0040, +0.0046] | -0.0014 [-0.0067, +0.0040] |

PRINT fp model - TOBIAS fp model: +0.0016 [-0.0026, +0.0062] (paired bootstrap, `4.2`, 10,342
pairs, 466 positives). Footprint coefficient +0.36, motif count -0.27, ARC.E2G.Score 1.69.
**Implications:** Two published footprinting methods give the same answer: erythroid footprints
carry ~0.1 partial rho of p300 beyond ATAC and motif presence (of which ~0.06-0.08 is above the
non-expressed-TF control), and none of it moves CRISPR AUPRC once ARC-E2G is in the model.
Feature B (element-level footprints as additive scE2G features) is closed. Reopen only with
(1) a footprint signal several times larger against p300 (e.g. seq2PRINT, trained on this
sample), or (2) footprints entering the activity term itself (ABC activity = ATAC x footprint
weight) rather than as a separate feature.
**Caveats:** Pretrained PRINT model, not seq2PRINT; motif set chosen with K562 knowledge; one
model family; the first full PRINT run lost three of four shards to scPrinter's in-memory
(backed=False) result handling and one Bus error, rerun in backed mode (identical scores on the
2,000-site check).
**Tags:** footprinting, scprinter, print, tobias, sce2g, crispr, p300, null, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-10-04 | 3.1 arrays 46428901 (shard 1), 46470149-51; 2.3 46480634; 2.4 46480635/6; scE2G 46481346; 4.2 46487342 | Xu K562 multiome, K562 CTCF/GATA1 ChIP, K562 p300, K562 CRISPR | 2026_1002_sce2g_tf_activity | PRINT fp - baseline -0.0020 [-0.0061, +0.0026]; gate 0.126 [0.119, 0.132] | null; feature B closed |

## F-052: Peak-TF Kendall (feature A) carries almost no p300 information beyond ATAC and generic peak covariation (partial rho 0.015-0.028); not taken to scE2G
**Status:** supported (one cell, one dataset, binarised ATAC)
**Claim:** `2026_1002_sce2g_tf_activity/scripts/5.1` over scE2G's Kendall ATAC matrix (157,587 peaks
x 7,821 Xu cells, binarised; exported by `5.0`) and log1p(CP10k) RNA of the 348 motif TFs detected
in >= 5% of cells. Kendall tau-b per peak x TF (rank-sum form, matches scipy to 1.6e-9); 4.44 FIMO
hits assign 10.9M peak-TF pairs. Element features = max over overlapping Kendall peaks (96.7% of
the 155,134 candidates). Gate as 2.4: partial Spearman with K562 p300 (ENCSR000EGE), 135,932
non-promoter elements, 200 element bootstraps.

| feature / test | Z = Xu ATAC | Z + ENCODE ATAC |
|---|---|---|
| kA_max (max tau over TFs with a motif) | 0.042 [0.037, 0.047] | 0.069 [0.064, 0.075] |
| kA_max_perm (motif -> TF labels permuted) | 0.029 [0.024, 0.035] | 0.050 [0.045, 0.056] |
| kA_any_max (max tau over all 348 TFs) | 0.048 [0.042, 0.053] | 0.074 [0.068, 0.079] |
| kA_ery_max (GATA1, TAL1, NFE2, KLF1) | 0.022 [0.016, 0.027] | 0.037 [0.032, 0.043] |
| TEST 2, kA_max beyond Z + n_tf | 0.035 [0.031, 0.040] | 0.060 [0.056, 0.066] |
| TEST 2, + perm + any (strict) | 0.015 [0.009, 0.019] | 0.028 [0.022, 0.034] |

Tau is small throughout (median 0.010, max 0.175 over 1M random peak-TF pairs).
**Implications:** Most of the peak-TF Kendall signal is generic: a peak that co-varies with any TF
(or with the wrong TF for its motifs) carries nearly as much p300 information as one that co-varies
with its own motif's TF. The motif-specific part, 0.015-0.028, is a fifth to a quarter of the
footprint signal (~0.1, F-049/F-051) that already failed to move CRISPR AUPRC, so feature A does not
get a scE2G run. K562 is a cell line with little TF expression variation across cells; the feature
may behave differently in a heterogeneous sample.
**Caveats:** Binarised ATAC (counts > 1 are 20% of non-zeros); one detection threshold (5%); the
perm control is one permutation; no metacell smoothing, which scE2G's Kendall also does not use.
**Tags:** kendall, multiome, tf-expression, p300, beyond-atac, gate, null, k562, sce2g

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-10-04 | 5.0 srun 46503422; 5.1 46503952 (Z = Xu ATAC), 46503953 (`--z-encode-atac`) | Xu K562 multiome ATAC + RNA, 4.44 FIMO hits, K562 p300 | 2026_1002_sce2g_tf_activity | strict gate 0.015 [0.009, 0.019] / 0.028 [0.022, 0.034] | weak; feature A not taken to scE2G |

## F-053: Footprints inside the ABC activity term do not improve scE2G CRISPR AUPRC; mild weighting is within noise (+0.004 to +0.005) and strong weighting hurts (-0.012 [-0.021, -0.004] vs mild)
**Status:** supported (one cell, two footprint methods, two weights)
**Claim:** `2026_1002_sce2g_tf_activity/scripts/4.3`, `4.4`: element activity reweighted as
ATAC x (1 + lam x q), q = percentile of the element's erythroid footprint sum (TOBIAS or PRINT,
F-049/F-051 element tables; 0 without an erythroid site, 61% of elements have one), ABC renormalised
per gene, and ARC-E2G recomputed (ARC = ABC x exp(21.5 x Kendall) in scE2G's output, R^2 1.0000).
Each model is the multiome v3 n6 set with ARC.E2G.Score swapped for the footprint-weighted score.
Run `$SCRATCH/sce2g_runs/2026_1004b_k562_xu_abcfp`; paired bootstrap (`4.2`, 2,000 resamples,
10,342 CRISPR pairs, 466 positives) of held-out-chromosome predictions:

| model | AUPRC [95% CI] | vs baseline |
|---|---|---|
| baseline n6 (ARC-E2G) | 0.632 [0.587, 0.677] | |
| TOBIAS, lam 1 | 0.635 [0.590, 0.681] | +0.0036 [-0.0052, +0.0123] |
| PRINT, lam 1 | 0.637 [0.591, 0.681] | +0.0052 [-0.0030, +0.0136] |
| TOBIAS, lam 4 | 0.623 [0.576, 0.670] | -0.0084 [-0.0243, +0.0069] |
| PRINT, lam 4 | 0.625 [0.578, 0.671] | -0.0062 [-0.0221, +0.0095] |

TOBIAS lam 4 - lam 1: -0.0120 [-0.0205, -0.0044]; PRINT lam 1 - TOBIAS lam 1: +0.0016 [-0.0022,
+0.0058]. The weighted score's coefficient falls from 1.73 (ARC-E2G) to 1.57-1.59 at lam 1 and
1.19-1.22 at lam 4.
**Implications:** Moving footprints from an extra feature (F-050/F-051) into the activity term
changes the sign of the point estimate but not the answer: no weight tested moves CRISPR AUPRC
beyond noise, and stronger weighting degrades the activity term. With feature A deferred (F-052),
element-level TF evidence from this multiome sample (footprints, peak-TF Kendall) does not improve
scE2G on the K562 benchmark. The lam 1 PRINT point estimate (+0.005) is the only lead, and it would
need replication (another sample or cell) before it counts.
**Caveats:** Two weights, chosen before the run, not tuned; the erythroid motif set uses K562
knowledge; the first run of this design (`2026_1004_k562_xu_abcfp`) is invalid because all four
models shared one feature name (see learnings 2026-10-04); one model family.
**Tags:** footprinting, abc, activity, tobias, print, sce2g, crispr, null, k562

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-10-04 | setup 46515688; driver 46515759; bootstrap 46515762 | Xu K562 multiome, TOBIAS/PRINT element footprints, K562 CRISPR | 2026_1002_sce2g_tf_activity | lam 1 +0.004 / +0.005 (CIs span 0); lam 4 -0.008 / -0.006; lam 4 - lam 1 -0.012 [-0.021, -0.004] | null; footprint line closed |
