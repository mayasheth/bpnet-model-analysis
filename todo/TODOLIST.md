# Todo List

Open work only. Completed items live in the git history and in `.living/`; the analysis
narrative is in `2026_0824_H3K27ac_model/h3k27ac_model_report.html`.

## Settled (do not redo)

Target = 5′ ends, ±500 bp window, `count_loss_weight = 10`. Panel is **ATAC-only** (DNase
and ATAC are not interchangeable inputs). PE H3K27ac targets use **read 1 only**.
Accessibility inputs should be ChromBPNet-style 5′ insertion counts; tracks are built and
validated but **no model uses them yet**. Residual-objective training helps only a
sequence-blind input and costs a multimodal one — replicated in K562 and GM12878, so it is
a property of the objective and needs no further per-cell-type testing.

## Highest value

- [ ] **W4. Motif syntax (SHAP / TF-MoDISco / FiNeMo).** Not started; the actual scientific
      goal. Run on the **residual-trained** model — its attributions are forced onto
      accessibility-independent signal, which multimodal attributions cannot separate. Use
      the ±500 bp window (zero neighbour contamination). Expect less signal than p300.

- [ ] **TeloHAEC training + transfer.** The only new cell type available under the ATAC-only
      rule. Tracks, elements and model-free coupling (0.33–0.37 top quintile) are all ready.
      3 modes × 5 folds, then the four-evaluation transfer set against K562 and GM12878.
- [ ] **TeloHAEC ±IL1b / ±TNFa / −VEGF.** Same genome and cell line, different regulatory
      state, no input-domain shift — a sharp and cheap test of whether the model tracks
      condition-specific change. Inference only if trained on ctrl.

## Interrogating the ABC negative result — the main thread

The CRISPR benchmark says predicted H3K27ac adds nothing detectable over ATAC alone (best
predicted arm 0.482 vs floor 0.457, CI overlapping; deployment-scenario arms at or below the
floor). `4.9.diagnose_abc_gap.py` localised why: the prediction is accurate (Spearman 0.82 on
CRISPR-tested regions, better than genome-wide) but **ranks the functional elements too low**,
and since ABC's qnorm removes scale by construction, rank is the only channel available. All
arms catch the same positives, so the deficit is in suppressing negatives.

`4.10`/`4.12` then localised the error to a specific population: every model that sees ATAC
over-predicts H3K27ac at accessible-but-unacetylated elements by 7–8× its own median, while
sequence-only elevates them 1.6× and observed H3K27ac not at all. Those elements are GC 0.59,
CpG o/e 0.55, 2.4× promoter-enriched, ATAC/K27ac ratio 13.25 vs 2.20 typical — a CpG-island /
CTCF phenotype. **The signal is in sequence; the additive trunk gives it no way to veto
accessibility.**

- [~] **RUNNING: gate × loss factorial** (`1.19`, 15 fold-jobs). Sequence-conditioned gate on
      the accessibility branch (`X_acc *= sigmoid(conv(X_seq))`, initialised open so the model
      starts as the ungated one), crossed with an asymmetric count loss weighting
      over-prediction 3×. `gate` alone is expected to do little because the gate has no
      gradient pressure while log1pMSE is ~400× more sensitive to missing signal than to
      inventing it — which is why it gets its own arm instead of being assumed.
      Backward compatibility for the shared `multimodal_bpnet.py` is gated by
      `scripts/test_asymmetric_loss.py`: weight 1.0 reproduces bpnetlite exactly, an open gate
      changes predictions by 0.39%, and pre-gate checkpoints still load.
- [ ] **Indicator-channel control.** GC and CpG-density tracks as extra accessibility
      channels, no gate, no loss change. If hand-supplied class information does as well as a
      learned gate, prefer it — simpler and interpretable. Both are sequence-derived, so
      unlike fragment channels they cost nothing in transferability. Needs `0.27` to build the
      tracks first.
- [ ] **CTCF/EP300 enrichment in the error strata** (`4.11`, submitted). Direct test against
      ENCODE peak calls rather than a PWM proxy. Prediction: over-predicted tail enriched for
      CTCF and depleted for EP300; under-predicted the reverse.
- [ ] **Characterise the 12 threshold-level false negatives** individually once `4.11` lands —
      observed H3K27ac 39× the genome median where the model predicts near-background.

## Sequence-gated multimodal for p300 — diagnose before training

Naming: the architecture is **sequence-gated multimodal**; `gate` is the flag and filename
token. The asymmetric loss is named separately on purpose, because the factorial exists to
attribute any gain to the architecture or to the loss.

- [ ] **Confirm p300 has the same failure mode BEFORE training p300 with the gate.** The
      premise is that p300 suffers the same accessibility domination, and it may not: p300's
      residual *r* is 0.654 against H3K27ac's ~0.55, so sequence adds MORE beyond accessibility
      there. That is consistent with the problem being milder, or with it being equally severe
      but more fixable — different expected gains, same number.
      Cheap test, reusing existing machinery: predict p300 on the ABC candidate regions with
      the existing `2026_0529_multimodal_p300_model` models via `4.1`, then run the `4.12`
      fold-elevation comparison against observed p300 (ENCSR000EGE,
      `ENCFF466WKF`/`ENCFF163FSR`). Confirmed if p300's sequence-only arm gets high-ATAC /
      low-p300 elements right while its multimodal arm does not.
- [ ] **Then train p300 with the gate**, 15 fold-jobs, only if the diagnostic confirms it.
      `multimodal_bpnet.py` is already shared, so no further code change is needed — the p300
      submit scripts just gain the two flags.

## Which model transfers best — do not assume the in-cell-type winner

- [~] **RUNNING: transfer matrix** (`2.24`). narrow/wide × flat/fragments, both directions,
      with each target's own models in the same table so the transfer drop is readable.
      The in-cell-type ranking is NOT the deployment ranking: wide gained +0.028/+0.016
      in-cell-type but only +0.012/+0.005 transferred, neither significant. Fragment channels
      have never been tested on transfer and are the arm most at risk, since fragment-size
      distributions are library properties.
- [ ] **TeloHAEC as a third cell type, with conditions.** The only new cell type under the
      ATAC-only rule. Tracks, elements and model-free coupling are ready. Caveats to carry:
      ATAC-derived elements (though derivation was shown not to matter, p = 0.83), 36 bp reads
      vs 95 bp, and shallower libraries — so read length and depth remain confounded even
      though element derivation does not. Fragment channels there would need PE BAMs.
- [ ] **Train on multiple cell types to optimise transferability.** The most promising route
      to a deployable model, since it directly optimises what deployment needs rather than
      in-cell-type fit. Do it after the transfer matrix says which architecture to carry.

## p300 as the activity term instead of H3K27ac

- [ ] **ABC + CRISPR with p300 predictions, and with observed p300.** p300 is a coactivator,
      one step closer to the sequence-specified event, and is arguably what ABC's activity
      term is really proxying. Observed p300 gives both a ceiling and — importantly — a qnorm
      REFERENCE, so the no-reference problem dissolves: build a p300 qnorm reference from the
      observed arm rather than running without qnorm and having to scale-match the geomean by
      hand. Data: EP300 ENCSR000EGE, `ENCFF466WKF/ENCFF163FSR.filtered.sorted.bam`, peaks
      `ENCFF702XPO.bed.gz`. p300 models already exist in `2026_0529_multimodal_p300_model`
      (residual r 0.654, higher than H3K27ac's).
- [ ] **Test removing qnorm** as its own arm (`use_qnorm: False`; `activity_base_no_qnorm` is
      already in every EnhancerList). Expectation is that it hurts — the model emits log1p
      counts over ±500 bp while observed H3K27ac is read counts over the element, so without
      qnorm the geomean multiplies incommensurate magnitudes, and ABC's thresholds are
      calibrated on qnorm'd values. Cheap enough to settle rather than argue.

## Downstream utility — decides whether the correlation metrics are the right target

- [ ] **Plug predicted H3K27ac into ABC and benchmark it.** The end-to-end test of whether
      this model is useful, and the concrete form of the standing question below. ABC's
      activity term is `geomean(accessibility, H3K27ac)`, so substitute the prediction for
      the observed mark, run ABC, and score against the CRISPR benchmark.
      **Three arms, or the result is uninterpretable:** observed H3K27ac (upper bound),
      predicted H3K27ac, accessibility alone (floor). If predicted lands near observed, the
      model is useful even at a top-quintile *r* of 0.69; if it lands near accessibility
      alone, it adds nothing downstream and the correlation gains we have been chasing do
      not matter.
      **Leakage trap:** ABC runs genome-wide but each fold model has seen four fifths of the
      genome. Assemble the genome-wide prediction from the 5 fold models, each applied only
      to its own held-out chromosomes, or the benchmark is contaminated.
      Repos: `~/Documents/ABC-Enhancer-Gene-Prediction`, CRISPR benchmark via
      `~/Documents/DC_TAP_Paper`.

- [ ] **ATAC → DNase converter.** DNase tracks H3K27ac and enhancer activity better than
      ATAC, so a converter would give DNase-like input in the cell types that only have ATAC.
      **Cheap gate first, before building anything:** measure model-free coupling of DNase
      vs ATAC against H3K27ac in K562 and GM12878 using the existing machinery
      (`0.12.atac_vs_h3k27ac.py`, `3.6`). Minutes of CPU, and it bounds what a perfect
      converter could buy. No DNase bigwigs exist in this project yet; the DNase used to
      call the element sets came from the rE2G runs (`reference/ELEMENT_DERIVATION.md`).
      **Confound in that gate:** the K562 and GM12878 element sets are DNase-derived, which
      favours DNase on element definition alone. Repeat on the ATAC-derived K562 set in
      `K562_ATAC_ChromBPNet/data/` before believing the gap.
      **If the gap is real:** train ATAC → DNase where both assays exist (K562, GM12878),
      then score H3K27ac prediction three ways — raw ATAC (floor), converted DNase, real
      DNase (ceiling). Without the real-DNase arm, a gain cannot be separated from the extra
      capacity the converter adds.
      This does not contradict the ATAC-only panel rule. That rule exists because ATAC and
      DNase are not interchangeable as *inputs*; a converter is the principled way to get a
      DNase-like input everywhere without mixing assays across cell types.
      Also unblocks the composite-metric item below, whose best form needs DHS.

## Adopt-or-not decisions, all measured

- [ ] **Turn on test-time RC averaging as the default?** Free (one extra forward pass, no
      retraining) and positive for every model: +0.0162 sequence only, +0.0081 sequence +
      ATAC, +0.0016 and non-significant for ATAC only, which is the control (report Fig. 13).
      Currently opt-in via `2.15 --rc-average` so existing numbers stay comparable. Adopting
      it means re-scoring the report's tables as a set.
- [ ] **Do the wider receptive field and the fragment channels combine?** `n_layers` 10 ×
      5 fragment channels, 5 fold-jobs. Both act on the accessibility side and may read the
      same neighbourhood structure, so +0.027 and +0.0135 may not sum. This decides what the
      deployed model is, so it should run before the ABC arms are treated as final.

## Open questions

- [ ] **Repeat the residual comparison for p300.** Resolved only for H3K27ac; p300's
      residual *r* is already 0.654 so the headroom may be smaller.
- [ ] **Predict the composite activity metric directly.** Downstream uses
      `geomean(accessibility, H3K27ac)`, best as `geomean(DHS, H3K27ac)`. Is the composite
      easier to predict than predicting H3K27ac and combining afterwards?
      **Design caution:** if the accessibility input is the same assay as the accessibility
      term in the target, half the target is readable off the input and the comparison is
      circular. Predict `geomean(DHS, H3K27ac)` from sequence + ATAC, and baseline against
      the two-step route scored on the *same* composite.
- [ ] **[Q for Maya] Is a mostly-accessibility model useful for the intended application,**
      or does the goal require the sequence component? The ABC benchmark above is the
      empirical form of this question and would settle it without needing an answer in
      advance. Also decides whether the +0.043 all-elements receptive-field gain matters:
      it is real, but it is entirely dead-vs-active separation, so it pays off only for an
      on/off task.

## Architecture, in expected order of value

- [~] **Wider receptive field — ADOPT for multimodal; the earlier "do not pursue" was wrong.**
      `n_layers` 10 (~4.2 kb vs ~1.1 kb), 5 folds, both cell types, paired within fold on the
      intersection of valid regions.
      Top quintile: sequence-only −0.006 (K562) and +0.010 (GM12878), both p=0.53;
      **multimodal +0.027 (p=0.006) and +0.014 (p=0.025)**, all five folds rising in both.
      Accessibility residual 0.502 → 0.547 and 0.397 → 0.469.
      Results in `wide_{k562,gm12878}_*.tsv`; report Fig. 10.
      - [ ] **ATAC-only wide arm is RUNNING** (10 fold-jobs, submitted 2026-09-03). If it
            reproduces the multimodal gain, the extra context is used purely as accessibility
            neighbourhood and sequence contributes nothing to it — which would also mean the
            deployed model should simply be widened.
      - [ ] Decide whether to re-run the transfer and deployment comparisons at
            `n_layers` 10, since those used the narrow multimodal model.
- [ ] **Decide what to do with the profile head — one decision, three options.** The 1 bp
      profile task is close to unlearnable: measured inter-replicate ceiling is 0.21 (K562)
      and 0.18 (GM12878) on the top quintile, rising to 0.72 and 0.70 at 50 bp binning
      (`profile_ceiling_binsize_*.tsv`, report Fig. 11). The options are (a) drop the profile
      loss entirely, (b) keep 1 bp, (c) bin to ~50 bp.

      **Test by hyperparameter before touching the architecture.** `loss = profile_loss +
      w · count_loss`, so a large `w` already approximates dropping the profile term, and the
      recorded sweep points AGAINST dropping it: `w` = 10 gave 0.496 while 100 gave 0.467 and
      1000 gave 0.464 (header of `1.11.submit_training_5prime_accs5p.sh`). Down-weighting the
      profile loss made the COUNTS worse. A head can be a poor predictor and still be a
      useful auxiliary task on the shared trunk, and those numbers are single-fold, so this
      is weak evidence pointing the opposite way from removal rather than a settled answer.
      The properly powered `count_loss_weight` sweep already listed under Statistical power
      answers it at zero architecture risk — run it jointly with bin size, since the two
      interact: binning changes how learnable the profile term is and therefore its optimal
      weight.

      **Gate on the real numbers.** Everything above about what the head achieves rests on a
      256-element smoke test (`profile_pearson` ~ 0.053 against a 0.21 ceiling). Job for the
      full 5-fold measurement is queued; use `prof_residual_grid_*.tsv` when it lands.

      **Cost of actual removal** (option a): it changes the model class, so every existing
      checkpoint becomes non-comparable and the whole grid needs retraining. It touches
      `multimodal_bpnet.py` and `train_multimodal_bpnet.py`, both shared with p300, so it
      needs a backward-compatibility regression — use the `2.18` pattern, exact on region
      counts and tolerant on metrics.
      **Payoff if it holds:** removes a term that is largely fitting Poisson noise, frees
      trunk capacity, and roughly halves the output tensor, buying a larger batch or a wider
      window at the same memory.
- [x] **Fragment-size ATAC channels — DONE and positive.** Top quintile 0.690 → 0.703,
      paired +0.0135 [+0.0075, +0.0195], p = 0.0034, every fold rising; report Fig. 12.
      Strict superset of the flat input (bins sum to it exactly), so the gain is added
      information. Results in `fragchan_k562_per_fold.tsv`.
      - [ ] **Do the two accessibility-side gains combine?** `n_layers` 10 × fragment
            channels. Both may be reading the same neighbourhood structure, so +0.027 and
            +0.0135 may not sum. 5 fold-jobs, and it decides what the deployed model is.
      - [ ] **GM12878 replication** needs its paired-end BAMs downloaded — fragment length
            lives in TLEN, which the per-read tagAligns discard.
      - [ ] **ATAC-only fragment arm** (`1.15 atac`) — asks whether fragment structure alone
            beats flat ATAC, with no sequence involved.
- [ ] **Switch the accessibility input to 5′ counts, as a set.** ~35 fold-jobs. Within a
      cell type this moves only the ATAC-only model; across cell types it removes a
      read-length confound (TeloHAEC 36 bp vs K562 95 bp). Mixing the two inputs is invalid.
- [ ] **Re-check ±500 vs ±1000** on the 5′ target once a retrain happens anyway.
- [ ] **Train on the ATAC-derived K562 element set** — easy (one path change in `1.11`),
      and now load-bearing for the ABC work rather than merely tidy.
      Verified path, 153,545 regions:
      `ENCODE_rE2G/results/2025_0226_ATAC_powerlaw_models/ATAC_H3K27ac_powerlaw/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed`
      (element set with class labels: the sibling `Neighborhoods/EnhancerList.bed`).
      **This file is byte-identical (md5 7d5995ce…) to the candidate regions of the July ABC
      run**, so the ABC experiment scores models on the ATAC-derived set while they were
      trained on the DNase-derived one (150,528 elements). The windows are convolutional so
      it is not fatal, but it is a genuine train/score element mismatch inside the ABC
      comparison, and training on this set removes it.
      Also the cheapest test of whether element derivation matters at all; if it does not,
      the panel-wide inconsistency stops being a caveat worth carrying.

## Panel data caveats

- [ ] **Element derivation is not uniform** (`reference/ELEMENT_DERIVATION.md`): K562 and
      GM12878 DNase-derived, TeloHAEC ATAC-derived. An ATAC-derived K562 set exists in
      `K562_ATAC_ChromBPNet/data/`; adopting it makes two of three consistent at the cost of
      re-deriving every K562 number.
- [ ] **HCT116** has ENCODE ATAC but its H3K27ac replicates differ in run type, so no
      inter-replicate ceiling is computable. Trainable, not normalisable.
- [ ] **H1, H9, Jurkat, THP-1 have no ENCODE ATAC at all** (all 559 released experiments
      checked). Reaching them means GEO/SRA fastqs through
      `Data/scripts/sra_paired_fastq_to_bam.sh` — a separate decision.
- [ ] **TeloHAEC_ctrl/ATAC holds 3 EA.hy926 files** (`SRR20809434/435/436`) under the same
      sample name. Always use explicit accession lists, never a directory glob.

## Statistical power

- [ ] **Re-run the `count_loss_weight` sweep with ≥3 folds and ≥2 seeds.** The current pick
      (10) came from single folds; 3/10/100 are within noise of each other.
- [ ] **Quantify run-to-run variance properly** — one config × 5 seeds.
- [ ] **Make submit scripts refuse to overwrite a completed fold directory.** A grid once
      silently overwrote a sweep result that shared an `OUT_DIR`.

## Housekeeping

- [ ] Delete the local git tag `backup-pre-msg-rewrite` on Oak once the rewritten history is
      confirmed good.
- [ ] `Data/ENCODE/K562/interm_ENCFF790GFL.se.filtered.sorted.bam/` is a stray directory.
- [ ] `scripts/0.3.make_training_bw.sh` does not sort its experiment bedGraphs before
      merging; works today only because inputs happen to be sorted.
