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
- [~] **Fragment-size ATAC channels** — tracks BUILDING, training not yet submitted.
      Five channels `[all, sub(≤139), mono(140–329), di(330–620), poly(≥621)]`, all
      single-base insertion counts, so they are a strict superset of the accs5p input and
      cannot regress. `0.23` builds them; the Tn5 shift was measured, not assumed (`0.22`,
      r = 1.0000 at +4/−5). Train with `1.15`, modes multimodal and atac.
      The old `atac_sub.bw`/`atac_mono.bw` are full-fragment coverage and are superseded.
      GM12878 needs its PE BAMs downloaded before the same thing is possible there.
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
