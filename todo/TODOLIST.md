# Todo List

Master list of future work items. Each item can have a detailed writeup
in a separate `.md` file in this directory.

## Items

### H3K27ac model — `h3k27ac-model.md`

Detailed writeup: [`h3k27ac-model.md`](h3k27ac-model.md). Ordered into waves by
dependency; everything inside a wave can run at the same time.

**Gate — DONE 2026-08-25.** 5′ and fragment-extended ceilings are *identical* on the
top quintile (0.760 vs 0.761 at ±500). Switched to 5′ on principle — MNLL well-posed, no
neighbour bleed — not for an accuracy gain. See `.living/learnings.md`.

- [x] **G1. Inter-replicate ceiling on the 5′ target** — decides whether 5′ ends or the
      250 bp fragment track is the model target. Every later training decision inherits
      this. ~10 min CPU. `0.3.replicate_ceiling_by_window.py` with the new
      `h3k27ac_rep{1,2}_5p_{plus,minus}.bw`. The 5′ track is sparse (max 27, ~6 reads per
      kb at background), so the ceiling could fall; if it does, test unextended 36 bp
      read coverage as a middle option.

**Wave 1 — independent of G1, all parallel**

- [x] **W1a. Residual EVALUATION** — done (jobs 40799753, 40803766). Residual metric
      is now the headline; it reverses the p300/H3K27ac ranking (F-002). Residual
      *training* is now item W2e below, promoted by that result.
- [ ] **W1a-bis. Residual / difference-from-ATAC model + evaluation** — the highest-value
      item and the one that matches the actual goal. Stage 1 already trained. ~4 GPU-h.
- [ ] **W1b. GM12878 cross-cell-type transfer** — separates model capacity from
      "H3K27ac is not sequence-determinable". Inference only on existing K562 models.
      ~0.5 GPU-h. Data all present (`ENCFF645BAL`, `ENCFF865OOP`).
- [x] **W1c-i. PE ATAC BAMs downloaded** to `$SCRATCH/atac_pe` (24.5 GB, 165M/153M/229M
      paired reads). Channels building: job 40808663.
- [ ] **W1c. Download paired-end ATAC BAMs → fragment-stratified channels** — unblocks
      the nucleosome-positioning input. `K562/log.sh:9` has the curl commands
      (commented). **Download to `$SCRATCH`, not Oak — Oak is at 97%.**
- [ ] **W1d. Implement nucleosome-resolution (binned) profile target** — code only, no
      compute. Touches `multimodal_bpnet.py` and `train_multimodal_bpnet.py`, both shared
      with the p300 models, so it needs the same backward-compatibility care and
      regression test as the unstranded change.

**Wave 2 — needs its Wave 1 / gate input**

- [~] **W2a. Re-sweep `count_loss_weight` on the 5′ target** — RUNNING (jobs 40808673/6/8,
      weights 10/100/1000). Expect the optimum to fall far below 1000 now that MNLL
      sees real read counts. Was: (needs G1) — expect it
      to fall a long way from 1000 once MNLL sees real read counts. 2–3 GPU jobs.
- [ ] **W2b. Train with fragment-size ATAC channels** (needs W1c) — 5–10 GPU jobs.
      **Revisit the split first.** The fragment distribution is a clean nucleosomal
      ladder (sub-nucleosomal mode ~42 bp, mono ~205 bp, di ~400 bp, tri ~600 bp; see
      `figures/atac_fragment_lengths.pdf`), but the current two channels capture only
      30% (sub 1–99) and 20% (mono 180–247) of fragments — **50% fall in neither**,
      including both higher-order peaks and the 100–180 bp trough.
      Preferred design: **[all, sub, mono, di]** — keep an all-fragments channel so the
      input is a strict superset of the current flat-ATAC baseline and cannot regress,
      plus a di-nucleosomal channel (~350–450 bp) for higher-order structure. Coverage
      fractions per replicate are in `results/atac_fragment_length_summary.tsv`.
      Architecture support landed 2026-08-25: `MultiModalBPNet(n_acc_channels=N)`; the
      forward pass already slices `X[:, 4:, :]` generically, so only the `acc_conv`
      in-channels needed changing. Still to do: multi-BigWig accessibility input in
      `extract_windows` / `--accessibility-bw`, which currently accepts one track.
- [ ] **W2c. Validate the binned profile target** (needs W1d) — 5–10 GPU jobs, ~4–8 GPU-h.
- [ ] **W2d. Re-check the ±500 vs ±1000 window on the 5′ target** (needs G1) — less
      smearing means less neighbour bleed, so the contamination penalty may differ.

- [ ] **W2e. Residual TRAINING — promoted, high value.** An independently-trained
      sequence model captures almost none of the ATAC residual (r = 0.100): it is
      *redundant* with accessibility, not complementary. So do not expect a sequence
      model to find the complement on its own — train it explicitly on
      `observed − atac_pred`. Needs a per-region count offset in the trainer, applied to
      the count loss only. ~4 GPU-h once implemented.
- [ ] **W2f. GM12878 transfer evaluation** — GM12878 fragment-extended target building
      (job 40808659); ATAC (`2026_0606_GM12878_transferability/data/atac.bw`) and
      elements already exist. Must match the target definition the models were trained
      on, so this runs against the fragment target until the 5′ retrain lands.

**Wave 3 — only once the target and weight are settled, so it is not redone**

- [ ] **W3a. Full retrain grid on the final target + weight** — 15–30 GPU jobs.
- [ ] **W3b. Wider receptive field (`--n-layers 10`)** — best single architecture bet, and
      the most expensive: `in_window` 2114 → 5186, ~3.1× FLOPs/epoch, 2–4 h per job,
      ~30–35 GB RSS. 30–60 GPU-h. Deliberately last of the training items.

**Wave 4**

- [ ] **W4. Motif syntax (SHAP / MoDISCo / FiNeMo)** on the final chosen model. Per F-001,
      expect less sequence signal to attribute here than for p300.

### Revisit with proper statistical power

- [ ] **Re-run the `count_loss_weight` sweep with >=3 folds and >=2 seeds.** The current
      choice of 10 rests on single-fold runs, and measured run-to-run variance is ~0.018
      (same config gave 0.4962 and 0.4785 on separate runs). Only the extremes are
      resolvable at that noise level: 10 clearly beats 1 and 1000, but 3 / 10 / 100 are
      indistinguishable. Weight 10 sits in a flat region so nothing is currently wrong —
      revisit if the weight is ever suspected of mattering, or after the profile target is
      binned (which changes MNLL's magnitude and so moves the optimum).
- [ ] **Quantify run-to-run variance properly.** One config x 5 seeds would give an error
      bar for every comparison in this project. Several conclusions so far rest on
      differences of 0.01-0.05 between single runs, which may or may not clear it.
- [ ] **Make submit scripts refuse to overwrite a completed fold directory** (or add a run
      tag to the path). The 5-prime grid silently destroyed the weight sweep's clw10 run
      because both wrote to `models/sequence5p_hw500_clw10/fold0`.

### Housekeeping

- [ ] Delete the local git tag `backup-pre-msg-rewrite` on Oak once the rewritten history
      is confirmed good.
- [ ] `$OAK/Users/sheth/Data/ENCODE/K562/interm_ENCFF790GFL.se.filtered.sorted.bam/` is a
      **1.3 GB leftover** intermediate from `bam_to_bigWig.sh` (its cleanup did not run).
      Oak is at 97%, so worth removing — owner's call, not deleted automatically.
- [ ] `scripts/0.3.make_training_bw.sh` does not sort its experiment bedGraphs before
      `bedGraphToBigWig`, which requires chrom+start order. It has worked so far because
      `genomecov` follows BAM header order, but that is not guaranteed to match
      `chrom.sizes`. `0.5.make_5prime_bigwigs.sh` sorts explicitly; consider backporting.

### Report (do not build yet)

- [ ] **Write up this analysis path** via `/engreitzlab-report` → `/analysis-report`
      (Quarto → self-contained HTML), once the story settles. Deliberately deferred: the
      narrative has already reversed twice (flanking-window idea rejected; residual metric
      flipping the p300/H3K27ac ranking), so wait for the residual-training and transfer
      results before fixing a storyline. Figures that would carry it, already generated:
      `h3k27ac_signal_vs_distance`, `profile_comparison`, `h3k27ac_model_comparison`,
      `h3k27ac_stratified_eval`, `residual_evaluation`, `h3k27ac_replicate_ceiling`.
