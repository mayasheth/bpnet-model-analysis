# Learnings

Append-only log of gotchas, surprises, and insights.

**Entry template:** copy from `skills/core/templates/learning-entry.md` (includes Category, What happened, Why it matters, Resolution, Tags fields). The `**Tags**:` line is consumed by `generate_index.py --summary-heuristic` to build the cluster summary in INDEX.md — use them.

### [2026-08-24] count_loss_weight must be calibrated to the actual loss magnitudes, not copied

**Category**: gotcha

**What happened**: The first H3K27ac model reached validation count Pearson of only 0.21 against a ceiling of ~0.87. bpnetlite's loss is `profile_loss + count_loss_weight * count_loss`. In this run profile MNLL was ~2800 while count MSE was ~3, so the p300 default of `count_loss_weight=1` (and even 10) gave the counts head under 1% of the gradient. Raising it to 1000 roughly doubled count Pearson to 0.41; 10000 gave no further gain.

**Why it matters**: The weight is not transferable between targets. Profile MNLL scales with both the output window width and the magnitude of the signal, so the same weight means something completely different for a 1000 bp vs 2000 bp window, or for a fragment-extended histone track vs a 5'-end TF track. Copying it from a working p300 config silently produces a model that barely optimizes the objective you care about.

**Resolution**: Read the reported `Training MNLL` and `Training Count MSE` columns from the first epochs and set the weight so `weight * count_loss` is at least comparable to `profile_loss`. Default in `2026_0824_H3K27ac_model/scripts/1.1.submit_training.sh` is now 1000, and the output directory includes `_clw{weight}` so sweeps do not overwrite each other.

**Tags**: bpnetlite, loss-weighting, hyperparameters, h3k27ac, training

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: A check at the end of epoch 1 in `train_multimodal_bpnet.py` that warns when `count_loss_weight * count_loss < 0.1 * profile_loss`, i.e. the counts head is receiving under ~10% of the gradient. Cheap, and would have caught this in 3 minutes instead of a full training run.

---

### [2026-08-24] Inter-replicate r is not a model performance ceiling without two corrections

**Category**: insight

**What happened**: Used `r(rep1, rep2)` on per-element counts as the upper bound on model performance. Models immediately exceeded it — multimodal H3K27ac scored 1.07x the supposed ceiling, which is impossible and revealed the bound was wrong.

**Why it matters**: An invalid ceiling silently invalidates every "fraction of achievable performance" claim, and in this case would have been reported as if models were near-perfect.

**Resolution**: Two corrections, both needed. (1) Spearman-Brown: the training target is the 2-replicate MERGE, which is less noisy than one replicate, so `reliability = 2r / (1 + r)`. (2) A model predicts the EXPECTED signal, not a noisy draw of it, and the maximum correlation between a perfect predictor of a true value and a noisy observation of it is `sqrt(reliability)`. Final ceiling `sqrt(2r / (1 + r))`: at +/-500 bp this turns a raw 0.604 into 0.868. Implemented and commented in `2026_0824_H3K27ac_model/scripts/2.1.aggregate_results.py`.

**Tags**: statistics, ceiling, reliability, spearman-brown, evaluation, replicates

**mitigation_type**: structural

**structural_mitigation_candidate**: Shipped — the correction lives in `2.1.aggregate_results.py` and `2.2.evaluate_stratified.py`, and `frac_of_ceiling > 1` is now the tripwire that reveals a mis-specified bound.

---

### [2026-08-24] bpnetlite's count target sums ALL channels, not one strand

**Category**: gotcha

**What happened**: While extending the stratified evaluator to the stranded p300 models, found that my observed-count computation used `signals[:, 0, :].sum()` — channel 0 only. bpnetlite's `_mixture_loss` does `y.reshape(n, -1).sum(-1)`, flattening channels and positions first, giving "a single count loss across all tracks".

**Why it matters**: Identical for an unstranded 1-channel target, so it passed unnoticed on H3K27ac. On a stranded 2-channel target it halves the observed counts, which would have made the p300 baseline look artificially bad in exactly the head-to-head comparison it was built for. A silent factor-of-two in a comparison is worse than a crash.

**Resolution**: Changed to `sigs.sum(axis=(1, 2))`, correct for both. Fixed before any p300 numbers were produced.

**Tags**: bpnetlite, stranded, count-target, evaluation, p300, off-by-factor

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: A test asserting that the evaluator's observed counts equal `bpnetlite.losses._mixture_loss`'s internal `y_` for both a 1-channel and a 2-channel target. That pins evaluation to training by construction rather than by reading the library source.

---

### [2026-08-24] Three Sherlock/SLURM traps that cost a job each

**Category**: gotcha

**What happened**: Three separate failures in one session, all mechanical rather than scientific. (1) `set -u` plus `"${ARR[@]}"` on an EMPTY array is an unbound-variable error in bash < 4.4 — killed a training job in 5 seconds on precisely the `sequence` path that skips accessibility. (2) Python block-buffers stdout when it is a file, so a running job showed an empty log for 35 minutes and looked hung. (3) A long-running python process started over `ssh` on the login node was killed when the calling command was backgrounded locally, exiting 0 with no output files written.

**Why it matters**: Each looks like a real failure and burns a debugging cycle. The buffering one is the most expensive: it makes every other failure undiagnosable while the job is alive.

**Resolution**: (1) `${ARR[@]+"${ARR[@]}"}`. (2) `export PYTHONUNBUFFERED=1` in every submit script. (3) Anything heavier than a few seconds goes through `sbatch`, never a login-node process over ssh.

**Tags**: slurm, bash, sherlock, buffering, submit-scripts, tooling

**mitigation_type**: convention

**structural_mitigation_candidate**: Promoted to `.living/conventions.md` as mandatory items for new submit scripts in this repo.

---

### [2026-08-24] Peak count scales the negative pool, which can OOM by 50x

**Category**: failure

**What happened**: `train_multimodal_bpnet.py` hardcoded `max_negs = len(train_peaks) * 10`. Fine for ~12k p300 peaks, but H3K27ac trains on ~120k DNase candidate elements, giving 1.2M negative windows — roughly 55 GB of extracted arrays against a 120 GB request, before counting the peaks themselves.

**Why it matters**: The 10x rule is invisible until the peak set gets large, and it fails during extraction after tens of minutes of work rather than immediately.

**Resolution**: Added `--max-negatives`, defaulting to the old behaviour so p300 runs are unchanged, set to 50000 for H3K27ac. The sampler only draws `--negative-ratio` (0.1) of each batch from negatives, so a large distinct pool was never doing much work. Observed peak RSS afterwards: 13.8 GB.

**Tags**: memory, oom, negatives, training, scaling

**mitigation_type**: structural

**structural_mitigation_candidate**: Shipped as `--max-negatives`. A stronger version would estimate bytes from `n_regions * in_window * channels * 4` before extraction and refuse to start if it exceeds the SLURM allocation.
