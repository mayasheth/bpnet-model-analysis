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

---

### [2026-08-25] Mycelium's decision-log template and its indexer disagree on heading level

**Category**: gotcha

**What happened**: Wrote three decision entries following `mycelium-upstream/skills/core/templates/decision-log-entry.md`, which specifies `## [YYYY-MM-DD] Title`. `generate_index.py` then reported "decisions.md | 0 entries". Its `extract_entries()` hardcodes `header_prefix = "### "` for both learnings AND decisions, while the decision template uses `## `. The learning template correctly uses `### `.

**Why it matters**: Decision entries written exactly as the upstream template instructs are invisible to `.living/INDEX.md`, which is the file the SessionStart hook surfaces and which future sessions are told to trust. The entries are on disk and greppable, but they silently do not exist as far as the knowledge index is concerned — the worst kind of failure, because nothing errors.

**Resolution**: Use `### [YYYY-MM-DD] Title` for decisions in this repo, matching the indexer rather than the template. Verified: index then reports the correct count. Upstream is a pinned submodule we do not edit, so this is a local convention.

**Tags**: mycelium, tooling, index, decisions, template-mismatch, silent-failure

**mitigation_type**: convention

**structural_mitigation_candidate**: A check in the repo that `grep -c '^### \[' .living/decisions.md` matches the count in `.living/INDEX.md`. Better still, report it upstream — either the template or `extract_entries` should change, since one of them is wrong.

---

### [2026-08-25] Excluding the nucleosome-free center makes the H3K27ac target worse, not better

**Category**: insight

**What happened**: Tested whether counting only the flanking windows (excluding +/-125 to +/-375 bp around the element center) gives a better H3K27ac target than a full symmetric window, on the reasoning that H3K27ac sits on flanking nucleosomes and the center is depleted. Inter-replicate ceiling on the top signal quintile fell monotonically with exclusion width: at +/-500 bp outer, 0.761 (no exclusion) -> 0.747 -> 0.729 -> 0.717. Same direction at +/-1000. Meanwhile the flanking-only count correlates 0.97-0.99 with the full-window count, so the two targets are barely different.

**Why it matters**: The premise looked sound from the meta-profile but was wrong quantitatively. The observed central dip is only ~9% below the shoulders, so the center still carries ~91% of peak signal — excluding it discards reads rather than noise, and fewer reads means more Poisson noise and a lower ceiling. Implementing it properly would have required decoupling the count target from the profile target in bpnetlite's loss, so testing the premise first avoided real work on a dead end.

**Resolution**: Dropped. Full symmetric window retained. `2026_0824_H3K27ac_model/scripts/0.4.flanking_vs_full.py`, results in `results/flanking_vs_full_window.tsv`.

**Tags**: h3k27ac, window-selection, ceiling, negative-result, nucleosome

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: Not a defect to prevent — the general lesson is to phrase target-definition questions as "does this raise the inter-replicate ceiling" and answer them with the replicate data before writing model code. That check is cheap and decisive.

---

### [2026-08-25] The IGV display track was the wrong target all along: use 5' ends

**Category**: gotcha

**What happened**: The H3K27ac target was `Data/share/IGV/ENCSR000AKP_coverage.bw`, built by `bam_to_bigWig.sh -r SINGLE`, which extends 36 bp reads to a fixed **250 bp fragment**. It was raw counts and unnormalized, so it looked valid as a training target — and it is, in the sense that nothing is scaled. But it was written for IGV display. Meanwhile the p300 pipeline (`scripts/0.3.make_training_bw.sh:83`) uses `bedtools genomecov -5 -dz`, i.e. single-base 5' ends.

**Why it matters**: Two concrete costs, both of which were visible as symptoms before the cause was identified. (1) The 250 bp extension smears signal across the nucleosome-free center, flattening the bimodality to a ~9% dip and blurring the boundary between adjacent elements — this is also why the flanking-window idea above looked plausible but failed. (2) It makes the profile loss ill-posed: bpnetlite's MNLL expects multinomial read counts, and fragment-extended coverage inflates per-window totals ~250x with near-identical adjacent positions. That is why profile MNLL sat at ~2800 versus ~500 for p300, which in turn is what made `count_loss_weight` need to be ~1000 instead of 1.

**Resolution**: Rebuilding the target from the source BAMs with `genomecov -5 -dz`, stranded, matching p300 — `scripts/0.5.make_5prime_bigwigs.sh`, merged plus per-replicate. Ceiling and models to be recomputed on the 5' target and compared against the fragment-extended results.

**Tags**: h3k27ac, target-definition, bigwig, fragment-extension, 5-prime, bpnetlite, mnll

**mitigation_type**: convention

**structural_mitigation_candidate**: Recorded in `.living/conventions.md`: a track under `Data/share/IGV/` is a DISPLAY artifact and must not be used as a model target without checking how it was built. Deriving training targets from the BAMs inside the analysis project makes the processing choices explicit and reviewable.

---

### [2026-08-25] 5'-end vs fragment-extended target: identical ceiling on active elements

**Category**: insight

**What happened**: Expected the 5'-end target to have a *lower* inter-replicate ceiling than the 250 bp fragment-extended track, because 5' counting gives far fewer counts per window (mean 17 vs 4271 at +/-500 bp) and so more Poisson noise. It does not. Top-quintile ceiling is 0.760 (5') vs 0.761 (fragment) at +/-500, and 0.785 vs 0.785 at +/-1000 — identical to three decimals in places. The all-elements ceiling *rises* sharply (0.844 vs 0.604), but that is the dead-vs-active contrast getting easier rather than a real gain: fragment extension gives every element hundreds of bleed-in counts, compressing log-scale dynamic range, whereas 5' leaves dead elements genuinely near-zero.

**Why it matters**: The switch to 5' is justified on principle — MNLL sees real multinomial read counts, no bleed across the nucleosome-free centre or into neighbouring elements — but it should NOT be expected to improve accuracy by itself. Predicting otherwise would have set up a false attribution when downstream numbers move.

**Resolution**: Adopt 5' as the target. Expect the benefit in the profile task and in `count_loss_weight` normalising, not in the counts ceiling.

**Tags**: h3k27ac, target-definition, ceiling, 5-prime, fragment-extension, negative-result

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: This is the third time an all-elements metric has been misleading in this project (ATAC-only skill, the ceiling comparison, and now this). The stratified number should be the DEFAULT output of every evaluation script here, with all-elements reported as secondary. `2.2.evaluate_stratified.py` and `2.4.evaluate_residual.py` already do this; the earlier scripts do not.

---

### [2026-08-25] An identical ceiling does not mean identical learnability

**Category**: insight

**What happened**: Predicted that switching from the 250 bp fragment-extended H3K27ac target to 5' ends would not improve accuracy, on the grounds that the inter-replicate ceiling was identical on the top quintile (0.760 vs 0.761 at +/-500). Wrong: sequence-only count Pearson went 0.4073 -> 0.4962 (+22% relative), same mode, same fold, same window.

**Why it matters**: The ceiling measures how REPRODUCIBLY a target can be measured. It says nothing about how PREDICTABLY it can be modelled from sequence. Fragment extension smears each read across 250 bp, so signal bleeds in from neighbouring elements — and that bleed is highly reproducible (both replicates see the same smear) while being unpredictable from the element's own sequence. Removing it raises predictability and leaves reproducibility untouched. Treating the ceiling as a proxy for "how much better can a model get" is therefore wrong in a specific, directional way: it is blind to reproducible-but-unattributable signal.

**Resolution**: Ceiling remains useful for ranking window widths and for sanity-checking that a model has not exceeded the possible. It must not be used to predict whether a target-definition change will help.

**Tags**: h3k27ac, ceiling, learnability, 5-prime, target-definition, prediction-was-wrong

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: When comparing two target definitions, run one cheap single-fold model on each rather than reasoning from the ceiling. One 30-minute job answered what the ceiling could not.

---

### [2026-08-25] count_loss_weight collapsed 1000 -> 10 on the 5' target, confirming the units diagnosis

**Category**: insight

**What happened**: On the fragment-extended target the optimal `count_loss_weight` was 1000 and profile MNLL sat at ~2800. On the 5' target, MNLL is ~119-195 and the optimum moved to 10 (weights 1 and 3 still testing). Best count Pearson by weight on 5': 10 -> 0.4962, 100 -> 0.4666, 1000 -> 0.4641.

**Why it matters**: Confirms MNLL magnitude was a units artifact of fragment extension, not a property of the data — the loss is a count-weighted sum, so inflating every count ~250x inflates it proportionally. It also reverses the earlier conclusion about the profile head: at weight 10 the profile term still outweighs counts roughly 12:1, and that is the BEST setting. On smeared coverage profile dominance was harmful because MNLL was fitting the smear; on real read counts it is helpful, acting as a useful auxiliary task. So "down-weight the profile head" was the right fix for the wrong target, and is not needed once the target is correct.

**Resolution**: Use `count_loss_weight` ~10 for the 5' target pending the 1/3 sweep. Do not carry the 1000 forward.

**Tags**: h3k27ac, loss-weighting, mnll, 5-prime, bpnetlite, profile-head

**mitigation_type**: convention

**structural_mitigation_candidate**: Already noted in `.living/conventions.md`: calibrate the weight from the first-epoch loss magnitudes rather than copying it between targets.

---

### [2026-08-25] The ATAC library is textbook-bimodal, but the 2-channel split covers only half of it

**Category**: gotcha

**What headline**: Fragment-length distribution across all three K562 ATAC replicates is a clean nucleosomal ladder: sub-nucleosomal mode at ~42 bp, mono-nucleosomal peak at ~205 bp, di- at ~400 bp, tri- at ~600 bp, replicates superimposed. Median 201 bp.

**What happened**: An earlier read off a truncated text histogram (first 16 bins only) suggested a short-fragment-dominated library with no mono-nucleosomal population, which would have undermined the fragment-channel idea. The full distribution shows the opposite. Separately, the chosen channels (sub 1-99, mono 180-247) capture only 30% and 20% of fragments — **50% fall in neither**, namely the 100-180 bp trough and everything above 247 bp including both higher-order nucleosomal peaks.

**Why it matters**: Two lessons. Never conclude a distribution's shape from a truncated view of it. And the current channel design silently discards half the data; a di-nucleosomal channel (~350-450 bp) or wider bounds would recover information about higher-order chromatin structure that is plainly present.

**Resolution**: Channels built as specified for now. `results/atac_fragment_length_summary.tsv` and `figures/atac_fragment_lengths.pdf` record the distribution and the coverage fractions so the split can be revisited.

**Tags**: atac, fragment-length, nucleosome, channels, data-inspection

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: `0.10.plot_fragment_lengths.py` now reports `frac_sub`, `frac_mono` and `frac_neither` explicitly, so the coverage gap cannot go unnoticed again.

---

### [2026-08-26] The 5' count_loss_weight optimum is 10, with a genuine interior peak

**Category**: insight

**What happened**: Full sweep on the 5' target, sequence-only, fold 0, +/-500: weight 1 -> 0.4366, 3 -> 0.4838, **10 -> 0.4962**, 100 -> 0.4666, 1000 -> 0.4641. A clean interior optimum, not a plateau or a monotone trend.

**Why it matters**: Confirms the profile head is genuinely useful on a well-posed target rather than merely harmless. Too little counts weight (1) is worse than balanced, and so is too much (100+). On the fragment-extended target the optimum was 1000 — two orders of magnitude higher — purely because MNLL was inflated ~15-24x by smearing. The weight is a property of the TARGET'S UNITS, not of the biology or the architecture, and must be re-swept whenever the target definition changes.

**Resolution**: Use `--count-loss-weight 10` for the 5' target. Re-sweep if the profile target is binned (planned) or the window changes, since both alter MNLL's magnitude.

**Tags**: h3k27ac, loss-weighting, mnll, 5-prime, hyperparameters

**mitigation_type**: convention

**structural_mitigation_candidate**: Already in `.living/conventions.md` — calibrate from first-epoch loss magnitudes. A 3-point sweep spanning two orders of magnitude costs ~90 min and finds it reliably.

---

### [2026-08-26] GM12878 is an EASIER H3K27ac target than K562, which makes the transfer failure worse

**Category**: insight

**What happened**: Computed the GM12878 inter-replicate ceiling to normalize the cross-cell-type transfer numbers, expecting some of the drop to be explained by GM12878 being a harder target. The opposite: GM12878's top-quintile ceiling is **0.8321** at +/-500 versus K562's **0.7601** — GM12878 is more reproducible (mean counts 21.1 vs 17.0 per window, i.e. deeper effective coverage).

**Why it matters**: The transfer drops cannot be excused as target difficulty. Normalized as a fraction of each cell type's own achievable ceiling, sequence-only goes from 47% of ceiling in K562 to just **18%** in GM12878; ATAC-only 71% -> 57%; multimodal 88% -> 63%. Normalizing makes the sequence-only collapse look *worse* than the raw correlations suggested, because it was measured against an easier target. This was the open question blocking F-003 from being read quantitatively, and it resolves against the sequence model.

**Caveat carried forward**: the transfer evaluation used the fragment-extended GM12878 target (matching what the models were trained on) while this ceiling is on the 5' GM12878 target. In K562 the two ceilings were identical to three decimals (0.760 vs 0.761), so the substitution is well-supported, but per-replicate GM12878 fragment tracks were never built and the equivalence is assumed rather than measured in GM12878.

**Resolution**: F-003 updated with normalized figures and promoted to supported. Architecture work should be scored on transfer, not held-out chromosomes.

**Tags**: h3k27ac, gm12878, ceiling, transferability, normalization, prediction-was-wrong

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: Compute the target ceiling in BOTH cell types before interpreting any transfer number. A raw cross-cell-type correlation is uninterpretable without it — the drop could be the model or the target, and here it was neither in the expected direction.

---

### [2026-08-26] The offset model must be trained on the same target as the residual model

**Category**: failure

**What happened**: Ran residual training (sequence model fitting `observed - atac_pred`) against the 5' target, but supplied `models/atac_hw500_clw1000` as the offset model — which was trained on the **fragment-extended** target. Fragment counts are ~250x larger, so the offsets arrived at mean 8.15 in log space while the 5' target sits near 2.9. Training converged (valCountPearson 0.813-0.854 across folds) but the numbers are not interpretable as a residual result.

**Why it matters**: Two separate problems, and the run looked healthy despite both. (1) The offset is not "what accessibility predicts for this target" but "what accessibility predicts for a different target, on a different scale". In log space a 250x factor is mostly a constant shift and correlation is shift-invariant, so nothing errored — but per-element neighbour bleed in the fragment track means the mismatch is more than a constant, so the learned component is residual PLUS a scale correction. (2) There was no ATAC-only baseline on the 5' target to compare against, and 5' all-elements correlations are inherently much higher than fragment ones (ceiling 0.844 vs 0.604), so 0.841 is not comparable to any previously reported figure.

**Resolution**: Retrain all three modes on the 5' target at weight 10 (15 jobs), which produces both the correct ATAC-only baseline and the correct offset model, then redo residual training against that. The first residual run was not wasted — it validated the offset plumbing end to end (offsets computed from the right stats, applied at the right place, negatives correctly lower than peaks at 5.21 vs 8.15) — but its numbers are discarded.

**Tags**: h3k27ac, residual, offset, target-mismatch, methodology, silent-failure

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: `train_multimodal_bpnet.py` should refuse to run when the offset model's target differs from the current one. The offset model directory has no record of which signal BigWig it trained on, so the check requires writing a small `training_target.json` next to `acc_normalization_stats.json` and comparing. That would have turned this into an immediate error instead of a plausible-looking result.

---

### [2026-08-26] Correction: the count_loss_weight "optimum" was within run-to-run noise

**Category**: gotcha

**Supersedes** the 2026-08-26 entry "The 5' count_loss_weight optimum is 10, with a genuine interior peak". That entry overstated the evidence.

**What happened**: The 5-prime grid retrained `sequence` fold 0 at weight 10 — the identical configuration the weight sweep had already run — and got **0.4785** where the sweep got **0.4962**. Same mode, fold, window, weight, target. That is ~0.018 of pure run-to-run variance (GPU nondeterminism plus early-stopping epoch choice; the sampler seed is fixed at 42).

**Why it matters**: The sweep was single-fold, and I called a "clean interior peak" off differences of ~0.01-0.03. With variance at 0.018, only the extremes are separable: weight 10 clearly beats 1 (0.437) and 1000 (0.464), gaps of 0.03-0.06. But **3 (0.484), 10 (0.478-0.496) and 100 (0.467) are not distinguishable from one run each**. Weight 10 is a defensible choice sitting in a flat region, not a measured optimum.

**Also lost a data point to a directory collision**: the sweep and the grid both write to `models/sequence5p_hw500_clw10/fold0`, so the grid overwrote the sweep run. The clw10 log preserved under `results/training_logs/` is the FRAGMENT-target run, not the 5-prime sweep, so the original 0.4962 is unrecoverable. Sweep logs for weights 1/3/100/1000 survive as `sweep5p_clw*.epochlog.tsv` only because the grid did not reuse those paths.

**Resolution**: Keep weight 10. Flagged in `todo/TODOLIST.md` as a point to revisit with a properly powered sweep (multi-fold, multi-seed) if the weight is ever suspected of mattering.

**Tags**: h3k27ac, loss-weighting, variance, reproducibility, overclaim, hyperparameters

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: Two things. (1) Hyperparameter comparisons should run >=3 folds before a winner is declared; a single fold cannot resolve differences below ~0.02 here. (2) The submit scripts should refuse to overwrite a fold directory that already contains a completed `multimodal_bpnet.log`, or include a run tag in the path, so a sweep result cannot be silently destroyed by a later grid.

---

### [2026-08-26] Correction: the 5' target's gain on all-elements is mostly the ceiling rising

**Category**: insight

**Refines** the 2026-08-26 entry "An identical ceiling does not mean identical learnability", which reported sequence-only going 0.4073 -> 0.4962 as a +22% gain from the 5-prime switch.

**What happened**: That comparison was single-fold and on ALL validation elements. Expressed as a fraction of each target's own corrected ceiling, the three modes barely move: sequence 47% -> 50%, ATAC-only 86% -> 85%, multimodal 93% -> 91%. The absolute all-elements numbers rise (sequence 0.410 -> 0.479, ATAC 0.746 -> 0.813, multimodal 0.809 -> 0.871) largely because the all-elements ceiling itself rises from 0.868 to 0.957 (corrected), 5-prime having far greater dynamic range once bleed-in stops filling dead elements.

**Why it matters**: The earlier entry's mechanism — bleed-through is reproducible but not predictable from the element's own sequence — is still sound, and the switch is still correct on principle (MNLL well-posed, no neighbour bleed). But the magnitude of the benefit was over-read from an all-elements single-fold number, which is exactly the metric this project has now been burned by three times. The honest test is the stratified top-quintile comparison, pending the full grid.

**Resolution**: Do not quote "+22% from 5-prime". Await the stratified comparison. Report fraction-of-ceiling alongside raw correlation whenever the target definition changes, since raw correlations are not comparable across targets.

**Tags**: h3k27ac, 5-prime, ceiling, all-elements-artifact, overclaim, stratification

**mitigation_type**: convention

**structural_mitigation_candidate**: Already the standing rule in `.living/conventions.md` (stratify every evaluation). Extend it: when the TARGET definition changes, raw correlations are incomparable and only fraction-of-corrected-ceiling should be compared across targets.

---

### [2026-08-26] The 5' switch is a real but modest gain — biggest for sequence-only

**Category**: insight

**What happened**: Full 5-fold stratified comparison, top quintile, where the two targets have essentially identical ceilings (fragment 0.7605 vs 5-prime 0.7601 raw; 0.9295 vs 0.9293 corrected), so the numbers ARE comparable:

| mode | fragment | 5-prime | delta |
|---|---|---|---|
| sequence | 0.357 | 0.3893 | +0.032 |
| atac-only | 0.543 | 0.5508 | +0.008 |
| multimodal | 0.668 | 0.6861 | +0.018 |

**Why it matters**: Settles two open questions. First, the switch genuinely helps rather than merely moving the ceiling — the all-elements gain was mostly ceiling, but the top-quintile gain is real. Second, the gain is **largest for sequence-only**, which fits the mechanism: fragment bleed-through is signal from neighbouring elements, unattributable to this element's own sequence, so removing it helps the model that has nothing but sequence the most.

Against the measured run-to-run variance of ~0.018 on a single fold, the standard error on a 5-fold pooled estimate is roughly 0.008. So sequence +0.032 is about 4 SE and solid; multimodal +0.018 about 2 SE and probably real; **atac-only +0.008 is about 1 SE and not distinguishable from noise**.

The earlier "+22% from 5-prime" claim remains wrong — the honest figure is **+9% relative** for sequence-only on the stratum that matters.

Critically, this does not move the core conclusion. Sequence's marginal contribution over accessibility is +0.135 on 5-prime versus +0.125 on fragment — essentially unchanged, and still roughly 2.2x smaller than p300's +0.301.

**Resolution**: 5-prime adopted as the target, with an honest account of what it bought. All superseded numbers now re-derived except GM12878 transfer and the p300 comparison's H3K27ac side.

**Tags**: h3k27ac, 5-prime, stratification, ceiling, variance, sequence-vs-accessibility

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: Report the pooled standard error next to every comparison. Roughly SE = single-fold sd / sqrt(n_folds); at ~0.018 sd, differences under ~0.016 across 5 folds are not resolvable and should not be described as improvements.

---

### [2026-08-26] A guard written to prevent a silent failure failed itself, because it was never tested

**Category**: failure

**What happened**: After a fragment-trained offset model was silently used against a 5-prime target (see the 2026-08-26 entry above), I added a guard to `train_multimodal_bpnet.py` that reads the offset model's `training_target.json` and refuses a mismatch. Then submitted 5 residual-training jobs against it. All 5 died within 17-288 seconds: `UnboundLocalError: local variable 'target_record' referenced before assignment`. The guard's comparison block sat early in `main()` among the argument validations, while `target_record` was built ~30 lines later next to `os.makedirs(args.output_dir)`.

**Why it matters**: The guard was written specifically to stop a silent failure, and shipped without a single execution. Syntax checks passed, which is exactly why it felt safe. A guard on a rarely-taken branch has no natural test coverage — the code path only runs when someone passes `--count-offset-model`, which nothing else in the project does. Cost was five queued GPU jobs and a full queue-wait cycle on a saturated partition.

**Resolution**: Moved the `target_record` build above the guard, then tested **both directions** before resubmitting: a fragment offset model is refused with the missing-record message, and the 5-prime offset model prints `Offset model target verified` and proceeds. Both verified by actually running the trainer, not by inspection.

Silver lining worth keeping: failing loudly in 17 seconds is the correct failure mode. The bug this guard exists to prevent trained cleanly for 40 minutes and produced a plausible wrong answer.

**Tags**: guards, validation, testing, silent-failure, methodology, slurm

**mitigation_type**: convention

**structural_mitigation_candidate**: Any new validation or guard must be exercised in both the pass and fail direction before it gates real work — a syntax check proves nothing about a branch that never runs. Cheapest form: invoke the entry point twice with deliberately good and bad inputs and confirm the message, which took ~2 minutes here versus five wasted jobs.

---

### [2026-08-26] A naive str.replace on a 2-character token corrupted a guard into always-fail

**Category**: failure

**What happened**: Patching the completion-marker guard into two evaluators, the script parameterized the model-path variable name with `CHECK.replace("mp", mvar)`. In `2.2.evaluate_stratified.py` the variable is `model_path`, so every occurrence of the substring "mp" became "model_path": `complete` -> `comodel_pathlete`, `preempted` -> `preemodel_pathted`, and critically the **filename literal** `training_complete.json` -> `training_comodel_pathlete.json`. The guard therefore looked for a file that can never exist and refused every fold unconditionally. `2.4.evaluate_residual.py` escaped only because its variable is literally named `mp`, making the replace a no-op there.

**Why it matters**: Both failure directions are bad. The guard would have blocked all evaluation, and the error message it printed was mangled gibberish that obscured the cause. The corruption is invisible to a syntax check — the file parsed fine, because only string contents changed. This is the second guard in two days shipped broken, and the first one also passed a syntax check.

**Resolution**: Repaired the block, then verified: the negative direction refuses with a clean message; the filename literal now appears 7 times consistently across the trainer and both evaluators; the `--allow-incomplete-folds` opt-out is defined in both. Testing is what caught it — the guard's own error message was the tell.

**Tags**: patching, str-replace, guards, testing, silent-failure, tooling

**mitigation_type**: convention

**structural_mitigation_candidate**: Never `str.replace` a short token when parameterizing patch text — use an unambiguous placeholder like `__MODELVAR__`, or better, write the block out per-file rather than templating it. And after any patch that rewrites string literals, grep for the literal across all touched files and confirm the count matches expectation; a syntax check cannot see this class of error.

---

### [2026-08-26] Jobs on the `owners` partition are preempted and restart from scratch

**Category**: gotcha

**What happened**: Residual fold 4 ran 15 epochs, was evicted from `owners`, and was requeued as PENDING with its SLURM start time reset. The trainer has no checkpoint resume, so it restarts from epoch 0 and the 15 epochs are lost. Meanwhile the fold directory still held a best-so-far checkpoint from the partial fit — indistinguishable from a finished model.

**Why it matters**: That partial checkpoint is a silent-failure trap: any evaluator pointed at the directory would have loaded an undertrained model and scored it without complaint. The lab's submit scripts all use `-p owners,gpu` because it schedules faster, so this is the default condition, not an edge case. It also means elapsed job times are not training times, and some earlier "slow" jobs in this project may have been preempted and restarted without my noticing.

**Resolution**: The trainer now writes `training_complete.json` as its final action and both evaluators refuse folds without it. `1.7.backfill_training_complete.py` infers completion from the epoch log for folds trained before the marker existed — validated against ground truth, correctly flagging the preempted fold and passing the four that SLURM confirmed complete.

**Tags**: slurm, owners, preemption, silent-failure, checkpoints, sherlock

**mitigation_type**: structural

**structural_mitigation_candidate**: Shipped — the marker plus the evaluator check. A stronger version would add checkpoint-resume to the trainer so preemption costs minutes rather than a whole run, which matters more as windows and receptive fields grow.

---

### [2026-08-27] The ATAC-only model transfers fine; only sequence collapses

**Category**: insight

**What happened**: The model-free ATAC-vs-H3K27ac correlation (no network involved, just summed ATAC against summed H3K27ac over the same element windows) is **lower in GM12878 than in K562**: top-quintile Pearson 0.409 vs 0.510, a 20% relative drop in how tightly the two assays are coupled. Maya predicted this.

**Why it matters**: It reassigns most of the ATAC-only model's apparent transfer loss. Comparing the model against what is available to it in each cell type:

| | K562 | GM12878 |
|---|---|---|
| model-free ATAC-H3K27ac coupling | 0.510 | 0.409 |
| ATAC-only model | 0.542 | 0.467 |
| model / coupling | 1.06 | 1.14 |

The ATAC-only model does not degrade across cell types at all — it actually extracts *more* relative to the raw coupling in GM12878 (1.14 vs 1.06). Its lower absolute score there is a property of the assay pair in that cell type, not a failure to generalize. Sequence-only, by contrast, falls 0.360 -> 0.146 with no comparable explanation available.

This sharpens F-003 rather than softening it. The earlier framing implied everything transfers somewhat worse; in fact **accessibility transfers essentially perfectly and only the sequence component collapses**, which is the more specific and more damning version of the claim.

**Resolution**: F-003 evidence ledger updated. Any future cross-cell-type comparison should report the model-free coupling for that cell type alongside the model, since the coupling is the denominator that makes accessibility-based numbers comparable.

**Tags**: h3k27ac, atac, transferability, gm12878, model-free-baseline, denominator

**mitigation_type**: convention

**structural_mitigation_candidate**: `scripts/0.12.atac_vs_h3k27ac.py` accumulates one row-set per cell type into a single TSV, so the coupling is cheap to compute for every new cell type before any model is trained. Make it a prerequisite for adding a cell type to the panel.

---

### [2026-08-27] Confidence intervals confirm the main claims and kill one of them

**Category**: insight

**What happened**: Switched every evaluation from a single pooled correlation to mean across the five chromosome-holdout folds with a t-based 95% CI. Results at +/-500 bp, top quintile:

| comparison | mean (95% CI) | verdict |
|---|---|---|
| sequence+ATAC vs ATAC-only, H3K27ac | 0.685 [0.667, 0.704] vs 0.548 [0.497, 0.599] | separated |
| ATAC-only vs sequence-only, H3K27ac | 0.548 [0.497, 0.599] vs 0.380 [0.323, 0.437] | separated |
| p300 sequence margin vs H3K27ac sequence margin | +0.304 vs +0.127, difference 0.177 +/- 0.027 | ~6.6 SE, decisive |

**Why it matters**: The headline claims survive with room to spare, so they were not artifacts of pooling. But the fold spread is large — `sd` is 0.041-0.046 for the sequence and ATAC-only models, giving CI half-widths near 0.05, which is **2-3x the run-to-run variance of 0.018 measured earlier**. Most of the uncertainty is between-fold (which chromosomes are held out), not between-run. That means: (a) differences under ~0.05 in this project are not resolvable by adding seeds, only by adding folds or elements; (b) the earlier `count_loss_weight` sweep, which compared single folds, could not possibly have resolved 3 vs 10 vs 100.

Note the multimodal model has a much tighter spread (sd 0.015) than the single-input models (0.041-0.046) — accessibility plus sequence is more stable across chromosome sets than either alone.

**Resolution**: Mean +/- CI with individual fold points is now the reporting standard; `2.2.evaluate_stratified.py` emits `*_stratified_per_fold.tsv` and `*_stratified_fold_summary.tsv` alongside the pooled table.

**Tags**: statistics, confidence-intervals, cross-validation, variance, reporting-standard

**mitigation_type**: structural

**structural_mitigation_candidate**: Shipped in the evaluator. The remaining gap is that folds are paired (same held-out chromosomes across models), so a paired test would be tighter than the independent-CI comparison used above; the current approach is conservative.

---

### [2026-08-27] Transfer is strongly asymmetric: the K562->GM12878 collapse is largely a property of the target cell type, not of sequence

**Category**: insight

**What happened**: Trained all three modalities on GM12878 with settings identical to the K562 runs and evaluated both directions on the 5-prime target, with each direction's in-cell-type reference. Top-quintile Pearson, mean over 5 folds, and as a percentage of that cell type's corrected ceiling:

| | sequence | ATAC only | seq+ATAC |
|---|---|---|---|
| K562-trained -> K562 | 0.380 (41%) | 0.548 (59%) | 0.685 (74%) |
| K562-trained -> GM12878 | 0.146 (15%) | 0.476 (50%) | 0.519 (54%) |
| GM-trained -> GM12878 | 0.221 (23%) | 0.493 (52%) | 0.579 (61%) |
| GM-trained -> K562 | 0.317 (34%) | 0.539 (58%) | 0.600 (65%) |

Sequence-only retention, transfer relative to that model's own in-cell-type score, both as a fraction of the respective ceiling: **K562->GM12878 38%, GM12878->K562 147%.**

**Why it matters**: F-003 was built on the forward direction alone and read the 38% retention as "the sequence component does not transfer". The reciprocal shows that is not the right reading. GM-trained models do *better* on K562 than in their own cell type — for every modality, and after normalising by each cell type's ceiling. So a large part of what looked like a transfer failure is that **GM12878 is simply a harder cell type to predict and K562 an easier one**, in a way the inter-replicate ceiling does not capture. The ceiling measures how reproducibly the target can be measured; it says nothing about how determined that target is by the available inputs, and those differ between cell types.

Two things do survive, and they are the parts worth keeping:
- **Sequence is weak in absolute terms in every direction** — 15-41% of ceiling, against 50-59% for accessibility alone and 54-74% for both. The core conclusion that accessibility dominates is unchanged.
- **The sequence margin shrinks on transfer in both directions**: +0.138 in-cell K562 -> +0.044 transferred (68% lost), +0.086 in-cell GM -> +0.061 transferred (29% lost). Some of the sequence contribution is genuinely cell-type-specific either way.

Also worth noting against the original worry that K562 might be unrepresentative: K562 is the *stronger* cell type both to train on (in-cell multimodal 74% of ceiling vs 61% for GM12878) and to predict. That does not make it representative — it may be atypically predictable — but it does rule out "K562 is a poor training cell type" as the explanation for weak sequence signal.

**Resolution**: F-003 must be restated. A single transfer direction cannot distinguish a model that fails to generalise from a target cell type that is harder to predict; both directions plus both in-cell-type references are the minimum for the claim.

**Tags**: h3k27ac, transferability, gm12878, k562, asymmetry, reciprocal, prediction-was-wrong, ceiling

**mitigation_type**: convention

**structural_mitigation_candidate**: Any cross-cell-type claim in this project requires four evaluations, not two: both transfer directions and both in-cell-type references, all on the same target processing. `scripts/2.6.submit_reciprocal.sh` runs exactly that set and should be the template when cell types are added.

### [2026-08-29] ATAC-only is the right rule, and it cuts the panel from five new cell types to one

**Category**: constraint

**What happened**: The multi-cell-type panel was planned around DNase, on the reasoning that
only 3 cell types had ATAC while 7 had DNase, and that swapping the accessibility bigwig
needed no architecture change. Maya rejected this: the project should consider ATAC only.
Eight panel tracks had already been built (5 DNase coverage + 3 H3K27ac 5' sets); all 23
bigwigs validated clean, and the 5 DNase files are now unused.

Scanning **all 559 released ENCODE ATAC-seq experiments** (matching cell types locally
rather than querying biosample term names one at a time, since a term-name miss returns
404 and is indistinguishable from "assay absent"):

| Cell type | ENCODE ATAC experiments | Usable now |
|---|---|---|
| K562 | 64 | yes — already trained |
| GM12878 | 2 | yes — already trained |
| HCT116 | 17 | input yes, target blocked |
| TeloHAEC (+IL1b/TNFa/noVEGF) | n/a, GEO | yes — ATAC BAMs already on Oak |
| H1, H9, Jurkat, THP-1 | **0** | no ENCODE ATAC under any biosample name |

**Why it matters**: DNase and ATAC are not interchangeable inputs. DNase I cutting and Tn5
insertion have different sequence biases and different footprint structure, and every
existing model here is ATAC-trained. A DNase-input panel would have confounded assay with
cell type in exactly the comparison the panel exists to make — cross-cell-type transfer —
and the confound would have been invisible in the results, appearing as a cell-type effect.

The cost is real and worth stating plainly: ATAC-only takes the expansion from five new
cell types (H1, H9, Jurkat, THP-1, HCT116) to **one** (TeloHAEC), plus three perturbation
conditions of that same line. HCT116 is recoverable on the input side but its H3K27ac
experiment ENCSR661KMA is `run_type = "se, pe"`, so under the one-(experiment, processing)-
pair-per-sample rule its replicates are not interchangeable and no inter-replicate ceiling
is computable — and every cross-cell-type number in this analysis is ceiling-normalised.

TeloHAEC is not a consolation prize: it is endothelial rather than a blood cancer line, so
it is a more distant lineage from K562/GM12878 than H1/H9 would have been, and its three
cytokine conditions give a within-cell-type perturbation axis that tests a different and
arguably more relevant kind of generalisation. Its elements are ATAC-derived, which under
the ATAC-only rule now makes it the *consistent* one and the DNase-derived element sets the
odd ones out.

**Resolution**: Panel restricted to ATAC. TeloHAEC is the one new cell type available
without new data processing, and its ctrl row pools GSE210489 + GSE210491, so it must be
split by experiment before use. H1/H9/Jurkat/THP-1 would require GEO/SRA fastqs through
`Data/scripts/sra_paired_fastq_to_bam.sh` — a separate decision, not a track rebuild.

**Tags**: atac, dnase, panel, generalization, encode, data-availability, confounding, scope

**mitigation_type**: convention

**structural_mitigation_candidate**: Accessibility input assay is part of the sample
definition, not a swappable file path. Any cell type entering the panel must supply ATAC;
record the assay alongside the bigwig path in `config/panel_bam_paths.tsv` so a DNase track
cannot be substituted silently.

### [2026-08-29] TeloHAEC's recorded blocker was the wrong blocker — the real one is a second cell line in the directory

**Category**: correction

**What happened**: The cell-type inventory carried `ceiling_computable = no` for all four
TeloHAEC rows, with the caveat "row pools 2 experiments (GSE210489 + GSE210491); split
before use". Under the ATAC-only rule TeloHAEC became the single available new cell type,
so the blocker was worth resolving rather than working around. Querying GEO for the series
of every GSM in `config/sample_metadata.tsv` shows the premise was wrong:

- **GSE210489 is the ATAC series; GSE210491 is the H3K27ac series**, both under SuperSeries
  GSE210523. One series per assay — not two experiments of the same assay pooled into one
  row. Every (condition, assay) group is a single experiment with a clean replicate set:
  ctrl H3K27ac n=4, ctrl ATAC n=6, and n=2 H3K27ac / n=3 ATAC for each of IL1b, TNFa and
  no-VEGF. **The inter-replicate ceiling is computable for all four conditions.**

The genuine defect was not recorded anywhere:

- **`TeloHAEC_ctrl/ATAC/` contains 9 BAMs, of which 3 are `cell_line = Eahy926`** —
  EA.hy926, a different endothelial line — not TeloHAEC replicates:
  `SRR20809434, SRR20809435, SRR20809436` (GSM6431138, GSM6431133, GSM6431132).
  They sit in the same directory under the same `sample_name`, so any glob of
  `TeloHAEC_ctrl/ATAC/*.bam` silently merges two cell lines into one accessibility track.
  The other three conditions and every H3K27ac group are pure TeloHAEC.

**Why it matters**: Two failure modes, opposite in direction. The recorded caveat would have
cost a real cell type — the only one ATAC-only leaves — for a reason that does not exist.
The unrecorded one would have produced a TeloHAEC ATAC input that is one-third a different
line, which no downstream check would catch: the track would be well-formed, the ceiling
would compute, and the number would simply be wrong. Directory layout encoded the condition
but not the cell line, and `sample_name` said `TeloHAEC_ctrl` for all nine.

**Resolution**: `results/celltype_inventory.tsv` corrected — all four TeloHAEC rows now
`ceiling_computable = yes`, with the Eahy926 accessions named in the caveat. Build TeloHAEC
tracks from an explicit BAM list, never a directory glob.

**Tags**: telohaec, data-provenance, cell-line-contamination, geo, sample-definition,
prediction-was-wrong, ceiling

**mitigation_type**: convention

**structural_mitigation_candidate**: A sample is (experiment, processing) *and cell line*.
Where a metadata table carries a `cell_line` column, group by it before merging BAMs, and
assert the group is single-valued — the directory name is not authoritative.

### [2026-08-30] The residual objective helps only when the input is blind to accessibility

**Category**: insight

**What happened**: Until now exactly ONE residual-objective model existed
(`residual5pFIXED_hw500_clw10`, mode=sequence), so "the residual is sequence-predictable"
and "residual training helps" were confounded. Trained the two missing cells (multimodal and
atac, 5 folds each, jobs 41326647-41326656 plus 41328487/41328488 after two `owners`
preemptions) and evaluated all five models together (job 41334970). Residual *r*, mean over
5 folds:

| input | trained on total signal | trained on residual |
|---|---|---|
| sequence | 0.149 [0.076, 0.222] | **0.459** [0.421, 0.496] |
| sequence + ATAC | **0.551** [0.510, 0.592] | 0.514 [0.473, 0.556] |
| ATAC | (baseline) | **-0.003** [-0.065, 0.060] |

Folds are paired across models, so paired differences are much tighter than the marginal CIs:
- `residual_sequence - sequence5p` = **+0.310** [0.224, 0.395], p = 5e-4
- `residual_multimodal - multimodal5p` = **-0.037** [-0.039, -0.034], p < 1e-4
- `multimodal5p - residual_sequence` = +0.093 [0.072, 0.113], p = 2e-4

**Why it matters**: The residual objective is not a general improvement — it is a fix for a
specific handicap. Where the input cannot see accessibility, forcing the target to be the
residual triples the score. Where the input already includes accessibility, it **costs**
0.037, and does so in every single fold (-0.0357, -0.0351, -0.0356, -0.0388, -0.0393). That
consistency is the point: the marginal CIs overlap and would have been read as "no
difference", but the paired comparison shows a small, highly reproducible harm. A multimodal
model can learn whatever transform of ATAC it needs; nailing it to a fixed offset it cannot
adjust only removes freedom. Joint training stays the best predictor.

The ATAC control is the load-bearing part. Asked to predict an ATAC model's own errors from
the same ATAC input, it scores -0.003 with incremental R^2 of -0.000 and a flat line across
all five |residual| quintiles. Had it come in clearly positive, the 0.459 would have been an
artifact of the metric rather than a finding.

**Correction to an earlier claim**: the residual metric's supposed free artifact channel --
that anything anti-correlated with `atac_pred` scores positive residual *r* -- is largely
CLOSED by construction, and this was overstated when first raised. `obs - atac_pred` is
approximately orthogonal to `atac_pred` by the least-squares residual property. The ATAC
control demonstrates it directly: `r(output, atac_pred) = -0.620 [-1.146, -0.094]`, i.e. the
channel is wide open, and it still scores -0.003. Partialling `atac_pred` out of `sequence5p`
*raised* its score (0.149 -> 0.205) rather than lowering it. The controls were worth running;
the alarm was disproportionate.

**Resolution**: Use joint (multimodal) training for prediction. Reserve residual training for
attribution/motif work, where forcing the model onto accessibility-independent signal is the
objective rather than a cost. Report Fig. 5 rebuilt as the full 2x2 plus control
(`figures/fig5_residual_grid.png`).

**Tags**: h3k27ac, residual, training-objective, multimodal, negative-control, paired-test,
prediction-was-wrong, interpretability

**mitigation_type**: convention

**structural_mitigation_candidate**: Compare models on PAIRED fold differences, not
overlapping marginal CIs -- the -0.037 multimodal effect is invisible in the marginals and
unambiguous when paired. And when a metric is defined against a model's own baseline, train
that baseline's own input on the residual as a negative control; it costs one model and
converts a suspected artifact into a measurement.

### [2026-08-31] The 5' ATAC rebuild works, and the read-length artifact is confirmed quantitatively

**Category**: verification

**What happened**: Built ChromBPNet-convention accessibility tracks (`genomecov -bg -5`,
single-base 5' insertion counts) for K562, GM12878 and the four TeloHAEC conditions as
`*_atac_5p.bw`, alongside the existing full-interval `atac.bw` (job 41326960, 6 tracks).
Validated with a check stronger than readability: `sum(atac.bw) / sum(atac_5p.bw)` must
equal the MEAN READ LENGTH, because full-interval coverage counts read_length bases per read
while 5'-end counting counts exactly one.

| | ratio measured | read length |
|---|---|---|
| K562 | 84.5 | 94-95 bp |
| GM12878 | 85.6 | 94-95 bp |
| TeloHAEC x4 | 30.9-35.5 | 35-36 bp |

The ratio lands at the read length in every case, so the smear being removed is exactly the
one diagnosed. With read length gone, the 5' totals are read COUNTS and the K562-vs-TeloHAEC
gap falls from 13-27x to **5.3-12.3x nuclear depth** -- the ~2.6x that was pure read-length
artifact is eliminated, and what remains is a real depth difference.

**Also found**: TeloHAEC tracks contain chrM (6-7% of reads) while K562/GM12878 do not --
their tagAligns were filtered upstream. This does NOT affect training or evaluation, because
windows are element-centred on nuclear chromosomes and `normalize_accessibility` derives its
statistics from those windows, so chrM never enters. It matters only if a genome-wide
statistic is ever computed from these tracks, and it means raw total-read comparisons across
these cell types are inflated ~6-7% on the TeloHAEC side.

**Why it matters**: The read-length dependence was diagnosed from tagAlign entry widths and
a coverage-scale mismatch; this confirms the diagnosis numerically rather than by argument.
It also bounds what the fix buys: it removes a technical 2.6x, not the whole gap. TeloHAEC
really is a shallower library, and that remains a genuine difference to handle when
transferring models to it.

**Not yet done**: no model uses `atac_5p.bw`. Switching requires retraining every ATAC-input
model as a set (~35 fold-jobs), which is deferred until the residual grids finish on the
current input so those comparisons stay internally consistent.

**Tags**: atac, chrombpnet, 5-prime, validation, read-length, telohaec, depth, chrM

**mitigation_type**: structural

**structural_mitigation_candidate**: When rebuilding a track to remove a suspected artifact,
verify with an arithmetic identity the fix implies (here: old/new ratio == read length), not
just that the new file opens. `scripts/0.21.validate_atac_5prime.py` encodes that check.

### [2026-08-31] The 5' ATAC input helps the ATAC-only model and does nothing for multimodal

**Category**: insight

**What happened**: Retrained every ATAC-input model on the ChromBPNet-convention 5' input
(`genomecov -bg -5`) and compared each against its original, full-interval counterpart. The
training scripts were byte-identical copies apart from the ATAC path and the output dir, and
both arms were scored in a single paired run (`2.2` reads `accessibility_bw` per config
entry), so the only difference is the input definition. Paired over folds, top signal
quintile:

| | ATAC-only | multimodal |
|---|---|---|
| K562 | **+0.032** [0.012, 0.053], p=0.011 | +0.003 [-0.009, 0.015], p=0.51 |
| GM12878 | **+0.023** [0.006, 0.039], p=0.019 | -0.001 [-0.013, 0.011], p=0.84 |
| p300 | +0.019 [-0.054, 0.092], p=0.51 | +0.009 [-0.052, 0.070], p=0.70 |

All-elements shows the same split: ATAC-only +0.016 (K562, p=0.001) and +0.015 (GM12878,
p<0.001); multimodal +0.001 in both, not significant.

**Why it matters**: The prediction going in was "no within-cell-type difference, because read
length is constant there and z-normalisation absorbs a scale change". That was right for
multimodal and **wrong for ATAC-only**. The missing step: within a cell type the smear is
not only a scale change. Scale is constant, but a 94-95 bp interval still BLURS positional
information, and that blur costs a model that has nothing else to go on.

The split is the informative part. Sharpening accessibility helps only the model that
depends on accessibility alone. The multimodal model already has sequence, which supplies
the fine positional detail the smeared ATAC lacked, so a sharper input is redundant with
what it already had. This is the same shape as the residual-grid result: joint training
already extracts what is available, and improving one input does not add to it.

**p300 is underpowered, not null.** 11,412 elements pooled against 57,144 for H3K27ac (2,283
vs 11,429 in the top quintile), giving CIs roughly 5x wider. Its point estimates trend the
same direction; it cannot resolve them.

**Log caveat**: the `ceiling` and fraction-of-ceiling values printed in
`log/accs5p_compare.*.txt` come from `2.2`'s DEFAULT `--ceiling`
(`results/replicate_ceiling_by_window.tsv`), which is the superseded fragment-target K562
ceiling. They are stale for K562 and the wrong cell type for GM12878 and p300. They are NOT
written to the saved TSVs, and no reported number depends on them, but do not quote them.

**Resolution**: The 5' switch is worth keeping — it is the field convention, it is
read-length independent, and it measurably helps the accessibility-only model. Its practical
effect on the multimodal model, which is the one used for prediction, is nil within a cell
type. The case for it remains strongest ACROSS cell types (TeloHAEC 36 bp vs K562/GM12878
94-95 bp), which this comparison does not test.

**Tags**: atac, chrombpnet, 5-prime, accessibility, paired-test, multimodal, p300,
prediction-was-wrong, input-definition

**mitigation_type**: convention

**structural_mitigation_candidate**: When a preprocessing change is claimed to be
"just a scale change that normalisation absorbs", check whether it also changes RESOLUTION.
Scale and blur are separable, normalisation only handles the first, and only models with a
single input feel the second.

### [2026-09-01] The residual-objective result reproduces in GM12878, so it is the objective and not K562

**Category**: insight

**What happened**: Repeated the full residual grid in GM12878 (Maya's suggestion, on the
grounds that it is the harder problem). All 15 folds trained; evaluated with the K562
evaluator copied and altered only in its four cell-type paths, asserted free of K562 paths,
so any difference is data rather than method.

Residual *r*, mean over 5 folds:

| input | trained on total | trained on residual | | K562 equivalents |
|---|---|---|---|---|
| sequence | 0.133 [0.090, 0.176] | **0.372** [0.334, 0.410] | | 0.149 -> 0.459 |
| sequence + ATAC | **0.449** [0.414, 0.485] | 0.414 [0.378, 0.451] | | 0.551 -> 0.514 |
| ATAC (control) | — | **0.010** [-0.019, 0.040] | | -0.003 |

Paired over folds, both cell types, all significant:

| difference | K562 | GM12878 |
|---|---|---|
| residual_sequence - sequence5p | +0.310, p=5e-4 | +0.239, p<1e-4 |
| residual_multimodal - multimodal5p | **-0.0369**, p<1e-4 | **-0.0350**, p=4e-4 |
| multimodal5p - residual_sequence | +0.093, p=2e-4 | +0.077, p<1e-4 |

**Why it matters**: The multimodal COST replicates to within 0.002 across two cell types with
different targets, different accessibility libraries and different element sets — -0.0369 vs
-0.0350, and -0.029 to -0.043 in every one of the ten folds. A number that stable across
that much variation is a property of the training objective, not of a dataset. Before this,
the K562 result was one cell type and the tight ATAC-H3K27ac coupling there (0.510 vs
GM12878's 0.409) was a live alternative explanation. It is now ruled out.

The negative control passes in both (CI includes zero), so neither result is an artifact of
the residual metric.

**The one difference is scale, not direction.** GM12878's sequence gain is smaller (+0.239 vs
+0.310) and every absolute level is lower (residual-trained sequence 0.372 vs 0.459;
multimodal 0.449 vs 0.551). That is consistent with F-003: GM12878 is the harder cell type to
predict. Ordering and sign of every effect are unchanged.

**Resolution**: Settled. Use joint (multimodal) training for prediction in both cell types —
residual training costs ~0.035 residual *r* when accessibility is already an input, reliably.
Reserve residual training for attribution, where forcing the model onto
accessibility-independent signal is the goal. The finding generalises across cell types and
should not need re-testing per cell type.

**Tags**: h3k27ac, residual, training-objective, gm12878, k562, replication, negative-control,
paired-test, generalization

**mitigation_type**: none

**structural_mitigation_candidate**: A result measured in one cell type here has twice been
partly a property of that cell type (F-003's transfer asymmetry; the ATAC-H3K27ac coupling
difference). Replicating in a second cell type before generalising is cheap relative to how
often it has changed the reading — 15 fold-jobs and one evaluation in this case.

### [2026-09-01] Residual and multimodal transfer equivalently; the significant difference is in the stratum we distrust

**Category**: insight

**What happened**: Tested cross-cell-type transfer of residual-objective vs total-target
models, both directions, in two variants that differ in where the offset comes from:

- **upper bound** — offset from the TARGET cell type's own ATAC model. Optimistic: that model
  is trained on target H3K27ac, which the deployment scenario lacks.
- **deployable** — offset from the SOURCE ATAC model, transferred. Only target ATAC is used,
  which is the real application (Maya: "we have ATAC-seq data for the target cell type but
  not H3K27ac").

Paired within fold, residual-multimodal minus multimodal:

| | all elements | top quintile |
|---|---|---|
| transfer K562->GM (local offset) | **+0.0051**, p=0.009 | +0.012, p=0.17 |
| deploy K562->GM | **-0.0054**, p=0.012 | -0.008, p=0.33 |
| transfer GM->K562 (local offset) | **+0.0111**, p=0.003 | +0.005, p=0.65 |
| deploy GM->K562 | **-0.0052**, p=0.022 | -0.002, p=0.76 |

**Why it matters**: The all-elements column tells a tidy story — the residual design wins
when the accessibility baseline is fitted locally and loses when it must be transferred, sign
flipping, all four significant. The top-quintile column says none of it is resolvable. The
all-elements stratum is the one this project has three times recorded as misleading, because
it is dominated by the dead-vs-active contrast that accessibility already resolves. **The
honest conclusion is that residual and multimodal transfer equivalently**, and the tidy story
was an artifact of reading the wrong stratum. Recording this because the temptation to lead
with the significant number was strong and nearly won.

What survives on the top quintile, deployable variant:

| | ATAC-only | multimodal | residual MM | sequence-only |
|---|---|---|---|---|
| K562 -> GM12878 | 0.477 | 0.520 | 0.512 | 0.147 |
| GM12878 -> K562 | 0.541 | 0.602 | 0.600 | 0.317 |

- Sequence genuinely adds on transfer: every sequence-containing model beats transferred
  ATAC-only by 0.04-0.06.
- Transferring the ATAC baseline rather than fitting it locally costs only ~0.017
  (0.494 -> 0.477), so the accessibility component is largely portable.
- Sequence-only transfer remains useless; accessibility must be measured in the target.

**Resolution**: For deployment into a cell type with ATAC but no H3K27ac, transfer the
MULTIMODAL model — not because it is better, but because it is equal on the stratum that
matters and simpler, with no offset model to ship and version alongside it. The residual
design's apparent transfer advantage requires a locally fitted accessibility baseline, which
requires the target H3K27ac that the scenario lacks.

**Tags**: h3k27ac, transferability, residual, deployment, gm12878, k562, stratification,
paired-test, prediction-was-wrong

**mitigation_type**: convention

**structural_mitigation_candidate**: When a difference is significant on all elements and not
on the top quintile, the default reading is that it lives in the dead-vs-active contrast, not
in the biology of interest. Report the top quintile first and the all-elements number second,
even when — especially when — the all-elements number is the significant one.

### [2026-09-01] A generated numbers manifest caught three wrong figures that proofreading missed

**Category**: tooling

**What happened**: Added `scripts/3.7.build_numbers_manifest.py`, which derives
`outputs/numbers.json` from the result TSVs at build time (896 values, 892 derived);
`render_report.py` then checks every number in the prose against it. On the first run it
flagged three claims that no result file supported:

| prose said | actual | source |
|---|---|---|
| centre carries ~91% of peak signal | **94%** | `center_over_peak = 0.945` |
| 49.8% of fragments in neither channel | **50.0%** | 30.3% + 19.7% sum to 50.0 |
| contamination rises to 42% | **41.5%** | `window_tradeoff.tsv` |

Two of the three contradicted a figure legend **in the same section** of the report —
Fig 2's legend already said 0.944 while the prose said 91%, and Fig 3's said 41.5% while the
prose said 42%.

**Why it matters**: Those two survived several deliberate editing passes, because
proofreading compares prose against prose. Only a check against the data catches them. The
manifest registers rounded variants alongside exact values, since prose legitimately writes
0.415 as "42%"; a genuine shift in an underlying number still matches nothing at any
precision and fails.

**Also**: `3.6` now caches per-element vectors keyed on input paths AND mtimes, taking a
replot from 3m33s to 3.5s. It previously saved only summary correlations, so changing a bar
colour cost a full recompute plus queue time. Keying on mtime means a rebuilt track
invalidates the cache instead of silently serving stale vectors.

**Resolution**: Re-run `3.7` after any pipeline change, then re-render; the linter reports
"Lint clean" only when every quoted number traces to a file.

**Tags**: tooling, reproducibility, report, numbers-manifest, caching, drift

**mitigation_type**: structural

**structural_mitigation_candidate**: Shipped. Generalises to any analysis with a Quarto
report — derive the manifest from result files rather than hand-registering values, or it
becomes another thing that drifts.

### [2026-09-01] Measure a coordinate shift by correlation instead of trusting the documented one

**Category**: verification

**What happened**: The fragment-size ATAC channels have to come from the paired-end BAMs,
because fragment length lives in TLEN and the per-read tagAligns discard it. But
`atac_5p.bw` was built from tagAligns that were already Tn5-shifted, while the BAMs are
unshifted alignments. The conventional answer is +4 on the plus strand and -5 on the minus
strand, and ENCODE's pipeline documents that. Rather than assume it,
`0.22.calibrate_tn5_shift.sh` dumps chr20 5' ends by strand from the BAMs and scores every
(plus_delta, minus_delta) in -8..+8 by Pearson r against `atac_5p.bw`.

| | plus | minus | r |
|---|---|---|---|
| best | +4 | -5 | **1.0000** |
| runner-up | +3 | -4 | 0.7705 |
| no shift | 0 | 0 | 0.3148 |

**Why it matters**: Two things fall out that an assumption could not have given. The offset
is confirmed against the actual file rather than against documentation, and r = 1.0000
rather than 0.99-something proves the BAM and tagAlign read sets are *identical* -- same
filtering, same duplicates removed. That is what licenses treating the flat channel as the
exact sum of the four length bins, which is the whole basis for calling the stratified input
a strict superset that cannot regress. Had filtering differed, `all` and `sub+mono+di+poly`
would disagree and the superset argument would quietly fail.

**Also**: a wrong shift here would not have crashed anything. The stratified channels would
have sat a few bp out of phase with the flat channel, and the model would have seen a phase
artifact that looks like nucleosome positioning.

**Tags**: atac, tn5, coordinates, validation, fragment-length, chrombpnet, superset

**mitigation_type**: structural

**structural_mitigation_candidate**: When deriving a new track from a different source file
than an existing one, calibrate the coordinate convention against the existing track over a
grid of candidate offsets before building genome-wide. It costs one chromosome of IO, it
fails loudly when the sources disagree, and the correlation at the argmax doubles as a check
that the two sources contain the same reads.

### [2026-09-02] Byte-identity is the wrong regression invariant for GPU inference

**Category**: process

**What happened**: `2.15` was generalised to score models with different receptive fields.
Because it is shared with the transfer and residual-grid results, the wide evaluation was
gated behind a regression: re-score an existing all-8-layer config and require the per-fold
table to be byte-identical to the stored one. The gate fired and blocked the run. The
differences were all last-digit at 4 dp -- 0.3284 vs 0.3283, 0.2625 vs 0.2626 -- and the
region counts matched exactly in all five folds. The config turned out to have exactly one
distinct accessibility input, so the new code path collapsed to the old one: a single
extraction pass and a centre-crop that is a no-op at in_window 2114. The code was
functionally identical, and cuDNN convolution is not bit-reproducible across GPU models, so
re-scoring on a different node moves the 4th decimal place on its own.

**Why it matters**: The gate worked in the sense of stopping the pipeline, but it stopped it
for a reason unrelated to the change, and a `diff` gate that cries wolf gets deleted rather
than fixed. The fix is to split the invariant by what is actually deterministic:

| | invariant | why |
|---|---|---|
| fold set, config labels, `n` per fold | EXACT | pure numpy over file contents; this is exactly what a geometry or region-set bug breaks |
| every metric | within 1e-3 | absorbs GPU float noise (~1e-4 observed) while staying far below the between-fold sd of 0.041-0.046 |

`2.18.compare_perfold_tables.py` implements this. Re-run: largest metric difference 1.0e-4,
region counts identical, PASS.

**Also**: the same generalisation fixed a latent bug. `2.15` used the BASELINE entry's
`accessibility_bw` for every entry and ignored per-entry values. The `accs5p_*` configs
already carry differing accessibility per entry, so any input-definition comparison routed
through `2.15` would have silently scored both arms on the same input. Entries can now
override it, with an assertion that distinct inputs yield identical valid-region masks.

**Tags**: testing, regression, gpu, determinism, tolerance, evaluation, latent-bug

**mitigation_type**: structural

**structural_mitigation_candidate**: When gating a refactor of an inference pipeline, assert
exact equality only on quantities computed deterministically (region counts, labels, shapes)
and compare model outputs within a tolerance chosen from the noise floor of the science, not
from float precision. State the tolerance and its justification in the comparator itself so
the next person does not tighten it back to zero.

### [2026-09-02] A 4.2 kb receptive field buys dead-vs-active separation and nothing within active elements

> **SUPERSEDED 2026-09-03 — this entry generalised from the sequence-only arm and the
> multimodal arm reverses it.** See the correction entry at the end of this file. The
> sequence-only numbers below are correct; the conclusion drawn from them is not.

**Category**: result

**What happened**: Trained `n_layers = 10` (trimming 2093, in-window 5186, ~4.2 kb receptive
field) against the matched 8-layer models (~1.1 kb), sequence mode, 5 folds, K562 and
GM12878. Scored on the intersection of valid regions so the receptive field is not
confounded with which regions near chromosome ends were scorable.

| | all elements | top quintile |
|---|---|---|
| K562 narrow / wide | 0.494 / 0.537 | 0.381 / 0.375 |
| K562 paired difference | **+0.0430** [+0.0245, +0.0616] p=0.0030 | −0.0063 [−0.0318, +0.0191] p=0.53 |
| GM12878 narrow / wide | 0.475 / 0.520 | 0.221 / 0.232 |
| GM12878 paired difference | **+0.0446** [+0.0222, +0.0670] p=0.0053 | +0.0101 [−0.0307, +0.0508] p=0.53 |

The all-elements gain replicates almost exactly across cell types (+0.043, +0.045) and is
significant in both. The top-quintile effect is null in both and FLIPS SIGN between them,
which is what a null looks like. The CIs exclude a top-quintile gain beyond about ±0.03, so
this is not a power problem.

**Interpretation**: 4.2 kb of context carries broad domain information. That is informative
about whether a region is active at all, and uninformative about how active an already-active
enhancer is. On the project reporting standard the wider receptive field does not help.

**Why it matters**: This is the FOURTH time an all-elements metric has pointed the opposite
way from the top quintile in this project, and the first time I generated the misleading
read myself. Ten folds of validation count Pearson showed +0.031 to +0.064, all positive,
and I reported it as promising. Validation count Pearson is computed over ALL validation
elements, so it is an all-elements metric; it reproduced as the +0.043 all-elements number
and evaporated under stratification. The gain was real, the attribution was wrong.

**Resolution**: Do not spend the 30-60 GPU-h on a wider receptive field for the graded
activity task. Revisit only if the application turns out to be on/off classification, where
+0.043 all-elements is a genuine improvement.

**Tags**: architecture, receptive-field, stratification, top-quintile, replication, k562, gm12878

**mitigation_type**: process

**structural_mitigation_candidate**: Never report a training-log metric as a preliminary
result. Validation count Pearson is unstratified by construction, so on this project it
carries the exact bias the reporting standard exists to remove. The only quotable number
comes from the stratified per-fold evaluator, even when the training logs look unanimous.


### [2026-09-03] The receptive-field conclusion inverted when the second input mode finished

**Category**: reversal

**What happened**: I reported the wide-receptive-field experiment as resolved and negative
on the basis of the sequence-only arm, which finished first: +0.043 and +0.045 on all
elements, null on the top quintile with the sign flipping between cell types. I wrote "do
not pursue" into the TODO, the report and a learning entry, and predicted the multimodal arm
would show even less because ATAC already supplies long-range context. The multimodal arm
shows the opposite.

| top quintile, paired within fold | K562 | GM12878 |
|---|---|---|
| sequence only | −0.006 (p=0.53) | +0.010 (p=0.53) |
| sequence + ATAC | **+0.027 (p=0.006)** | **+0.014 (p=0.025)** |

All five folds rise in both cell types for multimodal. The accessibility residual also rises,
0.502 → 0.547 in K562 and 0.397 → 0.469 in GM12878.

**Why it matters**: Two distinct errors, and the second is the more general one.

1. I generalised from one arm of a two-arm experiment while the other was still running, and
   picked the arm whose mechanism I had already reasoned about. The stated reason for
   expecting a small multimodal effect — accessibility already carries long-range context —
   is precisely the reason the effect is LARGE there: the wider window lets the model read
   accessibility over a wider neighbourhood. The same sentence supports both predictions,
   which means it was not evidence.
2. The project standing rule is "report the top quintile first". It caught the sequence arm's
   all-element artifact and it is what makes the multimodal result credible, but it does not
   protect against generalising across input modes. Fig. 10 now carries both arms for exactly
   this reason.

**Resolution**: The wide receptive field is worth adopting for the multimodal model, which is
the deployed one. An ATAC-only wide arm was submitted to test whether sequence contributes to
the gain at all; if it reproduces +0.027 the extra context is purely accessibility
neighbourhood.

**Also corrected**: the same report claimed the 1 bp profile head was "largely fitting
Poisson noise". First real numbers give top-quintile `profile_pearson` of 0.114 (ATAC), 0.144
(sequence) and 0.171 (multimodal) against an achievable ~0.21, so the head captures most of a
small reproducible signal. The case for binning is that it raises the ceiling, not that the
head fails. Note the two numbers are not exactly comparable — the ceiling excludes elements
flat in either replicate while `profile_pearson` scores them zero — so the ratio is
approximate and the model side is deflated.

**Tags**: reversal, receptive-field, architecture, multimodal, stratification, premature-conclusion, profile-head

**mitigation_type**: process

**structural_mitigation_candidate**: Do not record a conclusion for a multi-arm experiment
until every arm has finished, even when the finished arms agree with each other and with the
mechanism you expected. State the pending arms explicitly in the interim write-up, and put
all arms in the figure so a single-arm reading is not available to the next reader.


### [2026-09-03] Both architecture changes that work act on the accessibility input

**Category**: result

**What happened**: Two independent architecture changes were scored on the top signal
quintile in K562, paired within fold.

| change | top quintile | paired difference |
|---|---|---|
| receptive field 1.1 kb → 4.2 kb, sequence only | 0.381 → 0.375 | −0.006 (p=0.53) |
| receptive field 1.1 kb → 4.2 kb, sequence + ATAC | 0.690 → 0.717 | **+0.027 (p=0.006)** |
| flat ATAC → 5 fragment-size channels, sequence + ATAC | 0.690 → 0.703 | **+0.0135 (p=0.0034)** |

The receptive-field gain replicates in GM12878 (+0.014, p=0.025). The fragment-channel input
is a strict superset of the flat one — the four length bins sum to the flat 5′ track with
zero discrepancy across 545,661,218 insertions and zero largest single-base difference over
2 Mb — so its gain is added information rather than a changed input definition.

**Why it matters**: Everything that has moved the top quintile acts on the accessibility
side. Wider context helps only the model holding ATAC; enriching ATAC itself helps; widening
the window for sequence alone does nothing. Against that, the residual work showed the
sequence component is real but small, and the profile ceiling showed base-resolution shape is
nearly unreproducible. The consistent picture is that headroom on this target lives in how
accessibility is represented, and the sequence component is close to its practical limit at
this scale.

**Open and cheap**: whether the two gains are additive. Both plausibly read the same
neighbourhood structure — a wider window sees accessibility further out, and fragment length
sees nucleosome occupancy locally — so +0.027 and +0.0135 may not sum. One grid of
`n_layers` 10 × fragment channels answers it.

**Also**: the wide-vs-narrow comparison carries no region-coverage confound. K562 per-fold
element counts are identical at in-window 2114 and 5186 (10,143 / 14,392 / 12,801 / 12,021 /
9,298), so the intersection that the evaluator takes costs nothing. The earlier report text
implying the wider window drops chromosome-end elements was wrong and has been corrected.

**Tags**: architecture, accessibility, fragment-length, receptive-field, top-quintile, k562, gm12878

**mitigation_type**: none

**structural_mitigation_candidate**: When several architecture changes target the same input
stream, test the combination before adopting them additively — separately significant gains
on the same underlying signal do not compose.


### [2026-09-03] Test-time reverse-complement averaging is free, always positive, and carries its own control

**Category**: result

**What happened**: Training already augments with reverse complements, but nothing averaged
predictions over a sequence and its reverse complement at inference. Adding it (`2.15
--rc-average`, off by default) improves every model in the K562 residual grid, paired within
fold, and the size of the gain tracks reliance on sequence.

| model | top quintile Δ | p |
|---|---|---|
| sequence only | **+0.0162** [+0.0097, +0.0227] | 0.0023 |
| sequence + ATAC | **+0.0081** [+0.0057, +0.0105] | 0.0007 |
| sequence + ATAC, residual | +0.0073 [+0.0012, +0.0133] | 0.0295 |
| sequence, residual | +0.0054 [+0.0032, +0.0076] | 0.0023 |
| ATAC only | +0.0016 [−0.0004, +0.0036] | 0.086 |
| ATAC only, residual | +0.0018 [−0.0004, +0.0040] | 0.090 |

Profile shape improves for every model too, +0.0012 to +0.0060, all significant.

**Why it matters**: The ordering is a built-in control, and it is what makes the result
believable rather than just significant. RC averaging cancels residual strand asymmetry in
the sequence branch, so a model whose only input is unstranded accessibility coverage should
gain nothing — and the two ATAC-only rows are the only non-significant ones. A uniform gain
across all six models would have suggested something generic (numerical smoothing, an
evaluation artefact); the gradient points at the intended mechanism.

Cost is one extra forward pass and no retraining, which makes this the cheapest improvement
found on this project. For comparison, the fragment-size channels gained +0.0135 on the same
elements and cost a full retrain plus building four genome-wide tracks.

**Correctness check that mattered**: the transform was validated against a deliberately
wrong version that reverses positions without complementing the one-hot channels. The correct
transform agrees with the forward pass at r = 0.9918 and the wrong one at 0.9485, so the
complement demonstrably matters. Without that control a silently-wrong flip would have looked
like a small honest gain.

**Left off by default**, so every existing number stays comparable. Adopt it deliberately as
a set if the report's headline numbers are ever re-scored.

**Tags**: inference, reverse-complement, augmentation, free-win, negative-control, k562

**mitigation_type**: none

**structural_mitigation_candidate**: When adding a symmetry-based inference trick, include a
model that cannot benefit from the symmetry. A result that is uniform across every model is
usually measuring something other than the mechanism claimed.


### [2026-09-04] Reusing a Snakemake run's outputs needs the per-run files too, not just the obvious ones

**Category**: process

**What happened**: The nine ABC predicted-activity arms must share one candidate-region set,
so each arm's `Peaks/` was copied from a completed run and the dry run was gated to abort if
any region-calling rule got scheduled. It aborted three times.

1. `shutil.copytree` preserves mtimes, so the copies carried the original run's July dates
   and read as stale against September bigwigs.
2. Stamping every copied file with one timestamp fixed that but not the real requirement:
   Snakemake needs an output **strictly newer** than its input, and equal mtimes still
   scheduled the rule. Files had to be staggered along the chain
   `narrowPeak -> .sorted -> Counts.bed -> candidateRegions.bed`.
3. Still aborted. `snakemake -n -r` gave the actual reason, which was not an mtime at all:
   `sort_narrowpeaks` takes a SECOND input, the per-run
   `results/<run>/tmp/<chrom sizes>.bed` written by `generate_chrom_sizes_bed_file`. A fresh
   results directory lacks it, so Snakemake generated it and everything downstream re-ran
   with reason *"Input files updated by another job"*.

Copying `tmp/` as well, stamped earlier than the Peaks chain, produced a clean dry run:
9 x create_neighborhoods / create_predictions / filter_predictions / QC and nothing else.

**Why it matters**: I diagnosed twice by assumption and once by evidence, and only the
evidence was right. `-r` was available from the first abort and would have shown the
cascade immediately; instead two rounds went into an mtime theory that was true but not
sufficient. The general shape: when reusing part of a pipeline's output tree, the
per-biosample directories are the visible dependency and the per-RUN files are the invisible
one, and a missing per-run file re-runs everything through the "updated by another job" rule
rather than announcing itself.

**Why the gate was worth having**: without it the run would have completed, silently calling
regions per arm. All nine arms would then have had different candidate regions, and every
cross-arm comparison — the entire point of the experiment — would have been meaningless
while looking perfectly healthy.

**Tags**: snakemake, abc, reuse, mtime, dry-run, gating, diagnosis

**mitigation_type**: process

**structural_mitigation_candidate**: Run `snakemake -n -r` and read the stated reason before
theorising about why a rule is scheduled. When reusing outputs, diff the fresh results tree
against the source run at the top level, not just the directories you meant to copy.

### [2026-09-04] Nine ABC arms on one region set: what it took, and the two invariants that mattered

**Category**: process

**What happened**: Getting ABC to run nine predicted-activity arms took seven driver
submissions. Every failure but the last was caused by the fix for the previous one.

| attempt | failure | cause |
|---|---|---|
| 1 | region-calling rules scheduled | `copytree` preserves mtimes, so copied Peaks carried July dates |
| 2 | same | one uniform timestamp is not enough; Snakemake needs an output STRICTLY newer than its input |
| 3 | same | `sort_narrowpeaks` also takes the per-RUN `results/<run>/tmp/<chrom sizes>.bed`, absent in a fresh results dir, which cascades via input-files-updated-by-another-job |
| 4 | `CreateCondaEnvironmentException` | `--use-conda` needs mamba on PATH for the shell Snakemake spawns |
| 5 | `ModuleNotFoundError: pyranges` | the PATH fix for (4) put the wrapper env's python ahead of the rule env's |
| 6 | `LockException` | a driver cancelled in (5) was still alive and holding the lock |
| 7 | `IncompleteFilesException` | the dry-run gate omitted `--profile`, so it lacked the profile's `rerun-incomplete` while the real run had it |

Also found mid-run: two child jobs from a cancelled driver were still RUNNING 22 minutes
later. SLURM does not kill a driver's children, so they were about to write the same output
files as the new driver's children.

**Why it matters**: two invariants did all the protective work.

- **The dry-run gate.** It aborted whenever region-calling rules appeared. Without it the run
  would have completed while calling candidate regions independently per arm, and all nine
  arms would have had different region sets -- every cross-arm comparison meaningless, with
  nothing in the output looking wrong. Final verification: identical 153,546 rows and
  identical region md5 `7d5995ce` across all nine.
- **The guarded unlock.** Refusing to `--unlock` while another driver was RUNNING stopped two
  Snakemake instances from writing the same files.

**Why it took seven tries**: I diagnosed by assumption and resubmitted, twice, when
`snakemake -n -r` prints the reason directly. A ten-minute validation job
(`scripts/shimtest.sh`) that checked mamba resolution and python precedence ON A COMPUTE NODE
settled in seconds what four submissions had probed one hypothesis at a time.

**Tags**: abc, snakemake, gating, diagnosis, orphaned-jobs, mtime, conda, process

**mitigation_type**: process

**structural_mitigation_candidate**: For a cluster pipeline, write a seconds-long job that
asserts the environment preconditions on the node type the real work runs on, and run it
before the pipeline rather than after the third failure. Gate correctness invariants with a
dry run that uses the SAME flags as the real run. After cancelling a driver, always check for
orphaned children -- the scheduler leaves them running.

### [2026-09-05] Why better H3K27ac does not help ABC: rank displacement of the functional elements

**Category**: result

**What happened**: Observed H3K27ac beats ATAC alone by +0.062 AUPRC on the CRISPR benchmark
(0.519 vs 0.457). Predicted H3K27ac recovers none of it on precision-at-min-sensitivity, and
its AUPRC gain sits entirely inside the floor's CI. `4.9.diagnose_abc_gap.py` separates three
candidate explanations.

**Accuracy is not the problem.** Predicted vs observed H3K27ac per ABC region, Spearman:

| stratum | n | Spearman |
|---|---|---|
| all regions | 153,545 | 0.794 |
| CRISPR-tested regions | 3,016 | **0.824** |
| regions of regulated pairs | 339 | 0.663 |

Agreement on the tested regions is HIGHER than genome-wide, so the model is not failing where
the benchmark looks.

**The mechanism is rank displacement, and I had the wrong mechanism first.** I predicted
dynamic-range compression. That cannot be it: `normalized_atac` is identical across arms and
`normalized_h3k27ac` is qnorm'd onto the same K562 reference, so scale is mathematically
irrelevant — both arms draw values from the same distribution. The only way regulated
elements end up with lower activity is if the model RANKS THEM LOWER, which it does:

| stratum | p99/p50 observed | p99/p50 predicted |
|---|---|---|
| all regions | 32.4 | 25.4 |
| CRISPR-tested | 11.0 | 9.2 |
| regulated pairs | **6.00** | **4.03** |

Top-decile separation on regulated pairs: 6.46 observed against 4.23 predicted. The loss is
concentrated on exactly the elements CRISPR calls functional.

**The threshold behaviour is the surprise.** All three arms catch the same positives:
observed 322 of 429, predicted 325, ATAC-only 323, with 10-15 discordant either way. So the
observed-H3K27ac advantage is about RANKING (suppressing negatives), not about detecting more
true positives. Any effort aimed at "catching more enhancers" is aimed at the wrong quantity.

**One sharp failure class.** The 12 pairs observed catches and predicted misses sit at
elements with observed H3K27ac.RPM median 22.79 -- about 39x the genome-wide median of 0.59 --
where the model predicts 0.71. Very strong real signal treated as near-background.

**Also ruled out**: the train/score element mismatch. Training on the ATAC-derived set ABC
scores against changes nothing: -0.0011 [-0.0147, +0.0125], p = 0.83.

**Why it matters**: the quantity this project has been optimising is the wrong one. Overall
top-quintile r improved from 0.690 to 0.724 across the architecture work, and none of it
reached the endpoint, because the endpoint depends on ORDERING AMONG HIGH-SIGNAL ELEMENTS and
on separating them from negatives. A loss that weights high-signal elements, or a ranking
objective, targets that directly; more trunk capacity does not.

**Tags**: abc, crispr-benchmark, diagnosis, qnorm, ranking, negative-result, mechanism

**mitigation_type**: none

**structural_mitigation_candidate**: When a model improves on its training metric but not on
the downstream endpoint, diagnose which property of the prediction the endpoint consumes
before trying more architectures. Here the endpoint consumes rank order among high-signal
elements; the training metric averaged over everything. Also: check whether a hypothesised
mechanism is even reachable given the pipeline -- scale compression was impossible to blame
because qnorm removes scale by construction.

### [2026-09-05] Sequence already knows which elements are not acetylated; the additive trunk ignores it

**Category**: result

**What happened**: Every model that sees ATAC over-predicts H3K27ac at accessible-but-
unacetylated elements. Measuring each model's elevation in that stratum against ITS OWN
genome-wide median (the painted-bigwig scale differs from observed RPM, so cross-model raw
comparisons are meaningless):

| | in stratum | own median | fold |
|---|---|---|---|
| observed H3K27ac | 0.59 | 0.59 | **1.0x** |
| sequence only | 0.28 | 0.18 | **1.6x** |
| ATAC only | 1.00 | 0.12 | **8.3x** |
| sequence + ATAC | 0.80 | 0.11 | **7.3x** |

Those elements are GC 0.593, CpG o/e 0.551, 2.4x promoter-enriched, ATAC/K27ac ratio 13.25
against 2.20 typical -- a CpG-island / CTCF phenotype. The sequence-only model nearly gets
them right; the accessibility-using models do not, and multimodal inherits ATAC's error
almost intact.

**Why it matters**: the information is already in sequence and the architecture cannot use it.
`MultiModalBPNet` concatenates a 64-filter sequence branch with an 8-filter accessibility
branch and the trunk mixes them ADDITIVELY, so accessibility contributes equally everywhere.
There is no mechanism for sequence to veto accessibility at a CpG island. Adding
sequence-derived features (GC, CpG density, motif scores) would supply information the branch
already has -- the fix is an interaction, not another additive input.

**The loss is the other half.** log1pMSE is symmetric, so predicting 0.80 where truth is 0.59
costs ~0.015 squared log error while missing an 18.6 enhancer costs ~6.1 -- roughly 400x less
sensitive to inventing signal than to missing it. A gate has no gradient pressure to close
while that holds, so the architecture change and the loss change must go in together or the
experiment will wrongly conclude gating does not work.

**A mistake worth recording**: my diagnostic script printed an automated verdict of "sequence
cannot see it", the opposite of the truth. It applied an accessibility-residualised error to
the sequence-only model, which has no accessibility input -- charging it with failing to
predict a deviation from something it cannot observe. The raw numbers in the same output said
the opposite. Do not let a script print a conclusion that has not been checked against its
own inputs, and do not apply a control designed for one arm to an arm it does not fit.

**Tags**: architecture, gating, loss-design, cpg-island, ctcf, over-prediction, self-correction

**mitigation_type**: structural

**structural_mitigation_candidate**: Before adding a feature to fix a model error, check
whether a model that ALREADY has that feature makes the same error. If the sequence-only arm
gets these elements right, the problem is how the branches combine, not what they contain --
and no amount of extra input will fix it.
