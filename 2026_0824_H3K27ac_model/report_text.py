"""New prose for the three-way report split. Old sections are pulled by title."""

DATE = "2026-09-05"

# =====================================================================
# REPORT 1 -- data characterisation
# =====================================================================
R1_MAP = {8: 1, 2: 2, 3: 3, 11: 4, 9: 5, 1: 0, 4: 0, 5: 0, 6: 0, 7: 0, 10: 0,
          12: 0, 13: 0, 14: 0}

R1_HEAD = """## Headline figure

![](figures/fig8_coupling_across_celltypes.png)

**Figure 1 | How much of H3K27ac accessibility can explain is a property of the cell type, and the ordering inverts between strata.** **Headline figure.**

<details>
<summary>Full legend</summary>

See the coupling section below for the complete legend. The two points that make this the
headline: the model-free ATAC-H3K27ac correlation on the elements carrying the mark is
0.51 in K562 against 0.33-0.41 in every other cell type measured, so an accessibility-based
number is not comparable across cell types without this denominator; and on all elements the
ordering inverts, with TeloHAEC highest, because that stratum is dominated by the
dead-versus-active contrast every accessible-element set shares. Source:
`results/atac_vs_h3k27ac_by_celltype.tsv`, `results/coupling_panel_recomputed.tsv`.
</details>

## Summary

- **The target is not where the input is.** ATAC peaks on the element centre; H3K27ac is bimodal with shoulders at +/-250 bp and is still at 28% of maximum at +/-2 kb, so its extent exceeds the models' ~1.1 kb receptive field by an unmeasured margin (Fig. 2).
- **A 1 kb counting window is the widest choice with zero neighbour contamination**, and the achievable ceiling is nearly saturated there: top-quintile inter-replicate *r* rises only 0.760 to 0.798 from +/-500 to +/-2000 bp while the fraction of windows containing another element's centre rises from 0% to 41.5% (Fig. 3).
- **There is almost nothing at base resolution to predict.** A perfect model caps at *r* = 0.21 (K562) and 0.18 (GM12878) on the top quintile at 1 bp, rising to 0.72 and 0.70 at 50 bp binning (Fig. 4). Any profile-head number must be read against that cap, not against 1.
- **Accessibility-H3K27ac coupling varies two-fold across cell types and the strata disagree about the ordering** (Fig. 1), which is why every evaluation in this project is stratified and why a new cell type gets its coupling measured before any model is trained in it.
- **The ATAC library supports fragment-size stratification** -- a clean nucleosomal ladder with a sub-nucleosomal mode at 42 bp and a mono-nucleosomal peak near 205 bp (Fig. 5).

## Goals

This report records what the H3K27ac and ATAC data *are*, and what is predictable from them
in principle, separately from any model. It exists so that the other two reports can cite a
ceiling, a coupling or a track definition without re-deriving it, and so that a number
quoted anywhere in the project can be checked against the measurement it came from. Every
ceiling in the project lives here.

Three companion documents:

- **Report 2, evaluation methodology** -- how a comparison in this project is scored and what
  each metric can and cannot support. Standing cautions, each with the incident that
  motivated it.
- **Report 3, design decisions** -- the workbench. Which architecture and input changes were
  tried, what each did to the reporting metric, and the verdict.
- `reference/DATA_INVENTORY.md` -- the full assay inventory across all nine cell types,
  including which assays exist where and which files must not be used.
"""

R1_CONVENTIONS = """## Accessibility track conventions, and the read-length confound they hide

Two conventions are in play for the accessibility input and they are different quantities.

**Full-interval coverage** (`bedtools genomecov -bg` over the tagAlign interval) counts
`read_length` bases per read, so the track scales with read length and smears each Tn5
insertion across ~95 bases in our K562 and GM12878 libraries.

**5' insertion counts** (`genomecov -bg -5`) count exactly one base per read, which is the
ChromBPNet convention (`chrombpnet/helpers/preprocessing/reads_to_bigwig.py`: `-bg -5`,
unstranded, ATAC shift `4-plus_shift, -4-minus_shift`). Our tagAligns are already
Tn5-shifted, so the correct build reduces to adding `-5`; applying the shift again would
double-shift.

**The difference is a constant within a cell type and a confound across them.** Verified by
an arithmetic identity: `sum(full-interval) / sum(5')` must equal the mean read length.
Measured 84.5 and 85.6 for K562 and GM12878 (94-95 bp reads) and 30.9-35.5 for the four
TeloHAEC conditions (35-36 bp). So K562-versus-GM12878 comparisons were never affected,
but TeloHAEC's full-interval coverage is a different quantity from theirs. With read length
removed, the K562-versus-TeloHAEC gap falls from 13-27x to 5.3-12.3x nuclear depth, which
is genuine sequencing depth rather than an artifact of the convention.

The 5' tracks are read-length independent by construction and are the correct choice the
moment a third cell type enters. What the switch does to model performance is a design
decision, not a property of the data, and is reported in Report 3.

## Element derivation differs between cell types, which is a standing confound

Documented in full in `reference/ELEMENT_DERIVATION.md`. The K562 and GM12878 candidate
element sets are DNase-derived; TeloHAEC's are ATAC-derived. TeloHAEC is therefore the one
cell type whose elements match its own input assay, and the two cell types carrying most of
this project's results do not.

An ATAC-derived K562 set exists (`reference/K562_ATAC_candidate_elements.narrowPeak`, same
rE2G model type as TeloHAEC's), so the confound is testable rather than merely noted; the
test and its result are in Report 3. The ABC candidate regions used for the downstream
benchmark are ATAC-derived, from observed K562 ATAC.

## The training negative pool is not GC-matched, despite its filename

`reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed` is ChromBPNet's genome-wide
GC-*annotated* tiling: 3,088,298 bins at 1 kb stride with GC in column 4. That file is the
*input* to GC matching, not its output. Our training code reads only columns 1-3, discarding
the GC column, and the sampler draws uniformly at random from a 50,000-window random
subsample of it. So the training negatives are a uniform random sample of the genome, at mean
GC **0.389** against **0.466-0.593** for candidate elements, and no matching to either cell
type's positive set takes place. The same pool is used for K562 and GM12878, which would be
wrong if matching were happening -- since the positives differ per cell type -- and is moot
because it is not.

Three things follow, two reassuring and one not.

- **Negatives are not mislabelled.** They carry their *true* extracted H3K27ac signal rather
  than a hard zero, so the 1 kb bins that happen to overlap a real candidate element get
  their real signal. They are also restricted to each fold's training chromosomes.
- **No reported number is computed on them.** Evaluation runs on the candidate element set
  only (Report 2), so every correlation in these three reports is an element-only number.
- **But it changes what the model can learn.** A positive/negative contrast separable on GC
  content alone is a known way to turn a BPNet-style sequence branch into a coarse GC
  detector, and a large share of the sequence branch's training signal here is that contrast:
  negatives are 10% of every batch and differ from positives by 0.1 in mean GC. This is a
  live hypothesis for why the sequence branch has almost no dynamic range *within* candidate
  elements while still correlating 0.397 across them (Report 3), and it is the first item on
  that report's open list.

**Separately, the candidate element set itself spans active and inactive elements**, so most
of it carries no H3K27ac whatever the negatives are doing. That is the reason every
evaluation in this project is stratified, and it is independent of the issue above.
"""

R1_METHODS = """## Methods

All ChIP-seq targets are counted as **5' read ends** at single-base resolution
(`bedtools genomecov -5 -dz`, stranded), matching the convention already used for the p300
models in this repository. Accessibility inputs exist in both conventions described above;
which one a given result used is recorded per result in `results/TARGET_PROVENANCE.md`.

Ceilings are reported as the raw inter-replicate *r*. Converting one to a bound on model
performance requires two corrections, applied wherever a fraction-of-ceiling is quoted
anywhere in this project: Spearman-Brown for the fact that the model predicts the
two-replicate *merge* rather than one replicate (`rel = 2r/(1+r)`), then a square root,
since a model predicts the expected signal while a replicate is one noisy realisation of it.
At +/-500 bp this converts a raw 0.760 into a bound of 0.929 on the top quintile, and a raw
0.844 into 0.957 over all elements.

Windows are centred on candidate-element midpoints rather than summits, because the H3K27ac
summit sits on a flanking nucleosome (Fig. 2). The `summit` column of
`K562_DNase_candidate_elements.narrowPeak` equals `width/2`, so element-centred windows
require no change to the extraction code. Elements within `flank` of a contig end are
dropped. Signal is read with pyBigWig and `nan` treated as zero.

### Regenerate

```
Environment: pixi env `multimodal`, pixi.lock sha256:069960d36312
```

```bash
# Target and accessibility tracks
sbatch scripts/0.5.make_5prime_bigwigs.sh            # K562 H3K27ac 5' ends, merged + per-replicate
sbatch scripts/0.7.make_gm12878_5prime.sh            # GM12878 H3K27ac 5' ends
pixi run bash scripts/0.20.make_atac_5prime_bigwigs.sh
pixi run python scripts/0.21.validate_atac_5prime.py   # ratio == mean read length
sbatch scripts/0.22.calibrate_tn5_shift.sh
PLUS_DELTA=4 MINUS_DELTA=-5 sbatch scripts/0.23.make_atac_fragment_5p_channels.sh
pixi run python scripts/0.24.validate_atac_fragment_channels.py

# Characterisation
pixi run python scripts/0.6.compare_profiles.py   ...           # Fig 2
pixi run python scripts/0.3.replicate_ceiling_by_window.py ...  # Fig 3
pixi run python scripts/0.4.flanking_vs_full.py --outer 500 1000 --inner 0 125 250 375
sbatch scripts/0.25.submit.sh                                   # Fig 4
pixi run python scripts/0.10.plot_fragment_lengths.py ...       # Fig 5
pixi run python scripts/0.12.atac_vs_h3k27ac.py --label K562 ...  # Fig 1

# This report
python render_report.py report1_data_characterisation.qmd
```

Per-result provenance, including which processing of the target each result file used, is in
`results/TARGET_PROVENANCE.md`.
"""


# =====================================================================
# REPORT 2 -- evaluation methodology
# =====================================================================
R2_MAP = {1: 1}

R2_BODY = """## Headline figure

![](figures/fig1_three_mode_comparison.png)

**Figure 1 | The same five models, two strata, two different conclusions -- which is why the top quintile leads every comparison in this project.** **Headline figure.**

<details>
<summary>Full legend</summary>

**Figure 1 | Predicting K562 H3K27ac at candidate elements, by input modality.** Pearson *r*
between predicted and observed log counts. **a**, All elements. **b**, Top quintile by
observed H3K27ac -- the elements carrying the mark. Colour encodes input modality throughout
all three reports: **red sequence only, blue ATAC only, purple sequence + ATAC**. Dashed line
is the achievable ceiling; 0.957 all elements, 0.929 top quintile. The panels are not a
rescaling of one another: sequence-only falls from 0.494 to 0.380 while ATAC-only falls from
0.831 to 0.548, so the gap between the two modalities is 0.34 in **a** and 0.17 in **b**.
The multimodal model's fold-to-fold spread (sd 0.015) is a third of either single-input
model's (0.041-0.046). Source: `results/fiveprime_stratified_fold_summary.tsv`.
</details>

## Goals

This report is the standing set of rules for how a comparison in this project is scored, and
what each metric can and cannot support. It is deliberately short and each item is stated as
a caution together with the incident that motivated it, because every one of them was
written after a wrong reading was published internally, not in advance.

It covers: which stratum leads, how differences are tested, what the residual metric means
and the controls it needs, how a transfer evaluation must be *designed* to be interpretable,
and what the downstream ABC/CRISPR benchmark can decide. Report 1 holds the ceilings and
track definitions this report refers to; Report 3 holds the results.

## Report the top signal quintile first, and check every input mode

**The rule.** Every model comparison leads with the top quintile of elements by observed
signal. All-element numbers may be reported alongside but never alone, and a claim supported
only by the all-element stratum is not a finding.

**Why.** Evaluation runs on the candidate element set only -- training negatives are never
scored, so every number in these reports is an element-only number -- and that set spans
active and inactive elements, so most of it carries no H3K27ac. An unstratified correlation
over it is therefore dominated by the dead-versus-active contrast that accessibility already
resolves. The two panels of Fig. 1 are the same models and disagree about the size of every
effect.

**The incidents.** Four misleading readings so far, all of the same shape.

- The wider receptive field gains a significant +0.043 on all elements for the sequence-only
  model and nothing at all on the elements carrying the mark (-0.006, *p* = 0.53). Reported
  from the all-element stratum this reads as the architecture change working.
- The residual-transfer comparison flips sign between strata: residual-multimodal beats
  multimodal by +0.005 and +0.011 with a locally fitted baseline and loses by -0.005 when it
  is transferred, all four significant on all elements and none resolvable on the top
  quintile (*p* = 0.17-0.76).
- Cross-cell-type coupling inverts its ordering between strata: TeloHAEC is highest on all
  elements and lowest on the top quintile (Report 1, Fig.~1).
- The same architecture change points opposite ways in different input modes. In the wide
  receptive-field experiment the sequence and multimodal arms disagree on the top quintile,
  so a conclusion drawn from whichever arm finished first would have inverted.

**The corollary, which keeps being learned the hard way.** When an architecture change is
being tested in several input modes at once -- sequence only, ATAC only, sequence + ATAC --
do not draw a conclusion from whichever mode finishes training first. The receptive-field
question was closed as "resolved, do not pursue" on the strength of the sequence-only models
alone, then had to be reopened when the sequence + ATAC models landed and showed a
significant gain. Wait for the input mode the change is actually supposed to act on.

## Pair within fold; the between-fold spread dominates everything

**The rule.** Report the mean across the five chromosome-holdout folds with a t-based 95%
confidence interval, showing individual folds as points. For any comparison between two
models, difference them **within each fold and then average** -- a paired test.

**Why.** Fold `sd` is 0.041-0.046 for single-input models, giving CI half-widths near 0.05,
against a run-to-run variance of 0.018. Nearly every architecture effect in this project is
between 0.002 and 0.03, which is well inside that spread. Both models in a fold see the same
held-out chromosomes, so differencing within a fold removes the between-fold variation that
otherwise dominates. Unpaired, none of the fragment-channel, RC-averaging or
receptive-field results would be resolvable.

A single pooled correlation over concatenated folds is not an acceptable substitute: it gives
no error bar and can be biased when folds differ in mean or scale.

**The consequence for planning.** A change expected to move the metric by less than ~0.005
cannot be resolved with five folds, and either needs more folds or should not be run. State
the expected effect size before submitting.

## What the residual metric means, and the controls it needs

**Definition.** `residual_pearson = r(observed - atac_pred, model_pred - atac_pred)`, where
`atac_pred` is an **ATAC-only model's held-out prediction**, not the raw ATAC track. So the
residual is what accessibility genuinely cannot explain rather than what a linear function of
one bigwig cannot explain. Reported alongside incremental *R^2* against the observed signal,
and stratified by |true residual| so performance is scored where accessibility fails.

**The artifact to rule out.** Because `true_resid = obs - atac_pred`, anything merely
anti-correlated with `atac_pred` might earn positive residual *r* for free. It largely
cannot: `obs - atac_pred` is near-orthogonal to `atac_pred` by the least-squares residual
property. But the argument is not enough on its own, so three controls run with it.

1. **A negative control model.** An ATAC-input model trained on the residual -- asked to
   predict an ATAC model's own errors from the same ATAC input -- must score zero. It does:
   -0.003 [-0.065, 0.060], incremental *R^2* -0.000, flat across all |residual| quintiles,
   while its output correlates -0.620 with `atac_pred`, so the leak channel is wide open and
   still yields nothing. This control is load-bearing; a clearly positive value here would
   have made the whole residual result an artifact.
2. **Partialling.** Removing `atac_pred` from both sides leaves the residual-trained sequence
   model unchanged at 0.459, against 0.149 to 0.205 for the more entangled
   sequence-on-total model.
3. **Incremental *R^2* against observed signal**, which cannot be faked out of sample.

**A trap for anyone re-running this.** `MultiModalBPNet.forward()` does **not** add the count
offset -- that happens only inside `fit()` when scoring the loss -- so an offset-trained
model's raw output *is* the residual. The first evaluator subtracted `atac_pred` anyway,
double-counting the baseline and correlating a residual against observed signal for
`overall_pearson`. Both bugs produced plausible wrong numbers rather than failing. The fix
sits behind a `"residual": true` config key so earlier results stay reproducible, and
`scripts/2.7.diagnose_residual_offset.py` verifies the semantics empirically before any
scoring: a residual predictor's output must centre near 0 (measured -0.044) and a logcounts
predictor near the observed mean (2.835 against 2.786), with a plain sequence model as a
control proving the test discriminates.

**Residualisation is the wrong control for a stratum defined by accessibility disagreement.**
Later, and separately: the elements where the multimodal model over-predicts are by
construction "accessibility says high, H3K27ac says low", so the observed residual is
strongly negative there and **any** predictor that does not track accessibility downward --
including a constant -- scores as over-predicting. A residualised error column therefore
cannot distinguish "the sequence branch fails to see this" from "the sequence branch has no
dynamic range". It produced exactly that wrong verdict once (`scripts/4.12`), which is now
scored as fold elevation over each arm's own median, carrying the opposite tail as the arm's
own responsiveness scale. Report 3 has the corrected reading.

## Designing a transfer evaluation so that it is interpretable

A raw cross-cell-type correlation cannot distinguish a model that fails to generalise from a
target that is simply harder to predict. Four design requirements follow, and all four came
from readings that had to be withdrawn.

**Run both directions.** A one-directional analysis of K562 to GM12878 reads as a collapse of
the sequence component. The reciprocal shows GM12878-trained models scoring *higher* on K562
than in their own cell type for every modality, so much of the apparent failure is that
GM12878 is the harder cell type.

**Score against the ceiling of the cell type being predicted**, not the training cell type,
using the corrected ceiling from Report 1.

**Measure the model-free coupling in the target before training anything there.** It costs
minutes of CPU and a transfer number cannot be interpreted without it: the ATAC-only
"transfer drop" into GM12878 turned out to be entirely explained by weaker
accessibility-H3K27ac coupling there (Report 1, Fig.~1), not by any failure of the model.

**Include the target's own ATAC-only model as the floor.** The deployment question is not
"does the transferred model score well" but "does it beat what the target cell type's own
ATAC could have given for free". Without that arm a transfer table can look successful while
describing a model nobody should deploy.

**In-cell-type ranking is not the deployment ranking**, so any architecture change intended
for deployment is scored transferred as well as locally. The two rankings do in fact differ;
Report 3 has the matrix.

## ABC and the CRISPR benchmark as the downstream metric

The project's stated goal is a predicted activity track usable in ABC, so the terminal
metric is CRISPR-benchmark AUPRC rather than any correlation.

**Design.** ABC's `activity_base = sqrt(normalized_h3k27ac * normalized_atac)`. Model
predictions are painted into a bigwig as `predicted_counts / width` over the ABC candidate
regions and substituted for the H3K27ac term, with observed ATAC retained. Arms cover the
geometric mean with observed ATAC, predicted H3K27ac as the whole activity term, both
K562-trained and GM12878-trained models, and all three input modalities.

**Two scale questions, one of which turns out not to matter.** The pipeline counts reads from
BAMs while a painted prediction track has arbitrary units, so absolute scale is not
comparable. It does not need to be: `run_qnorm` is **rank-based**, so any monotone rescaling
of the activity term leaves the quantile-normalised value unchanged, and only the *ordering*
of the predictions matters. What does deserve a check is that the qnorm reference itself is
built from observed K562 signal, so a predicted track inherits the observed distribution's
shape -- which is why a qnorm-off arm is worth running and is listed as open.

**Two anchors are mandatory.** A floor of ATAC-only activity and a ceiling of ATAC with
observed H3K27ac. Every predicted arm is read as a position between them, and a predicted arm
that fails to clear the floor has not demonstrated anything, however good its correlation.

**What this benchmark can and cannot resolve.** Its CIs are roughly +/-0.05 wide while the
differences under discussion are ~0.02, so it cannot rank two similar predicted arms. It can
tell whether a predicted track reaches the observed-H3K27ac ceiling, and whether it clears
the floor. Present it as a forest plot with the CIs visible; bars hide the only thing that
matters.

**The mechanism check that goes with it.** A benchmark result of "no better than the floor"
is compatible with two very different causes -- inaccurate predictions, or accurate
predictions whose *dynamic range* is compressed so that ABC's ranking is displaced. These are
distinguished by comparing observed and predicted `activity_base` percentile ratios and
top-decile mean/median ratios on the tested regions, and by counting the pairs each arm
catches that the other misses. Run that diagnostic before concluding anything about the
model from an AUPRC.

## Regression gates, and what tolerance to set

**Do not byte-diff numerical output.** The first regression gate on the per-fold scoring
tables compared bytes and failed on 1e-4 differences from cuDNN nondeterminism.
`scripts/2.18.compare_perfold_tables.py` replaces it: exact on fold, label and *n* and on
reference columns; metrics compared within `--tol 1e-3`, a tolerance taken from the science's
noise floor rather than from what the hardware happens to reproduce. Added columns are
reported and skipped; missing columns fail.

**A regression test that cannot fail is worse than none.** The first version of the
asymmetric-loss test passed while verifying nothing: its over-prediction check had zero
over-predicted elements in the sample, and its gate comparison compared two models whose
random initialisations were never matched, because adding the gate consumed extra RNG draws.
Both now assert their own preconditions -- `3 <= n_over <= N-3`, and weight-matched
initialisation -- so a vacuous pass is impossible.

**Never quote a metric from a training log as a preliminary result.** Validation count Pearson
in the training log is unstratified, so it is the all-element number this report's first rule
warns about. A promising early gain reported from training logs vanished entirely on the top
quintile. Wait for `scripts/2.15.perfold_from_config.py`.

## Conventions

- **Reverse-complement averaging is on by default** in all scoring since 2026-09-03
  (`--no-rc-average` to disable). It is free, positive for every model, and applying it
  everywhere keeps tables comparable; the tables named `rc_*` are the re-scored set.
- **One scoring code path.** In-cell grids, transfer runs and input-definition swaps are all
  driven by the same config JSON through `2.15.perfold_from_config.py`, so a difference
  between two results cannot come from which evaluator ran. Model geometry is read per entry
  from `model.trimming` rather than hardcoded, after an earlier version silently could not
  score the wide models.
- **Report the change, the stratum and the *p* together**, or the number is not usable by
  anyone else.
- **Non-GPU work goes to the lab's `engreitz` partition**; see the repository `CLAUDE.md`.

## Methods

All correlations are on log counts, over element-centred windows of the half-width recorded
per result. Fold assignments come from `reference/hg38_five_folds.json` and are shared across
every model in the project so numbers are comparable. Per-fold tables are written to
`results/*_per_fold.tsv` and summaries to `results/*_fold_summary.tsv`; the paired
differences quoted anywhere in these reports come from the `--pair A B` output of
`2.15.perfold_from_config.py`, which is a t-test on the five within-fold differences.

Numeric claims in these reports are checked against `outputs/numbers.json` by
`render_report.py`, which flags any hand-typed number not registered at its compute site.

### Regenerate

```
Environment: pixi env `multimodal`, pixi.lock sha256:069960d36312
```

```bash
# the one scoring path, driven by a config JSON
pixi run python scripts/2.15.perfold_from_config.py config/<name>_configs.json \\
    <prefix>_ <elements.narrowPeak> --pair A B

# regression gates
pixi run python scripts/2.18.compare_perfold_tables.py results/old_per_fold.tsv \\
    results/new_per_fold.tsv --tol 1e-3
pixi run python scripts/test_asymmetric_loss.py

# residual-metric semantics check, before any residual scoring
pixi run python scripts/2.7.diagnose_residual_offset.py

# this report
python render_report.py report2_evaluation_methodology.qmd
```
"""


# =====================================================================
# REPORT 3 -- design decisions
# =====================================================================
# Old -> new figure numbers. New prose references figures as "Fig.~N" with the
# FINAL number, tilde-protected so the remapper leaves it alone; pulled prose
# uses old numbers and is remapped.
R3_MAP = {14: 1, 1: 2, 4: 3, 5: 4, 10: 5, 12: 6, 13: 7, 6: 8, 7: 9}
R3_FIXUPS = [
    (r"\(Fig\. 2\)", "(Report 1, Fig.~2)"),
    (r"\(Fig\. 9\)", "(Report 1, Fig.~5)"),
]

R3_HEAD = """## Headline figure

![](figures/fig14_crispr_benchmark.png)

**Figure 1 | No predicted-H3K27ac arm clears the ATAC-only floor on the CRISPR benchmark, and the deployment arms sit at or below it.** **Headline figure.** Superseded in part: switching the target to p300 does clear the floor, by +0.055 with a paired-bootstrap interval excluding zero -- see "Predicted p300 works where predicted H3K27ac does not" below. This figure remains the H3K27ac result.

<details>
<summary>Full legend</summary>

**Figure 1 | CRISPR-benchmark AUPRC for eleven predicted-activity arms, K562.** Forest plot
rather than bars because the confidence intervals are the point: they are roughly +/-0.05
wide while the differences under discussion are ~0.02. Vertical lines mark the floor
(ATAC-only activity, 0.457) and the ceiling (ATAC x observed H3K27ac, 0.519). Best predicted
arm is predicted H3K27ac used as the entire activity term from the K562 multimodal model, at
**0.482 [0.437, 0.530]** -- above the floor point estimate, interval overlapping both
anchors. The geometric mean of observed ATAC with predicted H3K27ac scores **0.469**, and the
two GM12878-trained deployment arms score **0.455** and **0.452**, at or below the floor.
Every sequence-only arm is far below distance-to-TSS (0.435): 0.393 for the geometric mean
and 0.275 for predicted-alone. Benchmark is
`EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz` with the scE2G intGENCODEv43 universes,
scored by `CRISPR_comparison`. Source:
`CRISPR_comparison_v3/.../results/2026_0904_predicted_activity/performance_summary.txt`.
</details>

## Summary

- **Predicted p300 clears the benchmark floor in the training cell type; predicted H3K27ac does not.** A K562 sequence + ATAC model predicting p300, substituted into ABC's activity term, beats the ATAC-only floor by **+0.055 AUPRC [+0.034, +0.074]** on a paired bootstrap, sign preserved in 100% of resamples, and is indistinguishable from *measured* H3K27ac (-0.007 [-0.032, +0.017]). No predicted-H3K27ac arm managed this (Fig. 1).
- **The benchmarked transfer direction fails, but transfer itself is strongly asymmetric and the benchmark can only test the failing direction.** GM12878-trained p300 applied to K562 scores +0.009 [-0.004, +0.022] over the floor on the benchmark, and on the correlation metric it lands exactly at the target's ATAC-only floor (0.312 against 0.306, *p*=0.72) with a *negative* accessibility residual (-0.088). But **K562-trained p300 applied to GM12878 retains 0.507 against a 0.608 local and a 0.299 floor** -- **+0.207 [+0.158, +0.257]** over the floor, 83% of the local advantage, residual *r* +0.139 [+0.067, +0.210]. CRISPR data exists only for K562, so the downstream benchmark is *structurally* unable to test the direction that works.
- **The model-quality confound is resolved, and it was not the explanation.** Both p300 models are equally strong at home -- 0.619 (K562) against 0.608 (GM12878), with matched ATAC-only floors of 0.306 and 0.299 -- so the GM12878 model is not weaker. What fails is transferring *out of* GM12878 specifically.
- **Observed p300 is simply a better ABC activity term than observed H3K27ac**, by +0.038 [+0.019, +0.058], independent of any model and replicated across two pipeline runs. That widens the available headroom over the floor from 0.062 to 0.099 and is a reason to prefer p300 as the target wherever it can be measured -- a conclusion about the assay, which stands regardless of the transfer result.
- **For p300 the sequence branch earns its place**, beating the ATAC-only p300 model by +0.051 [+0.032, +0.069]. For H3K27ac the sequence arms were catastrophic (0.393 and 0.275 against a 0.457 floor). This is the clearest evidence yet that target choice, not architecture, was the binding constraint.
- **For H3K27ac the mechanism is compressed dynamic range, not inaccuracy.** Predicted and observed H3K27ac correlate at Spearman 0.79 genome-wide, but the predicted activity term's p99/p50 ratio is 4.03 against 6.00 on the regions of regulated pairs, so ABC's ranking is displaced rather than wrong.
- **In-cell-type architecture gains do not travel.** Every transferred architecture is within 0.016 of the simplest one, and the richer models have the *larger* transfer drops. Transferring the best model buys +0.013 to +0.036 over simply using the target cell type's own ATAC-only model.
- **The failures have a clean signature, but it is not CTCF.** Over-predicted elements are accessible, GC-rich, CpG-rich, promoter-enriched and carry no more H3K27ac than a random element; CTCF-high elements are only 27% of that tail and removing them leaves the rest over-predicted more, not less, so the 2.6x pooled CTCF enrichment is a passenger. GC content correlates with the error 3-4x more strongly than CTCF does. Under-predicted elements are canonical active enhancers (EP300 3.4x, H3K4me1 3.3x, zero H3K27me3).
- **p300 as a second target is dropped.** Where the H3K27ac model fails, observed p300 is elevated only ~2.3x over its own input control and the p300 models predict it at 1.83x -- below that control -- so neither a multi-head model nor stacking can help.

## Goals

This is the workbench. It records every architecture and input change tried against the
H3K27ac target, what each did to the reporting metric, and the verdict, followed by the
downstream ABC/CRISPR result and the diagnosis of why it is negative.

It changes most often of the three reports. Report 1 holds the data characterisation and
every ceiling; Report 2 holds the scoring rules these results are read under, and should be
read first by anyone about to quote a number from here.

**Numbers in the decision table are reverse-complement-averaged**, which has been the
default in all scoring since 2026-09-03. The figures and the per-section prose below show the
single-pass values they were rendered from. The two differ by +0.002 to +0.016 and change no
verdict; the exact gap for every affected result is tabulated in Methods, so the discrepancy
is bounded rather than merely acknowledged. Figures are deliberately left uniformly
single-pass rather than partly migrated: two of the six affected figures have no `rc_*` table
in the format their plot script reads, and a half-migrated figure set would be harder to
describe correctly than a consistent one.

## Decision table

Effect is the paired within-fold change in top-quintile Pearson *r* on K562 H3K27ac unless
stated. Read every row against a between-fold `sd` of 0.041-0.046 (Report 2).

| Change | Effect on the top quintile | Verdict |
|---|---|---|
| **Inputs: sequence only** | 0.397 [0.344, 0.451] absolute | reference |
| **Inputs: ATAC only** | 0.551 [0.500, 0.601] absolute | reference |
| **Inputs: sequence + ATAC** | 0.695 [0.677, 0.713] absolute | **adopted** |
| **Target: p300 instead of H3K27ac** | sequence margin +0.304 against +0.127 | p300 is the better substrate for sequence; H3K27ac is the required target |
| **Objective: residual instead of total** | sequence +0.254, multimodal **-0.018** | adopt for attribution work only; joint training remains the better predictor |
| **Receptive field ~1.1 -> ~4.2 kb** | multimodal **+0.028** (*p*=0.004); GM12878 +0.016 (*p*=0.023); ATAC-only +0.017; sequence-only -0.007 (*p*=0.40) | **adopted for accessibility-input models**; the gain is accessibility neighbourhood, not sequence context |
| **Accessibility: 5 fragment-size channels vs flat** | **+0.016** (*p*=0.002); GM12878 +0.006 (*p*=0.031) | **adopted** |
| **Both: wide + fragments** | +0.034 vs narrow+flat; **+0.007 on top of wide alone** (*p*=0.053) | sub-additive at ~79%; take the receptive field, treat fragments as nearly redundant with it |
| **Inference: reverse-complement averaging** | +0.008 multimodal, +0.016 sequence-only, +0.002 ATAC-only (n.s.) | **adopted as the default everywhere** |
| **Accessibility: 5' insertion counts vs full-interval** | ATAC-only **+0.033** (*p*=0.010); multimodal +0.004 (n.s.) | **adopted**, for the cross-cell-type read-length confound rather than the gain (Report 1) |
| **Elements: ATAC-derived instead of DNase-derived** | -0.001 [-0.015, +0.013] (*p*=0.83) | irrelevant; the derivation confound is real but has no measurable effect |
| **Counting window +/-500 bp** | ceiling nearly saturated, zero neighbour contamination | keep (Report 1, Fig.~3) |
| **`count_loss_weight` = 10** | only the extremes 1 and 1000 are resolvable | keep; sits in a flat region |
| **Profile head at 1 bp** | ceiling 0.21; models reach 0.114-0.171 | not informative, but not failing either; removing it is untested |
| **Architecture choice on transfer** | no arm beats narrow+flat by a resolvable margin (best +0.016, *p*=0.066) | in-cell gains do not travel |
| **Predicted H3K27ac as the ABC activity term** | AUPRC 0.482 against a floor of 0.457 and a ceiling of 0.519 | does not clear the floor; deployment arms at or below it |
| **Predicted p300 as the ABC activity term, in-cell type** | **+0.055 [+0.034, +0.074]** over the floor, paired bootstrap | real; matches measured H3K27ac |
| **Predicted p300 as the ABC activity term, GM12878 -> K562** | **+0.009 [-0.004, +0.022]** over the floor | **null** in this direction |
| **Predicted p300, correlation metric, K562 -> GM12878** | **+0.207 [+0.158, +0.257]** over the target's ATAC-only model; 83% of local | **transfers**; the benchmark cannot test this direction |
| **Observed p300 instead of observed H3K27ac as the activity term** | **+0.038 [+0.019, +0.058]** | use p300 where it is measured; model-independent |
| **p300 as a second head alongside H3K27ac** | predicted p300 1.83x elevation against a 2.10x input control at the H3K27ac failure elements | **dropped**; p300 as the *primary* target is the version that works |
| **Sequence-gated accessibility, and asymmetric count loss** | in-cell null; **-0.013 (*p*=0.039) transferred when combined** | **closed** |
| **GC-matching the training negatives** | +0.001 (*p*=0.76) in-cell; -0.002 (*p*=0.33) transferred | **closed**; the GC-shortcut explanation is refuted |
"""

R3_INPUTS = """## Inputs: accessibility carries most of it, and sequence adds a real but small complement

**Q:** How much of H3K27ac at candidate elements is predictable from sequence, how much from
chromatin accessibility, and how much does sequence add beyond accessibility?

**A:** On the elements carrying the mark, ATAC alone reaches 0.551 [0.500, 0.601] and
sequence alone 0.397 [0.344, 0.451]; together 0.695 [0.677, 0.713], so sequence adds +0.144
over accessibility and accessibility adds +0.298 over sequence.

![](figures/fig1_three_mode_comparison.png)

<details>
<summary>Figure 2 legend</summary>

**Figure 2 | Predicting K562 H3K27ac at candidate elements, by input modality.** Pearson *r*
between predicted and observed log counts. **a**, All elements. **b**, Top quintile by
observed H3K27ac -- the elements carrying the mark. Colour encodes input modality throughout
all three reports: **red sequence only, blue ATAC only, purple sequence + ATAC**. Dashed line
is the achievable ceiling (Report 1); 0.957 all elements, 0.929 top quintile. Single-pass
top-quintile values as plotted: sequence 0.380 [0.323, 0.437], ATAC 0.548 [0.497, 0.599],
both 0.685 [0.667, 0.704] -- all three intervals disjoint. Reverse-complement-averaged, which
is the current default, the same three are 0.397, 0.551 and 0.695. The multimodal model's
fold-to-fold spread (sd 0.015) is a third of either single-input model's (0.041-0.046), so
combining inputs is also markedly more stable across chromosome sets. Target K562 H3K27ac
(ENCSR000AKP), 5' end counts, 1 kb element-centred window. Source:
`results/{,rc_}fiveprime_stratified_fold_summary.tsv`.
</details>

**Sequence + ATAC is the configuration everything else is measured against**, reaching 75% of
the corrected ceiling on the top quintile against 59% for accessibility alone and 43% for
sequence alone.

**The all-elements panel is not a rescaling of the top-quintile panel** and the two support
different conclusions about the size of every effect; Report 2 opens with why.

**Method.**

- Three input modes, five chromosome-holdout folds each, identical everything else.
- Reported on the top quintile by observed signal, with all elements alongside.

<details>
<summary>Full methods &amp; code</summary>

```bash
sbatch scripts/2.2.submit_evaluate.sh
sbatch scripts/2.22.submit_rc_rescore.sh    # the rc_ values quoted in the decision table
```
Source: `scripts/2.2.evaluate_stratified.py`, `scripts/2.15.perfold_from_config.py`
</details>
"""

R3_TRANSFER = """## Which architecture transfers best: none of them, by a resolvable margin

**Q:** The receptive-field and fragment-channel gains were measured in the training cell
type. Which architecture should actually be deployed to a new cell type?

**A:** The question has a null answer. Every transferred architecture lands within 0.016 of
the simplest one on the top quintile, in both directions, and the two gains that were
significant in-cell type are not significant transferred. The richer models also have the
*larger* transfer drops.

**Top-quintile Pearson, transferred against the target cell type's own ATAC-only model:**

| | K562 -> GM12878 | GM12878 -> K562 |
|---|---|---|
| target's own ATAC-only model (the floor) | 0.518 | 0.584 |
| narrow + flat, transferred | 0.531 | 0.614 |
| wide + flat, transferred | 0.542 | 0.620 |
| narrow + fragments, transferred | 0.530 | 0.611 |
| wide + fragments, transferred | 0.547 | not run |
| *best locally trained model, for scale* | *0.602* | *0.733* |

**Paired within-fold differences against narrow + flat, transferred:**

| | K562 -> GM12878 | GM12878 -> K562 |
|---|---|---|
| wide + flat | +0.0115 [-0.0059, +0.0288] *p*=0.14 | +0.0053 [-0.0161, +0.0268] *p*=0.53 |
| narrow + fragments | -0.0012 [-0.0185, +0.0162] *p*=0.86 | -0.0029 [-0.0267, +0.0208] *p*=0.75 |
| wide + fragments | +0.0159 [-0.0017, +0.0335] *p*=0.066 | not run |

**Fragment channels do not travel at all.** -0.0012 and -0.0029, both directions, both
comfortably null. Whatever the five channels add in-cell type -- and it is real, *p*=0.002
locally -- is either cell-type-specific or library-specific. That is the outcome the
comparison was designed to detect: fragment-length structure is a property of a particular
ATAC library as much as of the chromatin, and nothing in a paired in-cell-type test can tell
those apart.

**The wider receptive field survives directionally but not statistically.** +0.0115 and
+0.0053, and +0.0159 when combined with fragments, which is the largest transferred effect
in the matrix at *p*=0.066. It is the only arm worth carrying forward, and it is not
established.

**The better in-cell-type models lose more on transfer.** Paired local-minus-transferred, top
quintile: narrow + flat +0.0546 and +0.0847, wide + flat +0.0593 and +0.1073, narrow +
fragments +0.0613 and +0.1033 (all *p* < 0.006). Every architecture change that helped
locally increased the transfer drop, which is the signature of fitting cell-type-specific
structure rather than a portable rule.

**What transfer actually buys over doing nothing clever.** Between +0.013 and +0.036 on the
top quintile over the target cell type's own ATAC-only model, against a locally trained model
reaching 0.602 and 0.733. On the mechanistic metric the shortfall is larger: transferred
`residual_pearson` is 0.212-0.262 where a local model reaches 0.406-0.481, so roughly half
of what a model extracts beyond accessibility does not survive the move.

**Method.**

- One config JSON per direction through `2.15.perfold_from_config.py`, so every arm in the
  matrix is scored by the same code path with the same metric definitions.
- Each direction includes the target's own ATAC-only model as the floor and each
  architecture's locally trained counterpart, so the transfer drop is paired per architecture.
- All arms reverse-complement-averaged.

<details>
<summary>Full methods &amp; code</summary>

The four pairings answer four separate questions and are listed in the submit script's own
header so the next reader does not have to reconstruct them: transferred against
`atac_LOCAL` asks whether transferring beats the target's own accessibility; wide against
narrow and fragments against flat, both transferred, ask whether each gain travels; local
against transferred, per architecture, is the transfer drop.

The `WIDE_FRAG` transferred arm exists only in the K562 -> GM12878 direction because there is
no GM12878 wide + fragment model; GM12878 fragment channels required its paired-end BAMs,
which were downloaded later than the K562 set.

```bash
sbatch scripts/2.24.submit_transfer_matrix.sh
```
Source: `scripts/make_transfer_matrix.py`, `scripts/2.24.submit_transfer_matrix.sh`,
`results/txmatrix_{k562_to_gm12878,gm12878_to_k562}_{per_fold,fold_summary}.tsv`
</details>

## The downstream benchmark: predicted H3K27ac does not improve ABC

**Q:** Substituted for observed H3K27ac in ABC's activity term, do the model's predictions
recover CRISPR-benchmark performance?

**A:** No. The best predicted arm reaches AUPRC 0.482 [0.437, 0.530] against a floor of
0.457 (ATAC-only activity) and a ceiling of 0.519 (ATAC x observed H3K27ac), so it does not
clear the floor with a resolvable margin -- and the two deployment arms, GM12878-trained
models applied to K562, land at 0.455 and 0.452, at or below it.

See Fig.~1 for all eleven arms with their intervals.

**The ceiling is only 0.062 above the floor.** Observed H3K27ac itself buys ABC very little
here, so the experiment had limited room from the start. Any conclusion about the model has
to be read against that: a predicted track cannot demonstrate much when the observed track
it replaces demonstrates 0.062.

**Sequence-only arms are far below distance-to-TSS.** 0.393 for the geometric mean with
observed ATAC and 0.275 for predicted-alone, against 0.435 for distance to TSS. A
sequence-only activity track is worse than using no activity information at all.

**Method.**

- Model predictions painted over the ABC candidate regions as `predicted_counts / width`,
  substituted for the H3K27ac term, with observed ATAC retained.
- Floor and ceiling arms run in the same pipeline invocation so nothing but the activity
  track differs.
- Per-chromosome fold assembly from `val` union `test`, so no region is scored by a model
  that trained on its chromosome.

<details>
<summary>Full methods &amp; code</summary>

Nine biosample arms plus the two anchors. `run_qnorm` is rank-based, so the arbitrary scale
of a painted prediction track is irrelevant and only the ordering matters -- the scale
compatibility worry that motivated the check turned out to be moot. The qnorm *reference* is
built from observed K562 signal, so a predicted track still inherits the observed
distribution's shape; a qnorm-off arm is listed as open.

The pipeline is ABC on snakemake 7 with `--use-conda` and a cluster profile, so every rule
instance is a job. Setting up the arms required staggering the `Peaks` file mtimes along the
rule chain, because Snakemake requires outputs strictly newer than inputs and a copied
results tree has uniform timestamps; `4.3.setup_abc_arms.py` does that and then verifies the
chain is strictly increasing rather than assuming it.

```bash
pixi run python scripts/4.1.predict_h3k27ac_for_abc.py ...   # or 4.4 for the CPU sweep
pixi run python scripts/4.3.setup_abc_arms.py
bash scripts/4.5.submit_abc.sh
pixi run python scripts/4.6.setup_crispr_comparison.py
bash scripts/4.7.submit_crispr_comparison.sh
pixi run python scripts/4.8.plot_crispr_benchmark.py
```
Source: `scripts/4.1`, `4.3`-`4.8`
</details>

## Why it fails: the ranking is displaced, not wrong

**Q:** Is the benchmark result caused by inaccurate predictions, or by accurate predictions
whose dynamic range is too compressed for ABC's ranking?

**A:** The latter. Predicted and observed H3K27ac agree at Spearman 0.79 over all 153,545
ABC regions, but the predicted activity term's spread is compressed exactly where it matters:
on the regions carrying a regulated pair, p99/p50 is 4.03 against 6.00 observed, and the
top-decile mean/median ratio is 4.23 against 6.46.

| stratum | *n* | Spearman | Pearson (log) |
|---|---|---|---|
| all regions | 153,545 | 0.794 | 0.848 |
| CRISPR-tested regions | 3,016 | 0.824 | 0.785 |
| regions of regulated pairs | 339 | 0.663 | 0.666 |

**Agreement is worst exactly where the benchmark is decided.** 0.66 on the regions of
regulated pairs against 0.79 genome-wide. The regions that matter are the strongly acetylated
ones, which is the stratum the top-quintile reporting standard exists for.

**The arms disagree on only a handful of pairs.** Of 429 regulated pairs scored by all three,
the observed arm catches 322 and the predicted arm 325; 12 are caught by observed and missed
by predicted, 15 the other way. At those 12, observed H3K27ac RPM has a median of 22.79
against a predicted 0.71, where the genome-wide medians are 0.59 and 0.11 -- so observed
calls them 39x elevated and predicted calls them 6x. The model is not blind to them; it
under-states them by enough to move them out of the top of ABC's ranking.

**Method.**

- Join on ABC's `EnhancerList.txt` per arm, reconstructing coordinates rather than trusting
  the `name` column, whose format differs between the enhancer list and the prediction table.
- Every stratum asserted to hold more than 100 regions before any statistic is computed.

<details>
<summary>Full methods &amp; code</summary>

`EnhancerList` `name` is `class|chr:start-end` while `pred_elements` carries bare
coordinates, so an index built on `name` silently produces an empty join. The script now
indexes on reconstructed `chr:start-end` and asserts each stratum's size.

```bash
sbatch scripts/4.9.submit_diagnose.sh
```
Source: `scripts/4.9.diagnose_abc_gap.py`
</details>

## What the mis-predicted elements are: accessible but unacetylated, and CTCF is only a quarter of it

**Q:** Do the over- and under-predicted elements have a chromatin signature, and in
particular is the over-prediction the CTCF-site pattern -- open but not acetylated?

**A:** The under-predicted tail is unambiguous: canonical active enhancer, EP300 3.4x and
H3K4me1 3.3x above typical with zero H3K27me3. The over-predicted tail is **not** explained
by CTCF, despite a real 2.6x CTCF enrichment in pooled signal. Element by element, CTCF-high
elements are only **27%** of that tail, and removing them leaves the remaining 73% over-predicted
just as badly -- median error 0.82 against 0.76 for the CTCF-high subset, with the CTCF-*low*
subset over-predicted more (3.8x against 1.6x). CTCF is a passenger. The phenotype is
"accessible but unacetylated", of which CTCF sites are one instance.

![](figures/fig15_error_strata.png)

<details>
<summary>Figure 10 legend</summary>

**Figure 10 | The two error tails, what marks them, and why CTCF does not explain the
over-predicted one.** Strata defined by the multimodal model's accessibility-residualised
prediction error over 153,545 K562 ABC candidate regions; the 1% tails are n = 1,536 each and
the typical stratum is the middle 50%, n = 76,771. **a**, The phenotype in
(accessibility, H3K27ac) space. Grey is all regions; over-predicted elements are open and
unacetylated, under-predicted ones acetylated and comparatively closed. **b**, RPKM by mark as
fold change over the typical stratum, log axis, with IgG as the background control -- this is
a *pooled* per-stratum statistic and panel **c** shows why that matters. **c**, Peak overlap
against quantitative signal for the over-predicted tail: peaks report enrichment for every
mark, signal reports CTCF only. **d**, The attributable fraction. Splitting the over-predicted
tail at the 90th percentile of typical-stratum CTCF signal leaves 27% CTCF-high and 73%
CTCF-low, and both are over-predicted; the CTCF-low subset more so. Source:
`results/error_strata_{rpkm,peak_overlap,elements,ctcf_elements}.tsv`.
</details>

Strata are defined by the multimodal model's accessibility-residualised prediction error over
the 153,545 ABC regions.

**Composition (`results/prediction_error_strata.tsv`):**

| stratum | *n* | observed | predicted | ATAC | ATAC/H3K27ac | GC | CpG o/e | % promoter |
|---|---|---|---|---|---|---|---|---|
| under-predicted 1% | 1,536 | 18.56 | 0.65 | 3.90 | 0.21 | 0.466 | 0.289 | 15.6 |
| under-predicted 5% | 7,678 | 11.09 | 0.59 | 3.58 | 0.33 | 0.475 | 0.317 | 18.5 |
| typical (middle 50%) | 76,771 | 0.44 | 0.09 | 0.93 | 2.20 | 0.490 | 0.277 | 13.1 |
| over-predicted 5% | 7,678 | 0.52 | 0.35 | 5.23 | 9.09 | 0.574 | 0.475 | 25.0 |
| over-predicted 1% | 1,536 | 0.59 | 0.80 | 9.93 | 13.25 | 0.593 | 0.551 | 31.6 |

The `ATAC/H3K27ac` column is the accessible-but-unacetylated signature, and it spans a factor
of 60 across the strata. The over-predicted tail is GC-rich, CpG-rich and promoter-enriched;
the under-predicted tail is GC-*poor* and distal.

**Peak overlaps and quantitative signal disagree, and the quantitative version is the one to
believe.** Fraction of each stratum overlapping a peak:

| stratum | CTCF | EP300 | H3K4me1 | H3K27me3 | H3K27ac |
|---|---|---|---|---|---|
| over-predicted 1% | 23.4% | 16.5% | 30.4% | 8.1% | 33.3% |
| over-predicted 5% | 36.9% | 16.4% | 35.5% | 9.6% | 30.6% |
| typical | 17.4% | 8.8% | 28.3% | 11.5% | 15.8% |
| under-predicted 1% | 9.3% | 57.7% | 86.8% | 0.0% | 100.0% |
| under-predicted 5% | 10.8% | 52.3% | 83.2% | 0.0% | 99.6% |

Read as peaks, the over-predicted tail looks enriched for everything -- CTCF 2.1x, EP300
1.9x, H3K27ac 1.9x -- which supports no mechanism in particular. RPKM from the same BAMs,
with an IgG control, says something much sharper:

| stratum | CTCF | EP300 | H3K4me1 | H3K27me3 | H3K27ac | IgG |
|---|---|---|---|---|---|---|
| over-predicted 1% | 3.34 | 1.07 | 0.97 | 0.44 | 0.94 | 0.33 |
| over-predicted 5% | **4.69** | 1.16 | 1.24 | 1.17 | 1.17 | 0.39 |
| typical | 1.79 | 0.85 | 1.31 | 0.53 | 1.24 | 0.48 |
| under-predicted 1% | 1.20 | **2.90** | **4.35** | 0.31 | **17.52** | 0.68 |
| under-predicted 5% | 1.31 | 2.42 | 3.78 | 0.31 | 12.57 | 0.65 |

**The CTCF enrichment is real but it is a pooled statistic, and pooled statistics are how
this analysis went wrong twice.** RPKM per stratum is total reads over total kilobases, so a
minority of strongly bound elements sets the value. Element by element the over-predicted 1%
has CTCF **median 0.72 RPKM against 0.55 typical** -- a 1.3x difference, not 2.6x -- and only
**15.7%** of it exceeds the 90th percentile of typical-stratum CTCF (against 10.0% of typical
by construction). For the over-predicted 5% the figures are 26.5% CTCF-high and a median of
1.10. So the tail is not a set of CTCF sites; it contains rather more CTCF sites than average.

**The decisive test is what happens when the CTCF-high elements are removed, and CTCF fails
it.** Splitting the over-predicted 5% at that threshold:

| subset | *n* | % of tail | observed fold | predicted fold | median error | CTCF median |
|---|---|---|---|---|---|---|
| whole tail | 7,678 | 100 | 0.87 | 3.10 | 0.800 | 1.10 |
| CTCF-high | 2,038 | 26.5 | 0.74 | 1.63 | 0.759 | 13.46 |
| **CTCF-low** | **5,640** | **73.5** | 0.94 | **3.83** | **0.821** | 0.65 |
| typical | 76,771 | -- | 0.74 | 0.78 | 0.042 | 0.55 |

The CTCF-low three-quarters is over-predicted *more* than the CTCF-high quarter, on both the
error and the predicted-fold columns. Whatever causes the over-prediction is present with or
without CTCF binding.

**GC content is a much better correlate of the error than CTCF is.** Within accessibility
deciles -- necessary because CTCF sites are open and the model reads accessibility -- Spearman
rho(error, CTCF) is 0.02 to 0.23 with a median of 0.084, while rho(error, GC) is 0.24 to 0.37
in **every** decile. Partialling accessibility, GC, CpG and promoter class out leaves
rho(error, CTCF) at +0.191, above the IgG control's -0.126 in magnitude but with a control
that is not cleanly zero, so read that number as suggestive at best.

**What survives.** H3K27ac itself is at 0.9x in the over-predicted tail and H3K4me1 at 0.95x,
so these elements genuinely carry no more of the mark than a random element despite being
5-10x more accessible; and IgG is *below* typical there (0.39 against 0.48), so none of it is
an antibody-accessibility artifact. The description "accessible, GC-rich, CpG-rich,
promoter-enriched and unacetylated" is supported element by element. "CTCF sites" is not.

**The under-predicted tail is canonical active enhancer.** EP300 3.4x and H3K4me1 3.3x above
typical, against an IgG elevation of only 1.4x, so both are real; H3K27me3 is *depleted*
(0.31 against 0.53) and no element in the stratum overlaps an H3K27me3 peak; CTCF is depleted
(0.67x). These are the elements the model most needs to get right and most under-states.

**Method.**

- Strata residualised on accessibility before defining the error, so they are not simply
  "more open"; the composition table above shows the residualisation does not remove the ATAC
  difference, which is the point.
- Enrichment computed **both** ways: peak overlap with `bedtools`, and RPKM from the same
  BAMs with mapped-read totals as denominators and an IgG control alongside.

<details>
<summary>Full methods &amp; code</summary>

Peak overlaps alone were the first version and were misleading, for the reason above: a peak
call is a thresholded statement and a stratum can be enriched for weak peaks without carrying
more signal. Quantitative signal with a background control is the version to quote. Both are
kept because the disagreement between them is itself the finding.

Mapped-read totals used as RPKM denominators: CTCF 9,231,994; EP300 24,485,729; H3K4me1
8,476,072; H3K27me3 73,060,514; H3K27ac 6,335,377; IgG 24,933,797. Track sources are in
`reference/DATA_INVENTORY.md`.

```bash
sbatch scripts/4.10.characterize_prediction_errors.py --emit-beds
sbatch scripts/4.11.annotate_error_strata.sh    # peak overlaps
sbatch scripts/4.13.quantify_error_strata.sh    # RPKM + IgG
```
Source: `scripts/4.10`, `4.11`, `4.13`
</details>

## Does the sequence branch see the CTCF signature? It does not separate either tail

**Q:** CpG islands and CTCF motifs are strong local sequence signals and the first
convolution is 21 bp, wide enough for a 19 bp CTCF motif. Does the sequence branch already
know these elements are not acetylated, with the accessibility channel overriding it, or does
it fail to see the signature at all?

**A:** It fails to see it. On a scale-free statistic the sequence-only arm puts the
over-predicted and under-predicted tails at essentially the *same* level -- fold elevation
1.55x and 1.44x, ratio 1.07 -- where the truth separates them by a factor of 31 (ratio 0.032).
Its apparent near-correctness on the over-predicted tail is an artifact of having almost no
dynamic range.

**Fold elevation against each arm's own genome-wide median**, over the 153,545 ABC regions:

| stratum | *n* | observed | sequence | ATAC only | multimodal |
|---|---|---|---|---|---|
| under-predicted 1% | 1,536 | **31.33** | 1.44 | 5.01 | 5.81 |
| under-predicted 5% | 7,678 | 18.72 | 1.52 | 4.64 | 5.30 |
| typical (middle 50%) | 76,771 | 0.74 | 0.93 | 0.79 | 0.78 |
| over-predicted 5% | 7,678 | 0.87 | 1.17 | 4.19 | 3.10 |
| over-predicted 1% | 1,536 | **1.00** | 1.55 | 8.41 | 7.12 |

Responsiveness-normalised over-prediction, the over-1% elevation divided by the same arm's
under-1% elevation, which is scale-free within an arm: **observed 0.032, sequence 1.07,
multimodal 1.22, ATAC-only 1.68.**

**Reading the over-predicted column alone inverts the conclusion, and did.** Taken by itself,
sequence at 1.55x against a truth of 1.00x looks nearly right while multimodal is 7.12x, and
that reads as "the signal is in sequence and accessibility drowns it". The under-predicted
column shows why that is wrong: the sequence arm only reaches 1.44x where the truth is 31x,
so it is not correctly declining at the CTCF elements, it is flat everywhere. Both hypotheses
predict a low number in the over-predicted column and only one of them survives the control.

**This also explains the sequence arms' benchmark scores.** A predictor with almost no
dynamic range over ABC regions cannot rank them, which is what AUPRC 0.275 for
predicted-sequence-alone means.

**The indicated fix is representational.** Explicit GC and CpG-density input channels, or
motif-derived features, rather than rebalancing what the model already has. The gating and
loss-reweighting arms below were built on the opposite hypothesis and are worth scoring
anyway, since they were cheap and are already trained -- but this result lowers the prior on
them.

**Method.**

- Fold elevation against each arm's **own** genome-wide median, because the painted
  prediction tracks have medians of 0.11-0.18 against 0.59 for observed RPM and are not on a
  common scale.
- The opposite tail carried as each arm's own responsiveness scale, so a compressed predictor
  cannot pass by being flat.

<details>
<summary>Full methods &amp; code</summary>

The first version of this analysis used accessibility-residualised error and printed the
opposite verdict. It could not have worked: the over-predicted stratum is *defined* as
accessibility-says-high and H3K27ac-says-low, so the observed residual is strongly negative
there and any predictor that does not track accessibility downward -- including a constant --
scores as over-predicting. Report 2 records the general form of the mistake.

```bash
srun -p engreitz -t 20 --mem=24G python scripts/4.12.can_sequence_see_it.py
```
Source: `scripts/4.12.can_sequence_see_it.py`
</details>

## p300 as a second target: not learnable where H3K27ac fails

**Q:** The under-predicted elements are EP300-enriched. p300 is one step closer to the
sequence-specified events, so would a p300-target model -- or a multi-head model predicting
both -- recover them?

**A:** No. Where the H3K27ac model fails, observed p300 is elevated 4.89x but its own input
control is at 2.10x, so real enrichment is only about 2.3x; the p300 models predict 1.83x
(multimodal) and 1.55x (ATAC-only), *below* that control. The p300 models fail on the same
elements, so neither a second head nor stacking p300 predictions as an input can help.

**Fold elevation against each track's own genome-wide median:**

| stratum | *n* | observed H3K27ac | predicted H3K27ac | observed p300 | p300 input control | predicted p300 (mm) | predicted p300 (ATAC) |
|---|---|---|---|---|---|---|---|
| H3K27ac under-predicted 1% | 1,536 | 31.33 | 5.81 | **4.89** | **2.10** | **1.83** | 1.55 |
| H3K27ac under-predicted 5% | 7,678 | 18.72 | 5.30 | 3.88 | 1.80 | 1.76 | 1.53 |
| typical (middle 50%) | 76,771 | 0.74 | 0.78 | 0.79 | 0.90 | 0.80 | 0.85 |
| H3K27ac over-predicted 5% | 7,678 | 0.87 | 3.10 | 1.33 | 0.80 | 1.59 | 1.60 |
| H3K27ac over-predicted 1% | 1,536 | 1.00 | 7.12 | 1.36 | 0.50 | 1.92 | 1.79 |

**The input control changed the answer.** Peak overlap put EP300 at 6.6x enrichment in this
stratum, which made a p300 target look promising. Against its own input control the real
elevation is 2.3x -- the premise was substantially an accessibility artifact, which is why the
control was included before the expensive version was built.

**The p300 models fail the same way, in the same direction.** At the over-predicted tail they
predict 1.92x and 1.79x where observed p300 is 1.36x and its control is 0.50x. So p300
models over-predict the accessible-but-unmarked elements too. Whatever causes the H3K27ac
model's failure is not specific to the H3K27ac target.

**Consequence.** Both p300 ideas are dropped: the multi-head model, and using predicted p300
as an input feature. This closes a branch that would otherwise have cost roughly fifteen GPU
jobs. It does not affect the separate, already-reported finding that p300 is the better
substrate for *sequence* attribution (Fig.~3) -- that is about where motif work should be
done, not about predicting H3K27ac.

**Method.**

- Existing p300 models scored over the same ABC candidate regions through the same prediction
  path as the H3K27ac arms, with `4.4` generalised by environment overrides rather than
  forked, so the two targets cannot differ by which script ran.
- p300's own **input control** track carried through every stratum, since p300 ChIP without
  its control cannot distinguish binding from accessibility.
- The verdict threshold was written into the script before it ran: predicted p300 must reach
  at least half of observed p300's elevation for the multi-head idea to proceed.

<details>
<summary>Full methods &amp; code</summary>

p300 checkpoints predate the `mode` attribute, so loading them raises `AttributeError`;
`4.1` already handles that by setting `m.mode = mode` when the attribute is absent. Model
directories are `2026_0529_multimodal_p300_model/models/atac` (multimodal, despite the name)
and `models/atac_only`. Target and control tracks are
`2025_0703_retrain_p300_model/data/ENCSR000EGE_{plus,minus}.bigWig` and
`ENCSR000EGE_control_{plus,minus}.bigWig`.

```bash
MODEL_ROOT=... OUT_PREFIX=predp300 sbatch scripts/4.4.submit_predict_h3k27ac_cpu.sh atac 0
sbatch --dependency=afterok:... scripts/4.14.p300_at_h3k27ac_failures.py
```
Source: `scripts/4.4`, `scripts/4.14.p300_at_h3k27ac_failures.py`
</details>
"""



R3_P300_ABC = """## Predicted p300 works where predicted H3K27ac does not

**Q:** The H3K27ac arms failed the benchmark. p300 is one step closer to the
sequence-specified events -- does substituting *predicted p300* for the activity term do
better, is observed p300 even a better activity term to begin with, and does any of it
survive transfer to a cell type with only ATAC?

**A:** Yes to both, and the second explains part of the first. Observed p300 beats observed
H3K27ac as ABC's activity term by **+0.038 AUPRC [+0.019, +0.058]**. A sequence + ATAC model
predicting p300 beats the ATAC-only floor by **+0.055 [+0.034, +0.074]** with the sign
preserved in 100% of paired resamples, and is indistinguishable from *measured* H3K27ac
(-0.007 [-0.032, +0.017]). This is the first predicted-activity arm in the project to clear
the floor with a resolvable margin. **It does not transfer**: a GM12878-trained p300 model
applied to K562 lands at +0.009 [-0.004, +0.022] over the floor, an interval spanning zero.

![](figures/fig16_p300_benchmark.png)

<details>
<summary>Figure 11 legend</summary>

**Figure 11 | p300 as the ABC activity term: the gain is real in the training cell type and
absent on transfer.** Forest plot of AUPRC against the CRISPR benchmark for all nine arms plus
distance-to-TSS, on 10,342 element-gene pairs with 466 regulated, every arm scored on the
identical pair set. Error bars are the pipeline's **unpaired** per-predictor 95% intervals;
the dashed lines mark the ATAC-only floor (0.457) and observed H3K27ac (0.519). The paired
bootstrap delta against the floor is printed at the right of each row, because the unpaired
intervals drawn here overlap for every arm and cannot resolve any of the comparisons -- the
figure shows both so that the misleading layer and the informative one sit together. Observed
p300 (0.557) is the highest arm; the K562-trained predicted-p300 arm (0.512) reaches observed
H3K27ac; the GM12878-trained arm applied to K562 (0.466) sits at the floor with an interval
spanning zero. Geometry, floor and reference line match Fig.~1 so the two are directly
comparable, which the two shared anchor arms agreeing to 1e-4 across the runs licenses.
Source: `CRISPR_comparison_v3/.../results/2026_0906_p300_all/performance_summary.txt`,
`scripts/4.17` for the paired deltas.
</details>

**AUPRC on the CRISPR benchmark, 10,342 element-gene pairs, 466 regulated:**

| arm | AUPRC | unpaired 95% CI |
|---|---|---|
| ATAC x **observed p300** | **0.557** | [0.506, 0.605] |
| ATAC x observed H3K27ac | 0.519 | [0.468, 0.567] |
| ATAC x **predicted p300**, sequence + ATAC model | **0.512** | [0.463, 0.560] |
| predicted p300 alone, sequence + ATAC model | 0.502 | [0.454, 0.550] |
| **ATAC x predicted p300, GM12878 model -> K562** | **0.466** | [0.420, 0.516] |
| predicted p300 alone, GM12878 model -> K562 | 0.462 | [0.416, 0.512] |
| ATAC x predicted p300, ATAC-only model | 0.461 | [0.415, 0.511] |
| ATAC only (floor) | 0.457 | [0.408, 0.507] |
| predicted p300 alone, ATAC-only model | 0.443 | [0.394, 0.491] |
| distance to TSS | 0.435 | [0.386, 0.485] |

**The unpaired intervals resolve nothing and the paired ones resolve everything**, which is
the same lesson as the within-fold pairing in Report 2, applied downstream. Every interval
above overlaps every other; the arms are scored on the *identical* pair set, so resampling
pairs once and recomputing all arms on that resample removes the shared variance:

| comparison | delta | paired 95% CI | sign kept |
|---|---|---|---|
| observed p300 - observed H3K27ac | **+0.038** | [+0.019, +0.058] | 99.9% |
| observed p300 - floor | **+0.099** | [+0.079, +0.120] | 100% |
| **predicted p300 (seq + ATAC) - floor** | **+0.055** | **[+0.034, +0.074]** | **100%** |
| predicted p300 (seq + ATAC) - observed H3K27ac | -0.007 | [-0.032, +0.017] | 71% |
| predicted p300 alone (seq + ATAC) - floor | +0.045 | [+0.012, +0.077] | 99.6% |
| predicted p300 (ATAC-only model) - floor | +0.004 | [-0.012, +0.021] | 70% |
| predicted p300, seq + ATAC - ATAC-only model | **+0.051** | [+0.032, +0.069] | 100% |
| **transferred (GM12878 -> K562) - floor** | **+0.009** | **[-0.004, +0.022]** | 91% |
| transferred, predicted-alone - floor | +0.005 | [-0.019, +0.028] | 65% |
| **local - transferred (the transfer drop)** | **+0.046** | [+0.030, +0.061] | 100% |
| transferred - observed H3K27ac | -0.053 | [-0.074, -0.032] | 100% |

**Sequence is what makes it work, and that is the opposite of the H3K27ac result.** The
multimodal p300 model beats the ATAC-only p300 model by +0.051 with the sign preserved in
every resample, while the ATAC-only p300 model is at the floor (+0.004, interval spanning
zero). For H3K27ac the sequence-containing arms were the *worst* in the whole panel -- 0.393
for the geometric mean and 0.275 predicted-alone, both far below distance-to-TSS. Same
architecture, same accessibility input, same benchmark: the target was the binding constraint.

**Two independent things contribute and they should not be conflated.** Observed p300 is a
better activity term than observed H3K27ac (+0.038), so the *concept* has more headroom over
the floor -- 0.099 against 0.062. And the model captures more of the headroom it is given:
55% of p300's against a non-resolvable fraction of H3K27ac's. The first of these is a fact
about the assay and survives the transfer result; the second is a fact about a K562-trained
model on K562 and does not.

**The transfer result is the one that matters for the stated goal, and it is a null.** The
whole point is a model deployable in any cell type with only ATAC. Transferred, predicted p300
buys +0.009 [-0.004, +0.022] over that cell type's own ATAC-only activity -- nothing. The drop
of -0.046 [-0.061, -0.030] removes 83% of the in-cell-type gain, and the transferred arm is
-0.053 below measured H3K27ac. Deployment now fails for **both** targets: the H3K27ac
deployment arms scored 0.455 and 0.452 against the 0.457 floor, the p300 one scores 0.466. So
the barrier is transfer, not target choice, which is the same conclusion the architecture
programme reached (F-005) by a completely different route.

**The missing control has now run, and it changes the conclusion.** Scoring both p300 models
in both cell types on their own targets -- the local-versus-transferred 2x2 the H3K27ac matrix
used throughout -- separates "transfer is fatal" from "one model is weaker", and the answer is
neither. Top-quintile Pearson on each cell type's own EP300:

| | evaluated on K562 | evaluated on GM12878 |
|---|---|---|
| K562-trained p300 | **0.619** (local) | **0.507** (transferred) |
| GM12878-trained p300 | **0.312** (transferred) | **0.608** (local) |
| that cell type's ATAC-only p300 | 0.306 (floor) | 0.299 (floor) |

| paired within-fold comparison | delta | 95% CI | *p* |
|---|---|---|---|
| K562 local - K562 floor | +0.313 | [+0.286, +0.340] | <1e-4 |
| GM12878 local - GM12878 floor | +0.309 | [+0.263, +0.355] | 1e-4 |
| **K562 -> GM12878, transferred - GM12878 floor** | **+0.207** | **[+0.158, +0.257]** | 3e-4 |
| **GM12878 -> K562, transferred - K562 floor** | **+0.006** | **[-0.038, +0.050]** | **0.72** |
| K562 -> GM12878 drop (transferred - local) | -0.102 | [-0.121, -0.083] | 1e-4 |
| GM12878 -> K562 drop (transferred - local) | -0.307 | [-0.327, -0.286] | <1e-4 |

**Both models are equally strong at home**, +0.313 and +0.309 over matched floors of 0.306 and
0.299, so the GM12878 model is not the weaker model and that confound is closed.

**Transfer is strongly asymmetric, and the benchmark could only test the failing direction.**
K562-trained transferred to GM12878 keeps 0.507 of its 0.619, which is **+0.207 over the
target cell type's own accessibility model** and 83% of the local advantage, with a positive
accessibility residual (+0.139 [+0.067, +0.210]). GM12878-trained transferred to K562 lands
at 0.312 against a 0.306 floor -- indistinguishable, *p*=0.72 -- with a **negative** residual
(-0.088 [-0.140, -0.035]), meaning it is worse than the accessibility baseline at the thing
the baseline already does. CRISPR data exists only for K562, so the downstream benchmark can
*only* run the GM12878 -> K562 direction. Its null is real for that direction and says nothing
about the other.

**What this does and does not license.** It does not license "predicted p300 deploys": the
+0.207 is top-quintile Pearson on p300, not benchmark AUPRC, and no downstream test of the
working direction is possible with K562-only CRISPR data. It does retire two earlier readings
-- that the p300 gain simply does not transfer, and that the GM12878 model was weak. And it
makes the training cell type a first-class design variable: the same architecture and target
gives +0.207 or +0.006 over the floor depending only on which cell type it was trained in.

**Why the asymmetry: more training signal in K562, but not cleaner training signal.** Both
cell types are equally *predictable* locally, so the H3K27ac explanation for its own transfer
asymmetry -- that GM12878 is the harder cell type -- cannot apply here. Comparing the two
EP300 experiments directly:

| | K562, ENCSR000EGE | GM12878, ENCSR000DZG | ratio |
|---|---|---|---|
| peaks | 28,532 | 21,068 | 1.35x |
| mean peak width | 296 bp | 337 bp | 0.88x |
| mapped reads (pooled) | 51,127,010 | 30,001,681 | 1.70x |
| reads in peaks | 3,590,540 | 1,562,462 | **2.30x** |
| FRiP (pooled) | 0.070 | 0.052 | 1.35x |
| FRiP (per replicate) | 0.053, 0.086 | 0.067, 0.038 | -- |

**The quantity of training signal differs materially; its purity does not.** K562 carries
2.3x more reads inside peaks, from 1.7x the depth over 1.35x the peaks. But FRiP does **not**
separate the two experiments once replicate variance is accounted for: the within-experiment
spread is 1.6x in K562 (0.053 to 0.086) and 1.8x in GM12878 (0.038 to 0.067), both wider than
the 1.35x pooled difference, and GM12878's better replicate has a *higher* FRiP than K562's
worse one. So the specific "cleaner target" version of the hypothesis is not supported; the
"more signal to learn portable features from" version is consistent with the data and remains
untested.

This is diagnostic, not decisive -- two experiments cannot establish a relationship between
training-signal volume and portability. What it does is narrow the next step: a third cell
type is required, and it should be chosen for EP300 depth rather than convenience.

**This does not contradict the earlier p300 negative, and the distinction matters.** The
learnability test below shows p300 is *not* predictable at the specific elements where the
H3K27ac model fails, which correctly killed the multi-head and stacking ideas. It says
nothing about p300 as the *primary* target, where the relevant question is whether the
predicted track ranks elements well enough for ABC -- a different quantity, now measured, and
positive.

**Method.**

- Own ABC results directory and config, reusing the July floor and observed-H3K27ac ceiling
  so the numbers are directly comparable to the H3K27ac run.
- Observed p300 carried as this concept's own ceiling; a predicted arm with no observed anchor
  decides nothing.
- Differences tested by paired bootstrap over the shared pair set, 2,000 resamples.

<details>
<summary>Full methods &amp; code</summary>

**Cross-run comparability is verified, not assumed.** The two shared anchors reproduce across
the two independent pipeline runs: ATAC-only floor 0.4573 here against 0.4572 in the H3K27ac
run, and observed H3K27ac 0.5192 against 0.5192. Agreement to 1e-4 on both is what licenses
comparing 0.512 here against 0.482 there.

**A systematic offset between the two AUPRC estimators, which does not affect any delta.**
`4.17` computes average precision directly and reports 0.5673 / 0.5296 / 0.4680 for observed
p300 / observed H3K27ac / floor where the pipeline reports 0.5569 / 0.5192 / 0.4573 -- a
constant +0.0105 on all three, so it is a tie-handling or interpolation difference in the
estimator rather than a difference in the data. Deltas are unaffected, which is why the
paired table quotes deltas and the unpaired table quotes the pipeline's values.

The ABC driver for this run stalled at 18 of 21 steps with every child job COMPLETED and
nothing queued; the outstanding steps were a QC plot and aggregation, which nothing downstream
reads. The five `EnhancerPredictionsAllPutative.tsv.gz` files were verified as valid gzip at
10.1-10.3M rows each before the driver was cancelled and the benchmark run.

```bash
pixi run python scripts/4.3.setup_abc_arms.py --arms p300 \\
    --results-dir results/2026_0905_p300_activity
ABC_CFG=config/mine/config_p300_activity.yaml ABC_ARMS=p300 \\
    ABC_RESULTS_DIR=results/2026_0905_p300_activity sbatch scripts/4.5.submit_abc.sh
pixi run python scripts/4.6.setup_crispr_comparison.py --arms p300
CC_CFG=.../config/config_p300_activity.yml sbatch scripts/4.7.submit_crispr_comparison.sh
pixi run python scripts/4.17.paired_auprc_bootstrap.py <merged.txt.gz> --pair A B ...
```
Source: `scripts/4.3`, `4.5`-`4.7`, `4.17.paired_auprc_bootstrap.py`
</details>
"""


R3_GCMATCH = """## GC-matching the training negatives changes nothing

**Q:** Every model in this repository was trained with ChromBPNet's genome-wide GC-*annotated*
tiling passed straight to `--negatives`, which is the matching INPUT rather than its output
(Report 1). The pool sits at mean GC 0.389 against 0.466-0.593 for candidate elements, and it
is 10% of every batch, so a sequence branch could satisfy much of that contrast with a coarse
GC detector -- the leading explanation for why the H3K27ac sequence branch has almost no
dynamic range within candidate elements. Does actually GC-matching the pool help?

**A:** No, and not marginally: **+0.0011 [-0.0080, +0.0102]** in-cell type (*p*=0.76) and
**-0.0024 [-0.0083, +0.0036]** transferred to GM12878 (*p*=0.33). Absolute top-quintile
Pearson is 0.700 against 0.699 locally and 0.529 against 0.531 transferred. The hypothesis is
refuted.

| | matched negatives | unmatched (baseline) | paired difference |
|---|---|---|---|
| in-cell K562 | 0.700 [0.689, 0.711] | 0.699 [0.695, 0.703] | +0.0011 (*p*=0.76) |
| K562 -> GM12878 | 0.529 [0.501, 0.556] | 0.531 [0.500, 0.562] | -0.0024 (*p*=0.33) |
| transferred vs the target's own ATAC-only model | +0.0108 (*p*=0.28) | +0.0132 (*p*=0.25) | -- |

**The transfer arm was the one the hypothesis actually predicted**, since a GC shortcut should
cost most where accessibility does not generalise. It moved by -0.002. The transferred model
still beats the target's own ATAC-only model by the same non-significant ~0.011 it did before.

**The match was verified before the result was read**, which is what makes this a real null
rather than a failed intervention. `scripts/0.29` recovers GC by joining the matched pool's
coordinates back to the tiling's own annotation, so both sides are measured with one
estimator: the mean-GC gap to the elements fell from 0.0842 to 0.0128 (6.6x closer) in K562
and 0.0655 to 0.0141 in GM12878, and the largest per-bin discrepancy fell from 0.072 to 0.015.

**The one caveat, stated in advance.** The matched pool still undershoots the GC-rich tail --
p90 0.570 against the elements' 0.590 -- because a genome-background pool does not contain
enough very-high-GC windows, which is also why the matching tool hangs on the full element set
(Methods). So this does not exclude a shortcut confined to the extreme high-GC end. It does
exclude one large enough to be worth chasing: closing 85% of the mean gap and 80% of the
per-bin discrepancy moved the metric by 0.001.

**Consequence.** Training-element composition is closed as an explanation for the sequence
branch's flatness, and with it the cheapest remaining fix. The two hypotheses that survived
the session -- explicit GC/CpG input channels, and the accessible-but-unacetylated reweighting
-- are now the only composition-adjacent ideas left, and both are weaker a priori than the one
just refuted.

**Method.**

- `1.20` is `1.11` with exactly two changes, verified by diff on the comment-stripped scripts:
  the negatives file and the output directory. The baseline is `1.11`'s own output with
  identical settings, so the paired comparison isolates the negatives.
- Negatives built with `bpnet-gc-background`, the same tool used for the p300 v3 negatives in
  Oct 2025, rather than a reimplementation -- so a matched-versus-unmatched comparison cannot
  be confounded by two different matching algorithms.
- Scored in-cell type and transferred, paired within fold, RC-averaged.

<details>
<summary>Full methods &amp; code</summary>

`residual_pearson` is reported in these tables against the **unmatched multimodal model**,
because that is this config's baseline entry -- not against an ATAC-only model as everywhere
else in this report. It is therefore not comparable to any other `residual_pearson` here and
is not quoted.

Two failure modes of `bpnet-gc-background` are recorded in the repository learnings: it exits
0 having written a 0-byte file when `bedtools` is absent from its environment, and it hangs
rather than reporting a shortfall when the foreground is GC-rich enough to exhaust the high-GC
bins -- 93% in 15 seconds, then no progress for two hours. The fix for the second is to
subsample the foreground to 50,000 regions, since only the GC distribution defines the
matching target and 50,000 equals the trainer's `--max-negatives` cap.

```bash
sbatch scripts/0.28.make_gc_matched_negatives.sh
pixi run python scripts/0.29.verify_gc_match.py      # before reading the retrain
for F in 0 1 2 3 4; do sbatch scripts/1.20.submit_training_gcmatched.sh multimodal $F; done
sbatch scripts/2.26.submit_gcmatch_eval.sh
```
Source: `scripts/0.28`, `0.29`, `1.20`, `2.26`
</details>
"""

R3_GATE = """## Sequence-gated accessibility and an asymmetric count loss

**Q:** The multimodal trunk merges a sequence branch and an accessibility branch additively,
so it has no way to *withhold* accessibility at elements where accessibility is misleading.
Two cheap interventions target that: a gate that lets sequence modulate the accessibility
representation position by position, and a loss that penalises over-prediction more than
under-prediction. Do either help?

**A:** Neither helps, and together they hurt. In-cell type all three arms are
indistinguishable from the ungated symmetric baseline; transferred to GM12878 the gate is
null, the asymmetric loss trends negative, and the combination is **significantly worse**.

**Paired within-fold difference in top-quintile Pearson against the matched ungated baseline:**

| arm | in-cell K562 | K562 -> GM12878 |
|---|---|---|
| gate | -0.0001 [-0.0152, +0.0151] *p*=0.99 | -0.0015 [-0.0190, +0.0159] *p*=0.82 |
| asymmetric loss | +0.0001 [-0.0122, +0.0124] *p*=0.98 | -0.0078 [-0.0183, +0.0026] *p*=0.11 |
| both | -0.0032 [-0.0157, +0.0092] *p*=0.51 | **-0.0128 [-0.0246, -0.0010]** *p*=0.039 |

For scale, on the transfer arm the gated model beats the target cell type's own ATAC-only
model by +0.0117 (*p*=0.21) and the ungated baseline beats it by +0.0132 (*p*=0.25) -- so both
sit in the same place relative to the floor, and neither intervention moved anything.

**This is the outcome the corrected sequence analysis predicted.** Both interventions assume
the sequence branch already knows which elements are not acetylated and is being overridden
by accessibility. The section above shows it does not know: the sequence arm's elevation is
1.55x at the over-predicted tail and 1.44x at the under-predicted tail, where the truth
separates them 1.00x against 31.3x. A gate cannot surface information that is not in the
branch it gates, and reweighting the loss cannot either. The prediction was recorded before
the scoring ran.

**Why the combination is actively worse on transfer is worth one line of speculation and no
more.** Tripling the over-prediction penalty pushes predictions down; a gate that can also
suppress accessibility gives a second route to the same thing, and on transfer -- where the
accessibility term is the part that generalises -- suppressing it costs more than the
over-prediction it avoids. That is a hypothesis, not a result; the experiment was not designed
to test it and neither arm should be pursued.

**The design.** A 2x2 factorial against a matched ungated symmetric baseline -- same
accessibility bigwig, same `n_layers`, same `count_loss_weight`, differing only in the two
flags -- so the interaction is identifiable rather than inferred from two separate one-change
runs.

- **GATE.** `X_acc = X_acc * sigmoid(gate_conv(X_seq))`, a 21 bp convolution from the
  sequence representation to the accessibility channels, initialised open (`weight = 0`,
  `bias = 4.0`, so the gate starts at 0.982 and the model begins as its own baseline). This
  is a mixture-of-experts gate in the sense asked for: it learns, per position, how much to
  weight accessibility against sequence, so promoters and distal elements can be treated
  differently without that distinction being hand-specified.
- **ASYM.** `overprediction_weight = 3.0` in the count loss: squared error is tripled where
  the prediction exceeds the target. The default `log1pMSE` is already about 400x more
  sensitive to missing signal than to inventing it in the regime these elements sit in, so
  this pushes in the direction the error analysis says is needed.
- **GATE + ASYM**, and the ungated symmetric baseline.

**Scored in-cell type *and* transferred, because the two questions differ.** In-cell type asks
what the fix costs where accessibility is trustworthy; K562 -> GM12878 asks whether it helps
where over-reliance actually bites, which is the reason the gate exists.

**Verdict: both closed.** The gate and the asymmetric loss are the two cheapest
interventions against accessibility over-reliance and neither works, for a reason that is now
measured rather than guessed. Anything further in this direction has to make the sequence
branch discriminative first.

**Method.**

- Both flags live in `scripts/multimodal_bpnet.py` and
  `scripts/train_multimodal_bpnet.py` behind `--gate-accessibility` and
  `--overprediction-weight`, rather than in a forked model file, so the baseline and the
  variants are the same code.
- Guarded by `scripts/test_asymmetric_loss.py`, a four-check regression gate.

<details>
<summary>Full methods &amp; code</summary>

The asymmetric loss mirrors bpnetlite's `_mixture_loss` so that
`overprediction_weight = 1.0` reproduces it exactly; both the train and validation paths
route through the same function via `getattr(self, "overprediction_weight", 1.0)`, so a
pre-gate checkpoint still loads and scores identically.

The regression gate asserts its own preconditions after a first version that verified
nothing -- see Report 2. Its four checks: weight 1.0 matches bpnetlite to a relative
tolerance of 1e-5; a weight above 1 raises the loss by exactly the predicted amount, with
`3 <= n_over <= N-3` asserted so the check cannot be vacuous; a pre-gate checkpoint loads and
is deterministic; and a **weight-matched** gate comparison differs by 0.39%, the matching
being necessary because adding `gate_conv` consumes extra RNG draws and unmatched
initialisations differ by 108%.

```bash
pixi run python scripts/test_asymmetric_loss.py
pixi run python scripts/make_gate_scripts.py
for F in 0 1 2 3 4; do
    sbatch scripts/1.19.submit_training_gate.sh multimodal $F
    sbatch scripts/1.19.submit_training_asym.sh multimodal $F
    sbatch scripts/1.19.submit_training_gate_asym.sh multimodal $F
done
sbatch scripts/2.25.submit_gate_factorial.sh
```
Source: `scripts/multimodal_bpnet.py`, `scripts/train_multimodal_bpnet.py`,
`scripts/1.19.*`, `scripts/2.25.submit_gate_factorial.sh`
</details>

## Open questions / next steps

**The central question now has a positive answer, but for p300 rather than H3K27ac.**
Predicted H3K27ac does not improve ABC (Fig.~1); predicted p300 from a sequence + ATAC model
clears the floor by +0.055 [+0.034, +0.074] and matches measured H3K27ac. Two things follow
for planning: the p300 target is where the effort should go, and the H3K27ac diagnostics
(compressed dynamic range, a sequence branch that does not discriminate) explain why *that*
target underperformed rather than describing the approach as a whole.

**1. Explain the transfer asymmetry, because it is now the largest lever in the project.**
The same architecture and target gives +0.207 or +0.006 over the target's floor depending
only on which cell type it was trained in, and both models fit their own cell type equally
well. Compare the two EP300 experiments directly -- peak count, depth, FRiP,
signal-to-noise for ENCSR000EGE against ENCSR000DZG -- and add a third cell type to
distinguish "K562-trained models are portable in general" from "K562 happens to transfer to
GM12878". TeloHAEC is the cheapest third cell type, and it has no EP300, so this needs either
a DNase/H3K27ac proxy there or a fourth cell type with EP300.

**2. Get a downstream test of the working transfer direction, or accept that there is none.**
The CRISPR benchmark exists only for K562, so it can only ever run GM12878 -> K562, which is
the direction that fails. Every deployment claim therefore rests on correlation metrics unless
a benchmark in another cell type becomes available. This is a structural limitation worth
stating in any write-up rather than discovering later: the project's terminal metric cannot
evaluate its best transfer result.

Then, in order of expected value:

**3. Training-element composition, now largely closed.** GC-matching the negative pool was
the leading hypothesis for the sequence branch's flatness and it is refuted: +0.001 in-cell
and -0.002 transferred, with the match verified beforehand (see above). Two weaker ideas
remain in this family and neither is a priority: The negative pool is not GC-matched (Report 1): it is a uniform
random genome sample at mean GC 0.389 against 0.466-0.593 for candidate elements, and it is
10% of every batch. A sequence branch can satisfy that contrast with a coarse GC detector,
which is exactly what its behaviour looks like -- 0.397 correlation *across* candidate
elements and almost no dynamic range *within* them. Three experiments, cheapest first.

- **Harder negatives: accessible but unacetylated regions.** Sampling negatives from open
  chromatin rather than random genome removes accessibility as a discriminator too, forcing
  the sequence branch onto the distinction that actually matters downstream. Note these
  elements are already in the *positive* set with near-zero targets, so the honest framing is
  reweighting rather than relabelling -- up-weight the accessible-but-unacetylated quadrant
  and measure what happens to the over-predicted tail.
- **`negative_ratio` and the 50,000-window cap are untested.** Both were inherited from the
  p300 setup and neither has been swept. They set how much of the loss is spent on the
  positive-versus-background contrast rather than on grading elements against each other.

**4. Give the model the sequence features it is not learning.** Now the strongest surviving
idea in this family, though weakened: the GC-shortcut result says the branch's flatness is not
caused by what it was trained *against*, so explicit GC/CpG channels are a bet on
representation rather than on removing a confound. Explicit GC-content and CpG-density channels, computed by the same route as the
existing accessibility channels, with an indicator-channel control to confirm the gain comes
from the feature rather than from the extra capacity. Motif-derived channels are the heavier
version of the same idea. Run this after the composition experiments, since a GC channel and
a GC-matched negative set are two ways of attacking the same problem and the second is
cheaper.

**5. Attack the dynamic range directly.** The benchmark fails on compression, not on ranking
error, and no architecture change tried so far targets compression. Candidates: a loss on the
predicted *spread* over each batch, training on rank-transformed targets, or predicting the
composite `geomean(accessibility, H3K27ac)` that downstream consumers actually use. That last
one carries a design caution: if the accessibility input is the same assay as the
accessibility term in the target, half the target is readable off the input and the comparison
is circular, so it must be `geomean(DHS, H3K27ac)` predicted from sequence + ATAC and
baselined against the two-step route scored on the *same* composite.

**6. Motif syntax, which is the remaining scientific question and has never been started.**
SHAP / TF-MoDISco / FiNeMo on the **residual-trained** model, whose predictions are forced
onto signal accessibility cannot supply. Use the +/-500 bp window, which carries no neighbour
contamination. Expect less attributable signal than p300 (Fig.~3).

**Also open, smaller**

- **Run the dynamic-range diagnostic on the p300 arms.** `4.9` established that predicted
  H3K27ac fails through compression rather than inaccuracy. The same diagnostic on the p300
  arms would say whether p300 succeeds by being less compressed -- the mechanism behind the
  headline result is currently inferred, not measured.
- **Re-test the H3K27ac benchmark deltas with the paired bootstrap.** Those arms were read off
  overlapping unpaired intervals, which the p300 run demonstrates resolves nothing. The best
  H3K27ac arm was +0.025 over the floor unpaired, and +0.055 turned out to be resolvable for
  p300, so "no H3K27ac arm clears the floor" should be re-tested paired before it is treated
  as settled.
- **A qnorm-off ABC arm.** `run_qnorm` is rank-based so scale is irrelevant, but the reference
  is built from observed K562 signal, so a predicted track inherits the observed
  distribution's shape. Worth one arm.
- **Remove the profile head.** The base-resolution multinomial term is largely fitting Poisson
  noise (Report 1, Fig.~4). Removing it should be free and may help the counts head; untested.
- **Coarser profile resolution.** Binning the profile target to 20-50 bp raises its ceiling
  from 0.21 to about 0.72. The input stays single-base. Only worth doing if the profile head
  is kept.
- **Characterise the twelve threshold-level false negatives individually.** They are named in
  the diagnostic output; a per-element look may be more informative than the aggregate.
- **ATAC to DNase converter.** DNase tracks H3K27ac and enhancer activity better than ATAC
  does, so a converter would improve the input everywhere at once. Gate it on a cheap
  coupling comparison first, in the one cell type where both assays exist.
- **A third cell type and multi-cell-type training.** TeloHAEC is the only one available
  without new data processing, and its +/-IL1b / +/-TNFa / +/-VEGF conditions are a sharper
  and cheaper test than a new cell line: same genome, genuinely different regulatory state,
  no input-domain shift. Run the model-free coupling first (Report 2).
- **Repeat the residual-objective comparison for p300.** Resolved only for H3K27ac. p300's
  residual *r* is already 0.654, so the headroom may be smaller.
- **Train the p300 models with whichever multimodal architecture wins**, since the same
  accessibility over-reliance is likely present there.

**Closed, do not reopen without new evidence**

- **p300 as a second head or a stacked input** -- the target is not predictable where H3K27ac
  fails, and the premise was largely an input-control artifact.
- **Element derivation** -- ATAC-derived against DNase-derived training elements,
  *p* = 0.83.
- **Fragment channels for deployment** -- real in-cell type, exactly null transferred in both
  directions.
"""

R3_METHODS = """## Methods

Model is `MultiModalBPNet` (`scripts/multimodal_bpnet.py`), a BPNet variant with a
dilated-convolution trunk, a profile head trained with multinomial negative log-likelihood,
and a counts head trained with log1p mean-squared error; total loss is
`profile_loss + w * count_loss`. All results here concern the counts head unless a profile
metric is named. Baseline training uses `n_layers = 8` (a ~1.1 kb receptive field), a
+/-500 bp counting window with the required 2,114 bp input window, `w = 10`,
reverse-complement augmentation, chromosome-holdout 5-fold cross-validation, and a random
50,000-window subsample of the genome-wide negative pool, which is **not** GC-matched to the
positives (Report 1). Negatives are 10% of each batch and are never evaluated on. The wide arm uses `n_layers = 10`, whose trimming of
2,093 requires a 5,186 bp input window for the same output; the trimming is
`47 + sum(2^i for i in 1..n_layers)`.

`w = 10` was chosen from a single-fold sweep over 1, 3, 10, 100 and 1,000. Only the extremes
are resolvable at the measured noise level: 10 clearly beats 1 (0.437) and 1,000 (0.464),
while 3, 10 and 100 are within run-to-run variance of one another. With between-fold `sd` of
0.041-0.046, a single-fold sweep could not have resolved those differences at all, so this
choice rests on a flat region rather than on a measured optimum.

Scoring rules, metric definitions and the paired-test procedure are in Report 2; ceilings,
track definitions and the data inventory are in Report 1.

### Single-pass against reverse-complement-averaged, for every result where it matters

Test-time RC averaging became the scoring default on 2026-09-03. The decision table above
quotes the RC values; the figures and section prose quote single-pass. This is the complete
list of differences, so no reader has to guess whether a discrepancy is material.

| result, top quintile | single-pass | RC-averaged |
|---|---|---|
| sequence only, absolute | 0.380 | 0.397 |
| ATAC only, absolute | 0.548 | 0.551 |
| sequence + ATAC, absolute | 0.685 | 0.695 |
| wider receptive field, multimodal, K562 (paired) | +0.027 | +0.028 |
| wider receptive field, multimodal, GM12878 (paired) | +0.014 | +0.016 |
| fragment channels, K562 (paired) | +0.0135 | +0.0157 |
| residual *r*: sequence only | 0.149 | 0.148 |
| residual *r*: sequence + ATAC | 0.551 | 0.559 |
| residual *r*: residual-trained sequence | 0.459 | 0.470 |
| residual *r*: residual-trained multimodal | 0.514 | 0.529 |
| residual *r*: ATAC negative control | -0.003 | -0.002 |

Every RC value is equal or larger except the sequence-only residual *r*, and the negative
control stays at zero, which is what a change that removes strand-asymmetry noise rather than
adding signal should look like. No sign, ordering or significance verdict differs between the
two columns. Source: `results/{,rc_}{fiveprime,wide_*,fragchan_*,residual_grid_*}_*.tsv`.

**Why the figures were not migrated.** `3.1` and `3.4` read the `2.2`-evaluator table format,
for which no `rc_*` equivalent exists -- the RC re-scoring ran through `2.15`, which emits
different columns. Migrating them means either re-running the `2.2` evaluators with RC or
rewriting the plot scripts' data handling, both of which carry more risk of introducing an
error than the 0.002-0.016 they would correct. Recorded as open rather than done.

### Regenerate

```
Environment: pixi env `multimodal`, pixi.lock sha256:069960d36312
```

```bash
# Models: 3 modes x 5 folds, plus each architecture arm
for MODE in atac multimodal sequence; do for F in 0 1 2 3 4; do
    sbatch scripts/1.4.submit_training_5prime.sh $MODE $F        # K562 baseline
    sbatch scripts/1.8.submit_training_gm12878.sh $MODE $F       # GM12878 baseline
    sbatch scripts/1.9.submit_training_residual.sh $MODE $F      # residual objective
    sbatch scripts/1.13.submit_training_5prime_wide.sh $MODE $F  # wide receptive field
done; done
for F in 0 1 2 3 4; do
    sbatch scripts/1.15.submit_training_fragchan.sh multimodal $F   # fragment channels
    sbatch scripts/1.19.submit_training_gate.sh multimodal $F       # gate
done

# Evaluation: one code path, one config JSON per comparison
sbatch scripts/2.22.submit_rc_rescore.sh          # the rc_* tables the decision table quotes
sbatch scripts/2.23.submit_final_comparisons.sh   # factorial + element-derivation arms
sbatch scripts/2.24.submit_transfer_matrix.sh     # transfer matrix
sbatch scripts/2.25.submit_gate_factorial.sh      # gate 2x2

# Downstream: ABC + CRISPR benchmark, then the diagnostics
pixi run python scripts/4.3.setup_abc_arms.py && bash scripts/4.5.submit_abc.sh
pixi run python scripts/4.6.setup_crispr_comparison.py
bash scripts/4.7.submit_crispr_comparison.sh
sbatch scripts/4.9.submit_diagnose.sh
sbatch scripts/4.13.quantify_error_strata.sh

# Figures and this report
pixi run python scripts/3.7.build_numbers_manifest.py
python render_report.py report3_design_decisions.qmd
```

Per-result provenance, including which processing of the target each result file used, is in
`results/TARGET_PROVENANCE.md`.
"""


R1_FIXUPS = [
    # The source's Datasets row calls the negative pool GC-matched. It is not: the file is
    # ChromBPNet's GC-ANNOTATED tiling and the sampler discards the GC column. Corrected here
    # rather than in the superseded source, which is kept as the record of what was published.
    (r"\| \*\*GC-matched negatives\*\* \| Genome-wide GC-binned background, hg38 \|",
     "| **Training negative pool** | ChromBPNet genome-wide GC-*annotated* tiling, hg38; "
     "n = 3,088,298 bins, 1 kb stride, mean GC 0.389. **Not GC-matched to the positives** "
     "\u2014 see the negatives section below |"),
    (r"\| Training negatives \|", "| Training negatives only; never evaluated on |"),
]

R1_MAP2 = {8: 1, 2: 2, 3: 3, 11: 4, 9: 5}


def build(pull, write):
    p1 = lambda t, n=None: pull(t, R1_MAP2, n, R1_FIXUPS)
    p3 = lambda t, n=None: pull(t, R3_MAP, n, R3_FIXUPS)

    write("report1_data_characterisation.qmd",
          "What the H3K27ac and ATAC data are, and what is predictable from them",
          DATE,
          [R1_HEAD,
           p1("Datasets"),
           p1("Accessibility explains less of H3K27ac in GM12878, which accounts for the ATAC-only transfer drop",
              "Accessibility explains less of H3K27ac outside K562, and the strata disagree about the ordering"),
           p1("H3K27ac sits on flanking nucleosomes and is far broader than ATAC"),
           p1("A 1 kb counting window is the best trade-off between signal and neighbour contamination"),
           p1("H3K27ac has no reproducible base-resolution profile"),
           p1("The ATAC library supports fragment-size stratification"),
           R1_CONVENTIONS,
           R1_METHODS])

    write("report2_evaluation_methodology.qmd",
          "How a comparison in this project is scored, and what each metric can support",
          DATE, [R2_BODY])

    write("report3_design_decisions.qmd",
          "Architecture and input decisions: predicted p300 improves ABC in the training cell type, and not on transfer",
          DATE,
          [R3_HEAD,
           R3_INPUTS,
           p3("Sequence contributes less to H3K27ac than to p300, by a factor of 2.2"),
           p3("The residual objective helps only a sequence-blind input — replicated in two cell types"),
           p3("A wider receptive field helps only when accessibility is an input"),
           p3("Fragment-size-stratified accessibility adds real information"),
           p3("Test-time reverse-complement averaging is free and always helps"),
           p3("A sharper accessibility input helps the ATAC-only model and nothing else"),
           p3("Transfer is asymmetric: GM12878 is the harder cell type, in both directions"),
           p3("Deploying to a new cell type: transfer the multimodal model"),
           R3_TRANSFER,
           R3_P300_ABC,
           R3_GCMATCH,
           R3_GATE,
           R3_METHODS])
