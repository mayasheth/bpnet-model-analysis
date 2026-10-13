# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

### [2026-08-24] Center H3K27ac training windows on candidate elements, not ChIP peaks

**Context**: H3K27ac is deposited on the nucleosomes flanking a regulatory element, not on the nucleosome-free element itself. Measured meta-profile over 30k DNase candidate elements confirms this: shoulders at -275/+275 bp with a central dip, 3.48x over the distal plateau. The existing p300 pipeline centers windows on `start + summit` from a narrowPeak.

**Decision**: Train on `reference/K562_DNase_candidate_elements.narrowPeak`, centered on the element midpoint, with a counting window wider than the element.

**Alternatives considered**:
- Center on H3K27ac peak summits, rejected: the summit of a broad acetylation peak sits on a flanking nucleosome, so windows would be centered off the element and inconsistently so.
- Keep the element-sized window used for p300, rejected: a +/-250 bp window captures the shoulders but cuts deeply into real signal, which only reaches background around +/-2000 bp.

**Rationale**: The scientific question is about elements, so the element must define the coordinate system. Conveniently the candidate-elements file has its summit column set to `width/2`, so the trainer's existing `is_peak=True` path already centers on the midpoint with no code change.

**Consequences**: Negatives now mean "GC-matched genomic background" rather than "non-peak", and the element set spans active and inactive elements, so a large part of any reported correlation is the dead-vs-active contrast. This is why evaluation must be stratified by signal level.

**Tags**: h3k27ac, window-selection, element-centric, training-design

---

### [2026-08-24] Counting window is a trade-off between signal and neighbour contamination

**Context**: H3K27ac signal has no clean saturation point, it decays slowly and only reaches its distal plateau around +/-2000 bp. Meanwhile the 150,528 candidate elements cluster, so wide windows start containing other elements.

**Decision**: Test +/-500 and +/-1000 head to head rather than picking one. Expect +/-1000 to be the headline model and +/-500 the one trusted for attribution work.

**Alternatives considered**:
- A single wide window (+/-2000), rejected: the inter-replicate ceiling on active elements saturates by +/-1000 (0.785 vs 0.797 at +/-2000), so wider buys ~1 point of ceiling for 41% contamination.
- A single narrow window (+/-500), kept as an arm rather than the sole choice: zero contamination, but the ceiling is 2.4 points lower and it leaves real signal out.

**Rationale**: Neighbour contamination rises 0% -> 9.9% -> 19.9% -> 41.5% across +/-500/750/1000/2000, while the ceiling gains flatten after +/-1000. The two questions also want different things: prediction accuracy favours the wider window, while motif syntax is actively harmed by contamination, since a window containing a neighbour lets the model earn credit from the wrong element's motifs.

**Consequences**: Every window requires its own model geometry (`in_window = out_window + 2*557`), so this is a retraining sweep, not a re-evaluation. Attribution results from the +/-1000 model need to be read with the 20% contaminated subset in mind.

**Tags**: h3k27ac, window-selection, contamination, ceiling, trade-off

---

### [2026-08-24] Keep the profile head, down-weighted, rather than removing it

**Context**: H3K27ac counts are the quantity of interest, and the bimodal flanking profile is arguably not a meaningful "read profile" at all: the target is fragment-extended coverage, so MNLL (which expects multinomial read counts) sees ~250x inflated totals and near-identical adjacent positions. Profile MNLL sits at ~2800 here versus ~500 for the p300 runs.

**Decision**: Keep the profile head and control it with `count_loss_weight` rather than editing the loss to remove it.

**Alternatives considered**:
- Remove the profile head outright, deprioritized on evidence, not principle. At `count_loss_weight=10000` the counts term already takes ~91% of the gradient, making that run a close proxy for counts-only, and it scored 0.398 versus 0.407 at parity (weight 1000). Removal would land near 0.398, i.e. no better.
- Grid the profile-drop against mode and window, rejected: doubles the grid to 60 GPU jobs to test a lever the weight sweep already argues against.

**Rationale**: Reproducing the central dip is a useful diagnostic that the model learned H3K27ac structure rather than "accessible implies acetylated", and it costs nothing to keep once the weight is calibrated.

**Consequences**: `fconv` still consumes model capacity even with a small gradient, a second-order cost not measured. If ATAC and window width also plateau, revisit this alongside the more fundamental coverage-vs-read-counts mismatch.

**Tags**: h3k27ac, profile-head, loss-weighting, architecture, bpnetlite

**Update 2026-09-03, decision stands, and the reasoning was better than the evidence then
available.** Two measurements now exist. The inter-replicate ceiling on base-resolution
profile SHAPE is 0.21 on the top quintile at 1 bp, rising to 0.72 at 50 bp binning, so the
1 bp task is intrinsically near-unlearnable. Against that ceiling the head reaches
`profile_pearson` of 0.171 (multimodal), 0.144 (sequence) and 0.114 (ATAC only), so it
captures most of what is there. Separately, the `count_loss_weight` sweep still argues
against removal from the other direction: 10 gave 0.496 while 100 gave 0.467 and 1000 gave
0.464, so down-weighting the profile term made the COUNTS worse. Removal remains
deprioritized on evidence. The open question is now binning rather than removal.

---

### [2026-08-25] Residual correlation beyond ATAC becomes the headline metric

**Context**: ATAC alone predicts H3K27ac at 0.543 (top quintile), so a sequence+ATAC model at 0.668 looks strong while saying little about what sequence contributed. The deliverable is a model that takes ATAC + sequence and predicts activity, which makes the *departure from the ATAC expectation* the quantity of interest, not overall correlation.

**Decision**: Report `residual_pearson = r(observed - atac_pred, model_pred - atac_pred)` as the headline, with incremental R^2 and stratification by |true residual|. The baseline is the ATAC-only MODEL's held-out prediction, not the raw ATAC track, so the residual is what accessibility genuinely cannot explain. Implemented in `2026_0824_H3K27ac_model/scripts/2.4.evaluate_residual.py`.

**Alternatives considered**:
- Keep overall Pearson as the headline, rejected: it ranked the two targets backwards (see F-002).
- Partial correlation controlling for `atac_pred`, equivalent in spirit but less directly interpretable; the residual form states plainly "does the model's departure track the true departure".
- Use the raw ATAC track as the baseline, rejected: then the residual includes everything a linear read of the track misses, which flatters any model that merely learns a better ATAC transform.

**Rationale**: Overall correlation on this element set is dominated by the dead-vs-active contrast, which accessibility resolves on its own. The residual isolates the increment that motivates having a sequence model at all.

**Consequences**: Earlier reported numbers (all-elements and top-quintile Pearson) remain valid but are no longer the headline. Any future model comparison in this project must report the residual metric or it is not comparable to these results. Also implies a training change: see the residual-training item in `todo/TODOLIST.md`.

**Tags**: h3k27ac, evaluation, residual, metric-choice, atac

---

### [2026-08-29] Paired-end H3K27ac targets use read 1 only, not both mates

**Context**: TeloHAEC H3K27ac is paired-end; K562 and GM12878 are single-end. A plain
`bedtools genomecov -5` on PE data marks TWO 5' ends per fragment, one at each end, which is
not the quantity the SE tracks measure. Both variants were built and compared
(`scripts/0.16.compare_pe_5prime_variants.py`, job 41258528) across all four TeloHAEC
conditions, rather than assuming.

**Decision**: Build PE H3K27ac 5' targets from **read 1 only** (`samtools view -f 64`), so
each fragment contributes exactly one 5' end. Applies to TeloHAEC now and to any PE H3K27ac
entering the panel later. Both variants remain on disk (`*_r1_5p_*`, `*_both_5p_*`), so this
is reversible at no compute cost.

**Alternatives considered**:
- Count both mates, rejected on principle despite winning the top-quintile ceiling by
  0.003-0.005 in all four conditions. Two reasons. (i) It is a different quantity from the SE
  cell types, giving TeloHAEC a different count-to-molecule relationship. (ii) Each fragment
  contributes two *correlated* counts, violating the independent-read assumption behind
  bpnetlite's multinomial profile loss -- the same class of defect recorded as the reason for
  rejecting the 250 bp fragment-extended track.
- Decide from the ceiling numbers alone -- rejected: the differences are ~0.003-0.005
  corrected, far below the 0.041-0.046 between-fold sd that governs every model comparison
  here, so no model result could resolve them. The choice is measurably immaterial to
  performance and should therefore be made on assumptions, not on a tiebreak.

**Rationale**: The two variants are *not* distinguishable by signal quality. Normalised
meta-profiles agree at r = 0.9999 with max deviation below 0.9% of peak height, identical
shoulder positions and dip depth -- because R1 is equally likely to be the left or right end
of a fragment, so sampling one end recovers the same spatial distribution at half the
density. There is no smearing here, unlike fragment extension. With signal quality tied, the
deciding criteria are comparability with the existing cell types and not breaking the loss
function's assumptions, both of which favour r1.

**Consequences**: TeloHAEC counts will be ~half those of a both-mates track, which is
irrelevant to correlation but must be remembered if raw counts are ever compared across cell
types. The `*_both_5p_*` tracks stay on disk and must not be mixed into an analysis with r1
tracks. A `both_depthmatched` column exists in
`results/telohaec_*_pe_variant_ceiling.tsv` and is **not** a fair control -- matching total
counts halves the fragment count for `both`, because each fragment there contributes two
marks, so it is a half-depth library and loses for that reason alone. The script now says so.

**Tags**: h3k27ac, telohaec, paired-end, target-definition, profile-loss, panel, 5-prime

---

### [2026-08-30] Accessibility inputs should be 5'-end insertion counts (ChromBPNet convention), and ours are not

**Context**: The TeloHAEC coupling run exposed an accessibility-scale mismatch. K562 ATAC
mean coverage is 17.6 and GM12878's 18.3, but TeloHAEC's is 2.7-3.7 -- a 5-7x gap. Tracing
it: our `atac.bw` files are built by `bedtools genomecov -bg` over the FULL tagAlign
interval, so coverage is proportional to read length. K562 tagAlign entries are 94-95 bp;
TeloHAEC's are 35-36 bp. So roughly 2.6x of the gap is read length and ~2x is depth.
K562 and GM12878 happen to match each other closely, which is why this never surfaced
during their transfer work.

Checked against ChromBPNet (`chrombpnet/helpers/preprocessing/reads_to_bigwig.py`), verbatim:
```
awk (strand-specific shift) | bedtools genomecov -bg -5 -i stdin -g <genome>
ATAC:  plus_shift_delta, minus_shift_delta = 4-plus_shift, -4-minus_shift
DNase: plus_shift_delta, minus_shift_delta = -plus_shift, 1-minus_shift
```
Single-base 5' ends, unstranded, `-5` applied after the shift. Their pipeline auto-detects
pre-existing shifts, so an already-Tn5-shifted tagAlign gets zero additional shift.

**Decision**: **5'-end insertion counts are the standard for accessibility inputs going
forward** (Maya, 2026-08-30). Because our tagAligns are already Tn5-shifted
(`*.tn5.sorted.tagAlign.gz`), the correct build reduces to `bedtools genomecov -bg -5` with
no further shift. Build to NEW filenames (`atac_5p.bw`), never overwrite `atac.bw`.

**Alternatives considered**:
- Keep full-interval coverage and rely on z-normalisation -- rejected. Normalisation fixes
  mean and sd but not read-length-dependent spatial spread, and measurably not dynamic range:
  peak-to-background is ~9 for TeloHAEC against ~6 for K562/GM12878. Applying a K562-trained
  model's saved normalisation stats to TeloHAEC would put the input outside the training
  distribution for a purely technical reason, and it would present as a transfer failure.
- Truncate all tagAligns to a common width -- rejected: throws away real data and still
  depends on the arbitrary common width.
- Rebuild in place -- rejected: `2026_0529_multimodal_p300_model/data/atac.bw` is shared with
  the p300 models, so overwriting would silently invalidate those results too.

**Rationale**: This is the same defect, one level up, as the 250 bp fragment extension
rejected for the H3K27ac target: a read-length-dependent smear that breaks comparability
across samples. Matching the field-standard tool also means our accessibility input is
directly comparable to published ChromBPNet work.

**Consequences**: **Every ATAC-input model is affected and they must be redone as a set, not
piecemeal.** That is `atac5p_hw500_clw10`, `multimodal5p_hw500_clw10`,
`residual5pFIXED_hw500_clw10` (sequence input, but its offset comes from an ATAC model),
`gm12878_atac5p_*`, `gm12878_multimodal5p_*`, and the 10 residual runs launched 2026-08-30 --
roughly 35 fold-jobs plus re-evaluation. Sequence-only models are unaffected. Mixing a 5'
ATAC input with any existing number is invalid, so until the set is rebuilt, current results
stand as-is on the old input and must be quoted as such. The p300 models also use the shared
`atac.bw` and would need the same treatment before any p300/H3K27ac comparison is remade.

**Not yet decided**: whether to rebuild now (invalidating the in-flight residual grid) or
after the residual grid completes on the current input, so the grid comparison stays
internally consistent. Current course is the latter.

**Tags**: atac, accessibility, chrombpnet, 5-prime, input-definition, read-length,
comparability, transfer, telohaec

### [2026-09-03] Adopt the wider receptive field for the multimodal model, reversing the earlier verdict

**Context**: `n_layers` 8 gives a ~1.1 kb receptive field against an H3K27ac extent of several kb. The sequence-only arm of the `n_layers` 10 experiment finished first and was null on the top quintile (−0.006 K562, +0.010 GM12878, both p=0.53), which was recorded as "resolved, do not pursue".

**Decision**: Adopt ~4.2 kb for the multimodal model. The multimodal arm gains +0.027 (K562, p=0.006) and +0.014 (GM12878, p=0.025) on the top quintile, all five folds rising in both cell types, with the accessibility residual up from 0.502 to 0.547 and 0.397 to 0.469.

**Alternatives considered**:
- Keep `n_layers` 8, rejected: the gain is significant, replicated, and on the reporting standard.
- Adopt it for sequence-only too, rejected: sequence gains only on all elements, which is the dead-vs-active contrast.
- Wait for the ATAC-only arm before deciding, the arm is running and will refine the mechanism, but it cannot overturn a replicated multimodal gain.

**Rationale**: The useful long-range information is in the accessibility track. A wider window lets the model read accessibility over a larger neighbourhood; sequence-only cannot exploit 4 kb of sequence, so for it the extra window is noise.

**Consequences**: The deployed model geometry changes to trimming 2093 / in-window 5186, ~2-3x the training time per fold. Transfer and deployment comparisons were run at `n_layers` 8 and should be re-checked. Recorded reversal: the earlier decision generalised from one arm of a two-arm experiment.

**Tags**: architecture, receptive-field, multimodal, reversal, accessibility

---

### [2026-09-03] Fragment-size accessibility channels are 5' insertion counts and include the flat track

**Context**: Fragment length distinguishes nucleosome-free from nucleosome-occupied DNA and flat coverage discards it. An earlier pair of channels existed but used `genomecov -bg` over the full fragment interval, the read-length-dependent smear already removed from the flat track, and covered only half the fragments.

**Decision**: Five channels `[all, sub <=139, mono 140-329, di 330-620, poly >=621]`, every one a single-base Tn5 insertion count, with `all` being `atac_5p.bw` itself. Bin edges sit in the troughs of the measured fragment-length distribution.

**Alternatives considered**:
- Full-fragment coverage for the stratified channels, rejected: it would place two incompatible accessibility conventions in one input tensor, and a mono-nucleosomal interval marking occupancy directly is not worth that.
- The original `sub`/`mono` edges, rejected: they sat inside the modes and covered half the fragments.
- Stratified channels only, without the flat one, rejected: including `all` is what makes the input a provable superset, since the first convolution can zero the four bins and reproduce the baseline exactly.

**Rationale**: The bins partition the flat track exactly, zero discrepancy across 545,661,218 insertions, so any difference is added information rather than a changed input. The Tn5 shift was measured against the existing track (r = 1.0000 at +4/−5) rather than assumed, which also proved the PE BAMs and tagAligns hold the same reads.

**Consequences**: Gains +0.0135 on the top quintile (p = 0.0034). `atac_sub.bw` and `atac_mono.bw` are superseded and must not be used with a 5' model. GM12878 replication needs its PE BAMs downloaded, since fragment length lives in TLEN.

**Tags**: atac, fragment-length, accessibility, input-design, superset, tn5

---

### [2026-09-03] Inject predicted H3K27ac into ABC as a painted bigWig, with qnorm left on

**Context**: The ABC activity term is `geomean(accessibility, H3K27ac)` computed from read counts over candidate regions. Our model emits a per-element scalar from 5' end counts in a +/-500 bp window, which is not obviously commensurate with read counting over an element.

**Decision**: Write a bigWig with each candidate region painted at `predicted_counts / width` and pass it in the `H3K27ac` column. Keep `use_qnorm: True`. Assemble the genome-wide track from the five fold models, each applied only to the chromosomes it held out.

**Alternatives considered**:
- Patch ABC to accept a precomputed activity column, rejected as unnecessary: `neighborhoods.py:count_bigwig` already sums bigWig values per region, so a painted track is counted exactly as a real one.
- Turn qnorm off, rejected: it is what makes the injection scale-free, and the ABC score thresholds are calibrated on qnorm'd values.
- Paint the prediction itself, so the region sum scales with width as real read counts do, kept as `--paint density`; the default makes ABC's sum recover the model's prediction exactly, since the model predicts a fixed +/-500 bp window regardless of element width.
- A five-model ensemble for the cross-cell-type arms, where leakage is not a concern, rejected: it would give those arms an ensembling advantage the same-cell-type arms cannot have.

**Rationale**: `run_qnorm` is rank-based, mapping each region's within-sample quantile onto the K562 reference, so only the rank order of the injected values matters. That dissolves the units mismatch. Validated on chr22: predicted-vs-observed H3K27ac Spearman 0.716 all regions and 0.550 top quintile, against 0.633 and 0.469 for raw ATAC, so the prediction is a better proxy than the accessibility it would replace.

**Consequences**: ABC scores our models on the ATAC-derived candidate regions while they were trained on the DNase-derived set, a train/score element mismatch inside the comparison; a matched run is training. Every arm must share one region set, so `Peaks/` is pre-populated from the completed July run and must be copied after the bigwigs exist or Snakemake re-runs region calling.

**Tags**: abc, crispr-benchmark, activity, qnorm, injection, leakage

---

### [2026-09-03] Regression gates on inference code compare within a tolerance, never byte-identically

**Context**: `2.15` is shared with the transfer and residual-grid results, so generalising it needed a gate. The first gate required the per-fold table to be byte-identical to a stored one, and it failed on differences in the 4th decimal place while region counts matched exactly.

**Decision**: Assert exact equality only on deterministic quantities, fold set, config labels, `n` per fold, and the presence of every reference column, and compare metrics within 1e-3.

**Alternatives considered**:
- Keep byte-identity and pin the GPU model, rejected: it makes the gate depend on scheduling.
- Drop the gate, rejected: it is the only thing standing between a refactor and silently changed published numbers.

**Rationale**: cuDNN convolution is not bit-reproducible across GPU models, so re-scoring the same weights on a different node moves the 4th decimal place with the code untouched. The tolerance is chosen from the science's noise floor, between-fold sd is 0.041-0.046, rather than from float precision, so a real behavioural change moves numbers by far more than the gate absorbs.

**Consequences**: `2.18.compare_perfold_tables.py` implements it and states the tolerance and its justification in its own docstring, so the next reader does not tighten it back to zero. Added columns are reported and skipped; missing columns fail.

**Tags**: testing, regression, determinism, tolerance, evaluation

---

### [2026-09-03] Test-time reverse-complement averaging stays opt-in

**Context**: Averaging each prediction with its reverse complement at inference costs one extra forward pass and no retraining, and it improves every model: +0.0162 sequence only, +0.0081 multimodal, +0.0016 for ATAC only, which is the negative control and the only non-significant row.

**Decision**: Ship it behind `--rc-average`, off by default.

**Alternatives considered**:
- Make it the default immediately, rejected for now: it would shift every number in the report by a small amount, so past and future tables would not be comparable unless all are re-scored together.
- Leave it unimplemented, rejected: it is the cheapest measured gain on the project.

**Rationale**: The gain is real but small relative to the comparisons being made, and comparability across the report matters more than a few thousandths until the tables are re-scored as a set.

**Consequences**: Adopting it later means re-running every evaluation config in one batch. Listed as an explicit adopt-or-not decision in the TODO.

**Tags**: inference, reverse-complement, adoption, comparability

---

### [2026-09-03] Split the analysis into three reports along stability, not topic

**Context**: The single report has reached 17 sections and 13 figures, mixing stable characterisation with an actively churning set of architecture experiments.

**Decision**: Three documents, ordered by how often they change.

1. **Data characterisation**, what the data is and what is predictable in principle. H3K27ac position and width, the counting-window trade-off, both ceilings (counts and profile shape), ATAC fragment-length structure, ATAC-H3K27ac coupling across cell types, element derivation and panel caveats, ATAC input conventions. Every ceiling lives here and the other two reports cite it.
2. **Evaluation methodology**, short, concrete, framed as standing cautions with the incident that motivated each. Why the top quintile leads, what the residual metric means and its artifact controls, paired within-fold testing against a between-fold sd of 0.041-0.046, the transfer-versus-deployment evaluation DESIGN, and ABC/CRISPR as the downstream metric.
3. **Design decisions**, the workbench. Spine is a decision table of change / effect on the top quintile / verdict, with evidence below: input modality, residual versus total objective, receptive field, fragment channels, RC averaging, target definition, window size, `count_loss_weight`, the profile head, element derivation, and the transfer RESULTS with the deployment verdict.

**Alternatives considered**:
- Keep one report, rejected: the churn in category 3 forces re-reading stable material to find what changed.
- Split by audience, a short headline plus a technical appendix, rejected for a working project: it optimises for a reader who is not the one using this daily.
- Put transfer entirely in report 2 or entirely in report 3, rejected: transferability is both an evaluation axis and a deployment decision, so the design goes in 2 and the results in 3.

**Rationale**: Splitting on churn rate isolates the part that moves. Reports 1 and 2 become citable references; report 3 is expected to change every week.

**Consequences**: Figure numbering restarts per report, so every legend and prose cross-reference is touched, mechanical, and the numbers manifest lint catches any reference that stops resolving. `render_report.py` resolves `outputs/numbers.json` relative to the report, so all three sharing one directory keeps one manifest. Deferred until the runs in flight land, so the split happens once.

**Tags**: reporting, organisation, documentation, churn

---

### [2026-09-05] Drop both p300 ideas: no second head, no stacked p300 input

**Context**: The elements the H3K27ac model under-predicts overlap EP300 peaks at 57.7%
against 8.8% in the typical stratum, a 6.6x enrichment. That suggested either a multi-task
model with an H3K27ac head and a p300 head, or using predicted p300 as an input feature.
Option (b) -- do p300-target models predict the elements our H3K27ac model misses? -- is the
cheap prerequisite for both and was run first, with the pass threshold written into the
script before it ran.

**Decision**: Drop both. Observed p300 is elevated 4.89x at those elements but its own input
control is at 2.10x, so real enrichment is ~2.3x rather than 6.6x; the existing p300 models
predict 1.83x (multimodal) and 1.55x (ATAC-only), below the control. The p300 models also
over-predict the accessible-but-unmarked tail (1.92x where observed is 1.36x and its control
0.50x), so they fail in the same direction on the same elements.

**Alternatives considered**:
- Build the multi-head model anyway -- rejected: ~15 GPU jobs to learn a target that is not
  predictable where it would need to be.
- Use p300 as a model input -- rejected earlier on deployment grounds (the application must
  run on any cell type with only ATAC) and now on learnability grounds as well.

**Rationale**: A second target only helps if it is *predictable* exactly where the first one
fails. The peak-overlap premise was substantially an accessibility artifact, which is why the
input control was carried through every stratum.

**Consequences**: The separate finding that p300 is the better substrate for sequence
attribution is unaffected -- that governs where motif work should be done, not what to
predict. Recorded in report 3.

**Tags**: p300, multi-task, negative-result, controls

### [2026-09-05] Do not carry fragment channels into deployment

**Context**: Five fragment-size accessibility channels gain +0.016 on the top quintile in
K562 (*p*=0.002) and +0.006 in GM12878, replicated, with the channel partition verified exact
against the flat track. The transfer matrix scored them across cell types for the first time.

**Decision**: Keep them for in-cell-type work, do not carry them into any deployment
configuration. Transferred, they are -0.0012 (K562 -> GM12878) and -0.0029 (GM12878 -> K562)
against flat accessibility -- exactly null in both directions.

**Alternatives considered**:
- Adopt everywhere on the strength of the in-cell-type result -- rejected: in-cell-type
  ranking is not the deployment ranking, which is now measured rather than assumed.

**Rationale**: Fragment-length structure is a property of a particular ATAC library as much
as of the chromatin, and no paired in-cell-type test can separate those. The transfer test
can, and says the useful part does not travel.

**Consequences**: The wider receptive field is the only architecture change worth carrying
into deployment, and even it is only +0.0115 and +0.0053 transferred (*p*=0.14, 0.53). The
honest summary of the whole architecture programme is that its gains are in-cell-type gains.

**Tags**: fragment-channels, transfer, deployment, negative-result

### [2026-09-05] Close the gating / loss-reweighting direction; make training-element composition the next experiment

**Context**: The multimodal model over-predicts H3K27ac at accessible, CTCF-enriched,
CpG-rich elements. Two cheap interventions were built against that: a sequence gate on the
accessibility representation, and a 3x over-prediction penalty in the count loss. Both were
trained as a 2x2 factorial against a matched baseline. Separately, the corrected
responsiveness analysis showed the sequence branch does not separate the two error tails at
all (1.55x against 1.44x, where truth separates them 1.00x against 31.3x).

**Decision**: Close both. In-cell K562 all three arms are null; transferred to GM12878 the
gate is null, the loss reweighting trends negative, and the combination is significantly
worse (-0.0128, *p*=0.039). Make **training-element composition** the next experiment
instead, starting with actually GC-matching the negative pool -- which was never matched
(see learnings, same date), and which plausibly explains why the sequence branch behaves like
a coarse GC detector.

**Alternatives considered**:
- Tune the gate (different initialisation, wider kernel, per-element rather than per-position)
  -- rejected: a gate can only change how existing information is combined, and the branch it
  gates has been measured to carry no discriminative signal on these elements.
- Raise the over-prediction weight above 3.0 -- rejected: the transferred trend is already
  negative and the mechanism argument predicts it gets worse.
- Go straight to explicit GC/CpG input channels -- deferred, not rejected: a GC-matched
  negative set and a GC input channel attack the same problem and the former is cheaper and
  changes no architecture.

**Rationale**: The two interventions shared one premise, and that premise is now measured to
be false. Composition is upstream of all of it: if the positive/negative contrast is
separable on GC alone, the sequence branch is being trained to be exactly what we observe.

**Consequences**: `--gate-accessibility` and `--overprediction-weight` stay in the code
behind flags with their regression test, since they are harmless when unset and the negative
result is worth being able to reproduce. Report 3 records the numbers and the prediction that
preceded them.

**Tags**: architecture, gating, loss-design, training-composition, negatives, negative-result

### [2026-09-06] Make p300 the primary modelling target; train a GM12878 p300 model next

**Context**: Predicted H3K27ac never cleared the ABC benchmark floor. Substituting predicted
p300 does, by +0.055 AUPRC [+0.034, +0.074] on a paired bootstrap, matching *measured*
H3K27ac. Observed p300 also beats observed H3K27ac as the activity term by +0.038
[+0.019, +0.058], so the target is better independently of any model. The gain requires
sequence (+0.051 over the ATAC-only p300 model), which is the mirror image of H3K27ac, where
the sequence arms were the worst in the panel.

**Decision**: p300 becomes the primary target for the predicted-activity work. The next
experiment is a GM12878 p300 model and a transferred benchmark run; everything else drops
below it in priority. EP300 data exists in GM12878 (ENCSR000DZG), so this is five folds of
training on an existing scoring path.

**Alternatives considered**:
- Keep optimising the H3K27ac model -- rejected: same architecture, same accessibility input
  and same benchmark give +0.055 for p300 and nothing resolvable for H3K27ac, so the target
  was the binding constraint, not the architecture.
- Declare the transferability question answered from the K562 result -- rejected: this is
  K562-trained and K562-evaluated. F-005 showed in-cell-type architecture gains do not travel,
  and there is no reason to assume a target-choice gain travels either.
- Predict both and combine -- deferred: F-006 killed the multi-head design on learnability
  grounds and nothing here reopens it.

**Rationale**: p300 is one step closer to the sequence-specified events, which is also why
F-001 found sequence contributes 2.2x more to p300 than to H3K27ac. That earlier finding was
treated as guidance for *attribution* work only; it turns out to have predicted which target
would produce a usable predicted-activity track, which is worth remembering as a case where a
mechanistic result had a practical consequence we did not act on for six weeks.

**Consequences**: H3K27ac work is not wasted -- the ceilings, window choice, RC averaging,
receptive field and the whole evaluation methodology carry over unchanged, and the p300 models
share the negative pool so the GC-matching fix applies to both. The H3K27ac benchmark deltas
need re-testing with the paired bootstrap before "no arm clears the floor" is treated as
settled.

**Tags**: p300, target-choice, abc, crispr-benchmark, transferability, positive-result

### [2026-09-06] Use a paired bootstrap for every CRISPR-benchmark comparison

**Context**: `performance_summary.txt` reports an independent bootstrap CI per predictor.
Those are ~+/-0.05 wide while the differences of interest are 0.03-0.10, so every arm's
interval overlaps every other's and no comparison is possible from them. All arms are scored
on the identical element-gene pair set.

**Decision**: Compare arms with `scripts/4.17.paired_auprc_bootstrap.py` -- resample the pair
set once per iteration and recompute every arm on that same resample. Report the delta, a
percentile CI and the fraction of resamples preserving the sign. Never compare two arms by
eye from the unpaired intervals.

**Alternatives considered**:
- Read the unpaired intervals and describe overlapping arms as indistinguishable -- rejected,
  and this is not hypothetical: it is what produced the H3K27ac reading, and the p300 run
  shows a +0.055 difference sitting inside two overlapping unpaired intervals is in fact
  resolvable at 100% sign retention.

**Rationale**: Identical to the within-fold pairing already standard for the model metrics
(Report 2): shared sampling variance dominates, and pairing removes it. It was applied
upstream and not downstream purely by oversight.

**Consequences**: Report 3 quotes paired deltas for every benchmark comparison and keeps the
unpaired table only to show that it resolves nothing. `4.17` also revealed a constant +0.0105
offset between its own average-precision implementation and the pipeline's across all arms --
an estimator difference, not a data difference, which is why deltas rather than absolute
values are the quotable output.

**Tags**: statistics, paired-testing, crispr-benchmark, abc, bootstrap

### [2026-09-06] Close training-element composition: GC-matching the negatives changes nothing

**Context**: Every model in this repository -- H3K27ac and p300, both cell types -- was
trained with ChromBPNet's genome-wide GC-*annotated* tiling passed straight to `--negatives`.
That file is the matching input, not its output; the trainer reads three of its four columns
and samples uniformly, so the pool sat at mean GC 0.389 against 0.466-0.593 for candidate
elements while making up 10% of every batch. This was the leading hypothesis for why the
H3K27ac sequence branch has almost no dynamic range within candidate elements -- a coarse GC
detector would satisfy most of that contrast.

**Decision**: Closed. Retraining the multimodal model on GC-matched negatives built with
`bpnet-gc-background` gives **+0.0011 [-0.0080, +0.0102]** in-cell type (*p*=0.76) and
**-0.0024 [-0.0083, +0.0036]** transferred to GM12878 (*p*=0.33). Absolutes 0.700 against
0.699 and 0.529 against 0.531. The transferred model still beats the target's own ATAC-only
model by the same non-significant ~0.011.

**Alternatives considered**:
- Blame the residual mismatch -- rejected as the primary reading, though recorded: the matched
  pool undershoots the GC-rich tail (p90 0.570 against the elements' 0.590) because a
  genome-background pool lacks enough very-high-GC windows. But closing 85% of the mean-GC gap
  and 80% of the largest per-bin discrepancy moved the metric by 0.001, so any surviving
  shortcut is too small to chase.
- Reimplement matching to reach the high-GC tail -- rejected: it would require negatives drawn
  from accessible-but-unacetylated regions, which is the separate reweighting idea, and the
  effect size here gives no reason to expect it to matter.

**Rationale**: The match was verified before the result was read (`scripts/0.29`, joining the
matched pool back to the tiling's own GC annotation so both sides use one estimator), and the
transfer arm -- where a GC shortcut should cost most, since accessibility is the part that
generalises -- is the arm that moved least. A null with a verified manipulation and a
pre-committed failure condition is a real answer, not a failed intervention.

**Consequences**: The cheapest remaining fix for the sequence branch is gone, and the two
survivors in this family (explicit GC/CpG channels, accessible-but-unacetylated reweighting)
are both weaker a priori than the one just refuted -- the GC-shortcut result says the flatness
is not caused by what the branch was trained *against*, so a GC channel is now a bet on
representation rather than on removing a confound. `--gc-match-negatives` was NOT added to the
trainer: the matched file drops in through the existing `--negatives`, and a second matching
implementation would risk diverging from `bpnet-gc-background`.

**Tags**: negatives, gc-matching, training-composition, sequence-branch, negative-result


### [2026-09-10] Train the ATAC-to-DNase converter on ATAC-derived candidate elements

**Context**: Every other K562 model in this project trains on
`K562_DNase_candidate_elements.narrowPeak`. The converter predicts DNase from sequence and
ATAC, so that element set is chosen by the signal the model is trying to predict.

**Decision**: Train the converter on `K562_ATAC_candidate_elements.narrowPeak` (`1.24`), and
keep the DNase-derived set for the downstream H3K27ac arms (`1.26`) so those stay comparable
to their own anchors, `1.11` and `1.22`.

**Alternatives considered**:
- DNase-derived elements, matching every other model: rejected. At deployment the element set
  can only come from ATAC, so a DNase-derived set means the converter never trains on the
  regions ATAC nominates and DNase does not, which is exactly where a conversion does work
  rather than copying.
- Both, as a fifth arm: rejected as unjustified cost, given D-18's neighbour result.

**Rationale**: Element derivation, ATAC- against DNase-derived, made no difference on the
H3K27ac task (*p*=0.83), so the swap is free in measured terms and removes the circularity.

**Consequences**: The GM12878 converter (`1.25`) could NOT match this. `0.26` derives K562's
ATAC element set from a specific rE2G ATAC run with no GM12878 counterpart, so that arm
trains on the DNase-derived set and the two converters differ in element derivation. If the
GM12878 arm underperforms, that asymmetry is a live alternative to a transfer explanation and
has to be reported as one.

**Tags**: converter, elements, circularity, deployment, atac, dnase

### [2026-09-10] Score the converter's profile with the raw estimator, never `profile_pearson`

**Context**: `2.15`'s `profile_pearson` comes from bpnetlite's
`calculate_performance_measures` with `kernel_sigma=7, kernel_width=81`, so both profiles are
Gaussian-smoothed before correlating. The DNase ceiling that unblocked the converter (0.848 at
1 bp top quintile) is a RAW within-element correlation from `0.25`. The two are different
estimators on different scales.

**Decision**: `2.31` scores the converter with `0.25`'s `shape_corr` verbatim and computes the
inter-replicate ceiling on the SAME held-out windows in the same run. Report the bin-size
curve, not one number.

**Alternatives considered**:
- Quote `profile_pearson` against the `0.25` ceiling: rejected. On the same fold-0 model the
  smoothed metric reads 0.774 where the raw 1 bp value is 0.799 and the ceiling is 0.831, so
  the numbers are close enough to look comparable and are not.

**Rationale**: The converter exists to feed a branch that reads base resolution, so the
question is at what resolution the track is trustworthy. One number cannot answer that, and a
smoothed number answers a different question.

**Consequences**: F-011 rests on the raw estimator. `0.35` was added alongside, because
fraction-of-ceiling says how close the converter gets to DNase and not whether it beats the
ATAC track it would replace, which is the lower and more relevant bar.

**Tags**: converter, profile, metric-definition, ceiling, estimator

### [2026-09-11] Put every accessibility benchmark arm on one counting path

**Context**: The July run scored real ATAC (0.457) and real DHS (0.576) as ABC activity terms
from read files, which go through `count_bam`/`count_tagalign`. A painted prediction goes
through `count_bigwig`.

**Decision**: Run all four accessibility arms as 5-prime bigwigs, including the two real
tracks, and treat the July numbers as external reference points rather than arms.

**Rationale**: Comparing a converted track against a read-file arm would put a counting-path
difference inside the contrast the run exists to measure.

**Consequences**: Measured rather than assumed, and it mattered. Real ATAC as a bigwig scores
0.4572 against the July tagAlign arm's 0.4573, identical. Real DNase as a bigwig scores 0.5280
against the July DHS arm's 0.5763. The counting path leaves ATAC alone and moves DNase by
0.048, so quoting the July 0.576 beside these numbers would have inflated the real-DNase arm
and overstated the gap the converter had to close.

**Tags**: abc, counting, benchmark, confound, dnase, atac

### [2026-09-12] Zero painted accessibility below one read

**Context**: `4.20` paints `softmax(profile) x counts` per base. A softmax never emits zero, so
the painted track is dense where a real 5-prime track is 71.3% empty. That floor made the
genome-wide track 22 GB against real DNase's 828 MB and inflated the lowest observed-signal
quintile 2.7x.

**Decision**: `4.22 --min-value 0.5`. A real 5-prime track holds integer read counts, so a
painted base carrying less than one read is below what the real track can represent.

**Alternatives considered**:
- Leave it dense: rejected. The raised floor would have been a confound in the downstream arm,
  indistinguishable from the profile being useless.
- Threshold at 1.0: rejected on measurement. It overshoots real DNase's sparsity (86.7% zero)
  and starts costing shape.
- Quantile-map the painted track to a reference accessibility distribution: held in reserve.
  It is deployment-legal and would fix the residual distortion, but is only worth doing if the
  as-is arm comes out negative.

**Rationale**: Measured on chr22 before switching on. At 0.5 the sparsity reaches 74.5%
against real DNase's 71.3%, the 1 bp top-quintile shape correlation moves 0.7814 to 0.7793,
and the lowest-quintile inflation falls from 2.73x to 1.54x.

**Consequences**: Residual Q1 inflation of 1.54x remains and is the named alternative
explanation if the converted arm underperforms. The thresholded track's genome-wide total is
103.6M against real DNase's 301.1M; training computes its own accessibility normalisation, so
a uniform factor is absorbed, but the distortion is not uniform.

**Tags**: converter, painting, threshold, dynamic-range, confound


### [2026-09-14] Stop work on the ATAC-to-DNase converter as a model input

**Context**: F-011 established the converter reproduces the DNase profile at 95% of its
inter-replicate ceiling. F-012 showed it buys nothing as an ABC activity term, which scores
counts and discards the profile. F-013 showed it is resolvably WORSE than plain ATAC as a
model input, where the profile can matter, and diagnosed magnitude distortion as the likely
cause. F-014 tried two independent fixes for that and neither moved the number.

**Decision**: Stop. No third transform, no retraining of the converter, no further downstream
arm. The converted track is not usable as an accessibility input and the reason is not the
magnitude defect that was blamed.

**Alternatives considered**:
- A third post-hoc transform: rejected. Quantile mapping (0.37) fixed the marginal exactly and
  made dynamic range worse (18x -> 9x against real DNase's 39x); a power transform fitted
  out-of-cell (0.39) raised range to 24x at a shape cost of 0.006 and moved the downstream
  delta from -0.0099 to -0.0105. Two attacks on the same hypothesis, both null, is enough.
- Retraining the converter with a count-focused objective: rejected for now, and recorded
  rather than closed. Its counts are already top-quintile r 0.878, and the part of DNase
  counts carrying the benefit is by construction the part that differs from ATAC counts.
- Better shape: there is nowhere to go. 95% of the inter-replicate ceiling is the ceiling.

**Rationale**: The converter's measured strength is shape, and F-014 shows shape is not even
the operative property in every cell type: in-cell GM12878, smoothing DNase to destroy all
sub-250 bp structure costs +0.0011 (*p*=0.72). So the thing the converter is good at is worth
nothing where DNase is shallow, and where it is worth something the converter still loses.

**Consequences**: The DNase line keeps its positive results, which do not depend on the
converter: DNase as an input (F-010, F-014) and DNase as an ABC activity term (+0.0709, F-012).
Both require real DNase, so the deployment story is now "use DNase where it exists" rather than
"synthesise it where it does not". That makes the DNase-input panel across more cell types the
successor to this line, and THP-1 is built and ready for it. The painted tracks, `4.20`, `4.22`,
`0.36`, `0.37`, `0.39` and `0.40` are kept: the machinery is what a later converter attempt
would reuse, and the shape-destroyed control is reusable for any input-property question.

**Tags**: converter, dnase, atac, accessibility-input, stop-decision, negative-result

---

### [2026-10-01] Notes migrated from todo/TODOLIST.md

**Context**: `todo/TODOLIST.md` was migrated to the project-dashboard TODO standard (Now/Next/Later/Waiting/Inbox/Settled/Done; tasks only). Its narrative, results and reasoning are kept here verbatim so nothing is lost; open tasks were rewritten as one-line items in the new file, settled conclusions became `## Settled` lines.

**Decision**: Approved by Maya 2026-10-01, with two edits: no search for a second-cell-type benchmark, and the TF-activity follow-ups (finemo TF lists, leave-one-out Test 1, further channel variants and seed replicates) settled rather than scheduled, since F-034 to F-038 and F-041 put every variant within run noise.

**Consequences**: Snapshot of the old file at `todo/TODOLIST.pre-standard.md`. Section headings below are the old ones; item states are as of migration. The old header pointed to `2026_0824_H3K27ac_model/h3k27ac_model_report.html` for the analysis narrative.

**Tags**: todo, migration, housekeeping

#### Settled (do not redo)

Target = 5' ends, ±500 bp window, `count_loss_weight = 10`. Panel is **ATAC-only** (DNase
and ATAC are not interchangeable inputs). PE H3K27ac targets use **read 1 only**.
Accessibility inputs should be ChromBPNet-style 5' insertion counts; tracks are built and
validated but **no model uses them yet**. Residual-objective training helps only a
sequence-blind input and costs a multimodal one, replicated in K562 and GM12878, so it is
a property of the objective and needs no further per-cell-type testing.

#### What the ATAC-input models get wrong in transfer, F-021 (opened AND partly withdrawn 2026-09-17)

The best transferred ATAC-input arm is the **GM12878 p300 multimodal model into K562**, 0.4771
against a 0.4680 floor (+0.0091, ns). Transfer costs it +0.0458 [+0.0296, +0.0610], which is 83%
of its in-cell margin, and the cost is promoting accessible negatives (121 enter the top 10% of
pairs against in-cell's 102, while 25 positives leave and 5 arrive). Regulated pairs peak at ATAC
decile 9, not 10. `4.25.transfer_error_anatomy.py`.

**THE MECHANISM FIRST WRITTEN HERE WAS WRONG AND IS WITHDRAWN (same day).** The 0.64-to-2.16
accessibility slope was measured on ABC's post-qnorm activity over CRISPR-tested pairs only.
`4.28` and `4.29` decompose it: raw bigwigs over all candidate regions give x1.20, post-qnorm
over all elements x1.85, post-qnorm over CRISPR pairs x3.40. **The model contributes about 15%
of the swing; qnorm and the element subset manufacture the rest.** The two models' raw outputs
differ by a near-constant 0.48 factor. Do not describe this as the transferred model reading
accessibility too steeply.

- (done) **Accessibility-normalisation mismatch. DONE 2026-09-17. NOT the mechanism (F-022) and
      mildly HARMFUL downstream (F-023): -0.0058 [-0.0085, -0.0029] with 100% sign retention on
      the in-cell p300 geomean arm, null everywhere else. DO NOT ADOPT.** Re-predicted both p300 arms on region-derived statistics (`4.26` gives
      3.7232 / 1.1428 for K562 ATAC on the ABC regions, against the models' stored 4.155 /
      1.313 in-cell and 4.552 / 1.360 transferred). The tracks' scale moved x1.5 and x1.8 but
      the per-element ranking moved only 1.6% and 3.7%, and **ABC's rank qnorm absorbed the
      scale exactly**: post-qnorm activity medians 0.5372 -> 0.5374 and 0.5369 -> 0.5369.
      So the centring asymmetry was real and is now corrected, but it cannot be what costs the
      transferred arm 0.0458 AUPRC.
      Benchmark confirmation is pending in the `slopefix` comparison, but the structural result
      already settles it.
- (done) **Monotone recalibration of the predicted track. DROPPED 2026-09-17, the premise was
      wrong.** The plan was to fit a monotone correction on raw predictions in the source cell
      type. `4.28` shows the raw predictions carry only a x1.20 accessibility dependence, so
      there is almost nothing there to correct; the x3.40 came from qnorm and the CRISPR subset.
      A correction would have to act on post-qnorm activity, but ABC re-derives that per arm
      from the track it is given, so there is no stable target to fit against.
      **If this line is reopened, the object of study is the SHAPE of the predicted activity
      distribution**, since rank normalization against a shared reference is what converts a
      shape difference into an accessibility-dependent reordering. Note F-013 already failed at
      shape correction by quantile mapping, for a different track, and that failure is the
      nearest prior.
- (done) **Residual objective DOWNSTREAM on transfer. DONE 2026-09-17, NULL, F-023.**
      Transferred residual minus transferred plain is -0.0011 [-0.0114, +0.0096] (59.3%), in-cell
      -0.0006 [-0.0093, +0.0086] (54.9%), and the transferred residual arm is at or below the
      floor at -0.0059 [-0.0215, +0.0095]. The benchmark agreed with the upstream measurement
      this time. **Do not reopen without a new mechanism**; the objective relocates the
      accessibility channel into the offset model rather than removing it.
      Original reasoning, kept because the prediction was right:
      **CORRECTION 2026-09-17: transfer WAS measured upstream**, in
      `results/deploy_gm_to_k562_fold_summary.tsv` from `2.13.submit_residual_transfer.sh`,
      and I wrongly wrote earlier today that every residual measurement was in-cell.
      GM12878 -> K562: residual multimodal 0.823 overall / 0.600 top-quintile against plain
      multimodal 0.828 / 0.602, with an accessibility-only floor of 0.802 / 0.541. **The two
      objectives are indistinguishable on transfer upstream.**
      The reason to still run it downstream is that the benchmark and top-quintile Pearson
      have dissociated three times (F-004, F-018, F-019), twice in the direction of the
      benchmark rewarding something Pearson scored at zero. **Expect null; do not present it
      as promising.**
      **A caveat on the mechanism, which weakens the original rationale.** The residual
      objective does not remove the accessibility channel from the final prediction, it
      relocates it: prediction = residual + offset, and the offset IS an accessibility-only
      model, which when transferred carries its own accessibility-slope error. The in-cell
      grids show this cleanly, `r_out_vs_atac` 0.921 for plain multimodal against 0.020 for
      the residual output alone. So the test asks whether isolating the slope into a simpler
      model helps, not whether the slope is gone.
- (open; now in TODO or Settled) **Add an accessibility-stratified metric to the evaluation.** Top-quintile Pearson is
      stratified by SIGNAL quantile, which is why a 3.4x accessibility-slope error was invisible
      upstream. Report Pearson and the predicted/observed ratio by ATAC decile; `4.25` computes
      it on CRISPR elements and the same stratification belongs in the `2.x` genome-wide eval.
- (open; now in TODO or Settled) **Measure counts accuracy on CRISPR-tested elements, not only genome-wide.** For the p300
      arm, transfer loses 0.115 genome-wide overall Pearson but 0.236 on CRISPR-tested elements
      and 0.260 on regulated-pair elements, so the headline metric understates the damage where
      the benchmark is decided by about 2x. **The H3K27ac arm goes the other way**: its counts
      accuracy on CRISPR elements does not degrade on transfer at all (0.672 to 0.691, and 0.515
      to 0.621 on regulated-pair elements) while both arms sit at the floor. Do not assume either
      direction; measure it per arm.
- (done) **Stale comment in `4.3.setup_abc_arms.py` corrected 2026-09-17.** It claimed no GM12878
      p300 model exists. One does, at
      `2026_0606_GM12878_transferability/GM12878_multimodal_BPNet/models/atac`, and `P300_TX` has
      been using it all along. **It still has no `training_target.json`**, so what it predicts is
      recorded only in the June project's `config/input_data_gm12878_multimodal.json`; writing
      that file is worth doing so nobody has to correlate bigwigs to identify a model again.

#### The model under-uses TF information, 2026-09-28 (F-034, F-035)

- (done) Deployable channel `assoc_noctcf` (F-036): +0.010 ABC, +0.022 rE2G over MSE-only, best
      transferred model so far. Oracle channel + rank loss w=10: does not add.
- (done) Seed replicates (F-038): run sd ~0.002 ABC, ~0.005 rE2G; channel +0.006 ABC, +0.012 rE2G.
- (done) Frozen two-stage model (F-038): matches the end-to-end channel.
- (open; now in TODO or Settled) Corrected TF list (4.50 flagged + GC): needs 3 seeds to resolve its likely small gain.
- (done) **Teach the model to use the channel** (F-037), `1.38`, oracle weights: head / pair /
      headpair all within noise of the channel alone; pair gives the smallest GATA error.
  - (open; now in TODO or Settled) Paired loss on the deployable `assoc_noctcf` channel (untested).
- (open; now in TODO or Settled) Later candidates:
  - (open; now in TODO or Settled) Channel content: weights as departure from the training-cell mean truth; split positive
        and negative weights into two channels.
  - (open; now in TODO or Settled) Motif dropout / weight noise on the channel during training, and train on estimated
        (ATAC-derived) weights to match deployment.
- (open; now in TODO or Settled) Seed replicate of the oracle channel model (F-035 is one run; ~6 more epochs than comparator).

#### Nucleosome arrangement predicts H3K27ac, 2026-09-30 (F-039)

- (done) Fragment features vs H3K27ac (4.55): joint partial rho 0.48-0.52 beyond ATAC (p300
      0.16-0.22), opposite sign to p300, weights carry K562 <-> GM12878 losing 0.01.
- (done) Profiles (4.56, fig35): the signal is mostly WHERE accessibility lies (H3K27ac-high = broad
      domain, low = one sharp isolated peak), not nucleosome packing.
- (in progress at migration) H3K27ac correction (4.57, fitted in GM12878, applied to K562): element Spearman 0.732 ->
      0.787 (in-cell 0.791), mostly recalibration; ATAC +/- 1 kb matches fragments. ABC, CRISPR
      and rE2G done 2026-09-30; bootstraps running (k27_boot, re2g_boot). Write F-043.
- (in progress at migration) p300 x H3K27ac combined (4.58): observed ceiling and predicted (no training) arms in ABC and
      rE2G; same bootstraps. If observed combo beats both singles -> two-head model (see
      multi-head item below).
- (done) Flank fragment length (4.59, F-040): no +10-15 bp at H3K27ac-high or p300-high flanks, even
      at the adjacent nucleosome. Closed.

#### TF identity is the missing signal, 2026-09-27 (F-030 to F-033)

- (done) Ranking loss (F-030): weight 10 is best, about +0.01 AUPRC in ABC and rE2G.
- (open; now in TODO or Settled) Seed replicate of rank weight 10, to separate the gain from run-to-run noise (F-020).
- (done) Neighbour graph (F-031): null in ABC, same-class and cross-class. Not tried in rE2G.
- (done) Error anatomy and GATA oracle (F-032): GATA elements are about half the transfer gap.
- (done) TF-activity recovery (F-033): GATA and CEBP recovered from ATAC; whole ranking moderate.
- (done) Test 1 with oracle weights (F-035): rE2G +0.017, ABC +0.004; GATA error shrinks a third.
- (done) Deployable estimator (F-034): learned estimators lose the cell-specific part; ATAC
      association minus CTCF keeps it best.
- (open; now in TODO or Settled) finemo instead of FIMO, later: ChromBPNet (ATAC) vs BPNet (p300) finemo in the same cell
      separates "important for accessibility" from "important for p300 beyond accessibility".
- (open; now in TODO or Settled) Leave-one-out Test 1 for the other four cell types (upstream metrics only; no CRISPR).

#### Where the activity term loses, 2026-09-23 (F-028, F-029)

- (done) **Oracle decomposition** (F-028). 80% of the transfer gap is tested-element ranking.
- (done) **ENCODE-rE2G with p300 in the activity slot** (F-029). Observed p300 +0.12, predicted
      +0.005 to +0.008. Predictions are mostly ATAC re-encoded (partial rho 0.13-0.28).
- (open; now in TODO or Settled) **Ranking-loss models (1.36), weights 2 and 10.** Trained; ABC predictions running
      (44954760, 44954762). Then 4.3/4.5 ABC, 4.6/4.7 benchmark, 4.35 partial rho, and add to
      the rE2G run. Stopping epochs 36.4 / 34.0 against the MSE baseline's 34.6.
- (open; now in TODO or Settled) **rE2G promoter-feature ablation.** Tests F-029's leading hypothesis for why in-cell
      predicted p300 adds nothing while observed H3K27ac adds +0.083: drop
      normalizedH3K27ac_prom and combinedActivity_prom and retrain.
- (open; now in TODO or Settled) **Use partial rho given ATAC (4.35) as the upstream screening metric** for new models,
      over all candidate regions, never on the tested elements.
- **AlphaGenome: ON HOLD (Maya, 2026-09-23).** Every released AlphaGenome model was trained
  on K562 EP300 and H3K27ac, so no K562 evaluation is a transfer test. The only cell types
  it has not seen the target in (WTC11, Jurkat, HCT116 for EP300) have small, differently
  designed CRISPR sets that are not a reliable benchmark. K562 stays the benchmark.
  Weights and the PyTorch port are already on Oak if this reopens.

#### Highest value

- (open; now in TODO or Settled) **W4. Motif syntax (SHAP / TF-MoDISco / FiNeMo).** Not started; the actual scientific
      goal. Run on the **residual-trained** model, its attributions are forced onto
      accessibility-independent signal, which multimodal attributions cannot separate. Use
      the ±500 bp window (zero neighbour contamination). Expect less signal than p300.

- (open; now in TODO or Settled) **TeloHAEC training + transfer.** The only new cell type available under the ATAC-only
      rule. Tracks, elements and model-free coupling (0.33-0.37 top quintile) are all ready.
      3 modes x 5 folds, then the four-evaluation transfer set against K562 and GM12878.
- (open; now in TODO or Settled) **TeloHAEC ±IL1b / ±TNFa / −VEGF.** Same genome and cell line, different regulatory
      state, no input-domain shift, a sharp and cheap test of whether the model tracks
      condition-specific change. Inference only if trained on ctrl.

#### Interrogating the ABC negative result, the main thread

The CRISPR benchmark says predicted H3K27ac adds nothing detectable over ATAC alone (best
predicted arm 0.482 vs floor 0.457, CI overlapping; deployment-scenario arms at or below the
floor). `4.9.diagnose_abc_gap.py` localised why: the prediction is accurate (Spearman 0.82 on
CRISPR-tested regions, better than genome-wide) but **ranks the functional elements too low**,
and since ABC's qnorm removes scale by construction, rank is the only channel available. All
arms catch the same positives, so the deficit is in suppressing negatives.

`4.10`/`4.12` then localised the error to a specific population: every model that sees ATAC
over-predicts H3K27ac at accessible-but-unacetylated elements by 7-8x its own median, while
sequence-only elevates them 1.6x and observed H3K27ac not at all. Those elements are GC 0.59,
CpG o/e 0.55, 2.4x promoter-enriched, ATAC/K27ac ratio 13.25 vs 2.20 typical, a CpG-island /
CTCF phenotype. **The signal is in sequence; the additive trunk gives it no way to veto
accessibility.**

- (done) **Gate x loss factorial. CLOSED 2026-09-05** (`1.19`, decision 2026-09-05). All 15
      fold-jobs finished and were scored; this entry was left marked RUNNING by mistake and the
      marker was cleared 2026-09-14 after re-deriving the paired deltas from
      `gatefac_k562_per_fold.tsv` with `2.38`. In-cell K562, paired against the matched ungated
      baseline on top-quintile Pearson: GATE -0.0001 (*p*=0.99), ASYM +0.0001 (*p*=0.99),
      GATE_ASYM -0.0032 (*p*=0.51). Transferred to GM12878: GATE -0.0016 (*p*=0.82), ASYM
      -0.0078 (*p*=0.11), **GATE_ASYM -0.0128 (*p*=0.039), i.e. resolvably worse**. Both
      interventions shared one premise, that the sequence branch can identify the
      accessible-but-unacetylated elements, and that premise is measured false (1.55x against
      1.44x where truth separates the tails 1.00x against 31.3x). The flags and
      `test_asymmetric_loss.py` stay in the code so the negative result stays reproducible.
      **Not scored, and the one thing that would change the reading:** nobody ran the `4.12`
      fold-elevation comparison on the gate arms, so the closure rests on aggregate *r* rather
      than on the failure population the gate was built for. Reopen only if some other result
      revives the premise; on its own it is not worth two `4.1` prediction passes.
- (open; now in TODO or Settled) **Indicator-channel control.** GC and CpG-density tracks as extra accessibility
      channels, no gate, no loss change. If hand-supplied class information does as well as a
      learned gate, prefer it, simpler and interpretable. Both are sequence-derived, so
      unlike fragment channels they cost nothing in transferability. Needs `0.27` to build the
      tracks first.
- (done) **CTCF/EP300 enrichment in the error strata, PEAK OVERLAPS DONE** (`4.11`).
      Over-predicted 5%: CTCF 36.9% vs 17.4% typical (2.1x, confirming the CTCF hypothesis),
      but EP300 16.4% vs 8.8%, ENRICHED, contradicting my prediction of depletion, and
      H3K27ac peaks 30.6% vs 15.8% despite low H3K27ac RPM. Under-predicted 1%: EP300 57.7%
      (6.6x), H3K4me1 86.8% (3.1x), H3K27me3 0.0%, CTCF 9.3% (depleted), canonical active
      enhancers the model misses. CTCF enrichment is NON-MONOTONIC (36.9% at 5%, 23.4% at 1%),
      unexplained; the most extreme over-predictions may be a different population.
- (done) **Quantitative signal, not just peak overlaps. DONE 2026-09-05** (`4.13`, F-007);
      marker cleared 2026-09-14. Per-element RPKM from the local BAMs is in
      `error_strata_rpkm.tsv` and it changed the reading: the over-predicted 5% carries CTCF
      4.69 against 1.79 typical but H3K27ac 1.17 at 0.94x typical with **IgG 0.39 below
      background 0.48**, so the H3K27ac peaks called there at 2x background are threshold
      artefacts and the stratum really is unacetylated. The same control demoted the
      under-predicted stratum's p300 enrichment from 6.6x (peaks) to about 2.3x (signal minus
      input), which is what closed p300 as an auxiliary target (F-006). Peaks alone would have
      overstated both premises.
- (open; now in TODO or Settled) **Characterise the 12 threshold-level false negatives** individually once `4.11` lands -
      observed H3K27ac 39x the genome median where the model predicts near-background.

#### Sequence-gated multimodal for p300, diagnose before training

Naming: the architecture is **sequence-gated multimodal**; `gate` is the flag and filename
token. The asymmetric loss is named separately on purpose, because the factorial exists to
attribute any gain to the architecture or to the loss.

- (open; now in TODO or Settled) **Confirm p300 has the same failure mode BEFORE training p300 with the gate.** The
      premise is that p300 suffers the same accessibility domination, and it may not: p300's
      residual *r* is 0.654 against H3K27ac's ~0.55, so sequence adds MORE beyond accessibility
      there. That is consistent with the problem being milder, or with it being equally severe
      but more fixable, different expected gains, same number.
      Cheap test, reusing existing machinery: predict p300 on the ABC candidate regions with
      the existing `2026_0529_multimodal_p300_model` models via `4.1`, then run the `4.12`
      fold-elevation comparison against observed p300 (ENCSR000EGE,
      `ENCFF466WKF`/`ENCFF163FSR`). Confirmed if p300's sequence-only arm gets high-ATAC /
      low-p300 elements right while its multimodal arm does not.
- (open; now in TODO or Settled) **Then train p300 with the gate**, 15 fold-jobs, only if the diagnostic confirms it.
      `multimodal_bpnet.py` is already shared, so no further code change is needed, the p300
      submit scripts just gain the two flags.

#### Which model transfers best, do not assume the in-cell-type winner

- (done) **Transfer matrix. CLOSED 2026-09-05** (`2.24`, F-005). Finished and written up; this
      entry was left marked RUNNING by mistake and the marker was cleared 2026-09-14 after
      re-deriving every delta from `txmatrix_*_per_fold.tsv` with `2.38`, which reproduced
      F-005 exactly. Paired against narrow+flat transferred, top-quintile Pearson: wide
      +0.0114 (*p*=0.14) K562->GM12878 and +0.0053 (*p*=0.53) GM12878->K562; **fragment
      channels -0.0011 (*p*=0.86) and -0.0029 (*p*=0.75), exactly null in both directions**
      despite gaining +0.0055 (*p*=0.032) and +0.0156 (*p*=0.002) in-cell. Fragment-length
      structure is a property of a specific ATAC library, which no in-cell test can detect.
      Wide moves `overall_pearson` (+0.0079, *p*=0.009) and `residual_pearson` (+0.0495,
      *p*=0.0009) transferred while leaving top-quintile unresolved, so it survives transfer
      on the mechanistic metric only. A transferred multimodal beats the target's own
      ATAC-only floor by +0.0306 [+0.0130, +0.0482] (*p*=0.008) into K562 but only
      +0.0132 [-0.0141, +0.0405] (*p*=0.25) into GM12878.
- (done) **Re-test F-004's two transfer arms with the paired bootstrap. DONE 2026-09-17,
      `results/f004_transfer_paired.tsv`.** They are genuinely AT the ABC floor, not merely
      unresolvable on the unpaired CIs they were originally judged with. GM12878-trained
      applied to K562: ATAC-input **-0.0028 [-0.0132, +0.0078]** (sign kept 73.2%),
      multimodal **-0.0048 [-0.0216, +0.0110]** (72.7%). Pairing cut the intervals from
      about +/-0.05 to +/-0.01 and they still straddle zero.
      **Pairing DOES resolve the transfer penalty, which the unpaired intervals could not:**
      in-cell minus transferred is **+0.0169 [+0.0041, +0.0297]** (99.6%) for the multimodal
      model, +0.0053 [-0.0013, +0.0119] (94.5%) for the ATAC-input one. First paired
      DOWNSTREAM confirmation of the transfer failure F-005/F-016/F-017 found upstream.
      **Tension worth carrying:** the F-005 entry above has a transferred multimodal beating
      the target's own ATAC-only *model* floor by +0.0306 upstream into K562, while the same
      direction is at the ABC activity floor downstream. The two floors are not the same
      object (a model against a raw track in the activity slot), so this is not a
      contradiction, but do not quote the upstream number as evidence of downstream transfer.

- (open; now in TODO or Settled) **TeloHAEC as a third cell type, with conditions.** The only new cell type under the
      ATAC-only rule. Tracks, elements and model-free coupling are ready. Caveats to carry:
      ATAC-derived elements (though derivation was shown not to matter, p = 0.83), 36 bp reads
      vs 95 bp, and shallower libraries, so read length and depth remain confounded even
      though element derivation does not. Fragment channels there would need PE BAMs.
- (open; now in TODO or Settled) **Train on multiple cell types to optimise transferability.** The most promising route
      to a deployable model, since it directly optimises what deployment needs rather than
      in-cell-type fit. Do it after the transfer matrix says which architecture to carry.

#### Using p300 without needing p300 at deployment

The application target has ATAC and nothing else, so **p300 can never be a model input**, it
is unavailable in the target cell type, whatever we hold for training. It can be a target, an
auxiliary task, or an annotation. Ordered so the cheap prerequisite gates the expensive work:

- (open; now in TODO or Settled) **First: do p300-target models predict the elements the H3K27ac model misses?** The
      under-predicted stratum is 6.6x EP300-peak-enriched, 3.1x H3K4me1, CTCF-depleted, zero
      H3K27me3, canonical active enhancers the H3K27ac model treats as background. If the
      existing p300 models in `2026_0529_multimodal_p300_model` also miss them, neither idea
      below helps and both are dropped. Reuses `4.1` (predict p300 on the ABC regions) and
      `4.12` (fold-elevation vs observed p300, ENCSR000EGE).
- (open; now in TODO or Settled) **Then: multi-head model**, shared trunk, two counts heads, predicting H3K27ac and
      p300 from ATAC + sequence with a weighted loss. Deployment-valid because p300 is needed
      only at TRAINING time; inference reads the H3K27ac head. Directly supervises the failing
      population rather than adding capacity. Needs `extract_windows` to carry a second signal
      track and the loss to sum two count terms, more invasive than the gate, and
      `multimodal_bpnet.py` is shared with p300, so it needs the `test_asymmetric_loss.py`
      treatment.
      Trainable in K562 and GM12878, the same two cell types all the transfer work uses.
- (open; now in TODO or Settled) **Fallback only: stack predicted p300 as an input channel.** Information-theoretically
      redundant, a predicted p300 track is a deterministic function of ATAC + sequence, which
      the H3K27ac model already sees, so it can only act as an inductive bias, and it costs
      two training runs and two models at inference to get what multi-head gets in one. Try
      only if multi-head underperforms.

#### p300 as the activity term instead of H3K27ac

- (open; now in TODO or Settled) **ABC + CRISPR with p300 predictions, and with observed p300.** p300 is a coactivator,
      one step closer to the sequence-specified event, and is arguably what ABC's activity
      term is really proxying. Observed p300 gives both a ceiling and, importantly, a qnorm
      REFERENCE, so the no-reference problem dissolves: build a p300 qnorm reference from the
      observed arm rather than running without qnorm and having to scale-match the geomean by
      hand. Data: EP300 ENCSR000EGE, `ENCFF466WKF/ENCFF163FSR.filtered.sorted.bam`, peaks
      `ENCFF702XPO.bed.gz`. p300 models already exist in `2026_0529_multimodal_p300_model`
      (residual r 0.654, higher than H3K27ac's).
- (open; now in TODO or Settled) **Test removing qnorm** as its own arm (`use_qnorm: False`; `activity_base_no_qnorm` is
      already in every EnhancerList). Expectation is that it hurts, the model emits log1p
      counts over ±500 bp while observed H3K27ac is read counts over the element, so without
      qnorm the geomean multiplies incommensurate magnitudes, and ABC's thresholds are
      calibrated on qnorm'd values. Cheap enough to settle rather than argue.

#### Downstream utility, decides whether the correlation metrics are the right target

- (open; now in TODO or Settled) **Plug predicted H3K27ac into ABC and benchmark it.** The end-to-end test of whether
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

- (done) **ATAC to DNase converter. BUILT AND BENCHMARKED 2026-09-12.** The cheap gate passed
      and the whole line of work ran. What it settled, and what it did not:

      **The assay is settled.** DNase beats ATAC as a model input, +0.037 in-cell K562 and
      +0.086 in-cell GM12878 top-quintile Pearson, and it is the only input change in the
      project that survives transfer, clearing GM12878's own ATAC-only floor by
      +0.056 [+0.014, +0.099] (F-010). As the ABC activity term it beats ATAC by
      **+0.0709 [+0.0476, +0.0938]** (F-012), a larger downstream effect than anything the
      modelling has produced. DNase also has a base-resolution profile where H3K27ac has
      none, 1 bp top-quintile ceiling 0.848 against 0.21.

      **The converter reproduces the profile and nothing downstream yet.** 1 bp top-quintile
      shape against observed DNase is 0.792 against a ceiling of 0.834 in K562 (95%), where
      the raw ATAC track manages 0.292; sequence does the work, since the ATAC-only arm
      reaches 0.563 (F-011). Transferred to GM12878 it reaches 0.490 against that cell type's
      0.684 ceiling (72%) and keeps 91% of its sequence margin. But in the CRISPR benchmark
      neither converter clears real ATAC resolvably (+0.021 [-0.005, +0.044] K562-trained,
      +0.014 [-0.012, +0.036] GM12878-trained) and real DNase still beats the better one by
      +0.050 [+0.030, +0.071]. Training cell type does not matter, unlike p300 (F-009).

      **That benchmark could not test the converter's actual claim.** ABC sums the activity
      track over each region and discards the profile, so F-012 scored the converter's counts,
      which are its weak half. The confound the original gate worried about was handled the
      other way round from what was planned: the converter trains on the ATAC-derived element
      set (D-19) so it never trains on regions chosen by the signal it predicts.

- (done) **H3K27ac under four accessibility inputs. DONE 2026-09-12, F-013.** Top-quintile
      Pearson: real DNase 0.728, DNase smoothed at 250 bp 0.704, ATAC 0.690, converted DNase
      0.680. Paired against ATAC: real DNase +0.0379 [+0.0213, +0.0546], smoothed
      +0.0141 [-0.0057, +0.0338] (not resolvable), **converted -0.0099 [-0.0175, -0.0023],
      i.e. resolvably worse than the track it would replace**.

      **The control answered: shape, not magnitude.** Destroying structure finer than 250 bp
      costs real DNase +0.0239 [+0.0027, +0.0451] of its +0.0379, and the remainder no longer
      separates from ATAC. So base-resolution structure is where the advantage lives, which is
      the branch that keeps a converter conceptually alive.

      **But this converter is worse than raw ATAC.** Reproducing the DNase profile at 95% of
      ceiling (F-011) while hurting downstream means the painted track carries a defect that
      outweighs its shape fidelity. Named candidate: the 1.54x magnitude inflation of the
      lowest observed-signal quintile (D-22), plus a 34% genome-wide total. Shape fidelity is
      necessary and demonstrably not sufficient.

      **Also measured:** `profile_pearson` is 0.063-0.064 across all four arms, so the
      H3K27ac profile head learns nothing whatever accessibility it is given.

- (done) **Magnitude fixes for the painted track. BOTH FAILED 2026-09-14, F-014.** Quantile
      mapping to another cell type's DNase distribution (`0.37`) matched the marginal exactly
      and made dynamic range WORSE, 18x to 9x against real DNase's 39x, costing 0.12 of 1 bp
      top-quintile shape. A power transform with gamma fitted out-of-cell (`0.39`, `0.40`)
      raised range to 24x for a shape cost of 0.006, and moved the downstream delta from
      -0.0099 to -0.0105, i.e. nothing. Magnitude was not what was costing the converter.

- (done) **The converter is finished as a model input. STOPPED 2026-09-14** (decision
      2026-09-14, F-014). Resolvably worse than plain ATAC in all four directions of the 2x2:
      -0.0105 in-cell K562, -0.0579 K562->GM12878, -0.0206 in-cell GM12878. Do not attempt a
      third transform. The machinery is kept (`4.20`, `4.22`, `0.36`, `0.37`, `0.39`, `0.40`)
      because a later attempt would reuse it and the shape-destroyed control generalises.

- (done) **Filter-width sweep for the shape control. SUPERSEDED 2026-09-14, F-014.** The 2x2
      answered the question the sweep was meant to tighten, and answered it differently than
      expected: whether DNase wins on shape or magnitude is CELL-TYPE DEPENDENT. Smoothing at
      250 bp costs real DNase +0.0239 (*p*=0.035) in K562 and +0.0677 (*p*=0.0002) transferred
      to GM12878, but **+0.0011 (*p*=0.72) in-cell GM12878**. Shape matters where DNase is deep
      enough to have reliable shape: GM12878's DNase is 53.8M reads against K562's 301.1M, with
      a profile ceiling of 0.686 against 0.848. A width sweep within K562 would have tightened
      a claim that does not generalise; the informative axis is depth, not width.

- (done) **DNase-input panel, three cell types. DONE 2026-09-14, F-016.** THP-1 trained
      (`1.30`), the two missing DNase-only floors trained (`1.22`/`1.23` at `MODE atac`),
      panel scored (`2.40` + `2.41`). **The answer is negative and it retires the deployment
      thesis for H3K27ac: all six transferred arms land at or below the target's own
      DNase-only floor**, five resolvably. Top-quintile, transferred minus the target's
      assay-matched floor: K562<-GM12878 -0.133, K562<-THP-1 -0.178, GM12878<-K562 -0.045,
      GM12878<-THP-1 -0.016 (*p*=0.11), THP-1<-GM12878 -0.070, THP-1<-K562 -0.388. A LOCAL
      sequence+DNase model beats that floor in all three cell types (+0.121, +0.048, +0.057,
      all *p*<0.006), so sequence does add real information wherever it can be trained
      locally.

      **F-010's transfer claim is withdrawn.** It compared a DNase-input model against an
      ATAC-only floor, and that floor is much the weaker comparator (GM12878 0.518 against
      DNase-only 0.615). **General rule now recorded: a floor must use the same assay as the
      model it bounds**, or the comparison credits an assay swap to the architecture.
      F-010's in-cell gains (+0.037 K562, +0.086 GM12878) are unaffected.

      **Neither floor is a true deployment baseline**, since both need target-cell H3K27ac to
      train. The honest deployment comparator is the raw accessibility track used directly as
      the activity proxy, which is what ABC does and what F-012 measured: real DNase over real
      ATAC, +0.0709 [+0.0476, +0.0938]. That is still the only DNase result here with a
      downstream effect, and it needs no model at all.

      THP-1's own caveats stand: one DNase replicate, so no DNase shape ceiling and no use as
      a converter target; no ATAC, so it cannot join an ATAC-input comparison. Tracks and the
      H3K27ac ceiling (0.974 raw, 0.993 corrected) came from `0.38`; DNase alignments at
      `Users/sheth/Data/ENCODE/THP1/DNase/AG81591.filtered.bam`.

- (done) **Sanity-check the GM12878->K562 collapse. ANSWERED 2026-09-14 by F-016, and it was
      not GM12878.** K562 is the hardest TARGET, taking -0.133 from GM12878 and -0.178 from
      THP-1, while GM12878 as a source transfers into THP-1 at only -0.070. So the asymmetry
      is a property of K562 as a destination rather than of the GM12878 models being weak,
      which is what the in-cell GM12878 arms (0.663) had suggested but could not establish.
      The same question for p300 (F-009) is NOT closed by this: it was never re-scored against
      an assay-matched floor and p300 has no THP-1 arm.

- (done) **Depth-match the DNase inputs and re-score the panel. DONE 2026-09-15, F-017.**
      Thinned K562 and THP-1 to GM12878's 50.3M reads by binomial subsampling (`0.41`,
      within 0.01%), retrained both cells' arms (`1.26`/`1.30` with `ACC_BW`), re-scored
      (`2.41 DEPTH_MATCHED=1`). Per-window spread 7.1x -> 1.4x. Element sets unchanged, so
      the two panels score identical regions.

      **Depth WAS the asymmetry F-009 and F-014 both ended on.** GM12878->K562, the severe
      collapse at -0.133, reaches parity with K562's own floor once matched (+0.003,
      *p*=0.68). A model trained shallow and handed a 7x deeper library is out of
      distribution on its own input. The mild direction K562->GM12878 was NOT depth and is
      unchanged (-0.045 -> -0.039).

      **Transfer still fails, so F-016 stands.** Five of six directions remain resolvably
      below the target's floor, the sixth only reaches parity, and a local model wins
      everywhere (+0.136, +0.048, +0.064).

      **Two things to carry forward.** (1) Depth-match before ANY cross-cell-type comparison
      in this project, or a collapse gets misattributed to the cell types; this generalises
      F-005 from fragment-length structure to read depth. (2) Thinning costs the deep cell
      type real performance (K562 floor 0.607 -> 0.570, local 0.728 -> 0.707), so it is the
      right control for a transfer question and the wrong choice for a production model.

      **My own reasoning to not repeat:** F-016 argued depth was unlikely because the
      mismatch ratios do not order with transfer quality. They still do not, and the effect
      was real anyway. Ratio-ordering was the wrong test of the hypothesis.

- (open; now in TODO or Settled) **Re-run the thinning at a second seed.** One random draw at seed 0 stands behind
      F-017. Cheap (`0.41 --seed 1`, then retrain and re-score) and it bounds the sampling
      noise on a finding that now carries the transfer conclusion. Only worth doing if F-017
      is going into a figure or the manuscript.

- (done) **Multi-task arm: DNase profile head, H3K27ac counts head. DONE 2026-09-14, F-015
      established.** Kept, but it is a small architectural win and not progress on the real
      gap. Against an **epoch-matched** baseline: overall Pearson **+0.0026 [+0.0006,
      +0.0046]** (*p*=0.022), `incremental_r2` 0.005 [0.002, 0.007] and `incremental_r2_topq`
      0.015 [0.001, 0.030] both clearing zero. **Top-quintile Pearson does NOT clear against
      that control: +0.0076 [-0.0078, +0.0230], *p*=0.24** (0.701 against 0.693 and the
      early-stopped baseline's 0.690). The +0.0109 (*p*=0.043) first reported was against the
      early-stopped baseline and about a third of it was the epoch difference.

      **The auxiliary task itself is decisively learnable:** validation profile *r*
      0.530-0.561 in all five folds against 0.060-0.068, and the head drops H3K27ac shape to
      0.004. Needs no DNase at inference.

      **The epoch confound was real and is excluded** (`1.29`, jobs 43440391-96). Given the
      same 100-epoch budget the baseline moves +0.0000 [-0.0008, +0.0009] (*p*=0.92) and its
      best checkpoint still lands at epoch 33-47. Longer training was a CONSEQUENCE of the
      different loss, not the cause of the gain. **Carry this forward as a method rule: any
      intervention that changes the loss also changes the stopping rule, so an epoch-matched
      arm is the only way to attribute the result.** Three earlier arms in this project
      changed the loss (`1.19` asym, the residual objective, every `clw` sweep point) and
      none of them had one.

      **Only open follow-up, and it is optional:** `--profile-loss-weight 0.0561` is the
      measured depth ratio and the smallest defensible weight, so the ceiling on this effect
      is unknown. A sweep upward is five fold-jobs per weight. Not worth it unless something
      else revives interest, since even a doubling of the effect stays an order of magnitude
      below swapping the accessibility INPUT to real DNase (+0.038, F-010).

      Original scoping, for the record: Same trunk, same counts objective, profile head's target swapped from
      H3K27ac to DNase. Trainer change is in (`44e0e3c`): `extract_windows` takes an optional
      second signal pair, the dataset carries it as a fourth/fifth slot, and the loss takes
      `y_profile`. Gated by `scripts/test_profile_target.py`;
      `scripts/test_asymmetric_loss.py` still passes unchanged, and all ten existing callers
      of `extract_windows` keep the 4-tuple. Submit with `1.28`, score with `2.39`
      (`config/multitask_k562_configs.json`), paired against
      `multimodal5p_accs5p_hw500_clw10`.

      **`--profile-loss-weight` was NOT in the original plan and the arm is uninterpretable
      without it.** MNLL is `-sum(y * log_softmax)` plus a y-only term, so it scales with the
      target's read depth, and `loss = profile + clw * count`. K562 DNase carries **17.8x**
      the reads of K562 H3K27ac over these 1 kb windows (mean per-window total 910.2 against
      51.1). Swapping the target in unweighted would have grown the profile term about 18x
      and divided the effective count weight by the same factor, so a loss on counts would
      have measured "we trained it to care less about counts" rather than anything about the
      auxiliary task. The arm runs at the measured ratio, 0.0561, which starts the profile
      term at the baseline's magnitude. Default stays 1.0, so no earlier run moves.

      **Read only the counts columns.** `profile_*` in `2.39`'s table scores the profile head
      against H3K27ac for both arms, so for this arm it measures how well a DNase-trained head
      happens to predict H3K27ac shape. A drop there is expected and says nothing.

      **If it comes out null, the weight is the first thing to question, not the last.** 0.0561
      matches gradient magnitude to the baseline, which is the clean isolation of the change,
      but it is also the smallest defensible weight. A sweep upward is the follow-up, and it
      is cheap: no new tracks, no painting, five fold-jobs per weight.

      **Its premise held up and is now measured twice.** `profile_pearson` sits at 0.063-0.064
      across ALL FOUR accessibility inputs in K562 (F-013), so the H3K27ac profile head learns
      nothing whatever it is fed. It is dead weight with a target whose 1 bp inter-replicate
      ceiling is 0.21.

      **Run it in K562, not GM12878, and F-014 is the reason.** A DNase profile head is only a
      learnable auxiliary task where DNase has reliable base-resolution shape. K562's DNase
      profile ceiling is 0.848; GM12878's is 0.686 on a 5.6x shallower library, and F-014
      showed that in GM12878 destroying all sub-250 bp DNase structure costs the model nothing
      (*p*=0.72). Training the auxiliary head against GM12878 DNase would be fitting the same
      noise the H3K27ac head already fits.

      **Note what it does NOT inherit from the converter's failure.** This uses DNase as a
      TARGET at training time and needs no DNase at inference, so it is unaffected by the
      painted track being unusable (F-014) and is deployment-legal by the same argument that
      makes a multi-head model preferable to feeding p300 in (DATA_INVENTORY, deployment
      constraint).

      **Why.** The H3K27ac profile head currently trains against a target whose 1 bp
      inter-replicate ceiling is 0.21, i.e. mostly noise. D-3 kept it anyway, down-weighted,
      and that entry's own consequences note flags that `fconv` still consumes capacity for a
      near-zero gradient. DNase's 1 bp ceiling is 0.848, so this gives the head a learnable
      task for the first time. Hypothesis: a learnable base-resolution task makes the shared
      trunk better at H3K27ac COUNTS, which is the quantity of interest.

      **Why not joint prediction of both signals.** A second counts objective competes for the
      same gradient, and the p300 auxiliary head is the precedent: dropped 2026-09-05 after
      the p300 models predicted 1.83x enrichment at H3K27ac's failure elements against a
      2.10x INPUT control, i.e. below the control. What makes DNase different is specific.
      p300 failed because it could not identify the right elements; DNase would supply
      base-resolution structure that H3K27ac simply does not have.

      **Why it waits.** If DNase's input advantage is magnitude rather than shape the case
      weakens, though not fatally: DNase as an auxiliary TARGET teaches the trunk a
      sequence-to-structure mapping, a different channel from DNase as an INPUT. The control
      informs the prior rather than deciding it.

      **Implementation.** Needs a trainer change, not just a config: the profile and counts
      heads currently read one target. Cheapest correct route is a second pair of signal
      bigwig arguments consumed only by the profile loss. No painting and no second model, so
      materially cheaper than the converter path.

- (open; now in TODO or Settled) **Converter count accuracy.** Considered and DEFERRED 2026-09-12, not closed. The
      converter already predicts DNase counts at top-quintile r 0.878 and still bought nothing
      in ABC, and the part of DNase counts carrying the downstream benefit is by construction
      the part that differs from ATAC counts. Revisit only if the H3K27ac arms say counts are
      the binding constraint. If they do, quantile-map the painted track to a reference
      accessibility distribution first, which is deployment-legal and fixes the residual
      magnitude distortion.

      **Caveat to carry into reading those arms.** The painted track compresses dynamic range:
      painted/observed ratio by observed-DNase quintile runs 2.73, 0.88, 0.73, 0.71, 0.73.
      Q2 to Q5 is near-uniform and training's own normalisation absorbs it, but the 1.54x
      residual inflation of the lowest quintile after thresholding (D-22) is the named
      alternative explanation if the converted arm underperforms.

      This does not contradict the ATAC-only panel rule. That rule exists because ATAC and
      DNase are not interchangeable as *inputs*; a converter is the principled way to get a
      DNase-like input everywhere without mixing assays across cell types.
      Also unblocks the composite-metric item below, whose best form needs DHS.

#### Training longer is NOT resolved and should not be adopted, F-019

- (done) **Is the +0.0086 the epoch budget or run variance? MEASURED 2026-09-17, and the
      run pair CANNOT attribute it. `1.31.epoch_budget_vs_run_variance.py`.** The proposed
      test (score the same run's epoch-35 and epoch-100 checkpoints) is impossible: training
      saves only two checkpoints per fold, best-validation (`multimodal_bpnet.torch`) and
      final (`multimodal_bpnet.final.torch`). What the logs do allow is comparing the two
      runs AT THE SAME EPOCH, where the budget cannot yet have had any effect, so every
      difference there is run variance.
      **Run variance at matched epochs 20/25/30 is sd 0.1753 MNLL and sd 0.00180 count
      Pearson. The scored best-checkpoint difference between the two runs is -0.0636 MNLL
      and +0.00003 count Pearson.** The claimed effect is 2.8x smaller than the noise on
      MNLL and 60x smaller on count Pearson. Worst single matched-epoch draw, fold 1 at
      epoch 25, is -0.4905 MNLL, about 8x the whole effect.
      **The budget was not even binding in every fold.** Baseline fold 1 ran to epoch 55
      with its save at 45, LATER than the 100-epoch run's save at 40. Patience was 10 in
      every baseline fold, so 4 of 5 folds saved later under the longer budget, but not
      fold 1.
      **DO NOT adopt 100 epochs as the default**, and do not re-run the adopt-or-not arms
      at 100 epochs as a replication; that repeats the one-run-per-condition design that
      caused this.

- (done) **Does training past the validation minimum help? NO, 5 of 5 folds, zero run
      variance. 2026-09-17.** Within the 100-epoch run, the final epoch is worse than its
      own best-validation checkpoint on every fold: MNLL +0.018 to +0.096, count Pearson
      -0.00046 to -0.00154. This is the one epoch-budget question a single run CAN answer,
      because both checkpoints come from the same weights lineage. Checkpoint selection is
      doing its job; the gain is not "more gradient steps are better".

- (open; now in TODO or Settled) **If the +0.0086 still matters, it needs seeds, not another single run.** Upstream
      metrics cannot see the effect at all, so the only way to resolve it is downstream at
      2-3 seeds per condition, paired via `4.17`. Cost is 2-3 extra trainings plus scoring
      per condition, which is why it is no longer cheap and no longer "free".
      **The benchmark has disagreed with Pearson three times (F-004, F-018, F-019), so a
      +0.0000 upstream effect with a real downstream one is not impossible** and this is
      not a refutation of the downstream number. It is a refutation of attributing that
      number to the epoch budget on this evidence.

- (done) **Multi-task model through ABC. DONE 2026-09-16, F-019.** It clears the floor
      (+0.0255 [+0.0068, +0.0439]) and is the only arm to do so while reading ATAC ALONE at
      inference, which is the real deployment constraint. **But most of that is the longer
      training, not the DNase supervision.** Against the epoch-matched baseline the
      auxiliary task is +0.0047 [-0.0035, +0.0122] (ns) in the geomean form and only
      +0.0129 [+0.0022, +0.0240] in the prediction-alone form.
      Do not quote the +0.0133 against the early-stopped baseline; that comparison carries
      the epoch confound and I made it before running the control.

#### The DNase-input model is downstream-useful, and that reopens the line

- (done) **Real-DNase-input model through ABC and CRISPR. DONE 2026-09-15, F-018. IT WORKS.**
      Predicted H3K27ac from the DNase-input model clears the benchmark floor by
      **+0.0839 [+0.0613, +0.1048]** (DNase in the activity slot) and **+0.0403 [+0.0230,
      +0.0563]** (ATAC in it), and beats real DNase alone by **+0.0238 [+0.0051, +0.0429]**.
      First predicted-activity arm in the project to clear the floor, and the first direct
      evidence that the model contributes information accessibility does not carry.

      **The input assay is what did it, not the architecture and not the statistics.** The
      identical architecture reading ATAC cannot separate from the floor even under the
      paired test (+0.0121 [-0.0053, +0.0285]), so F-004 survives, rescoped to ATAC-input
      models.

      **Do not quote "beats the ceiling".** The observed-H3K27ac arm counts reads through
      `count_bam`/`count_tagalign` while every bigwig arm uses `count_bigwig`, so the
      +0.0222 over it sits inside a counting-path confound. `4.6` records which pairs are
      legitimate; all the numbers above are matched on counting path.

- (open; now in TODO or Settled) **Re-measure the dynamic range of the DNase-input prediction.** F-004 named compressed
      spread on regulated-pair regions as the mechanism (predicted p99/p50 4.03 against 6.00
      observed, agreement falling to Spearman 0.663 there). F-018 does not test whether the
      DNase-input prediction fixed that, so the mechanism of the only positive downstream
      result in the project is currently unknown. Reuses `4.9`/`4.10` against
      `predk27ac_k562_dnase.bw`. Cheap, and it is the difference between "we found something
      that works" and "we know why it works".

- (open; now in TODO or Settled) **Does the DNase-input ABC arm transfer?** Expected to fail: F-016 and F-017 showed no
      H3K27ac model beats an assay-matched floor in any of six directions. Worth one arm
      anyway, because the downstream benchmark and the top-quintile metric have disagreed
      before (F-004 against every architecture gain), so a transfer failure on Pearson does
      not strictly imply one on AUPRC. Needs the GM12878 DNase-input model predicted onto
      K562 regions, which is one `4.4` call.

- (open; now in TODO or Settled) **THP-1 and GM12878 ABC arms.** The benchmark is K562-only for want of CRISPR data
      elsewhere, so this is blocked on benchmark availability rather than on compute. Check
      whether the DC-TAP or ENCODE-rE2G benchmark sets cover either cell type before
      planning it.

#### Adopt-or-not decisions, all measured

- (open; now in TODO or Settled) **Turn on test-time RC averaging as the default?** Free (one extra forward pass, no
      retraining) and positive for every model: +0.0162 sequence only, +0.0081 sequence +
      ATAC, +0.0016 and non-significant for ATAC only, which is the control (report Fig. 13).
      Currently opt-in via `2.15 --rc-average` so existing numbers stay comparable. Adopting
      it means re-scoring the report's tables as a set.
- (done) **Do the wider receptive field and the fragment channels combine? ANSWERED, they are
      SUB-ADDITIVE, and it does not matter for deployment.** Closed 2026-09-14; the model
      (`multimodal5p_fragchan_wide_hw500_clw10`) was trained and scored earlier and the answer
      was sitting unread in `txmatrix_gm12878_to_k562_per_fold.tsv` and
      `rc_fragwide_k562_per_fold.tsv`. In-cell K562 top-quintile, paired: wide+frag over
      wide-alone **+0.0065 (*p*=0.053)** against fragment channels' +0.0156 on their own, and
      over frag-alone +0.0188 (*p*=0.015) against wide's +0.0280 on its own. So each addition
      keeps roughly half its solo gain and they do not sum. **Deployment is decided by transfer,
      not by this:** wide+frag over wide-alone transferred is +0.0045 (*p*=0.52), and F-005
      already showed fragment channels are null transferred in both directions, so the
      combination adds nothing a deployed model could use. Prefer wide+flat.

#### Open questions

- (open; now in TODO or Settled) **Repeat the residual comparison for p300.** Resolved only for H3K27ac; p300's
      residual *r* is already 0.654 so the headroom may be smaller.
- (open; now in TODO or Settled) **Predict the composite activity metric directly.** Downstream uses
      `geomean(accessibility, H3K27ac)`, best as `geomean(DHS, H3K27ac)`. Is the composite
      easier to predict than predicting H3K27ac and combining afterwards?
      **Design caution:** if the accessibility input is the same assay as the accessibility
      term in the target, half the target is readable off the input and the comparison is
      circular. Predict `geomean(DHS, H3K27ac)` from sequence + ATAC, and baseline against
      the two-step route scored on the *same* composite.
- (open; now in TODO or Settled) **[Q for Maya] Is a mostly-accessibility model useful for the intended application,**
      or does the goal require the sequence component? The ABC benchmark above is the
      empirical form of this question and would settle it without needing an answer in
      advance. Also decides whether the +0.043 all-elements receptive-field gain matters:
      it is real, but it is entirely dead-vs-active separation, so it pays off only for an
      on/off task.

#### Architecture, in expected order of value

- (in progress at migration) **Wider receptive field, ADOPT for multimodal; the earlier "do not pursue" was wrong.**
      `n_layers` 10 (~4.2 kb vs ~1.1 kb), 5 folds, both cell types, paired within fold on the
      intersection of valid regions.
      Top quintile: sequence-only −0.006 (K562) and +0.010 (GM12878), both p=0.53;
      **multimodal +0.027 (p=0.006) and +0.014 (p=0.025)**, all five folds rising in both.
      Accessibility residual 0.502 -> 0.547 and 0.397 -> 0.469.
      Results in `wide_{k562,gm12878}_*.tsv`; report Fig. 10.
      - (done) **ATAC-only wide arm. DONE, and it reproduces the multimodal gain.** Closed
            2026-09-14; the 10 fold-jobs finished and were scored into
            `rc_wide_{k562,gm12878}_per_fold.tsv`, where the marker was left stale. Top-quintile,
            paired within fold: ATAC-only widening gains **+0.0168 [+0.0048, +0.0288]**
            (*p*=0.018) in K562 and **+0.0144 [+0.0071, +0.0217]** (*p*=0.005) in GM12878,
            against the multimodal gains of +0.0280 and +0.0162 in the same tables. In GM12878
            the ATAC-only arm captures essentially the whole effect, so the extra context is
            being read as accessibility neighbourhood rather than as sequence context; K562
            leaves room for a sequence contribution (+0.0280 against +0.0168) but the intervals
            overlap. Conclusion as pre-registered: widen the deployed model.
            **Caveat:** both arms come from the RC-averaged tables, which is internally
            consistent since every arm in those tables is RC-averaged, but the numbers are not
            directly comparable to the non-RC `wide_*` tables.
      - (open; now in TODO or Settled) Decide whether to re-run the transfer and deployment comparisons at
            `n_layers` 10, since those used the narrow multimodal model.
- (open; now in TODO or Settled) **Decide what to do with the profile head, one decision, three options.** The 1 bp
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
      answers it at zero architecture risk, run it jointly with bin size, since the two
      interact: binning changes how learnable the profile term is and therefore its optimal
      weight.

      **Gate on the real numbers.** Everything above about what the head achieves rests on a
      256-element smoke test (`profile_pearson` ~ 0.053 against a 0.21 ceiling). Job for the
      full 5-fold measurement is queued; use `prof_residual_grid_*.tsv` when it lands.

      **Cost of actual removal** (option a): it changes the model class, so every existing
      checkpoint becomes non-comparable and the whole grid needs retraining. It touches
      `multimodal_bpnet.py` and `train_multimodal_bpnet.py`, both shared with p300, so it
      needs a backward-compatibility regression, use the `2.18` pattern, exact on region
      counts and tolerant on metrics.
      **Payoff if it holds:** removes a term that is largely fitting Poisson noise, frees
      trunk capacity, and roughly halves the output tensor, buying a larger batch or a wider
      window at the same memory.
- (done) **Fragment-size ATAC channels, DONE and positive.** Top quintile 0.690 -> 0.703,
      paired +0.0135 [+0.0075, +0.0195], p = 0.0034, every fold rising; report Fig. 12.
      Strict superset of the flat input (bins sum to it exactly), so the gain is added
      information. Results in `fragchan_k562_per_fold.tsv`.
      - (open; now in TODO or Settled) **Do the two accessibility-side gains combine?** `n_layers` 10 x fragment
            channels. Both may be reading the same neighbourhood structure, so +0.027 and
            +0.0135 may not sum. 5 fold-jobs, and it decides what the deployed model is.
      - (open; now in TODO or Settled) **GM12878 replication** needs its paired-end BAMs downloaded, fragment length
            lives in TLEN, which the per-read tagAligns discard.
      - (open; now in TODO or Settled) **ATAC-only fragment arm** (`1.15 atac`), asks whether fragment structure alone
            beats flat ATAC, with no sequence involved.
- (open; now in TODO or Settled) **Switch the accessibility input to 5' counts, as a set.** ~35 fold-jobs. Within a
      cell type this moves only the ATAC-only model; across cell types it removes a
      read-length confound (TeloHAEC 36 bp vs K562 95 bp). Mixing the two inputs is invalid.
- (open; now in TODO or Settled) **Re-check ±500 vs ±1000** on the 5' target once a retrain happens anyway.
- (open; now in TODO or Settled) **Train on the ATAC-derived K562 element set**, easy (one path change in `1.11`),
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

#### Panel data caveats

- (open; now in TODO or Settled) **Element derivation is not uniform** (`reference/ELEMENT_DERIVATION.md`): K562 and
      GM12878 DNase-derived, TeloHAEC ATAC-derived. An ATAC-derived K562 set exists in
      `K562_ATAC_ChromBPNet/data/`; adopting it makes two of three consistent at the cost of
      re-deriving every K562 number.
- (open; now in TODO or Settled) **HCT116** has ENCODE ATAC but its H3K27ac replicates differ in run type, so no
      inter-replicate ceiling is computable. Trainable, not normalisable.
- (open; now in TODO or Settled) **H1, H9, Jurkat, THP-1 have no ENCODE ATAC at all** (all 559 released experiments
      checked). Reaching them means GEO/SRA fastqs through
      `Data/scripts/sra_paired_fastq_to_bam.sh`, a separate decision.
- (open; now in TODO or Settled) **TeloHAEC_ctrl/ATAC holds 3 EA.hy926 files** (`SRR20809434/435/436`) under the same
      sample name. Always use explicit accession lists, never a directory glob.

#### Statistical power

- (open; now in TODO or Settled) **Re-run the `count_loss_weight` sweep with ≥3 folds and ≥2 seeds.** The current pick
      (10) came from single folds; 3/10/100 are within noise of each other.
- (open; now in TODO or Settled) **Quantify run-to-run variance properly**, one config x 5 seeds.
- (open; now in TODO or Settled) **Make submit scripts refuse to overwrite a completed fold directory.** A grid once
      silently overwrote a sweep result that shared an `OUT_DIR`.

#### Housekeeping

- (open; now in TODO or Settled) Delete the local git tag `backup-pre-msg-rewrite` on Oak once the rewritten history is
      confirmed good.
- (open; now in TODO or Settled) `Data/ENCODE/K562/interm_ENCFF790GFL.se.filtered.sorted.bam/` is a stray directory.
- (open; now in TODO or Settled) `scripts/0.3.make_training_bw.sh` does not sort its experiment bedGraphs before
      merging; works today only because inputs happen to be sorted.

#### Carried over from the old root TODO.md (2026-09-12)

Open items that were stranded in `TODO.md` and `todo/h3k27ac-model.md`. Figure status and the
p300-model CV key values moved to `reference/PUBLICATION_FIGURES.md`; completed sections were
dropped, since git history and `.living/` hold them.

- (open; now in TODO or Settled) **Revisit the p300+ definition across all figures.** `EP300_peak_overlap` (from
      `finemo_peaks_all_chr.chromatin_annotations.tsv`) flags any 1 kb window overlapping a
      p300 peak call, and **~40% of "p300+" elements by that definition have
      `true_logcounts = 0`**, most likely because the peak overlaps only the window edge,
      outside the 500 bp window where reads are counted. Consider replacing it with the top
      20% of elements by observed p300 counts. Figures that would change: `model_comparison`
      (1d), `k562_multimodal_comparison` (2 K562), `training_region_comparison` (S1b),
      `transferability_bar` (2d), every `*_p300plus.pdf`, and the `mean_predictions` scatters.
- (open; now in TODO or Settled) **Two BioRender architecture schematics never started**: Fig 1a (BPNet) and Fig 2a
      (multimodal).
- (open; now in TODO or Settled) **ATAC ChromBPNet true CV performance.** Fig 1d still carries the manuscript's
      Pearson r = 0.70 as a placeholder, with no p300+ number at all.
- (open; now in TODO or Settled) **Section 3 and later publication figures.** MoDISCo motif logos, FiNeMo hits, and the
      motif spacing/pair experiments. Two pieces are already done: the individual-motif
      insertion violin plot (`scripts/plot_individual_motif_insertions.py`) and the motif-pair
      heatmaps for max log2FC and synergy (`scripts/plot_motif_pair_heatmaps.py`). The logos
      and hit panels are not.
- (open; now in TODO or Settled) **Compare the GM12878 MoDISCo motifs to the K562 set.** Raised in the GM12878
      handover and never done. 26 GM12878 motifs exist
      (`modisco/max_seqlets_250k_30_10_0/logos/`). A shared motif vocabulary across the two
      cell types is the sequence-side counterpart to the transfer results.
- (open; now in TODO or Settled) **Confirm the identity of `REPEAT_G`, `NF2L_NFE` and `ELF`** in the v2 FiNeMo motif set.
      Raised in the FiNeMo TF-analysis handover and unresolved.
- (open; now in TODO or Settled) **hashFrag chromosome splits** to reduce train/test sequence-similarity leakage. Never
      assessed; every result in the project rests on plain chromosome holdout.
- (open; now in TODO or Settled) **Update `P300_INTERACTORS`** in `plot_finemo_composite_figure.py` after a
      literature/BioGRID review. The composite figure's p300-interaction star depends on it.
- (open; now in TODO or Settled) **SHAP on the K562 ATAC multimodal model.** Never submitted. Overlaps W4 above, but W4
      specifies the residual-trained model, so this is a separate run.
- (open; now in TODO or Settled) **Reweight the accessible-but-unacetylated quadrant.** Those elements are already
      positives with near-zero targets, so up-weight rather than relabel. Aims directly at
      the over-predicted tail. Survivor of the training-composition family after GC-matching
      was refuted (D-18).
- (open; now in TODO or Settled) **Sweep `negative_ratio` (0.1) and `--max-negatives` (50,000).** Both inherited from the
      p300 setup and never swept for H3K27ac.
- (open; now in TODO or Settled) **Bin the profile target to nucleosome scale (~50-150 bp).** H3K27ac has no meaningful
      1 bp structure, so the profile head is fitting noise. Cheap to run and expensive to
      develop: touches `multimodal_bpnet.py` and `train_multimodal_bpnet.py`, which the p300
      models share. **Read alongside the DNase multi-task arm above, which attacks the same
      dead profile head from the other direction**, by giving it a learnable target instead of
      coarsening an unlearnable one. Doing both would confound them.
- (open; now in TODO or Settled) **Re-check the counting-window choice on the 5-prime target.** Less smearing means less
      neighbour bleed-through, so the contamination penalty at +/-1000 may differ from what
      the fragment-extended target showed.

#### Where the project goes next (moved from TODO.md 2026-09-12; written 2026-09-07)

Every model so far fits one cell type's accessibility-to-signal mapping. Accessibility is the
part that generalises, so a transferred model decays toward the ATAC-only floor by
construction. The one arm that transferred well (K562-trained p300, +0.207 over the target's
floor, 83% of its local advantage) is also the one trained on 2.3x more reads in peaks, so
portable sequence grammar is learnable here and the models that failed to learn it were the
thinner ones. Order below reflects that.

##### Running now

- (done) **Depth-subsampling causal test. REFUTED 2026-09-08**: the depth-matched K562 model
      still transfers at +0.202 against the full-depth +0.207 (paired -0.005, *p*=0.40), and
      peak-set geometry is excluded too. The asymmetry survives equalising volume, so a third
      cell type is now required rather than merely desirable. Original plan: Subsample K562 EP300 to GM12878's 30.0M mapped reads
      and 21,068 peaks, retrain the multimodal p300 model, re-run the transfer 2x2 (`2.27`).
      If the subsampled K562 model stops transferring, training-signal volume explains the
      asymmetry and the answer is more sequencing. If it still transfers, something else about
      K562 as a training cell type does, and a third cell type becomes urgent.
      Five folds, no new tooling beyond the subsampling and bigwig rebuild.

##### Top section

- (open; now in TODO or Settled) **Multi-cell-type joint training, leave-one-cell-type-out.** K562 + GM12878 + TeloHAEC,
      shared sequence trunk, per-cell-type accessibility input, held-out cell type as the
      evaluation. Measures the deployment case directly instead of inferring it from a two-way
      transfer. Needs a cell-type-invariant target scale, which the quantile item supplies.
- (open; now in TODO or Settled) **Quantile-normalised targets within cell type, plus a ranking loss** in place of log1p
      MSE on counts. ABC consumes an ordering; the measured downstream failure is compressed
      dynamic range on the strongly acetylated elements; nothing tried so far targets spread.
      Also makes targets commensurable across cell types, which joint training requires.
- (open; now in TODO or Settled) **Pretrained sequence trunk** (Borzoi or Enformer embeddings, frozen first, then
      fine-tuned) replacing the 8-layer dilated stack. Those representations were fit across
      hundreds of cell types, so cross-cell-type grammar is imported rather than learned from
      two experiments.
- (open; now in TODO or Settled) **DNase-input panel across six or more cell types** instead of ATAC across three. The
      confound that ruled DNase out was about mixing assays inside one comparison; a
      DNase-only panel is internally consistent and doubles the available cell types.

##### Lower confidence, cheap once the above exist

- (open; now in TODO or Settled) Gradient-reversal on a cell-type classifier reading the sequence embedding, penalising
      anything cell-type-identifiable in the sequence representation. Meaningless before joint
      training exists.
- (open; now in TODO or Settled) Multiplicative two-tower factorisation: sequence predicts regulatory potential,
      accessibility predicts availability, output is their product. The additive trunk lets
      accessibility dominate. The gate tested sequence modulating accessibility and failed;
      the reverse conditioning is untested.
- (open; now in TODO or Settled) Multi-task across marks available in six cell types (H3K4me1, H3K27me3, CTCF) rather
      than p300, which exists in two. The p300 multi-head died on p300's learnability at
      specific elements, not on the multi-task idea.

##### Do not reopen without new evidence

- Further single-cell-type architecture variants. Wide receptive field, fragment channels,
  the sequence gate, the asymmetric loss and GC-matched negatives have all either failed
  outright or failed to travel.

##### Structural blocker

- (open; now in TODO or Settled) **A downstream benchmark in a second cell type.** CRISPR data exists only for K562, so
      the terminal metric can only test the transfer direction that already fails. MPRA or
      eQTL data in a second cell type would be worth more right now than another model.

### [2026-10-05] Pre-register the ENCODE-rE2G test of the GC and CTCF corrections (F-055)

**Context**: F-055 found element GC (+0.0139 [+0.0069, +0.0206]) and CTCF motif count
(+0.0083 [+0.0013, +0.0153]) improve scE2G on the K562 CRISPR benchmark, but F-054 chose both
from residuals on those same pairs. Before believing them, test them in a model that played no
part in the choice: ENCODE-rE2G with ENCODE bulk K562 ATAC on the July candidate regions
(`2026_0824_H3K27ac_model/scripts/4.62`, `4.32`, launcher `log/gcctcf_launch.sh`). This entry is
written before the run.

**Decision**: Same feature definitions as scE2G (GC of midpoint +/- 250 bp; JASPAR MA0139.2 sites at
p < 1e-4 by TOBIAS's scan), added to the 10 ATAC-only rE2G features, `_all` feature tables, rE2G's
own chromosome CV, paired bootstrap `4.34` (2,000 resamples).
- **Primary**: `ATAC_GCCTCF` - `ATAC`. PASS if the 95% interval excludes 0 on the positive side
  AND both coefficients are negative in the full model, as in scE2G. Anything else is a FAIL.
- **Secondary, descriptive only**: `ATAC_GC` - `ATAC`; `ATAC_CTCF` - `ATAC`;
  `ATAC_H3K27ac_GCCTCF` - `ATAC_H3K27ac` (expected smaller: if GC and CTCF stand in for "open but
  inactive", observed H3K27ac already carries it); `ATAC_GCCTCF` - `ATAC_H3K27ac` (share of the
  gap to observed H3K27ac).

**Alternatives considered**:
- Another CRISPR set: none in hand (K562 is the only benchmark, Settled 2026-10-01); still the
  only test free of the selection, and the question stays open for Maya.
- Re-pick features in rE2G by its own error analysis: rejected, it repeats the selection.

**Rationale**: A different model family, ATAC sample and candidate set, with the comparison and
pass rule fixed in advance, removes the freedom to choose features or comparisons after seeing the
result.

**Consequences**: The test shares about 470 positives with F-055, so a PASS shows the gain
transfers across model families, not that it is free of selection on these pairs. On PASS, keep
GC + CTCF as candidate scE2G defaults pending another CRISPR set. On FAIL, F-055 stays a
scE2G-on-Xu result and is not adopted.

**Tags**: pre-registration, re2g, gc, ctcf, crispr-benchmark, k562, sce2g

### [2026-10-05] Pre-register the scE2G confirmation of the p300 feature screen (7.2)

**Context**: The 7.2 screen (`A2/plans/7_feature_screen.md`) ranked ~1,800 deployable element
features against observed K562 p300 beyond Xu ATAC, distance, promoter, mappability, GC and CpG,
over 155,134 Xu elements, without CRISPR labels. One feature passed the shortlist rule: the share of
Xu insertions from fragments < 100 bp (partial rho 0.19 [0.18, 0.20] on replication chromosomes).
Out-of-fold gradient-boosted trees recover 0.25-0.26 of the residual rank variance (elastic net
0.10-0.12). Written before the scE2G run.

**Decision**: One scE2G run (`scripts/7.3`, `config/spec_screen.tsv`, `log/screen_launch.sh`),
paired bootstrap `4.2` (2,000 resamples). Arms: n6 + GC (reference, the F-055 arm); n6 + GC +
short-fragment share; n6 + GC + composite (1 + percentile of the out-of-fold tree prediction).
- PASS for an arm if its paired interval vs n6 + GC lies above 0. Two arms, no multiplicity
  correction; both are reported whatever the outcome.
- Descriptive: each arm vs n6; observed p300 (F-056 run) vs the composite arm, the share of the
  ceiling reached.

**Alternatives considered**:
- Add the 5 strongest non-passing features (GATA motifs ~0.12): rejected, the rule was fixed
  before reading the output and the erythroid footprint at 0.116 was null on CRISPR.

**Rationale**: Features chosen without CRISPR labels and tested once against a fixed rule carry no
selection on the benchmark; the GC reference carries F-055's selection, so gains are read on top
of it.

**Consequences**: On PASS, the arm's feature becomes a candidate scE2G default pending another
cell. The composite is fitted on K562 p300, so outside K562 it is a transferred model; a PASS for it
is a K562 result. On FAIL, p300-beyond-ATAC is not a usable proxy for this benchmark at the
recovered strength, and the screen line stops.

**Tags**: pre-registration, sce2g, feature-screen, fragments, crispr-benchmark, k562

### [2026-10-07] Pre-register the scE2G confirmation of the nested CRISPR feature discovery (plan 8, Phase 4)

**Context**: `A2/plans/8_crispr_feature_discovery.md`. Feature selection inside nested chromosome CV
(`scripts/8.3`; 5 outer chromosome groups; base = n6 with ARC split into ABC + exp(Kendall), 0.627)
gained +0.027 [+0.011, +0.042] (forward selection) and +0.023 [+0.011, +0.034] (CpG-island family) on
held-out chromosomes. Both beat all 100 runs of a null that shuffles candidate vectors within TSS-distance
bins (null max +0.002 / +0.004; `8.4`, `8.5`). The within-gene label shuffle was rejected as a null
(Maya: 218 of 429 positives in genes it cannot change; it moves positives distal). With log TSS distance
in the base the gains held (+0.032, +0.025) with the same picks. Picked in >= 4 of 5 outer folds:
log distance to the nearest CpG island (5/5), GATA motif cluster GATA6 / Gata3 / TRPS1 (5/5), IKZF2 motif
(4/5). Written before the run.

**Decision**: One scE2G run (`scripts/8.6`, `config/spec_discovery.tsv`, `log/discovery_launch.sh`),
features as 1 + percentile over all Xu elements, paired bootstrap `4.2` (2,000 resamples).
- **Primary**: n6 + the three features vs n6. PASS if the paired interval lies above 0. Predicted gain
  from the nested estimate: about +0.02 over n6 (the split base sits 0.006 below n6).
- Secondary, descriptive: n6 + CpG-island distance alone; n6 + GC (F-055 arm); n6 + GC + the three, vs
  n6 and vs n6 + GC (are GC and CpG-island distance the same signal?).

**Alternatives considered**:
- Adopt the nested estimate without a pipeline run: rejected, the official run checks the features
  survive scE2G's own joins and encodings.
- Include the L1 model (~100 terms): rejected, its terms are not stable across folds.

**Rationale**: Features were chosen without seeing the held-out chromosomes they were scored on, and
the null measures what searching the library yields by chance; the official run is the deployable form.

**Consequences**: The run reuses the K562 CRISPR pairs, so a PASS confirms the pipeline form of an honest
internal estimate, not an independent replication; a second CRISPR set or cell is still the external test.
On PASS, the three features become candidate scE2G defaults. On FAIL, investigate the difference between
the replica and pipeline joins before concluding.

**Tags**: pre-registration, sce2g, feature-discovery, nested-cv, cpg-island, gata, crispr-benchmark, k562

### [2026-10-12] Pre-register the eQTL transfer test of the F-059 features in unseen immune cell types (plan 9)

**Context**: F-059's three features (log distance to the nearest CGI, GATA cluster, IKZF2 motif) were
chosen and confirmed on the K562 CRISPR benchmark only; F-060 to F-062 suggest CGI distance may
transfer only in part. Maya (2026-10-12): use the eQTL benchmark (eQTLEnrichment, eQTL Catalogue v7
fine-mapped, fine-grained cell-type matching, every other cell type dropped), all clusters, both
readouts; no features computed across co-analysed cell types. Written before any cluster is scored.

**Decision**: `2026_1002_sce2g_tf_activity/scripts/9.0` -> `9.1` (array, one task per biosample) -> `9.2`.
- Arms, each scE2G's logistic regression fitted on all K562 CRISPR pairs (2026_1007 discovery run
  table): n6; n6 + CGI distance; n6 + the three. Applied to the 31 biosamples of
  `config_predictions_hd_fine.tsv` (scE2G multiome_powerlaw_v3 cluster predictions), new features
  as 1 + percentile over each cluster's own elements. GATA and IKZF2 counts from the genome-wide FIMO
  scan (4.44) for training and application alike (not BINDetect, which clusters lack).
- Scores on scE2G's quantile scale with the TPM >= 1 rule, as run_e2g_cv.py; eQTLs PIP >= 0.5,
  distal noncoding, gene universe, as eQTLEnrichment; counts pooled over matched pairs.
- Readouts: R1 enrichment at linking recall 0.10; R2 mean enrichment over recall 0.02-0.20 (0.9 x
  the n6 maximum if n6 does not reach 0.20). Paired bootstrap over eGenes, 1,000 resamples.
- **Primary**: n6 + the three vs n6, pooled over all matched pairs except K562. PASS if both
  intervals lie above 0; PARTIAL if one does; WORSE if either lies below 0; otherwise FAIL.
- Secondary: n6 + CGI distance vs n6, same rule. Descriptive: per sample (PBMC, BMMC22, GM12878),
  per eQTL cell type, K562 alone (in-sample), the published v3 scores, K562 LOCO AUPRC of the FIMO
  version.

**Alternatives considered**:
- Rerun scE2G per cluster: rejected; the saved predictions carry every n6 feature, so rescoring
  offline is exact for these arms and costs minutes.
- Normalise the new features across clusters: rejected (Maya), predictions must not depend on
  which other cell types are analysed.

**Rationale**: eQTLs in cell types the model was not trained on are the only available test
outside K562; the pipeline's own definitions keep the numbers comparable to the scE2G paper's.

**Consequences**: PASS -> the three (or CGI distance, if only it passes) become candidate scE2G
defaults, pending a CRISPR set in a second cell. FAIL or WORSE -> the F-059 gain is K562-specific;
record it and stop adding K562-selected features without a transfer check.

**Tags**: pre-registration, sce2g, eqtl, transfer, cpg-island, gata, ikzf2, immune

**Amendment [2026-10-12], before any cluster was scored**: the FIMO scan has no IKZF2 site (MA2326.1
is 6 bp wide; a perfect match has p ~ 2.4e-4, above the scan's 1e-4 cutoff), so IKZF2 is counted as
exact consensus matches (AGGAAG / CTTCCT) within the element midpoint +/- 250 bp, for training and
application alike. GATA stays on FIMO (Spearman 0.957 with the BINDetect counts over Xu elements).
The arms are refitted with this definition; everything else as above.
