# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

### [2026-08-24] Center H3K27ac training windows on candidate elements, not ChIP peaks

**Context**: H3K27ac is deposited on the nucleosomes flanking a regulatory element, not on the nucleosome-free element itself. Measured meta-profile over 30k DNase candidate elements confirms this: shoulders at -275/+275 bp with a central dip, 3.48x over the distal plateau. The existing p300 pipeline centers windows on `start + summit` from a narrowPeak.

**Decision**: Train on `reference/K562_DNase_candidate_elements.narrowPeak`, centered on the element midpoint, with a counting window wider than the element.

**Alternatives considered**:
- Center on H3K27ac peak summits — rejected: the summit of a broad acetylation peak sits on a flanking nucleosome, so windows would be centered off the element and inconsistently so.
- Keep the element-sized window used for p300 — rejected: a +/-250 bp window captures the shoulders but cuts deeply into real signal, which only reaches background around +/-2000 bp.

**Rationale**: The scientific question is about elements, so the element must define the coordinate system. Conveniently the candidate-elements file has its summit column set to `width/2`, so the trainer's existing `is_peak=True` path already centers on the midpoint with no code change.

**Consequences**: Negatives now mean "GC-matched genomic background" rather than "non-peak", and the element set spans active and inactive elements, so a large part of any reported correlation is the dead-vs-active contrast. This is why evaluation must be stratified by signal level.

**Tags**: h3k27ac, window-selection, element-centric, training-design

---

### [2026-08-24] Counting window is a trade-off between signal and neighbour contamination

**Context**: H3K27ac signal has no clean saturation point — it decays slowly and only reaches its distal plateau around +/-2000 bp. Meanwhile the 150,528 candidate elements cluster, so wide windows start containing other elements.

**Decision**: Test +/-500 and +/-1000 head to head rather than picking one. Expect +/-1000 to be the headline model and +/-500 the one trusted for attribution work.

**Alternatives considered**:
- A single wide window (+/-2000) — rejected: the inter-replicate ceiling on active elements saturates by +/-1000 (0.785 vs 0.797 at +/-2000), so wider buys ~1 point of ceiling for 41% contamination.
- A single narrow window (+/-500) — kept as an arm rather than the sole choice: zero contamination, but the ceiling is 2.4 points lower and it leaves real signal out.

**Rationale**: Neighbour contamination rises 0% -> 9.9% -> 19.9% -> 41.5% across +/-500/750/1000/2000, while the ceiling gains flatten after +/-1000. The two questions also want different things: prediction accuracy favours the wider window, while motif syntax is actively harmed by contamination, since a window containing a neighbour lets the model earn credit from the wrong element's motifs.

**Consequences**: Every window requires its own model geometry (`in_window = out_window + 2*557`), so this is a retraining sweep, not a re-evaluation. Attribution results from the +/-1000 model need to be read with the 20% contaminated subset in mind.

**Tags**: h3k27ac, window-selection, contamination, ceiling, trade-off

---

### [2026-08-24] Keep the profile head, down-weighted, rather than removing it

**Context**: H3K27ac counts are the quantity of interest, and the bimodal flanking profile is arguably not a meaningful "read profile" at all: the target is fragment-extended coverage, so MNLL (which expects multinomial read counts) sees ~250x inflated totals and near-identical adjacent positions. Profile MNLL sits at ~2800 here versus ~500 for the p300 runs.

**Decision**: Keep the profile head and control it with `count_loss_weight` rather than editing the loss to remove it.

**Alternatives considered**:
- Remove the profile head outright — deprioritized on evidence, not principle. At `count_loss_weight=10000` the counts term already takes ~91% of the gradient, making that run a close proxy for counts-only, and it scored 0.398 versus 0.407 at parity (weight 1000). Removal would land near 0.398, i.e. no better.
- Grid the profile-drop against mode and window — rejected: doubles the grid to 60 GPU jobs to test a lever the weight sweep already argues against.

**Rationale**: Reproducing the central dip is a useful diagnostic that the model learned H3K27ac structure rather than "accessible implies acetylated", and it costs nothing to keep once the weight is calibrated.

**Consequences**: `fconv` still consumes model capacity even with a small gradient, a second-order cost not measured. If ATAC and window width also plateau, revisit this alongside the more fundamental coverage-vs-read-counts mismatch.

**Tags**: h3k27ac, profile-head, loss-weighting, architecture, bpnetlite

**Update 2026-09-03 — decision stands, and the reasoning was better than the evidence then
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
- Keep overall Pearson as the headline — rejected: it ranked the two targets backwards (see F-002).
- Partial correlation controlling for `atac_pred` — equivalent in spirit but less directly interpretable; the residual form states plainly "does the model's departure track the true departure".
- Use the raw ATAC track as the baseline — rejected: then the residual includes everything a linear read of the track misses, which flatters any model that merely learns a better ATAC transform.

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
- Count both mates — rejected on principle despite winning the top-quintile ceiling by
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
- Keep `n_layers` 8 — rejected: the gain is significant, replicated, and on the reporting standard.
- Adopt it for sequence-only too — rejected: sequence gains only on all elements, which is the dead-vs-active contrast.
- Wait for the ATAC-only arm before deciding — the arm is running and will refine the mechanism, but it cannot overturn a replicated multimodal gain.

**Rationale**: The useful long-range information is in the accessibility track. A wider window lets the model read accessibility over a larger neighbourhood; sequence-only cannot exploit 4 kb of sequence, so for it the extra window is noise.

**Consequences**: The deployed model geometry changes to trimming 2093 / in-window 5186, ~2-3x the training time per fold. Transfer and deployment comparisons were run at `n_layers` 8 and should be re-checked. Recorded reversal: the earlier decision generalised from one arm of a two-arm experiment.

**Tags**: architecture, receptive-field, multimodal, reversal, accessibility

---

### [2026-09-03] Fragment-size accessibility channels are 5' insertion counts and include the flat track

**Context**: Fragment length distinguishes nucleosome-free from nucleosome-occupied DNA and flat coverage discards it. An earlier pair of channels existed but used `genomecov -bg` over the full fragment interval, the read-length-dependent smear already removed from the flat track, and covered only half the fragments.

**Decision**: Five channels `[all, sub ≤139, mono 140–329, di 330–620, poly ≥621]`, every one a single-base Tn5 insertion count, with `all` being `atac_5p.bw` itself. Bin edges sit in the troughs of the measured fragment-length distribution.

**Alternatives considered**:
- Full-fragment coverage for the stratified channels — rejected: it would place two incompatible accessibility conventions in one input tensor, and a mono-nucleosomal interval marking occupancy directly is not worth that.
- The original `sub`/`mono` edges — rejected: they sat inside the modes and covered half the fragments.
- Stratified channels only, without the flat one — rejected: including `all` is what makes the input a provable superset, since the first convolution can zero the four bins and reproduce the baseline exactly.

**Rationale**: The bins partition the flat track exactly — zero discrepancy across 545,661,218 insertions — so any difference is added information rather than a changed input. The Tn5 shift was measured against the existing track (r = 1.0000 at +4/−5) rather than assumed, which also proved the PE BAMs and tagAligns hold the same reads.

**Consequences**: Gains +0.0135 on the top quintile (p = 0.0034). `atac_sub.bw` and `atac_mono.bw` are superseded and must not be used with a 5' model. GM12878 replication needs its PE BAMs downloaded, since fragment length lives in TLEN.

**Tags**: atac, fragment-length, accessibility, input-design, superset, tn5

---

### [2026-09-03] Inject predicted H3K27ac into ABC as a painted bigWig, with qnorm left on

**Context**: The ABC activity term is `geomean(accessibility, H3K27ac)` computed from read counts over candidate regions. Our model emits a per-element scalar from 5' end counts in a ±500 bp window, which is not obviously commensurate with read counting over an element.

**Decision**: Write a bigWig with each candidate region painted at `predicted_counts / width` and pass it in the `H3K27ac` column. Keep `use_qnorm: True`. Assemble the genome-wide track from the five fold models, each applied only to the chromosomes it held out.

**Alternatives considered**:
- Patch ABC to accept a precomputed activity column — rejected as unnecessary: `neighborhoods.py:count_bigwig` already sums bigWig values per region, so a painted track is counted exactly as a real one.
- Turn qnorm off — rejected: it is what makes the injection scale-free, and the ABC score thresholds are calibrated on qnorm'd values.
- Paint the prediction itself, so the region sum scales with width as real read counts do — kept as `--paint density`; the default makes ABC's sum recover the model's prediction exactly, since the model predicts a fixed ±500 bp window regardless of element width.
- A five-model ensemble for the cross-cell-type arms, where leakage is not a concern — rejected: it would give those arms an ensembling advantage the same-cell-type arms cannot have.

**Rationale**: `run_qnorm` is rank-based, mapping each region's within-sample quantile onto the K562 reference, so only the rank order of the injected values matters. That dissolves the units mismatch. Validated on chr22: predicted-vs-observed H3K27ac Spearman 0.716 all regions and 0.550 top quintile, against 0.633 and 0.469 for raw ATAC, so the prediction is a better proxy than the accessibility it would replace.

**Consequences**: ABC scores our models on the ATAC-derived candidate regions while they were trained on the DNase-derived set, a train/score element mismatch inside the comparison; a matched run is training. Every arm must share one region set, so `Peaks/` is pre-populated from the completed July run and must be copied after the bigwigs exist or Snakemake re-runs region calling.

**Tags**: abc, crispr-benchmark, activity, qnorm, injection, leakage

---

### [2026-09-03] Regression gates on inference code compare within a tolerance, never byte-identically

**Context**: `2.15` is shared with the transfer and residual-grid results, so generalising it needed a gate. The first gate required the per-fold table to be byte-identical to a stored one, and it failed on differences in the 4th decimal place while region counts matched exactly.

**Decision**: Assert exact equality only on deterministic quantities — fold set, config labels, `n` per fold, and the presence of every reference column — and compare metrics within 1e-3.

**Alternatives considered**:
- Keep byte-identity and pin the GPU model — rejected: it makes the gate depend on scheduling.
- Drop the gate — rejected: it is the only thing standing between a refactor and silently changed published numbers.

**Rationale**: cuDNN convolution is not bit-reproducible across GPU models, so re-scoring the same weights on a different node moves the 4th decimal place with the code untouched. The tolerance is chosen from the science's noise floor — between-fold sd is 0.041-0.046 — rather than from float precision, so a real behavioural change moves numbers by far more than the gate absorbs.

**Consequences**: `2.18.compare_perfold_tables.py` implements it and states the tolerance and its justification in its own docstring, so the next reader does not tighten it back to zero. Added columns are reported and skipped; missing columns fail.

**Tags**: testing, regression, determinism, tolerance, evaluation

---

### [2026-09-03] Test-time reverse-complement averaging stays opt-in

**Context**: Averaging each prediction with its reverse complement at inference costs one extra forward pass and no retraining, and it improves every model: +0.0162 sequence only, +0.0081 multimodal, +0.0016 for ATAC only, which is the negative control and the only non-significant row.

**Decision**: Ship it behind `--rc-average`, off by default.

**Alternatives considered**:
- Make it the default immediately — rejected for now: it would shift every number in the report by a small amount, so past and future tables would not be comparable unless all are re-scored together.
- Leave it unimplemented — rejected: it is the cheapest measured gain on the project.

**Rationale**: The gain is real but small relative to the comparisons being made, and comparability across the report matters more than a few thousandths until the tables are re-scored as a set.

**Consequences**: Adopting it later means re-running every evaluation config in one batch. Listed as an explicit adopt-or-not decision in the TODO.

**Tags**: inference, reverse-complement, adoption, comparability

---

### [2026-09-03] Split the analysis into three reports along stability, not topic

**Context**: The single report has reached 17 sections and 13 figures, mixing stable characterisation with an actively churning set of architecture experiments.

**Decision**: Three documents, ordered by how often they change.

1. **Data characterisation** — what the data is and what is predictable in principle. H3K27ac position and width, the counting-window trade-off, both ceilings (counts and profile shape), ATAC fragment-length structure, ATAC-H3K27ac coupling across cell types, element derivation and panel caveats, ATAC input conventions. Every ceiling lives here and the other two reports cite it.
2. **Evaluation methodology** — short, concrete, framed as standing cautions with the incident that motivated each. Why the top quintile leads, what the residual metric means and its artifact controls, paired within-fold testing against a between-fold sd of 0.041-0.046, the transfer-versus-deployment evaluation DESIGN, and ABC/CRISPR as the downstream metric.
3. **Design decisions** — the workbench. Spine is a decision table of change / effect on the top quintile / verdict, with evidence below: input modality, residual versus total objective, receptive field, fragment channels, RC averaging, target definition, window size, `count_loss_weight`, the profile head, element derivation, and the transfer RESULTS with the deployment verdict.

**Alternatives considered**:
- Keep one report — rejected: the churn in category 3 forces re-reading stable material to find what changed.
- Split by audience, a short headline plus a technical appendix — rejected for a working project: it optimises for a reader who is not the one using this daily.
- Put transfer entirely in report 2 or entirely in report 3 — rejected: transferability is both an evaluation axis and a deployment decision, so the design goes in 2 and the results in 3.

**Rationale**: Splitting on churn rate isolates the part that moves. Reports 1 and 2 become citable references; report 3 is expected to change every week.

**Consequences**: Figure numbering restarts per report, so every legend and prose cross-reference is touched — mechanical, and the numbers manifest lint catches any reference that stops resolving. `render_report.py` resolves `outputs/numbers.json` relative to the report, so all three sharing one directory keeps one manifest. Deferred until the runs in flight land, so the split happens once.

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
