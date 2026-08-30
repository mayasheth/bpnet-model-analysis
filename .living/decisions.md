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
