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
