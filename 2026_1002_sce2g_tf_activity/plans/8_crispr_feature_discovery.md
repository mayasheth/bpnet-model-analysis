# Direct feature discovery for the scE2G activity term (draft, 2026-10-06)

Select features against CRISPR itself, with every selection step inside cross-validation, so the
gain we report is the gain of the selection procedure on chromosomes it never saw. p300 as an
intermediate target failed both ways: footprints moved p300 and not CRISPR (F-050), the 7.2
composite recovered 0.25 of p300 beyond ATAC and moved CRISPR by -0.001 (F-058). A feature that is
silent on p300 can still help CRISPR, so p300 drops out as a gate.

## Constraints that shape the design

- About 470 regulated pairs, clustered by gene. A paired-bootstrap interval on a single AUPRC
  difference is roughly +/-0.007, so gains under ~0.007 are not resolvable on this benchmark.
- Any feature chosen by looking at these pairs (F-054 -> F-055 GC) carries selection bias. The fix
  is to make selection part of what is cross-validated, then measure how much a procedure gains by
  chance with a permutation null.
- Each official scE2G run costs ~2.5 h, so screening needs a fast offline copy of scE2G's model.

## Phase 0: offline scE2G replica (half a day)

scE2G's model is a logistic regression on log(|x| + 0.01) features with leave-one-chromosome-out
(LOCO) predictions. Rebuild it in Python from the run's CRISPR-pair feature table and check it
reproduces the baseline n6 AUPRC (0.632) and n6 + GC (0.646) to within 0.002. Everything below
runs on this replica in seconds per model.

## Phase 1: candidate library, grouped before any label is read

- Element features: the 7.1 library (~1,800: composition, 879 motifs, phyloP, repeats, CpG
  islands, fragment lengths, insertion profile, cell-level accessibility, neighbourhood).
- Collapse redundancy first, label-free: motifs into ~100 archetype clusters (max over a
  cluster), dinucleotides into a few composition axes. Target ~150-250 candidates in ~10 families.
- Add activity-adjacent composites, each trained without CRISPR labels and out-of-fold by chromosome:
  the 7.2 p300 composite, plus the same on H3K27ac, and on K562 lentiMPRA / STARR-seq activity if
  ENCODE has them for these elements (sequence-intrinsic activity, a different notion from p300).
- Interaction terms with distance only for families that pass Phase 2 (keeps the search small).

## Phase 2: nested LOCO selection

Outer loop: hold out one chromosome group (5 groups, balanced for positives). Inner loop, on the
training groups only: run a selection procedure, refit, predict the held-out group. Pool held-out
predictions; compare to the baseline replica with the paired bootstrap.

Procedures, all pre-specified:
1. L1-penalised logistic regression over all candidates on top of the n6 features (penalty from
   inner CV).
2. Forward selection, at most 5 features, inner-CV AUPRC as the criterion.
3. Family tests: add one whole family at a time (group lasso or ridge within the family), so the
   unit of discovery is "fragment features" or "conservation", not one motif.
4. Gradient-boosted trees on the baseline logit's residual (shallow, monotone where sensible), to
   ask whether the activity term is missing non-linear structure rather than new inputs.

Outputs per procedure: held-out AUPRC gain with interval; how often each feature or family is
chosen across the 5 outer folds (the discovery readout); rank shifts by distance bin and by the
F-056 groups (accessible elements with and without erythroid sites).

## Phase 3: how much is chance

Rerun each whole procedure on 100 label permutations (shuffle labels among pairs within each gene,
keeping each gene's positive count and the distance structure). A procedure's real gain counts only
if it beats the 95th percentile of its permuted gains. This also gives the expected optimism of
the earlier GC selection, which the F-055 interval does not include.

## Phase 4: confirmation

Take at most 5 features or 2 families that are chosen in >= 4 of 5 outer folds and beat the
permutation null. Pre-register them with the nested-CV gain as the predicted effect, then run one
official scE2G training run (via `4.8`). An independent CRISPR set or a second cell is still the
only test free of this benchmark; the nested estimate is the best internal substitute.

## Decisions needed

- Permutation scheme (within-gene label shuffle proposed).
- Whether to fetch K562 lentiMPRA / STARR-seq for the extra composites (download).
- Compute: Phase 3 is ~100 x 4 procedures on the replica, a few hours on one node.

## Reused

7.1 library and UCSC files (`$SCRATCH/feature_screen/`, `$SCRATCH/ucsc_hg38/`), 7.2 composite, 4.2
paired bootstrap, 4.8 setup for the final run.
