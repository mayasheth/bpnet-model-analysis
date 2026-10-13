# EP300_BPNet: scE2G activity evals (CRISPR feature discovery, eQTL transfer)

**Goal:** improve scE2G's activity term with deployable features that hold in unseen cell types. Selection on K562 CRISPR, transfer checked on fine-mapped eQTLs in matched immune cell types (Maya, 2026-10-12: no features computed across co-analysed cell types).

## State (verified 2026-10-13, evening)
- Repo `$OAK/Users/sheth/EP300_BPNet`, main at the plan 10 commit, 28 commits UNPUSHED (never push unasked). Committed today: `97e61dd` (plans 7-9), `4c7577d` (report 7 Part C, F-057 to F-065), then plan 10. Left untracked on purpose: `2026_0824_H3K27ac_model/data/*` (GB), September per-region `results/*.regions.tsv.gz` and other >5 MB tables, `*.new` report drafts, `2026_1004_sce2g_tf_error_analysis.pairs.tsv.gz`.
- No jobs running.
- Results:
  - F-059: three K562-CRISPR-selected features, ~+0.021 out of sample on K562 CRISPR.
  - F-063 to F-065: plan 9 enrichment readouts are mostly RNA-driven; F-063 inconclusive for activity effects.
  - F-066 (plan 10, pre-registered decisions [2026-10-13]): within expressed genes (non-promoter candidates <= 1 Mb, E1 pooled AUPRC, E2 within-gene AP), the ABC swap moves nothing: E1 -0.0002 [-0.0023, +0.0021], E2 +0.0007 [-0.0027, +0.0040] -> NOT calibrated. n6 and ABC rank a gene's candidates alike (Spearman 0.994); eQTL positives are the proximal top element (median 22 kb). n6_three / n6_cgi never unblinded.
- Conclusion: eQTLs cannot test activity features of F-059's size. Only a CRISPR set in a second cell can, which needs the Settled "K562 CRISPR only" rule lifted.

## Next steps
1. Maya: lift the "K562 CRISPR is the only benchmark" rule or not (TODOLIST Waiting). If lifted, find a second-cell CRISPR set and score the three on it.
2. If not lifted: pick the next direction from TODOLIST Next (Corgi ideas 1-2, H3K4me1 models).
3. Report 7 Part C does not yet cover F-066 (one table; add with the next report update).

## Open questions for Maya
- Lift the K562-only benchmark rule?
- EngreitzLabVault location for the report 7 pointer note?
- [EngreitzLab/scE2G#118](https://github.com/EngreitzLab/scE2G/pull/118) review pending.

## Pointers
- Plan 9 machinery: `A2/scripts/9.0`-`9.3`; plan 10: `10.1` (per-cluster candidates, `$SCRATCH/eqtl_transfer/withingene/`), `10.2` (readouts, `--blind`, calibration), `10.3` (diagnostics, unblinded arms only).
- Score conventions: quantile = rank of nonzero scores over all pairs; multiome zeroes TPM < 1.
- Replica: `A2/scripts/8.0`; CRISPR pair -> Xu element mapping: `.living/learnings.md` 2026-10-09.

## Read on resume
1. `$OAK/Users/sheth/EP300_BPNet/CLAUDE.md`
2. `todo/TODOLIST.md`, Waiting: the K562-only rule question
3. `.living/findings/predicting-regulatory-element-function-at-scale.md`, F-063 to F-066 at the end
4. `.living/decisions.md`, the [2026-10-13] plan 10 entry
