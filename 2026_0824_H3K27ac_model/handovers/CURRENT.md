# EP300_BPNet: scE2G activity evals (CRISPR feature discovery, eQTL transfer)

**Goal:** improve scE2G's activity term with deployable features that hold in unseen cell types. Selection on K562 CRISPR, transfer checked on fine-mapped eQTLs in matched immune cell types (Maya, 2026-10-12: no features computed across co-analysed cell types).

## State (verified 2026-10-13)
- Repo `$OAK/Users/sheth/EP300_BPNet`, main at `045a26f`, 25 commits UNPUSHED (never push unasked). Uncommitted: findings F-057 to F-065, decisions (plan 9 pre-registration [2026-10-12] + IKZF2 amendment), learnings, TODOLIST, `analysis/ANALYSIS_MANIFEST.md`, A2 scripts `8.0`-`8.10`, `9.0`-`9.3` + submit scripts, `log/eqtl_transfer_launch.sh`, results; H3K27ac dir `scripts/5.5`, `5.6`, report 7, `figures/fig40`-`fig49`.
- No jobs running.
- Results:
  - F-059: three K562-CRISPR-selected features, +0.036 in-sample, ~+0.021 out of sample on K562 CRISPR.
  - F-063: they do not transfer. eQTLs, 62 matched immune pairs: R1 (enrichment at recall 0.10) -0.34 [-1.15, +0.11], R2 (mean enrichment, recall 0.02-0.20) +0.11 [-0.18, +0.43] -> FAIL; CGI alone WORSE. No gain on K562's own eQTLs.
  - F-064: scATAC - multiome R1 -2.63 [-3.48, -1.79], R2 -1.72 [-2.27, -1.12] for -0.049 AUPRC; mostly the RNA side (revised by F-065).
  - F-065 (activity-only control, arm `n6_abc`, `results/2026_1013_eqtl_activity_control.*` = all arms): ARC-E2G -> ABC (-0.050 AUPRC) moves R2 only -0.35 [-0.47, -0.20], R1 -0.40 [-0.78, +0.01]. eQTL enrichment barely sees activity changes; F-063 is inconclusive (predicted +0.15 R2 inside its interval), not a refutation. F-063/F-064 carry dated revision notes.
- Report 7 updated 2026-10-13: Part C Fig. 16 (`figures/fig50`, `scripts/5.6`) covers F-063 to F-065; title, summary, open questions, datasets, methods updated.

## Next steps
1. Maya to decide the transfer benchmark for activity features (TODOLIST P1): a CRISPR set in a second cell, or an eQTL readout restricted to expressed genes where the activity term decides the ranking, calibrated with the `n6_abc` arm.
2. Commit the plan 8/9 and report 7 work when Maya asks.

## Open questions for Maya
- Which transfer benchmark for activity features (see Next steps 1)?
- EngreitzLabVault location for the report 7 pointer note?
- [EngreitzLab/scE2G#118](https://github.com/EngreitzLab/scE2G/pull/118) review pending.

## Pointers
- Plan 9 machinery: `A2/scripts/9.0` (models, features, variants; `$SCRATCH/eqtl_transfer/`), `9.1` (score one cluster per arm), `9.2` (curves, readouts, eGene bootstrap), `9.3` (published prediction file as an extra arm). Minutes per arm.
- Score conventions: quantile = rank of nonzero scores over all pairs; multiome zeroes TPM < 1, scATAC has no TPM rule (checked against qnorm per cluster).
- Replica: `A2/scripts/8.0`; CRISPR pair -> Xu element mapping: `.living/learnings.md` 2026-10-09.
- Lessons sent to `~/Documents/my-scripts-for-me/learnings/INBOX.md` (2026-10-08 outage, 2026-10-13 short motifs and FIMO).

## Read on resume
1. `$OAK/Users/sheth/EP300_BPNet/CLAUDE.md`
2. `todo/TODOLIST.md`, Now: the eQTL positive-control and transfer-check items
3. `.living/findings/predicting-regulatory-element-function-at-scale.md`, F-063 to F-065 at the end
4. `.living/decisions.md`, the [2026-10-12] plan 9 entry and amendment
5. `A2/scripts/9.2.eqtl_transfer_gather.py` docstring (readout definitions)
