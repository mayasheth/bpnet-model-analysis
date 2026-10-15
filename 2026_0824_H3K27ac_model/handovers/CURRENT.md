# EP300_BPNet: Corgi base done; review of reports 1-7 awaiting scope

**Goal:** improve the transferred (multi-cell, K562 held out) p300 activity track for ABC/rE2G on K562 CRISPR, and keep the eight reports trustworthy.

## State (verified 2026-10-14)
- Repo `$OAK/Users/sheth/EP300_BPNet`, main at `f03c503`, 32 commits UNPUSHED (never push unasked). Working tree clean except long-standing untracked data/large tables (deliberately out of git).
- No jobs running.
- F-067 (Corgi idea 1, `--panel-mean`, `A/scripts/1.42`): pre-registered PASS. ABC Corgi - MSE-only +0.0104 [+0.0010, +0.0204] pairs, [+0.0008, +0.0201] gene-level. rE2G +0.0117 [+0.0002, +0.0235] pairs, NOT resolved gene-level [-0.0050, +0.0275]. Corgi is now the multi-cell base: `A/models/p300_multicell_panelmean_holdK562_hw500_clw10{,_seed1,_seed2}`.
- New: `A/scripts/4.17` and `4.34` take `--cluster gene` (genes resampled with replacement). The default `pair` mode reproduces earlier tables exactly. A guard refuses one-gene-per-pair parses; tested both ways.
- Review 2026-10-14 (`/engreitzlab-review`, 9 reviewers): `.living/outputs/reviews/2026-10-14-all-reports.md` (committed `d477d25`; Mac copy `reviews/`), 82 findings, 23 Major. Reports NOT yet edited. Maya asked "did you fix the reports"; I offered: do text/doc fixes first, figure fixes second, and hold the text the reruns could change. Awaiting her answer.

## Next steps
1. Get Maya's scope for the review fixes (TODOLIST "Act on the 2026-10-14 review").
2. Urgent regardless of scope: copy the K562 PE ATAC BAMs (`$SCRATCH/atac_pe/*.pe.bam`, mtime 2026-08-25) and the baseline scE2G run (`$SCRATCH/sce2g_runs/2026_1002_k562_xu`) to Oak before ~2026-11-20 (review F57). Ask before moving large data.
3. Corgi idea 2 (z-score TF list) on the Corgi base (TODOLIST Next).

## Open questions for Maya
- Scope of review fixes; rescore upstream on test chromosomes (F8) or relabel as validation numbers?
- Lift the "K562 CRISPR is the only benchmark" rule (TODOLIST Waiting)?
- EngreitzLabVault pointer note: register the project (review F70)?

## Pointers
- Findings: F-063 to F-067 at the end of `.living/findings/predicting-regulatory-element-function-at-scale.md`.
- Corgi tables: `A/results/2026_1013_panelmean_{auprc,re2g}_pairs{,.gene}.tsv`.
- Review merged report validated with a scratch copy of mycelium's validator plus the 3 lab categories (upstream validator hard-codes 6).
- Handover is edited here on the Mac; `update-dashboard.sh` pushes it to Oak.

## Read on resume
1. `$OAK/Users/sheth/EP300_BPNet/CLAUDE.md`
2. `todo/TODOLIST.md`: Next (review item, Corgi idea 2) and Waiting
3. `.living/outputs/reviews/2026-10-14-all-reports.md`: Key decisions, Questions, Finding tally
4. `.living/findings/...`: F-067
