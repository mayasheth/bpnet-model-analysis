# EP300_BPNet: review text pass done; figure pass next

**Goal:** improve the transferred (multi-cell, K562 held out) p300 activity track for ABC/rE2G on K562 CRISPR, and keep the eight reports trustworthy.

## State (verified 2026-10-14)
- Repo `$OAK/Users/sheth/EP300_BPNet`, main at `5c57ecb`, 35 commits UNPUSHED (never push unasked). Working tree clean except long-standing untracked data/large tables (deliberately out of git).
- No jobs running.
- Review 2026-10-14 (`.living/outputs/reviews/2026-10-14-all-reports.md`, 82 findings): Maya chose text fixes first, figure fixes second, holding text that reruns could change.
- **Text pass done and committed (`5c57ecb`).** All eight reports edited, linted (no new errors or warnings vs HEAD; pre-existing slide-structure and untraceable-number warnings remain) and re-rendered (image counts match HEAD). Report 5 now passes the required-section lint. Also: decisions.md logs RC-averaging default-on (F30); F-063/F-064 retitled inconclusive, F-035 wording (F36, F27); ANALYSIS_MANIFEST lists reports 1-7 (F37); DATA_MANIFEST has the scratch-only sets (F55); env docs and CLAUDE.md describe the pixi envs (F62); INDEX.md regenerated (F41, 67 findings).
- Fixed in passing: report 3 RPKM table H3K27me3 over-predicted 5% was 1.17, table says 0.50; report 7 TOBIAS interval lower bound 0.05 -> 0.04; report 6 fig26a neighbour baseline is 0.215, not the review's 0.268.
- Scratch data (F55): Maya said track only, don't copy; decide by 2026-11-14 (TODOLIST Waiting). BAMs (24.5 GB) purge ~2026-11-23, scE2G runs (106 GB) from ~2027-01-01.
- F-067 Corgi is the multi-cell base (`A/models/p300_multicell_panelmean_holdK562_hw500_clw10{,_seed1,_seed2}`).

## Next steps
1. Figure pass (TODOLIST review item): F64-F82 plus the figure parts of F10, F17, F44, F65, F68, F69, F71, F75. F78 (route 5.x/4.56/4.59 through `nature_style`) is the largest. Re-render and recheck captions after.
2. Reruns held for text (TODOLIST review item): F3, F7, F13, F19, F20, F21, rank-loss w10 seeds (F4).
3. Corgi idea 2 (z-score TF list) on the Corgi base (TODOLIST Next).

## Open questions for Maya
- Upstream tables: rescore on `test` chromosomes or relabel as validation numbers (F13)?
- Lift the "K562 CRISPR is the only benchmark" rule (TODOLIST Waiting)?
- Register an EngreitzLabVault pointer note (F63)?

## Pointers
- `build_reports.py` regenerates reports 1-4 from `report_text.py` and would overwrite today's edits; do not run it (TODOLIST code item).
- Edits to Oak files: stage in `~/Documents/ep300-h3k27ac/pending-upload/`, diff, rsync with `-rt` (macOS rsync has no `--chmod`), then empty it.
- Lint only: `python3 render_report.py <r>.qmd --no-render` (module python/3.12.1); render in a small `srun -p engreitz` job.
- Handover is edited here on the Mac; `update-dashboard.sh` pushes it to Oak.

## Read on resume
1. `$OAK/Users/sheth/EP300_BPNet/CLAUDE.md`
2. `todo/TODOLIST.md`: the review item (Next) and Waiting
3. `.living/outputs/reviews/2026-10-14-all-reports.md`: figure findings F64-F82
