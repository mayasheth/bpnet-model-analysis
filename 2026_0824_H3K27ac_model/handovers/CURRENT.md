# EP300_BPNet: Bolero idea 2 null (F-069); idea 1 waits on Maya

**Goal:** improve the transferred (multi-cell, K562 held out) p300 activity track for ABC/rE2G on K562 CRISPR. Current line: Bolero preprint ideas (doi:10.64898/2026.10.19.760026), started 2026-10-20 at Maya's request.

## State (verified 2026-11-04)
- Repo `$OAK/Users/sheth/EP300_BPNet`, main at `aaa2542` (F-069 `7732d64`, `abcfloor` launcher stage `0c3d0bb`), tree clean; 46 commits UNPUSHED (never push unasked).
- **Idea 2 done, pre-registered NULL (F-069):** ABC FILM - BASE -0.0056 [-0.0197, +0.0098] gene-level. rE2G +0.0115 [-0.0045, +0.0283]. FILM - CLOSEST +0.0113 [+0.0069, +0.0170] (K562's own vector beats GM12878's). BASE - Corgi5 +0.0139 [-0.0004, +0.0294] (panel expansion alone). Tables `A/results/2026_1021_film13_{auprc,re2g}_pairs{,.gene}.tsv`.
- **Floor context** (each arm vs ATAC only, not pre-registered): the original call failed (4.17 KeyError, it was passed groups it does not load); split into stage `abcfloor`, job 51095622 (queued, Priority) -> `A/results/2026_1021_film13_auprc_floor{,.gene}.tsv`. Add to F-069 when it lands.
- **Idea 3 / descriptive 3:** `6.1` + `6.2` on FILM seed 0 done. GATA1 WGATAA - E-box K562 - HepG2 up to +0.0259 (Corgi +0.0006); max context difference 0.0285 vs Corgi 0.0041; but most pairs are more cooperative in K562 than HepG2, MCF-7 also has the GATA term, FOXA1 not HepG2-favored. Tables `A/results/2026_1022_motif_pair_{spacing,seam_control}_p300_multicell_panelmean_film_..._condcontext*`.
- **Idea 1:** H3K27ac atlas + p300 head, TODOLIST Next; Maya decides the data fetch now that idea 2 is in.
- 89 GB of panel BAMs on `$SCRATCH/p300_panel14/` (purge clock runs from 2026-10-21).

## Next steps
1. When 51095622 finishes: add the floor rows to F-069's table and commit.
2. Report F-069 to Maya; she decides idea 1 (H3K27ac fetch) and whether to run `--cond-cell HepG2` on the K562 context to separate ATAC from the vector.

## Open questions for Maya
- Idea 1 go / no-go after the idea 2 null.
- Delete the 89 GB of BAMs now?
- Report node sh03-15n10 (CUDA driver init failure, 2026-10-22) to SRCC?
- Unanswered: F13 (upstream test chromosomes vs relabel), F63 Vault note, lift "K562 CRISPR only".

## Pointers
- FiLM in the motif assay: `A/scripts/6.1 --model-dir --cond-cell`, `6.2` same; submit env MODEL_DIR / COND_CELL; test `A/scripts/6.1.test_cond.sh`.
- Lessons routed 2026-11-04: `.living/learnings.md` "2026-10-23: bootstrap calls" (4.17 groups); `~/Documents/my-scripts-for-me/learnings/INBOX.md` (afterok on purged jobs, held "user env retrieval failed" tasks). `owners` preemption without resume is already in `learnings/sherlock.md`; use `-p gpu` for long folds.
- Review figure pass (2026-10-14 review) still pending in TODOLIST Next.
- Handover edited here on the Mac; `update-dashboard.sh` pushes it to Oak.

## Read on resume
1. `$OAK/Users/sheth/EP300_BPNet/CLAUDE.md`
2. `.living/findings/predicting-regulatory-element-function-at-scale.md`: F-069
3. `todo/TODOLIST.md`: Now
