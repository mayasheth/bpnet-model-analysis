# EP300_BPNet: Bolero idea 2 training running; idea 3 seam control done

**Goal:** improve the transferred (multi-cell, K562 held out) p300 activity track for ABC/rE2G on K562 CRISPR. Current line: Bolero preprint ideas (doi:10.64898/2026.10.19.760026), started 2026-10-20 at Maya's request.

## State (verified 2026-10-21)
- Repo `$OAK/Users/sheth/EP300_BPNet`, main at `6be88ab`, 41 commits UNPUSHED (never push unasked). Tree clean apart from the handover sync and long-standing untracked data.
- **Idea 2 pre-registered** (decisions [2026-10-21]): FILM (Corgi + FiLM on `config/cell_cond_chromvar13pc.tsv`) vs BASE (Corgi) on panel 13 (`config/p300_panel_13.json`: 12 training cells, K562 held out; SK-N-SH_RA in, lung out per Maya). Primary FILM - BASE, ABC, seed-averaged, gene-level bootstrap, PASS if interval > 0. Descriptive: CLOSEST (FILM with GM12878's row), BASE - Corgi5, `6.1` on FILM.
- **Running:** training FILM 48843057/58/59, BASE 48843060/61/62 (`A/scripts/1.43`, 192 GB; smoke 48824689 peaked 43 GB at 2,000 negatives/cell). Scoring launcher 48844639 (`A/log/film13_score_launch.sh`, afterok on all six, kill-on-invalid-dep): 4.26 -> 9 x 4.4 -> 4.5 ABC (region gate) -> 4.6/4.7 CRISPR -> 4.17; rE2G 4.32 -> 4.33 -> 4.34. Child IDs in `A/log/film13_score_launch.48844639.txt`.
- **Outputs to read:** `A/results/2026_1021_film13_auprc_pairs{,.gene}.tsv` (FILM-BASE primary, FILM-CLOSEST, BASE-Corgi5), `..._auprc_floor{,.gene}.tsv`, `..._re2g_pairs{,.gene}.tsv`.
- **As of 18:51:** FILM seed-0 folds 0-1 (48843057_0/_1) running, the other 28 tasks pending (fairshare 0.087).
- **Idea 3:** F-068 restated with the seam control (`A/scripts/6.2`, decisions [2026-10-21] "Seam control for F-068" + amendment): 21/31 seam rows syntax, 9 junction k-mers (AP1-CTCF, AP1-CEBPB d = 0, ...), 1 not partner-specific, same in K562 and HepG2; CEBPB replicate QC passes. The three leading pairs (AP1-IRF4 +/-2, AP1-ETS1 0) pass, but where a junction base changed they keep only 0.33-0.47 of the spike (descriptive). Jobs 48896361 smoke, 48899901 array, 48904583 combine -> `A/results/2026_1021_motif_pair_seam_control_corgi_k5_cs1.tsv`. 6.1 refactored (`load_models`, `build_acc`), gate shows identical output.
- **Idea 1:** H3K27ac atlas + p300 head, in TODOLIST Next; data fetch waits on idea 2's readout (Maya).
- 89 GB of panel BAMs on `$SCRATCH/p300_panel14/`; tracks verified; Maya to delete or let purge.

## Next steps
1. Check jobs (`squeue -u $USER`); if a training fold fails, the scoring launcher is cancelled: fix, resubmit that array, resubmit `log/film13_score_launch.sh` with the new afterok IDs.
2. When the bootstraps land: read the primary against the pre-registration (gene-level table decides), write F-069 with both resampling units and the descriptive contrasts, update TODOLIST and this handover.
3. Rerun `A/scripts/6.1` and `6.2` on the FILM base run (needs FiLM conditioning passed into 6.1's forward) (descriptive 3): GATA1 WGATAA - E-box, K562 minus HepG2 context, vs F-068.
4. Report to Maya; she decides the H3K27ac fetch for idea 1.

## Open questions for Maya
- Delete the 89 GB of BAMs now?
- Uncommitted on Oak: 6.1 refactor, 6.2 + submit, decisions, F-068, TODOLIST (commit when Maya asks).
- Unanswered: F13 (upstream test chromosomes vs relabel), F63 Vault note, lift "K562 CRISPR only".

## Pointers
- Code: FiLM `scripts/multimodal_bpnet.py attach_cell_conditioning`, `train_multimodal_bpnet.py --cell-cond-tsv`; `A/scripts/4.1|4.4 --cond-cell/COND_CELL`; `4.63 --panel-json`, `4.64`; `1.33 --exclude`; arm family `film13` in 4.3/4.6/4.32.
- Lessons: `.living/learnings.md` 2026-10-20. Pre-existing bug logged in TODOLIST: `4.6 --arms dnaseinput|multitask` assert.
- Review figure pass (2026-10-14 review) still pending in TODOLIST Next.
- Handover edited here on the Mac; `update-dashboard.sh` pushes it to Oak.

## Read on resume
1. `$OAK/Users/sheth/EP300_BPNet/CLAUDE.md`
2. `.living/decisions.md`: "[2026-10-21] Pre-register Bolero idea 2" and "Seam control for F-068"
3. `todo/TODOLIST.md`: Now (Bolero ideas 2, 3)
4. `A/log/film13_score_launch.48844639.txt` and the result TSVs above, once present
