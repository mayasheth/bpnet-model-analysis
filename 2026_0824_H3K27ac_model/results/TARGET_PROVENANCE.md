# Which target definition each result used

**Current target: 5-prime ends.** `data/h3k27ac_5p_{plus,minus}.bw`, built by
`scripts/0.5.make_5prime_bigwigs.sh` with `bedtools genomecov -5 -dz`, matching the p300
convention in `scripts/0.3.make_training_bw.sh`.

**Do not train or evaluate anything new on the 250 bp fragment-extended track**
(`$OAK/Users/sheth/Data/share/IGV/ENCSR000AKP_coverage.bw`). It smears each 36 bp read
across 250 bp, which breaks the multinomial read-count assumption behind bpnetlite's
profile loss and bleeds signal between neighbouring elements. It is retained only so
earlier results remain reproducible. `scripts/train_multimodal_bpnet.py` now writes
`training_target.json` into every model directory and refuses a `--count-offset-model`
whose target does not match.

## Superseded — on the fragment-extended target, do not quote

| File | What it is |
|---|---|
| `replicate_ceiling_by_window.tsv` | K562 ceiling. Superseded by `fiveprime_replicate_ceiling_by_window.tsv`. Top-quintile figures were near-identical (0.7605 vs 0.7601 at +/-500). |
| `signal_vs_distance.tsv`, `window_tradeoff.tsv` | Meta-profile and window trade-off. Conclusion (+/-500) still holds: `profile_comparison.tsv` shows the 5-prime and fragment profile shapes overlay almost exactly. |
| `flanking_vs_full_window.tsv` | Flanking-only test. Conclusion (rejected) unaffected — it was a within-target comparison. |
| `stratified_evaluation.tsv`, `training_results_{per_fold,summary}.tsv` | The 3-mode x 2-window grid. **Being replaced** by the 5-prime grid. |
| `residual_evaluation.tsv`, `extreme_residual_elements.tsv` | Residual metric on fragment models. Metric design stands. **Superseded 2026-08-29** by the `residual5p_*` files below. |
| `p300_vs_h3k27ac_stratified_evaluation.tsv` | **Mixed** — the p300 side used p300's own 5-prime target (`ENCSR000EGE_{plus,minus}.bigWig`), the H3K27ac side used fragment. Re-derive the H3K27ac side. |
| `gm12878_transfer_stratified_evaluation.tsv` | Cross-cell-type transfer. Underpins F-003; to be re-derived on 5-prime (GM12878 5-prime tracks already built). |
| `gm12878_frag250_replicate_ceiling_by_window.tsv` | GM12878 fragment ceiling. Kept because it *closes* a caveat: 0.8339 top-quintile vs 0.8321 on 5-prime, so the two are equivalent and F-003 normalization holds either way. |

## Current — on the 5-prime target

| File | What it is |
|---|---|
| `fiveprime_replicate_ceiling_by_window.tsv` | K562 ceiling. Top quintile 0.7601 at +/-500. |
| `gm12878_fiveprime_replicate_ceiling_by_window.tsv` | GM12878 ceiling. Top quintile 0.8321 at +/-500 — GM12878 is the *easier* target. |
| `profile_comparison{,_summary}.tsv` | ATAC vs both H3K27ac targets, deliberately side by side. |
| `residual5p_residual_evaluation.tsv`, `residual5p_extreme_residual_elements.tsv` | Residual metric, pooled, on the 5-prime target. Baseline `atac5p_hw500_clw10`, compared against `sequence5p`, `multimodal5p` and the offset-trained `residual5pFIXED_hw500_clw10`. |
| `residual5p_per_fold.tsv`, `residual5p_fold_summary.tsv` | The same comparison per fold with t-based 95% CIs, plus the two artifact controls (`r_out_vs_atac`, `partial_r`) and top-signal-quintile columns. **Quote these, not the pooled table.** |

## Target-independent

| File | What it is |
|---|---|
| `atac_fragment_length_summary.tsv` | ATAC fragment-length distribution. Concerns the ATAC *input*, not the H3K27ac target. |
| `training_logs/sequence_hw500_clw10.epochlog.tsv` | Preserved epoch log from the mis-weighted first run. |


## A note on offset-trained models

`models/residual5pFIXED_hw500_clw10` was trained with `--count-offset-model
models/atac5p_hw500_clw10`. `MultiModalBPNet.forward()` does **not** add that offset back —
it is applied only inside `fit()` when scoring the loss — so the model's raw output is the
RESIDUAL `observed - atac_pred`, not logcounts. Any evaluator must add `atac_pred` back to
recover a comparable prediction. `2.4.evaluate_residual.py` handles this behind a
`"residual": true` key in the config JSON; without that key it would subtract `atac_pred`
twice and report a plausible wrong number rather than failing.
`2.7.diagnose_residual_offset.py` checks the semantics empirically (a residual predictor's
output centres near 0, a logcounts predictor near the mean observed logcount) with the plain
sequence model as a control, and should be re-run if the trainer's offset handling changes.

## Added since the 5-prime switch

All on the 5-prime H3K27ac target and the ORIGINAL full-interval `atac.bw` input.

| File | What |
|---|---|
| `residual_grid_{residual_evaluation,per_fold,fold_summary}.tsv` | K562 residual grid: 3 standard + 3 residual-objective models against the ATAC-only baseline. Supersedes `residual5p_*`. |
| `gm12878_residual_grid_*.tsv` | The same grid in GM12878. |
| `{transfer,deploy}_{k562_to_gm,gm_to_k562}_*.tsv` | Cross-cell-type transfer. `transfer_*` uses the TARGET cell type's own ATAC model as the offset (optimistic — needs target H3K27ac); `deploy_*` transfers the source ATAC model (the real application). |
| `accs5p_{k562,gm12878,p300}_stratified_*.tsv` | Paired comparison of the ChromBPNet 5-prime accessibility INPUT against full-interval coverage. The only results here whose accessibility input is `atac_5p.bw`. |
| `atac_vs_h3k27ac_by_celltype.tsv` | Model-free coupling, now 6 cell types/conditions. |
| `coupling_panel_recomputed.tsv` | Coupling recomputed for Fig 8; agrees with the above to 0.002. |

**Two accessibility inputs now exist and must never be mixed.** `atac.bw` is
`genomecov -bg` over the full tagAlign interval, so coverage scales with read length;
`atac_5p.bw` is `genomecov -bg -5`, single-base insertion counts, read-length independent
(the ChromBPNet convention). Every result except `accs5p_*` uses `atac.bw`.

**Number traceability.** `outputs/numbers.json` is generated from these files by
`scripts/3.7.build_numbers_manifest.py`; `render_report.py` checks every quoted number
against it. Re-run the generator after any pipeline change, then re-render — prose that no
longer matches its source is flagged.
