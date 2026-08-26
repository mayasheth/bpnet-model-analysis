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
| `residual_evaluation.tsv`, `extreme_residual_elements.tsv` | Residual metric on fragment models. Metric design stands; numbers to be re-derived. |
| `p300_vs_h3k27ac_stratified_evaluation.tsv` | **Mixed** — the p300 side used p300's own 5-prime target (`ENCSR000EGE_{plus,minus}.bigWig`), the H3K27ac side used fragment. Re-derive the H3K27ac side. |
| `gm12878_transfer_stratified_evaluation.tsv` | Cross-cell-type transfer. Underpins F-003; to be re-derived on 5-prime (GM12878 5-prime tracks already built). |
| `gm12878_frag250_replicate_ceiling_by_window.tsv` | GM12878 fragment ceiling. Kept because it *closes* a caveat: 0.8339 top-quintile vs 0.8321 on 5-prime, so the two are equivalent and F-003 normalization holds either way. |

## Current — on the 5-prime target

| File | What it is |
|---|---|
| `fiveprime_replicate_ceiling_by_window.tsv` | K562 ceiling. Top quintile 0.7601 at +/-500. |
| `gm12878_fiveprime_replicate_ceiling_by_window.tsv` | GM12878 ceiling. Top quintile 0.8321 at +/-500 — GM12878 is the *easier* target. |
| `profile_comparison{,_summary}.tsv` | ATAC vs both H3K27ac targets, deliberately side by side. |

## Target-independent

| File | What it is |
|---|---|
| `atac_fragment_length_summary.tsv` | ATAC fragment-length distribution. Concerns the ATAC *input*, not the H3K27ac target. |
| `training_logs/sequence_hw500_clw10.epochlog.tsv` | Preserved epoch log from the mis-weighted first run. |
