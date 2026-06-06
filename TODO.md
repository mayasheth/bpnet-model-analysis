# TODO

## Model improvements

- [ ] Consider switching chromosome fold splits to reduce train/test leakage: evaluate hashFrag (https://github.com/de-Boer-Lab/hashFrag) for sequence-similarity-aware splits
- [ ] Multimodal BPNet (`2026_0529_multimodal_p300_model/`) — see `HANDOVER.md` there for full detail
  - [x] ATAC variant trained (all 5 folds) and predicted; CV Pearson = 0.785 vs v1 = 0.651
  - [ ] Run SHAP on ATAC multimodal model
  - [ ] Train DNase variant (BigWig not yet generated)
  - [ ] Cross-cell-type transferability test (GM12878)

## Evaluation

- [x] Compute inter-replicate p300 ChIP-seq Pearson/Spearman correlation
  - All elements: Pearson = 0.876, Spearman = 0.822
  - p300+ elements: Pearson = 0.746, Spearman = 0.673
  - p300− elements: Pearson = 0.796, Spearman = 0.735
  - Results: `2025_0517_official_EP300_K562_model/replicate_correlations.tsv`
  - Script: `scripts/compute_replicate_correlations.py`

- [x] Fix prediction performance metric — `scripts/2.3.compute_prediction_performance.py` rewritten:
  - Correct metric: `log1p(expm1(+strand) + expm1(−strand))` applied per-strand then summed
  - Reads from h5 files (have coordinates); uses direct interval overlap for peak annotation
  - Bug fixed: mean h5 files have inconsistent element ordering across folds — now sorts by (chrom, start) before averaging

- [x] p300 BPNet v1 performance computed:
  - CV Pearson = **0.651** (all), **0.521** (p300+), 0.555 (p300−)
  - Results: `2025_0517_official_EP300_K562_model/predictions_cv/all_folds/prediction_accuracy.tsv`
  - Mean predictions Pearson pending (job 27917452 running)

- [x] Multimodal ATAC BPNet performance computed:
  - CV Pearson = **0.785** (all), **0.663** (p300+), 0.727 (p300−)
  - Results: `2026_0529_multimodal_p300_model/predictions/atac/prediction_accuracy.tsv`

- [ ] GATA1 BPNet performance — job 27917041 running (fold 2 CV h5 truncated; re-prediction job 27913526 submitted)
  - After fold 2 finishes, re-run: `sbatch` with `2.3.compute_prediction_performance.py` and `--h5-name model_split000_predictions.h5`
  - Results will go to `K562_GATA1_BPNet/predictions_cv/all_folds/prediction_accuracy.tsv`

- [ ] ATAC ChromBPNet performance — predictions exist in `K562_ATAC_ChromBPNet/predictions/` but h5 only has predicted counts (no true counts); using manuscript Pearson r = 0.70 as placeholder in bar chart

- [ ] Compute inter-replicate p300 ChIP-seq correlation for GM12878 BAM files (at `$OAK/Users/sheth/Data/ENCODE/GM12878/EP300`) — use `scripts/compute_replicate_correlations.py`; needed to establish ceiling for GM12878 transferability evaluation

## Publication figures (section 1: BPNet predicts p300 binding)

Scripts: `scripts/plot_model_comparison.py` (bar chart), `scripts/2.3.compute_prediction_performance.py`

- [ ] **Fig 1a** — Training schematic (BioRender/Illustrator): BPNet architecture + p300 peaks vs. GC-matched negatives, 5-fold CV
- [ ] **Fig 1b** — Predicted vs. observed log p300 counts, all 150k accessible elements (CV) — regenerate from `2025_0517.../predictions_cv/all_folds/cv_predictions.tsv.gz` with corrected metric
- [ ] **Fig 1c** — Same as 1b, restricted to p300+ elements — panel from same script output
- [ ] **Fig 1d** — Model comparison bar chart — `figures/model_comparison.pdf`; re-run `plot_model_comparison.py` once GATA1 results complete
- [ ] **Fig S1a** — Per-fold CV scatter (5 panels) — `cv_predictions_by_fold.pdf`; needs formatting
- [ ] **Fig S1b** — Training region comparison: v1 (GC negatives) vs. v2/v3 — exists as separate PDFs; needs combined figure
