# TODO

## Model improvements

- [ ] Consider switching chromosome fold splits to reduce train/test leakage: evaluate hashFrag (https://github.com/de-Boer-Lab/hashFrag) for sequence-similarity-aware splits
- [ ] Train a p300 BPNet-like model that takes DNA sequence AND base-pair-level chromatin accessibility as input (multimodal architecture)
  - Architecture implemented: `scripts/multimodal_bpnet.py` (MultiModalBPNet), `scripts/train_multimodal_bpnet.py`, `scripts/shap_multimodal_bpnet.py`
  - Model dir: `2026_0529_multimodal_p300_model/`; pixi env: `multimodal` (see `pixi.toml`)
  - ATAC variant trained (all 5 folds); predictions running (job 27879460)
  - Future improvement: bias-correct accessibility input signal; train DNase variant

## Evaluation

- [x] Compute inter-replicate p300 ChIP-seq Pearson/Spearman correlation — done
  - Results: `$OAK/Users/sheth/TF_analysis/2026_0603_EP300_reps/EnhancerList.rpm_correlations.tsv`
  - R1 vs R2: Pearson = 0.821, Spearman = 0.822 (across all ~150k accessible elements)
  - [ ] Also compute separately for p300-peak-overlapping vs. non-overlapping elements (data already in `EnhancerList.chromatin_annotations.tsv`, just filter and recalculate)

- [ ] Fix prediction performance metric for consistency across models
  - Current `2.3.compute_prediction_performance.py` (CV mode) uses `log1p(+strand) + log1p(-strand)` — incorrect; should use `log1p(+counts + -counts)`
  - Multimodal model (`2.1.predict_multimodal.py`) correctly uses `log1p(total)` — use this as the standard
  - Need to re-run or correct performance numbers for primary BPNet (v1) before comparing across models
  - CV predictions for primary BPNet lack coordinates/peak_overlap — need to regenerate in the clean format from the multimodal script

- [ ] Compute performance on the same 150k elements for comparison across models (all using `log1p(total counts)` metric, split by p300 peak overlap):
  - [ ] Primary p300 BPNet v1 (sequence only) — re-run with corrected metric
  - [ ] ATAC ChromBPNet — predictions may already exist in `K562_ATAC_ChromBPNet/predictions/`
  - [ ] GATA1 BPNet — predictions may already exist in `K562_GATA1_BPNet/predictions_mean/`

## Publication figures (section 1: BPNet predicts p300 binding)

- [ ] **Fig 1a** — Training schematic (BioRender/Illustrator): BPNet architecture + p300 peaks vs. GC-matched negatives, 5-fold CV
- [ ] **Fig 1b** — Predicted vs. observed log p300 counts, all 150k accessible elements (CV, all folds combined) — exists as `2025_0517.../predictions_cv/all_folds/pred_vs_observed.all_folds.pdf`; needs formatting + corrected metric
- [ ] **Fig 1c** — Same as 1b, restricted to p300-peak-overlapping elements — currently panel 2 of mean predictions PDF; needs to be regenerated for CV predictions
- [ ] **Fig 1d** — Model comparison bar/dot chart: p300 BPNet v1, ATAC ChromBPNet, GATA1 BPNet, inter-replicate ceiling (all on same elements, same metric)
- [ ] **Fig S1a** — Per-fold CV scatter (5 panels) — exists as `pred_vs_observed_by_fold.pdf`; needs formatting
- [ ] **Fig S1b** — Training region comparison: v1 (GC negatives) vs. v2/v3 — exists as separate PDFs; needs combined figure
