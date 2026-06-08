# TODO

## Model improvements

- [ ] Consider switching chromosome fold splits to reduce train/test leakage: evaluate hashFrag (https://github.com/de-Boer-Lab/hashFrag) for sequence-similarity-aware splits
- [ ] Multimodal BPNet (`2026_0529_multimodal_p300_model/`) — see `HANDOVER.md` there for full detail
  - [x] ATAC variant trained (all 5 folds) and predicted; CV Pearson = 0.785 vs v1 = 0.651
  - [ ] Run SHAP on ATAC multimodal model
  - [ ] Train DNase variant (BigWig not yet generated)
  - [x] Cross-cell-type transferability test (GM12878) — eval complete (jobs 28116257/346/354)

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

- [x] GM12878 cross-cell-type transferability computed:
  - Results: `2026_0606_GM12878_transferability/predictions/*/mean/all_folds/prediction_accuracy.tsv`

  | Model | Pearson — all | Pearson — p300+ |
  |-------|--------------|-----------------|
  | GM12878 BPNet (in-cell-type ceiling) | 0.432 | 0.328 |
  | K562 v1 BPNet (seq only) | 0.277 | 0.114 |
  | K562 multimodal ATAC BPNet | **0.793** | **0.628** |
  | Inter-replicate ceiling | 0.881 | 0.835 |

- [ ] GATA1 BPNet performance — resubmitted 2026-06-08 (jobs 28329432 / 28329433):
  - **Root cause:** multiprocessing generator deadlock in `bpnet-refactor/bpnet/generators/generators.py`
    - `_stealer` called `mpq.get()` with no timeout → hung forever when worker process died
    - **Fixed 2026-06-08:** `_proc_target` now puts `None` sentinel on exception; `_stealer` handles sentinel; `gen()` raises `RuntimeError` — fails fast instead of hanging
    - Mean predictions: chrM excluded (`chrM:15952-16568` 2114bp window exceeds BigWig size)
    - CV fold 2: unknown trigger in chr4/11/12/15/Y; generator fix will now surface exact error
  - CV fold 2: job 28329432; mean predictions folds 0–4: job array 28329433
  - After completion: run `mean_predictions.py` then `2.3.compute_prediction_performance.py --h5-name model_split000_predictions.h5`
  - See `K562_GATA1_BPNet/HANDOVER.md` for full details

- [ ] ATAC ChromBPNet performance — predictions exist in `K562_ATAC_ChromBPNet/predictions/` but h5 only has predicted counts (no true counts); using manuscript Pearson r = 0.70 as placeholder in bar chart

- [ ] Compute prediction performance for p300 BPNet v2 (`2025_0703_retrain_p300_model/`) and v3 (`2025_1016_p300_model_v3/`) using `scripts/2.3.compute_prediction_performance.py` — for supplemental figure comparing training region strategies

- [x] Compute inter-replicate p300 ChIP-seq correlation for GM12878 — results: `2026_0606_GM12878_transferability/GM12878_replicate_correlations.tsv`
  - All: Pearson = 0.881, p300+: 0.835, p300−: 0.785

## Publication figures (section 1: BPNet predicts p300 binding)

Scripts: `scripts/plot_model_comparison.py` (bar chart), `scripts/2.3.compute_prediction_performance.py`

- [ ] **Fig 1a** — Training schematic (BioRender/Illustrator): BPNet architecture + p300 peaks vs. GC-matched negatives, 5-fold CV
- [ ] **Fig 1b** — Predicted vs. observed log p300 counts, all 150k accessible elements (CV) — regenerate from `2025_0517.../predictions_cv/all_folds/cv_predictions.tsv.gz` with corrected metric
- [ ] **Fig 1c** — Same as 1b, restricted to p300+ elements — panel from same script output
- [ ] **Fig 1d** — Model comparison bar chart — `figures/model_comparison.pdf`; re-run `plot_model_comparison.py` once GATA1 results complete
- [ ] **Fig S1a** — Per-fold CV scatter (5 panels) — `cv_predictions_by_fold.pdf`; needs formatting
- [ ] **Fig S1b** — Training region comparison: v1 (GC negatives) vs. v2/v3 — exists as separate PDFs; needs combined figure

## Publication figures (section 2: multimodal model and cross-cell-type transferability)

Scripts: `scripts/plot_transferability.py`, data in `2026_0606_GM12878_transferability/predictions/`

- [ ] **Fig 2a** — Schematic: multimodal architecture (DNA + ATAC → p300), middle fusion design
- [ ] **Fig 2b** — Predicted vs. observed log p300 counts on K562 elements (multimodal CV) — use `2026_0529_multimodal_p300_model/predictions/atac/cv_predictions.tsv.gz`; parallel to Fig 1b
- [ ] **Fig 2c** — Same as 2b, restricted to p300+ elements — parallel to Fig 1c
- [x] **Fig 2d** — Transferability bar chart — `figures/transferability_bar.pdf`; `scripts/plot_transferability.py`
- [x] **Fig S2a** — 3-panel scatter, all GM12878 elements — `figures/transferability_scatter_all.pdf`
- [x] **Fig S2b** — 3-panel scatter, p300+ GM12878 elements — `figures/transferability_scatter_peaks.pdf`

## GM12878 SHAP and MoDISCo

- [ ] **SHAP on GM12878 BPNet (peaks)** — job array 28324270 (folds 0–4, 8h limit); logs: `2026_0606_GM12878_transferability/log/shap_fold{N}.28324270.txt`
  - Previous job 28132251 timed out at 4h; resubmitted with 8h limit; DONE.txt guard skips completed chroms
  - After completion: merge per-chr h5s, average across folds, run MoDISCo — all commands in `2026_0606_GM12878_transferability/0.0.log.sh` stages 5–6
