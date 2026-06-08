# Handover: K562 GATA1 BPNet

**Date:** 2026-06-08 (updated)
**Directory:** `/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/K562_GATA1_BPNet/`

---

## Goal

Train a BPNet model on GATA1 ChIP-seq in K562, compute predictions and SHAP values, and run motif spacing experiments to compare GATA1 motif syntax to p300 syntax.

---

## Current status (as of 2026-06-08)

### In progress
- [ ] **Regenerate `../figures/model_comparison.pdf`** — run `python ../scripts/plot_model_comparison.py` now that GATA1 metrics are available

### Completed
- [x] **Prediction performance metrics** — completed 2026-06-08:
  - `predictions_cv/all_folds/prediction_accuracy.tsv` — CV metrics (all 5 folds)
  - `predictions_mean/all_folds/prediction_accuracy.tsv` — mean predictions metrics
  - **Command:**
    ```bash
    cd /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/K562_GATA1_BPNet
    module load devel pixi/0.53.0
    pixi run -e ism python ../scripts/2.3.compute_prediction_performance.py \
      --cv-pred-dir predictions_cv \
      --mean-pred-dir predictions_mean \
      --peaks data/peaks_inliers.bed.gz \
      --overlap-col GATA1_peak_overlap \
      --h5-name model_split000_predictions.h5
    ```
  - Peaks: `data/peaks_inliers.bed.gz` — ENCODE ENCFF509ZLE (4,603 GATA1 K562 peaks)
  - Results (`predictions_cv/all_folds/prediction_accuracy.tsv`):

  | Subset | Pearson r | Spearman rho |
  |--------|-----------|--------------|
  | CV — all elements | 0.597 | 0.556 |
  | CV — GATA1+ | 0.544 | 0.508 |
  | CV — GATA1− | 0.548 | 0.524 |
  | Mean — all elements | 0.612 | 0.571 |
  | Mean — GATA1+ | 0.620 | 0.587 |
  | Mean — GATA1− | 0.562 | 0.539 |

### Completed
- [x] Models downloaded from TFAtlas (ENCSR000EWM) — `models/fold{0-4}/model_split000/`
- [x] SHAP (all elements) — all 5 folds, all 24 chromosomes — `shap/fold{0-4}/{chr}/`
- [x] CV predictions — all 5 folds — `predictions_cv/fold{0-4}/model_split000_predictions.h5`
  - Fold 2 completed 2026-06-08 (job 28338346); corrupt `chr15:69972829-69973829` BigWig block handled by fill-with-zeros in generator fix
- [x] CV TSV/plots — `predictions_cv/all_folds/cv_predictions.tsv.gz`, `cv_predictions.pdf`, `cv_predictions_by_fold.pdf`
  - TSV regenerated 2026-06-08 to include all 5 folds (150,519 regions)
- [x] Mean predictions — all 5 folds — `predictions_mean/fold{0-4}/model_split000_predictions.h5` (~1.3 GB each)
  - Merged: `predictions_mean/all_folds/mean_predictions.h5` (Jun 8 09:15)
- [x] Motif spacing experiment — `motif_spacing/GATA_50bp_n234/`
  - Config: GATAA motif, 2/3/4 copies, 50 bp spacing
  - Outputs: `raw_results.tsv`, `results_summary.tsv`, `all_spacings_effects.pdf`, `all_spacings_heatmap.pdf`, `model_variability.pdf`

---

## Key file paths

| File | Description |
|------|-------------|
| `models/fold{0-4}/model_split000/` | TF SavedModel format (from TFAtlas ENCSR000EWM) |
| `data/plus.bigWig`, `data/minus.bigWig` | GATA1 ChIP-seq signal (stranded) |
| `data/peaks_inliers.bed.gz` | GATA1 ChIP-seq peaks — **use for `--peaks` in metrics** |
| `config/input_data_predict.json` | Prediction config (BigWig paths) |
| `predictions_cv/all_folds/cv_predictions.tsv.gz` | CV predictions (all 5 folds, 150,519 regions) |
| `predictions_cv/all_folds/prediction_accuracy.tsv` | CV Pearson/Spearman metrics (pending) |
| `predictions_mean/all_folds/mean_predictions.h5` | Merged mean predictions h5 |
| `predictions_mean/all_folds/prediction_accuracy.tsv` | Mean Pearson/Spearman metrics (pending) |
| `shap/fold{0-4}/` | SHAP scores, all elements, all chroms |
| `motif_spacing/GATA_50bp_n234/` | GATA motif spacing experiment results |
| `resubmit_cv_fold2.sh` | CV fold 2 resubmission script (final working version) |
| `resubmit_mean_predict.sh` | Mean predictions resubmission script (final working version) |
| `log.sh` | Full workflow log with all commands |

---

## Prediction failure history

| Job | Date | Issue | Fix |
|-----|------|-------|-----|
| original jobs | 2026-06-06 | Jobs hung for 8–24h (appeared slow, actually deadlocked) | See below |
| 28137653/28137654 | 2026-06-07 | Root cause 1: no GPU constraint → ran on CPU; still deadlocked after chrM fix | GPU constraint added |
| 28137654 (round 2) | 2026-06-07 | All 5 mean h5s = 96 bytes after chrM exclusion; OAK I/O error or BigWig block issue | Generator patched |
| 28329432/28329433 | 2026-06-08 | Generator deadlock fixed but element-level BigWig errors still crashed workers | Fill-with-zeros patch added |
| 28338346/28338347 | 2026-06-08 | **Final run — all predictions complete** | See below |

### Root cause: multiprocessing generator deadlock

`_stealer` in `bpnet-refactor/bpnet/generators/generators.py` calls `mpq.get()` with no timeout. When a worker process crashes (pyBigWig `RuntimeError`), the stealer blocks indefinitely. Jobs appeared to run for hours but were actually frozen.

**Round 1 fix (2026-06-07):** GPU constraint added; chrM excluded from mean predictions.

**Round 2 fix (2026-06-08):** Sentinel-based error propagation:
- `_proc_target`: wrapped `_generate_batch` in try/except; sends `None` sentinel on any exception
- `_stealer`: checks for `None` sentinel; propagates to regular queue and exits
- `gen()`: raises `RuntimeError` with a clear message if it receives `None`
- Also: explicit BigWig handle close after each batch to prevent fd accumulation

**Round 3 fix (2026-06-08):** Element-level error tolerance:
- Signal and bias BigWig reads inside `_generate_batch` wrapped in try/except
- Any pyBigWig error fills the affected window with zeros and logs a WARNING
- This handles corrupt blocks (e.g., `chr15:69972829-69973829`) without crashing the worker

**Identified triggers:**
- Mean predictions: `chrM:15952-16568` — 2114 bp window exceeds BigWig chrM size (16569 bp). Fixed by excluding chrM from `--chroms` in `resubmit_mean_predict.sh`.
- CV fold 2: `chr15:69972829-69973829` — corrupt BigWig block. Fixed by fill-with-zeros in round 3. Final job: 28338346.

---

## Notes

- GATA1 model uses `model_split000` (not `ENCSR000EWM_split000`) — h5-name must be `model_split000_predictions.h5` for `2.3.compute_prediction_performance.py`
- Environment: `conda activate bpnet_37` + `module load cuda/11.1.1 cudnn/8.1.1.33` for BPNet predict/SHAP; `pixi run -e ism` for metrics
- CV fold 2 test chroms: chr4, chr11, chr12, chr15, chrY
