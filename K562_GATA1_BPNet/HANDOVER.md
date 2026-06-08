# Handover: K562 GATA1 BPNet

**Date:** 2026-06-08
**Directory:** `/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/K562_GATA1_BPNet/`

---

## Goal

Train a BPNet model on GATA1 ChIP-seq in K562, compute predictions and SHAP values, and run motif spacing experiments to compare GATA1 motif syntax to p300 syntax.

---

## Current status (as of 2026-06-08)

### In progress
- [ ] **CV fold 2 predictions** — job 28329432 (2h limit, `--threads 1`)
  - Script: `resubmit_cv_fold2.sh`; log: `log/cvp_f2.28329432.txt`
  - Folds 0, 1, 3, 4 complete (good h5 sizes); fold 2 = 96 bytes (failed previously)
- [ ] **Mean predictions folds 0–4** — job array 28329433 (4h limit, `--threads 2`)
  - Script: `resubmit_mean_predict.sh`; logs: `log/mp_f{0-4}.28329433.txt`
  - All 5 folds currently 96 bytes (failed previously)
  - After completion: run `mean_predictions.py` then `2.3.compute_prediction_performance.py` — see `log.sh` section 2D

### Completed
- [x] Models downloaded from TFAtlas (ENCSR000EWM) — `models/fold{0-4}/model_split000/`
- [x] SHAP (all elements) — all 5 folds, all 24 chromosomes — `shap/fold{0-4}/{chr}/`
- [x] CV predictions — folds 0, 1, 3, 4 — `predictions_cv/fold{0,1,3,4}/model_split000_predictions.h5`
- [x] CV TSV/plots — `predictions_cv/all_folds/cv_predictions.tsv.gz`, `cv_predictions.pdf`, `cv_predictions_by_fold.pdf`
  - Pearson r numbers will be in the TSV; full metrics pending mean predictions completion
- [x] Motif spacing experiment — `motif_spacing/GATA_50bp_n234/`
  - Config: GATAA motif, 2/3/4 copies, 50 bp spacing
  - Outputs: `raw_results.tsv`, `results_summary.tsv`, `all_spacings_effects.pdf`, `all_spacings_heatmap.pdf`, `model_variability.pdf`

---

## Key file paths

| File | Description |
|------|-------------|
| `models/fold{0-4}/model_split000/` | TF SavedModel format (from TFAtlas ENCSR000EWM) |
| `data/plus.bigWig`, `data/minus.bigWig` | GATA1 ChIP-seq signal (stranded) |
| `data/peaks_inliers.bed.gz` | GATA1 ChIP-seq peaks |
| `config/input_data_predict.json` | Prediction config (BigWig paths) |
| `predictions_cv/all_folds/cv_predictions.tsv.gz` | CV predictions (folds 0,1,3,4) |
| `predictions_mean/fold{0-4}/model_split000_predictions.h5` | Per-fold mean prediction h5s (pending) |
| `shap/fold{0-4}/` | SHAP scores, all elements, all chroms |
| `motif_spacing/GATA_50bp_n234/` | GATA motif spacing experiment results |
| `resubmit_cv_fold2.sh` | CV fold 2 resubmission script |
| `resubmit_mean_predict.sh` | Mean predictions resubmission script (array 0–4) |
| `log.sh` | Full workflow log with all commands |

---

## After predictions complete

```bash
# 1. Merge per-fold h5s into mean prediction
conda activate bpnet_37
MEAN_PRED_DIR=predictions_mean
PRED_H5=model_split000_predictions.h5
PRED_H5_LIST="$MEAN_PRED_DIR/fold0/$PRED_H5,$MEAN_PRED_DIR/fold1/$PRED_H5,$MEAN_PRED_DIR/fold2/$PRED_H5,$MEAN_PRED_DIR/fold3/$PRED_H5,$MEAN_PRED_DIR/fold4/$PRED_H5"
python $BPNET_DIR/utils/mean_predictions.py \
  --prediction_h5s $PRED_H5_LIST \
  --chrom_sizes $CHR_SIZES \
  --output_dir $MEAN_PRED_DIR/all_folds

# 2. Compute Pearson/Spearman metrics
module load devel pixi/0.53.0
pixi run -e ism python scripts/2.3.compute_prediction_performance.py \
  --mean-pred-dir predictions_mean \
  --cv-pred-dir predictions_cv \
  --peaks reference/K562_DNase_candidate_elements.narrowPeak \
  --overlap-col "GATA1_peak_overlap" \
  --h5-name model_split000_predictions.h5
```

Full commands in `log.sh` section 2D.

---

## Prediction failure history

| Job | Date | Issue | Fix |
|-----|------|-------|-----|
| original jobs | 2026-06-06 | Jobs hung for 8–24h (appeared slow, actually deadlocked) | See below |
| 28137653/28137654 | 2026-06-07 | Root cause 1: no GPU constraint → ran on CPU; still deadlocked after chrM fix | GPU constraint added |
| 28137654 (round 2) | 2026-06-07 | All 5 mean h5s = 96 bytes after chrM exclusion; OAK I/O error or BigWig block issue | Generator patched |
| 28329432/28329433 | 2026-06-08 | Current attempt — generator deadlock fixed | See below |

### Root cause: multiprocessing generator deadlock

`_stealer` in `bpnet-refactor/bpnet/generators/generators.py` calls `mpq.get()` with no timeout. When a worker process crashes (pyBigWig `RuntimeError`), the stealer blocks indefinitely. Jobs appeared to run for hours but were actually frozen.

**Fix applied 2026-06-08** to `bpnet-refactor/bpnet/generators/generators.py`:
- `_proc_target`: wrapped `_generate_batch` in try/except; sends `None` sentinel on any exception and returns
- `_stealer`: checks for `None` sentinel; propagates to regular queue and exits (no more hang)
- `gen()`: raises `RuntimeError` with a clear message if it receives `None` — job fails in seconds with a log entry identifying the problematic batch

**Identified triggers:**
- Mean predictions: `chrM:15952-16568` — 2114 bp window ends at 17317 but GATA1 BigWig chrM size is 16569 bp → out-of-bounds pyBigWig read. Fixed by excluding chrM from `--chroms` in `resubmit_mean_predict.sh`.
- CV fold 2: cause not yet identified (chr4/11/12/15/chrY have no out-of-bounds elements); likely intermittent OAK I/O error. With the generator patched, the job will now fail fast with the exact error instead of hanging.

---

## Notes

- GATA1 model uses `model_split000` (not `ENCSR000EWM_split000`) — h5-name must be `model_split000_predictions.h5` for `2.3.compute_prediction_performance.py`
- Environment: `conda activate bpnet_37` + `module load cuda/11.1.1 cudnn/8.1.1.33` for BPNet predict/SHAP; `pixi run -e ism` for metrics
- CV fold 2 test chroms: chr4, chr11, chr12, chr15, chrY
