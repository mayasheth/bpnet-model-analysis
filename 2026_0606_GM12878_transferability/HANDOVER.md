# Handover: GM12878 cross-cell-type transferability

**Date:** 2026-06-08 (updated)
**Directory:** `/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0606_GM12878_transferability/`

---

## Goal

Test whether K562-trained p300 models generalise to GM12878, and compare against a GM12878-trained model as an in-cell-type ceiling.

| Model | Trained on | Evaluated on | Purpose |
|-------|-----------|-------------|---------|
| GM12878 EP300 BPNet (TFAtlas/ENCODE) | GM12878 | GM12878 | In-cell-type ceiling |
| K562 p300 BPNet v1 (sequence-only) | K562 | GM12878 | Cross-cell-type transfer, seq only |
| K562 multimodal ATAC BPNet | K562 | GM12878 | Cross-cell-type transfer + chromatin |

---

## Results summary

| Model | Pearson — all | Pearson — p300+ | Spearman — all | Spearman — p300+ |
|-------|--------------|-----------------|----------------|------------------|
| GM12878 BPNet (in-cell-type ceiling) | 0.432 | 0.328 | 0.394 | 0.320 |
| K562 v1 BPNet (seq only, cross-cell-type) | 0.277 | 0.114 | 0.267 | 0.114 |
| K562 multimodal ATAC BPNet (cross-cell-type) | **0.793** | **0.628** | **0.819** | **0.636** |
| Inter-replicate ceiling (GM12878) | 0.881 | 0.835 | — | — |

Key finding: the K562-trained multimodal model substantially outperforms the GM12878-trained sequence-only BPNet on GM12878 elements. Chromatin accessibility generalises well across cell types and more than compensates for the cross-cell-type training gap.

---

## Current status (as of 2026-06-08)

### In progress
- [ ] **SHAP on GM12878 BPNet** — job array 28324270 (folds 0–4); 8h time limit
  - Script: `scripts/3.1.submit_gm12878_shap.sh`; logs: `log/shap_fold{0-4}.28324270.txt`
  - Output: `shap_peaks/fold{0-4}/{chr}/counts_scores.h5 + DONE.txt`
  - DONE.txt guard: already-completed chromosomes are skipped on resubmission
  - Previous attempt 28132251 timed out at 4h (6–19/24 chroms done per fold); limit raised to 8h
  - After completion: run merge (step 5.2) then mean (step 5.3) then MoDISCo (step 6.1) — see `0.0.log.sh`

### Completed
- [x] **Transferability plots** — script: `scripts/plot_transferability.py`; outputs in `figures/`:
  - `figures/transferability_bar.pdf` (Fig 2d) — 4 bar groups: GM12878 BPNet ceiling, K562 BPNet v1, K562 multimodal ATAC, GM12878 inter-replicate ceiling; dual bars (all / p300+ elements)
  - `figures/transferability_scatter_all.pdf` (S2a)
  - `figures/transferability_scatter_peaks.pdf` (S2b)

### Completed
- [x] GM12878 candidate elements narrowPeak — `reference/GM12878_candidate_elements.narrowPeak`
  - 154,224 regions (500 bp windows), derived from H3K27ac megamap candidate regions
- [x] GM12878 EP300 BAM files filtered and sorted (`$OAK/Users/sheth/Data/ENCODE/GM12878/EP300/`)
  - Rep 1: `ENCFF515HYM.filtered.sorted.bam` — 14.4M reads; Rep 2: `ENCFF215GSQ.filtered.sorted.bam` — 15.6M reads
- [x] GM12878 EP300 peaks (21,068 peaks) — `ENCFF926AKK.bed.gz` (same directory as BAMs)
- [x] EP300 BigWigs from BAMs — job 27925774; `data/EP300_plus.bw`, `data/EP300_minus.bw`
- [x] ENCODE BPNet signal BigWigs — job 27925868; `data/ENCFF960OFK_plus.bw`, `data/ENCFF941MGK_minus.bw`
  - Same files used to train the GM12878 BPNet model; used as signal in SHAP config
- [x] GM12878 EP300 BPNet model — job 27925831; `GM12878_EP300_BPNet/models/fold_{0-4}/ENCSR000DZG_split000/`
- [x] GM12878 inter-replicate p300 correlation — `GM12878_replicate_correlations.tsv`
  - All: Pearson = 0.881; p300+: Pearson = 0.835; p300−: Pearson = 0.785
- [x] GM12878 ATAC BigWig — job 27969397; `data/atac.bw`
- [x] GM12878 BPNet predictions on GM12878 elements — job 27974981; `predictions/gm12878_bpnet/mean/fold{0-4}/`
- [x] K562 v1 BPNet predictions on GM12878 elements — job 28041042; `predictions/k562_bpnet_v1/mean/fold{0-4}/`
- [x] K562 multimodal ATAC BPNet predictions on GM12878 elements — job 28041421; `predictions/k562_multimodal_atac/`
- [x] Evaluate all 3 models — jobs 28116257/346/354; `predictions/*/mean/all_folds/prediction_accuracy.tsv`

### Next steps (after SHAP completes)
1. Merge per-chr h5 per fold: `merge_shap_across_chrom.py` (step 5.2 in `0.0.log.sh`)
2. Mean across folds: `mean_shap_plus_peaks.py` → `shap_peaks/all_folds/counts_mean_shap_scores.h5`
3. Run TF-MoDISCo: `sbatch scripts/4.1.submit_counts_modisco.sh ...` (step 6.1)
4. Inspect motif hits and compare to K562 MoDISCo results

---

## Key file paths

| File | Description |
|------|-------------|
| `reference/GM12878_candidate_elements.narrowPeak` | 154,224 GM12878 accessible elements |
| `data/ENCFF960OFK_plus.bw` / `ENCFF941MGK_minus.bw` | ENCODE BPNet signal BigWigs (used for SHAP) |
| `data/EP300_plus.bw` / `EP300_minus.bw` | BAM-derived signal (5' end, merged reps) |
| `data/atac.bw` | GM12878 ATAC BigWig (merged 3 reps) |
| `GM12878_EP300_BPNet/models/fold_{0-4}/ENCSR000DZG_split000/` | GM12878-trained BPNet (TF SavedModel) |
| `config/input_data_gm12878_shap.json` | SHAP config (ENCODE BigWigs + K562 controls as bias proxy) |
| `shap_peaks/` | SHAP output directory (in progress) |
| `modisco/` | MoDISCo output directory (pending) |
| `scripts/3.1.submit_gm12878_shap.sh` | Single-fold SHAP script (loops all chroms internally) |
| `$OAK/Users/sheth/Data/ENCODE/GM12878/EP300/ENCFF926AKK.bed.gz` | GM12878 p300 peaks (21,068) |
| `GM12878_replicate_correlations.tsv` | Inter-replicate p300 Pearson/Spearman |

---

## Notes

- SHAP uses `conda activate bpnet_37 + module load cuda/11.1.1 cudnn/8.1.1.33` (TF BPNet requires this; no pixi equivalent)
- Bias BigWigs: K562 controls used as proxy (`ENCSR000EGE_control_{plus,minus}.bigWig`) — GM12878 controls not downloaded
- ENCODE plus strand BigWig (ENCFF960OFK, 150M) is larger than our BAM-derived version (85M) — difference worth investigating but not blocking
- Model path uses `fold_` (underscore) not `fold` — pattern `fold_{fold}/ENCSR000DZG_split000`
- SHAP approach: one SLURM job per fold (array 0–4), looping chromosomes internally — faster than per-chr arrays given only 21k peaks (vs 150k for K562 all-elements SHAP)

---

## Failed attempts (reference)

- **Eval jobs 28079524/525/529** (2026-06-06): conda activation fails in SLURM `--wrap`; invalid `--output-dir` arg; `--peaks` passed to wrong script. Fixed by switching to `pixi run -e ism`, making `--cv-pred-dir` optional, removing invalid args.
- **Per-chr SHAP arrays 28131714–28131722** (2026-06-07): cancelled — per-chromosome parallelization was overkill for 21k peaks. Replaced with single-job-per-fold approach (28132251).
- **SHAP job 28132251** (2026-06-07): timed out at 4h with 6–19/24 chromosomes done per fold. Resubmitted as 28324270 with 8h limit.
- **Transferability figures originally in wrong directory** (2026-06-07): `plot_transferability.py` defaulted to `EP300_BPNet/figures/` (project root) instead of subproject. Fixed `--output-dir` default; moved PDFs to `figures/` here.
