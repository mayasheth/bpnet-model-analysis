# Handover: multimodal p300 BPNet

**Date:** 2026-06-08 (updated)
**Directory:** `/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0529_multimodal_p300_model/`

---

## What this project is

A multimodal BPNet variant that takes **DNA sequence + base-pair-resolution chromatin accessibility** as input and predicts stranded p300 ChIP-seq signal. The goal is to test whether providing accessibility as an explicit input improves p300 prediction beyond sequence alone, and to understand what chromatin context the model uses.

Architecture is defined in `scripts/multimodal_bpnet.py`. See `multimodal_bpnet_architecture.html` for a full interactive description of every parameter and design decision.

---

## Current status (as of 2026-06-07)

### Completed
- [x] Model architecture (`scripts/multimodal_bpnet.py`) — 5-channel input, middle fusion, n_outputs=2 (stranded), PyTorch
- [x] Training script (`scripts/train_multimodal_bpnet.py`) — stranded signal inputs, log1p z-score accessibility normalization, RC augmentation, 5-fold CV
- [x] ATAC BigWig generated: `data/atac.bw` (1.9G, merged from 3 tagAlign replicates)
- [x] Training scripts updated with real paths for ATAC variant (`scripts/1.1.submit_training_atac.sh`)
- [x] pixi environment (`multimodal`) working — bedtools/samtools/bedGraphToBigWig, leidenalg/igraph via conda to avoid source compilation issues
- [x] **ATAC model trained — all 5 folds complete.** Models saved to `models/atac/fold{0-4}/multimodal_bpnet.torch`. Training summary: `models/atac/training_summary.tsv`.

  | Fold | Best epoch | Val MNLL | Val profile Pearson | Val count Pearson | Val count MSE |
  |------|-----------|----------|---------------------|-------------------|---------------|
  | 0 | 19 | 494.9 | 0.268 | 0.665 | 0.232 |
  | 1 | 28 | 512.1 | 0.265 | 0.662 | 0.226 |
  | 2 | 32 | 540.1 | 0.255 | 0.703 | 0.215 |
  | 3 | 23 | 561.3 | 0.259 | 0.686 | 0.214 |
  | 4 | 23 | 520.2 | 0.272 | 0.690 | 0.236 |
  | mean | 25 | 525.7 | 0.264 | 0.681 | 0.225 |

- [x] Prediction scripts written: `scripts/2.1.predict_multimodal.py`, `scripts/2.2.plot_prediction_accuracy.py`, `scripts/2.1.submit_predict_multimodal.sh`
- [x] **ATAC predictions complete** (job 27879460). Outputs in `predictions/atac/`:
  - `mean_predictions.tsv.gz`, `cv_predictions.tsv.gz`, `mean_predictions.pdf`, `cv_predictions.pdf`, `prediction_accuracy.tsv`
  - CV Pearson r: **0.785** (all elements), **0.663** (p300+ only), 0.727 (p300−)
  - Mean Pearson r: 0.798 (all), 0.699 (p300+), 0.739 (p300−)
  - Compares favourably to v1 BPNet CV Pearson = 0.651 (all), 0.521 (p300+)

### Immediate next steps

- [x] **Submit ATAC training** (5 folds):
  ```bash
  cd /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0529_multimodal_p300_model
  for FOLD in 0 1 2 3 4; do sbatch scripts/1.1.submit_training_atac.sh $FOLD; done
  ```
  - Jobs 27789169–27789173 (folds 0–4); logs written to `log/train_atac.27789169.txt` and similar
  - Previous attempts failed and were fixed:
    - **27619407–27619411**: SLURM copies sbatch scripts to `/var/spool/slurmd/`, causing `${BASH_SOURCE[0]}` to resolve incorrectly. Fixed by hardcoding `SCRIPT_DIR`/`PROJECT_DIR` as absolute paths in `1.1.submit_training_atac.sh`.
    - **27621027–27621033**: Two bugs: (1) fold JSON was per-fold file (`fold{N}.json`, key `'0'` only) instead of `hg38_five_folds.json` — fixed in submit script. (2) OOM: all ~2.5M genome-wide negatives were loaded into RAM before subsampling. Fixed in `scripts/train_multimodal_bpnet.py` to subsample to `10 × n_train_peaks` (~220k) before extracting windows. Unused `val_negs` load also removed.
    - **27655596–27655600**: `AssertionError: Torch not compiled with CUDA enabled`. Conda-forge strict channel priority excluded `nvidia` channel packages, so pixi installed CPU-only pytorch. Fixed by switching `torch` and `captum` from conda to PyPI in `pixi.toml`, using the PyTorch cu118 wheel index (`https://download.pytorch.org/whl/cu118`). Verified `torch 2.6.0+cu118` installs correctly.
    - **27707218–27707222**: `RuntimeError: Input type (torch.cuda.FloatTensor) and weight type (torch.FloatTensor) should be the same`. `fit()` moved inputs to CUDA but never called `self.to(device)` on the model. Fixed by adding `self.to(device)` at the start of `fit()` in `scripts/multimodal_bpnet.py`.
    - **27758489–27758493**: `TypeError: 'NoneType' object is not subscriptable` — `peak_ordering` initialized to `None` and only set when `idx==0` in `__getitem__`. DataLoader workers receive random idx values and may never see `idx=0` first. Fixed by initializing `peak_ordering = self.rng.permutation(self.n_peaks)` in `__init__` in `scripts/train_multimodal_bpnet.py`.
    - **27787994–27787999**: Full code review before next submit revealed two more bugs: (1) `labels` tensor not moved to device in `fit()` — `y_hat_logits[labels == 1]` in `_mixture_loss` would crash indexing a GPU tensor with a CPU mask. (2) `num_workers=4` with a stateful dataset breaks epoch-level peak shuffling — each worker has an isolated copy of `n_peaks_seen`/`rng` state and workers other than worker-0 never trigger the `idx==0` reshuffle. Fixed: added `labels = labels.to(device)` in `fit()`; set `num_workers=0` in DataLoader (all data is pre-loaded in RAM, no I/O to parallelize). Also removed dead `y_valid_counts` line and extracted `device_type = device.split(':')[0]` for `torch.autocast`.

- [x] **Submit ATAC predictions on 150k candidate elements** (job 27879460):
  ```bash
  sbatch scripts/2.1.submit_predict_multimodal.sh atac
  ```
  - Outputs to `predictions/atac/`: `mean_predictions.tsv.gz`, `cv_predictions.tsv.gz`, `mean_predictions.pdf`, `cv_predictions.pdf`, `prediction_accuracy.tsv`
  - Prediction script bugs fixed during submission:
    - `RuntimeError: Invalid interval bounds!` — ATAC BigWig has different chrom sizes than FASTA. Fixed by computing min chrom size across all BigWigs + FASTA before querying.
    - `ModuleNotFoundError: No module named 'multimodal_bpnet'` — `torch.load` requires the model class to be importable; fixed by inserting `EP300_BPNet/scripts/` (not this subproject's scripts/) into sys.path.

- [x] **Check prediction results and compare to standard BPNet counts Pearson**
  - Multimodal ATAC CV Pearson = 0.785 vs v1 BPNet CV Pearson = 0.651 — clear improvement
  - `2.2.plot_prediction_accuracy.py` updated: `--max-counts` arg (default 10 for zoom), CV now shows all/p300+/p300- panels; `--from-tsv` flag for fast replotting without loading h5 files
  - `scripts/2.3.compute_prediction_performance.py` (shared script) also updated with `--from-tsv` and `--max-counts` args
  - p300 BigWig files (`2025_0703_retrain_p300_model/data/ENCSR000EGE_{plus,minus}.bigWig`) verified correct: these are the Vivek/TFAtlas files; peak max ~2–4 per base at 1bp resolution is expected for this ChIP-seq depth

- [ ] **Run SHAP on ATAC multimodal model** — use `scripts/shap_multimodal_bpnet.py`; submit similar to v1 BPNet SHAP jobs

- [ ] **Generate DNase BigWig** when ready to run DNase variant:
  ```bash
  sbatch --partition=owners --time=8:00:00 --mem=64G --cpus-per-task=4 \
    --job-name=dnase_bw \
    --output=log/dnase_bw.%j.txt \
    --error=log/dnase_bw.%j.txt \
    --wrap="module load devel pixi/0.53.0 && pixi run -e multimodal bash \
      /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0529_multimodal_p300_model/scripts/0.1.make_accessibility_bigwig.sh \
      --input \
        /oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562/ENCFF205FNC.filtered.sorted.bam \
        /oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562/ENCFF860XAE.filtered.sorted.bam \
      --output /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0529_multimodal_p300_model/data/dnase.bw \
      --chrom-sizes /oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes \
      --type dnase"
  ```
- [x] **Cross-cell-type transferability (GM12878)** — complete; see `../2026_0606_GM12878_transferability/HANDOVER.md`
  - K562 multimodal on GM12878: Pearson = **0.793** (all), **0.628** (p300+)
  - Outperforms GM12878-trained seq-only BPNet (0.432 / 0.328) — chromatin accessibility generalises cross-cell-type
  - [x] Plots written — `../2026_0606_GM12878_transferability/figures/transferability_{bar,scatter_all,scatter_peaks}.pdf`
    - Bar chart (Fig 2d): 4 bar groups with dual all/p300+ bars; inter-replicate ceiling added as 4th group

- [x] **GM12878 multimodal BPNet training** — all 5 folds complete 2026-06-08 (job 28359131; ~6 min/fold on GPU)
  - Same architecture as K562 ATAC model (64 filters, 8 layers, middle fusion)
  - Models: `../2026_0606_GM12878_transferability/GM12878_multimodal_BPNet/models/atac/fold{0-4}/`
  - Purpose: in-cell-type multimodal ceiling to compare against K562 → GM12878 transfer
- [ ] **GM12878 multimodal predictions** — job 28387044 submitted 2026-06-08
  - Script: `../2026_0606_GM12878_transferability/scripts/2.1.submit_predict_gm12878_multimodal.sh`
  - Output: `../2026_0606_GM12878_transferability/predictions/gm12878_multimodal_atac/`
  - After completion: evaluate with `2.3.compute_prediction_performance.py`, add to Fig 2d

- [ ] Submit DNase training (5 folds) after `data/dnase.bw` exists:
  ```bash
  for FOLD in 0 1 2 3 4; do sbatch scripts/1.2.submit_training_dnase.sh $FOLD; done
  ```

---

## Key file paths

| File | Description |
|------|-------------|
| `data/atac.bw` | ATAC-seq accessibility BigWig (1 bp, merged 3 replicates) |
| `data/dnase.bw` | DNase-seq accessibility BigWig — **not yet generated** |
| `scripts/1.1.submit_training_atac.sh` | SLURM training script (ATAC variant, takes fold number as arg) |
| `scripts/1.2.submit_training_dnase.sh` | SLURM training script (DNase variant) |
| `scripts/0.1.make_accessibility_bigwig.sh` | BigWig generation for ATAC (tagAlign) or DNase (BAM) |
| `scripts/0.0.log.sh` | Full workflow log with commands |
| `multimodal_bpnet_architecture.html` | Detailed architecture documentation |
| `config/input_data_atac.json` | All input file paths for ATAC variant |
| `config/input_data_dnase.json` | All input file paths for DNase variant |

### Signal BigWigs (p300 ChIP-seq, stranded)
- Plus: `/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig`
- Minus: `/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model/data/ENCSR000EGE_minus.bigWig`

### Input data
| Modality | Files | Location |
|----------|-------|----------|
| ATAC tagAlign (3 reps) | `ENCFF077FBI/ENCFF128WZG/ENCFF534DCE.tn5.sorted.tagAlign.gz` | `/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562/` |
| DNase BAM (2 reps) | `ENCFF205FNC/ENCFF860XAE.filtered.sorted.bam` | `/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE/K562/` |
| p300 peaks | `ENCSR000EGE_peaks_inliers.narrowPeak` | `reference/` |
| GC-matched negatives | `genomewide_gc_stride_1000_flank_size_1057.gc.bed` | `reference/` |
| Genome | `hg38.fa` | `/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/` |
| Chrom sizes | `GRCh38.main.chrom.sizes` | `/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/` |
| Fold JSONs | `fold{0-4}.json` | `reference/` |

---

## Architecture summary

- **Input:** (N, 5, 2114) — 4-channel one-hot DNA + 1-channel log1p z-scored accessibility
- **Fusion:** Middle fusion — separate initial Conv1d for each modality (kernel=21), then concatenate → (N, 72, 2114)
- **Core:** 8 dilated residual Conv1d layers (dilation 2→256), receptive field 1067 bp
- **Profile head:** Conv1d(72→2, kernel=75) → trim to (N, 2, 1000) — plus and minus strand logits
- **Counts head:** Global mean pool → Linear(72→1) → log counts
- **Loss:** MNLL (profile) + log1pMSE (counts), weight 1.0
- **Parameters:** ~141k

---

## Environment

```bash
module load devel pixi/0.53.0
pixi run -e multimodal python scripts/train_multimodal_bpnet.py ...
```

The `multimodal` pixi environment is defined in `pixi.toml` at the project root. Key packages: pytorch, pybigwig, pyfaidx, bedtools, samtools, ucsc-bedgraphtobigwig, leidenalg (conda, not PyPI — avoids igraph C compilation), pysam (conda).

**Known pixi issue:** deeptools cannot be installed in this env due to channel priority conflicts with the old glibc/kernel requirements. The BigWig generation scripts use bedtools + bedGraphToBigWig instead (no RPGC normalization, but training normalizes accessibility anyway).

---

## Design decisions made this session

1. **n_outputs=2 (stranded)** — switched from n_outputs=1 to match standard BPNet convention and preserve strand-asymmetric read pileup in the profile target.
2. **Middle fusion** — separate initial convolutions per modality, merged before dilated stack.
3. **No deeptools** — replaced bamCoverage with bedtools genomecov to avoid env conflicts; RPGC normalization dropped since training z-scores accessibility anyway.
4. **sort temp to /tmp** — sort in BigWig generation uses `-T /tmp` (node-local disk) rather than OAK to avoid network filesystem slowness.
5. **Accessibility normalization** — log1p then z-score using training-peak statistics saved to `acc_normalization_stats.json` per fold.
