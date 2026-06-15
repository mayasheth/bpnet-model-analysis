#!/bin/bash
# Log of commands for the multimodal p300 BPNet project
#
# Three model variants sharing the same architecture (scripts/multimodal_bpnet.py):
#   multimodal  5-channel input (DNA seq + ATAC profile), middle fusion
#   sequence    4-channel input (DNA seq only)  — standard BPNet equivalent
#   atac        1-channel input (ATAC profile only), 64 filters
#
# Training script: scripts/train_multimodal_bpnet.py (--mode flag selects variant)
# Environment: pixi run -e multimodal
#   module load devel pixi/0.53.0

OAK=/oak/stanford/groups/engreitz
PROJECT_DIR=$OAK/Users/sheth/EP300_BPNet
THIS_DIR=$PROJECT_DIR/2026_0529_multimodal_p300_model

# ============================================================
# Stage 0: Preprocess accessibility BigWigs
# ============================================================
# Signal BigWigs (stranded, n_outputs=2) already exist from retrained p300 model:
#   plus:  2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig
#   minus: 2025_0703_retrain_p300_model/data/ENCSR000EGE_minus.bigWig

CHROM_SIZES=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes

# ATAC: tagAlign files (Tn5-shift corrected) -> merged BigWig — DONE
# Files: ENCFF077FBI, ENCFF128WZG, ENCFF534DCE (K562 ATAC-seq, ENCODE)
# module load devel pixi/0.53.0
# pixi run -e multimodal bash $THIS_DIR/scripts/0.1.make_accessibility_bigwig.sh \
#     --input \
#         $OAK/Users/sheth/Data/ENCODE/K562/ENCFF077FBI.tn5.sorted.tagAlign.gz \
#         $OAK/Users/sheth/Data/ENCODE/K562/ENCFF128WZG.tn5.sorted.tagAlign.gz \
#         $OAK/Users/sheth/Data/ENCODE/K562/ENCFF534DCE.tn5.sorted.tagAlign.gz \
#     --output $THIS_DIR/data/atac.bw \
#     --chrom-sizes $CHROM_SIZES \
#     --type atac

# DNase: BAM files -> raw-count BigWig — NOT YET GENERATED
# Files: ENCFF205FNC, ENCFF860XAE (K562 DNase-seq, ENCODE)
# pixi run -e multimodal bash $THIS_DIR/scripts/0.1.make_accessibility_bigwig.sh \
#     --input \
#         $OAK/Users/sheth/Data/ENCODE/K562/ENCFF205FNC.filtered.sorted.bam \
#         $OAK/Users/sheth/Data/ENCODE/K562/ENCFF860XAE.filtered.sorted.bam \
#     --output $THIS_DIR/data/dnase.bw \
#     --chrom-sizes $CHROM_SIZES \
#     --type dnase


# ============================================================
# Stage 1: Train models (all variants, 5 folds each)
# ============================================================
cd $THIS_DIR

## [1.1] Multimodal ATAC (DNA + ATAC profile) — COMPLETE; jobs 27789169-27789173
# for FOLD in 0 1 2 3 4; do sbatch scripts/1.1.submit_training_atac.sh $FOLD; done
# Models: models/atac/fold{0-4}/multimodal_bpnet.torch

## [1.2] Multimodal DNase (DNA + DNase profile) — BLOCKED (data/dnase.bw not generated)
# for FOLD in 0 1 2 3 4; do sbatch scripts/1.2.submit_training_dnase.sh $FOLD; done

## [1.3] ATAC-only (ATAC profile only, no DNA) — SUBMITTED job array 28413076 2026-06-08
# sbatch --array=0-4 scripts/1.3.submit_training_atac_only.sh
# Logs: log/train_atac_only.28413076_{0-4}.txt
# Models: models/atac_only/fold{0-4}/multimodal_bpnet.torch


# ============================================================
# Stage 2: Predict on K562 candidate elements
# ============================================================

## [2.1] Multimodal ATAC predictions — COMPLETE; job 27879460
# sbatch scripts/2.1.submit_predict_multimodal.sh atac
# Output: predictions/atac/  (cv_predictions.tsv.gz, mean_predictions.tsv.gz)
# CV Pearson r: 0.785 (all), 0.663 (p300+)

## [2.2] ATAC-only predictions on K562 elements — COMPLETE 2026-06-08 (job 28442190); CV Pearson = 0.606 (all), 0.435 (p300+)
# sbatch --partition=owners,engreitz,normal --time=8:00:00 --mem=64G \
#   --job-name=atac_only_pred \
#   --output=$THIS_DIR/log/predict_atac_only.%j.txt \
#   --error=$THIS_DIR/log/predict_atac_only.%j.txt \
#   --wrap="module load devel pixi/0.53.0 && pixi run -e multimodal python \
#     $THIS_DIR/scripts/2.1.predict_multimodal.py \
#     --mode atac \
#     --elements $PROJECT_DIR/reference/K562_DNase_candidate_elements.narrowPeak \
#     --signal-plus-bw $PROJECT_DIR/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig \
#     --signal-minus-bw $PROJECT_DIR/2025_0703_retrain_p300_model/data/ENCSR000EGE_minus.bigWig \
#     --accessibility-bw $THIS_DIR/data/atac.bw \
#     --model-dir $THIS_DIR/models/atac_only \
#     --fold-json $PROJECT_DIR/reference/hg38_five_folds.json \
#     --peaks $PROJECT_DIR/reference/ENCSR000EGE_peaks_inliers.narrowPeak \
#     --output-dir $THIS_DIR/predictions/atac_only \
#     --batch-size 512 --device cpu"


# ============================================================
# Stage 3: ATAC vs p300 correlation
# ============================================================

## [3.1] Pearson/Spearman correlation of observed ATAC and p300 logcounts — job 28409190 2026-06-08
# sbatch --partition=owners,engreitz,normal --time=1:00:00 --mem=32G \
#   --job-name=atac_p300_corr \
#   --output=$PROJECT_DIR/log/atac_p300_corr.%j.txt \
#   --error=$PROJECT_DIR/log/atac_p300_corr.%j.txt \
#   --wrap="module load devel pixi/0.53.0 && pixi run -e ism python $PROJECT_DIR/scripts/compute_atac_p300_correlation.py \
#     --elements $PROJECT_DIR/reference/K562_DNase_candidate_elements.narrowPeak \
#     --atac-bw $THIS_DIR/data/atac.bw \
#     --cv-predictions $THIS_DIR/predictions/atac/cv_predictions.tsv.gz \
#     --output $THIS_DIR/predictions/atac/atac_p300_correlation.tsv"


# ============================================================
# Stage 4: SHAP attributions (multimodal ATAC model)
# ============================================================
# Run after predictions complete.
# sbatch scripts/2.1.submit_shap.sh atac fold0
