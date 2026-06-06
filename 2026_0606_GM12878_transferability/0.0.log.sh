#!/bin/bash
# GM12878 cross-cell-type transferability analysis
# Goal: evaluate whether K562-trained p300 models generalize to GM12878
# Models tested:
#   (1) p300 BPNet v1 — sequence-only (K562-trained)
#   (2) Multimodal ATAC BPNet — DNA + ATAC accessibility (K562-trained)
# Date: 2026-06-06

# --- FILE PATHS --- #

OAK=/oak/stanford/groups/engreitz
PROJECT_DIR=$OAK/Users/sheth/EP300_BPNet
SCRIPTS_DIR=$PROJECT_DIR/scripts
THIS_DIR=$PROJECT_DIR/2026_0606_GM12878_transferability
REF_DIR=$THIS_DIR/reference
DATA_DIR=$THIS_DIR/data
LOG_DIR=$THIS_DIR/log
PRED_DIR=$THIS_DIR/predictions

# GM12878 candidate elements (154,224 regions; derived from H3K27ac megamap)
# Source: /oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_main/ENCODE_rE2G/results/
#         2025_0122_validation_CTCF_H3K27me3/GM12878_H3K27ac_megamap/Peaks/
#         macs2_peaks.narrowPeak.sorted.candidateRegions.bed
GM12878_ELEMENTS=$REF_DIR/GM12878_candidate_elements.narrowPeak  # DONE (2026-06-06)

# GM12878 ENCODE data being downloaded
ATAC_DIR=$OAK/Users/sheth/Data/ENCODE/GM12878/ATAC   # tagAlign files (in progress)
EP300_DIR=$OAK/Users/sheth/Data/ENCODE/GM12878/EP300  # BAM files (in progress)

# Reference genome
CHR_SIZES=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
GENOME=$OAK/Users/sheth/hg38_resources/hg38.fa

# K562-trained models
V1_MODEL_PATTERN="$PROJECT_DIR/2025_0517_official_EP300_K562_model/models/release_run_1/fold{fold}/ENCSR000EGE"
MULTIMODAL_MODEL_DIR=$PROJECT_DIR/2026_0529_multimodal_p300_model

# K562 p300 signal BigWigs (for computing K562 training performance reference)
K562_PLUS=$PROJECT_DIR/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig
K562_MINUS=$PROJECT_DIR/2025_0703_retrain_p300_model/data/ENCSR000EGE_minus.bigWig


# ============================================================
# STAGE 0: PREPARE INPUT DATA
# ============================================================

## [0.1] GM12878 candidate elements narrowPeak — DONE 2026-06-06
# Converted from candidateRegions.bed using 0.2.bed_to_narrowPeak.sh (flank_size=1057)
# bash $SCRIPTS_DIR/0.2.bed_to_narrowPeak.sh \
#   .../macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
#   $GM12878_ELEMENTS \
#   1057


## [0.2] Generate GM12878 ATAC BigWig from tagAlign files
# (once tagAlign files are available at $ATAC_DIR)
ATAC_BW=$DATA_DIR/atac.bw
# sbatch --partition=owners --time=8:00:00 --mem=64G --cpus-per-task=4 \
#   --job-name=gm_atac_bw \
#   --output=$LOG_DIR/atac_bw.%j.txt \
#   --error=$LOG_DIR/atac_bw.%j.txt \
#   --wrap="module load devel pixi/0.53.0 && pixi run -e multimodal bash \
#     $PROJECT_DIR/2026_0529_multimodal_p300_model/scripts/0.1.make_accessibility_bigwig.sh \
#     --input <tagAlign1> <tagAlign2> ... \
#     --output $ATAC_BW \
#     --chrom-sizes $CHR_SIZES \
#     --type atac"


## [0.3] Generate GM12878 EP300 BigWigs (plus/minus strand) from BAM
# (once BAM files are available at $EP300_DIR)
EP300_PLUS=$DATA_DIR/EP300_plus.bigWig
EP300_MINUS=$DATA_DIR/EP300_minus.bigWig
# Use same approach as 2025_0703_retrain_p300_model/scripts/0.3.make_training_bw.sh
# or the make_accessibility_bigwig.sh script with --type chip


# ============================================================
# STAGE 1: PREDICT WITH p300 BPNET V1 (SEQUENCE-ONLY)
# ============================================================
# Model: 2025_0517_official_EP300_K562_model (5-fold, TF SavedModel format)
# Loci: GM12878 candidate elements
# Config: needs input_data_predict.json pointing to GM12878 EP300 signal + candidate elements

V1_PRED_DIR=$PRED_DIR/bpnet_v1
V1_CONFIG=$THIS_DIR/config/input_data_predict.bpnet_v1.json

## [1.1] Create prediction config for v1 BPNet on GM12878
# (copy from K562 config, update signal paths to GM12878 EP300 BigWigs and loci to GM12878 elements)
# Template: $PROJECT_DIR/2025_0517_official_EP300_K562_model/config/input_data_predict.json

## [1.2] Submit mean predictions (all genome, all folds)
# (run from THIS_DIR so log paths resolve correctly)
# sbatch \
#   --array=0-4 \
#   --output=$LOG_DIR/v1_mean_pred_f%a.%A.txt \
#   --error=$LOG_DIR/v1_mean_pred_f%a.%A.txt \
#   $SCRIPTS_DIR/2.2.submit_mean_predict.sh \
#     "$V1_MODEL_PATTERN" \
#     $CHR_SIZES \
#     $GENOME \
#     $V1_PRED_DIR/mean \
#     $V1_CONFIG


# ============================================================
# STAGE 2: PREDICT WITH MULTIMODAL ATAC BPNET
# ============================================================
# Model: 2026_0529_multimodal_p300_model (5-fold, PyTorch)
# Loci: GM12878 candidate elements
# Requires: GM12878 ATAC BigWig ($ATAC_BW)

MULTIMODAL_PRED_DIR=$PRED_DIR/multimodal_atac

## [2.1] Submit multimodal predictions on GM12878
# (once ATAC BigWig is ready)
# sbatch --partition=owners,engreitz --time=6:00:00 --mem=64G --gpus=1 \
#   --job-name=gm_mm_pred \
#   --output=$LOG_DIR/multimodal_pred.%j.txt \
#   --error=$LOG_DIR/multimodal_pred.%j.txt \
#   $MULTIMODAL_MODEL_DIR/scripts/2.1.submit_predict_multimodal.sh atac \
#     --elements $GM12878_ELEMENTS \
#     --atac-bw $ATAC_BW \
#     --out-dir $MULTIMODAL_PRED_DIR


# ============================================================
# STAGE 3: EVALUATE TRANSFERABILITY
# ============================================================
# Compare predicted vs. observed GM12878 p300 log counts (Pearson r)
# Subsets: all elements, p300+ (GM12878 peaks), p300-

## [3.1] Compute prediction performance for v1 BPNet on GM12878
# conda activate tfmodisco
# python $SCRIPTS_DIR/2.3.compute_prediction_performance.py \
#   --mean-pred-dir $V1_PRED_DIR/mean \
#   --peaks <GM12878_EP300_peaks.narrowPeak> \
#   --overlap-col EP300_peak_overlap \
#   --h5-name ENCSR000EGE_split000_predictions.h5

## [3.2] Compute prediction performance for multimodal ATAC BPNet on GM12878
# conda activate tfmodisco
# python $MULTIMODAL_MODEL_DIR/scripts/2.2.plot_prediction_accuracy.py \
#   --pred-dir $MULTIMODAL_PRED_DIR \
#   --peaks <GM12878_EP300_peaks.narrowPeak> \
#   --overlap-col EP300_peak_overlap
