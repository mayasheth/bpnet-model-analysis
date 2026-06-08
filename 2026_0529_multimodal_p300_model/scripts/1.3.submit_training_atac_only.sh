#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 24:00:00
#SBATCH --mem=120G
#SBATCH -o /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0529_multimodal_p300_model/log/train_atac_only.%A_%a.txt
#SBATCH -e /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0529_multimodal_p300_model/log/train_atac_only.%A_%a.txt
#SBATCH --job-name=k562_atac_only
#SBATCH -n 1
#SBATCH --ntasks 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'

# Train ATAC-only BPNet on K562 p300 ChIP-seq peaks.
# Input: 1-channel base-pair resolution ATAC profile (no DNA sequence).
# Usage: sbatch --array=0-4 scripts/1.3.submit_training_atac_only.sh

set -euo pipefail

FOLD=$SLURM_ARRAY_TASK_ID

OAK=/oak/stanford/groups/engreitz
PROJECT_DIR=$OAK/Users/sheth/EP300_BPNet
THIS_DIR=$PROJECT_DIR/2026_0529_multimodal_p300_model

SIGNAL_PLUS_BW=$PROJECT_DIR/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig
SIGNAL_MINUS_BW=$PROJECT_DIR/2025_0703_retrain_p300_model/data/ENCSR000EGE_minus.bigWig
ATAC_BW=$THIS_DIR/data/atac.bw

PEAKS=$PROJECT_DIR/reference/ENCSR000EGE_peaks_inliers.narrowPeak
NEGATIVES=$PROJECT_DIR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed
FOLD_JSON=$PROJECT_DIR/reference/hg38_five_folds.json
OUTPUT_DIR=$THIS_DIR/models/atac_only/fold${FOLD}

mkdir -p $OUTPUT_DIR

echo "$(date): Starting K562 ATAC-only training fold $FOLD"
echo "GPU: $CUDA_VISIBLE_DEVICES"

module load devel pixi/0.53.0

pixi run -e multimodal python $PROJECT_DIR/scripts/train_multimodal_bpnet.py \
    --mode atac \
    --peaks "$PEAKS" \
    --negatives "$NEGATIVES" \
    --accessibility-bw "$ATAC_BW" \
    --signal-plus-bw "$SIGNAL_PLUS_BW" \
    --signal-minus-bw "$SIGNAL_MINUS_BW" \
    --fold "$FOLD_JSON" \
    --fold-key "$FOLD" \
    --output-dir "$OUTPUT_DIR" \
    --n-filters 64 \
    --n-layers 8 \
    --count-loss-weight 1.0 \
    --batch-size 64 \
    --max-epochs 100 \
    --early-stopping 10 \
    --lr 1e-3 \
    --negative-ratio 0.1 \
    --max-jitter 50 \
    --device cuda

echo "$(date): Training complete for fold $FOLD"
