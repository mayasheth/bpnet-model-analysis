#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 24:00:00
#SBATCH --mem=120G
#SBATCH -o log/train_atac.%j.txt
#SBATCH -e log/train_atac.%j.txt
#SBATCH --job-name=mm_train_atac
#SBATCH -n 1
#SBATCH --ntasks 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'

set -euo pipefail

FOLD=${1:?Usage: sbatch 1.1.submit_training_atac.sh FOLD_NUM}

SCRIPT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0529_multimodal_p300_model/scripts"
PROJECT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"

# Paths
GENOME="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
SIGNAL_PLUS_BW="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig"
SIGNAL_MINUS_BW="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model/data/ENCSR000EGE_minus.bigWig"
ATAC_BW="$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac.bw"

PEAKS="$PROJECT_DIR/reference/ENCSR000EGE_peaks_inliers.narrowPeak"
NEGATIVES="$PROJECT_DIR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed"
FOLD_JSON="$PROJECT_DIR/reference/hg38_five_folds.json"
OUTPUT_DIR="$SCRIPT_DIR/../models/atac/fold${FOLD}"

mkdir -p "$SCRIPT_DIR/../log" "$OUTPUT_DIR"

module load devel pixi/0.53.0

pixi run -e multimodal python "$PROJECT_DIR/scripts/train_multimodal_bpnet.py" \
    --peaks "$PEAKS" \
    --negatives "$NEGATIVES" \
    --genome "$GENOME" \
    --signal-plus-bw "$SIGNAL_PLUS_BW" \
    --signal-minus-bw "$SIGNAL_MINUS_BW" \
    --accessibility-bw "$ATAC_BW" \
    --fold "$FOLD_JSON" \
    --fold-key "$FOLD" \
    --output-dir "$OUTPUT_DIR" \
    --n-filters 64 \
    --n-acc-filters 8 \
    --n-layers 8 \
    --count-loss-weight 1.0 \
    --batch-size 64 \
    --max-epochs 100 \
    --early-stopping 10 \
    --lr 1e-3 \
    --negative-ratio 0.1 \
    --max-jitter 50 \
    --device cuda

echo "Training complete for fold $FOLD"
