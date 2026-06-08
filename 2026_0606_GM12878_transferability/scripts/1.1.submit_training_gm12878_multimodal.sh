#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 24:00:00
#SBATCH --mem=120G
#SBATCH -o /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0606_GM12878_transferability/log/train_gm12878_multimodal.%j.txt
#SBATCH -e /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0606_GM12878_transferability/log/train_gm12878_multimodal.%j.txt
#SBATCH --job-name=gm12878_mm_f%a
#SBATCH -n 1
#SBATCH --ntasks 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'

# Usage: sbatch --array=0-4 scripts/1.1.submit_training_gm12878_multimodal.sh

set -euo pipefail

FOLD=$SLURM_ARRAY_TASK_ID

OAK=/oak/stanford/groups/engreitz
PROJECT_DIR=$OAK/Users/sheth/EP300_BPNet
THIS_DIR=$PROJECT_DIR/2026_0606_GM12878_transferability

GENOME=$OAK/Users/sheth/hg38_resources/hg38.fa
SIGNAL_PLUS_BW=$THIS_DIR/data/ENCFF960OFK_plus.bw
SIGNAL_MINUS_BW=$THIS_DIR/data/ENCFF941MGK_minus.bw
ATAC_BW=$THIS_DIR/data/atac.bw

PEAKS=$OAK/Users/sheth/Data/ENCODE/GM12878/EP300/ENCFF926AKK.bed.gz
NEGATIVES=$PROJECT_DIR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed
FOLD_JSON=$PROJECT_DIR/reference/hg38_five_folds.json
OUTPUT_DIR=$THIS_DIR/GM12878_multimodal_BPNet/models/atac/fold${FOLD}

mkdir -p $OUTPUT_DIR

echo "$(date): Starting GM12878 multimodal training fold $FOLD"
echo "GPU: $CUDA_VISIBLE_DEVICES"

module load devel pixi/0.53.0

pixi run -e multimodal python $PROJECT_DIR/scripts/train_multimodal_bpnet.py \
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

echo "$(date): Training complete for fold $FOLD"
