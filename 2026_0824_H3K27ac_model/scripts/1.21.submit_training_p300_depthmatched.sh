#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 24:00:00
#SBATCH --mem=120G
#SBATCH -o /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model/log/p300depth.%j.txt
#SBATCH -e /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model/log/p300depth.%j.txt
#SBATCH --job-name=p300depth
#SBATCH -n 1
#SBATCH --ntasks 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'

#
# EXACTLY THREE changes vs 2026_0529/scripts/1.1.submit_training_atac.sh: the signal bigwigs,
# the peak set, and the output directory. Everything else, including --max-negatives being
# unset and --count-loss-weight 1.0, is inherited so the only difference from the full-depth
# model is the training data volume.
#
# Tests whether the p300 transfer asymmetry is explained by training-signal volume. K562
# transferred to GM12878 keeping 83% of its local advantage; GM12878 transferred to K562
# keeping none, while both fit their own cell type equally well. K562 carries 2.30x the reads
# in peaks. This arm gives K562 GM12878's budget: 30.0M mapped reads (from 51.1M) and 21,068
# peaks (from 28,532), built by scripts/0.31.
#
# The full-depth arm is the EXISTING 2026_0529/models/atac. That is only valid because 0.31
# rebuilt the full-depth track by the same code path and it matched the existing target at
# r = 1.0000 on both strands with a sum ratio of 0.9984, so construction is not a confound.
set -euo pipefail

FOLD=${1:?Usage: sbatch 1.1.submit_training_atac.sh FOLD_NUM}

SCRIPT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0529_multimodal_p300_model/scripts"
PROJECT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"

# Paths
GENOME="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
SIGNAL_PLUS_BW="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model/data/p300_depthmatched/ep300_depthmatched_plus.bw"
SIGNAL_MINUS_BW="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model/data/p300_depthmatched/ep300_depthmatched_minus.bw"
ATAC_BW="$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac.bw"

PEAKS="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model/data/p300_depthmatched/ep300_peaks_depthmatched.narrowPeak"
NEGATIVES="$PROJECT_DIR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed"
FOLD_JSON="$PROJECT_DIR/reference/hg38_five_folds.json"
OUTPUT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model/models/p300_depthmatched/fold${FOLD}"

mkdir -p "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model/log" "$OUTPUT_DIR"

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
