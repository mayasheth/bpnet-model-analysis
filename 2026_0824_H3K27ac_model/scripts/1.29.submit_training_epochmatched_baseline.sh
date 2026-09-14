#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 24:00:00
#SBATCH --mem=120G
#SBATCH -o log/epmatch_%x.%j.txt
#SBATCH -e log/epmatch_%x.%j.txt
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# EPOCH-MATCHED CONTROL FOR THE MULTI-TASK ARM. 1.11 with early stopping effectively off.
#
# WHY THIS EXISTS. The multi-task arm (1.28) beat this baseline by +0.0109 [+0.0006,
# +0.0211] top-quintile, p=0.043. But early stopping watches the TOTAL validation loss, and
# for the multi-task arm that total includes the reweighted DNase profile term, so the two
# arms stop on different objectives and trained for different lengths: 53-99 epochs against
# 32-55, every fold longer. Two explanations fit the win equally well, and one of them has
# nothing to do with the auxiliary task:
#   (a) a learnable base-resolution task improves the shared trunk
#   (b) the multi-task arm simply trained longer
# This control removes (b) by letting the baseline run the same 100 epochs. Then either the
# baseline catches up, and the multi-task result is an artefact of the stopping rule, or it
# does not, and the auxiliary task is doing real work.
#
# EARLY STOPPING IS SET TO THE EPOCH BUDGET, NOT REMOVED. The checkpoint saved is still
# best-validation-loss, so this does not hand the baseline an overfitted model; it only
# stops the run from ending before epoch 100. That keeps the two arms' checkpoint-selection
# rule identical, which is the point of a control.
#
# Usage: sbatch 1.29.submit_training_epochmatched_baseline.sh FOLD
# Env:   MAX_EPOCHS (default 100, matching 1.28's budget), COUNT_LOSS_WEIGHT (default 10),
#        HALF_WINDOW (default 500)
set -euo pipefail
export PYTHONUNBUFFERED=1

FOLD=${1:?Usage: sbatch 1.29.submit_training_epochmatched_baseline.sh FOLD}

PROJECT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
PROJ="$PROJECT_DIR/2026_0824_H3K27ac_model"
PY="$PROJECT_DIR/.pixi/envs/multimodal/bin/python"

N_LAYERS=8
TRIMMING=$(( 47 + 2 + 4 + 8 + 16 + 32 + 64 + 128 + 256 ))
HALF_WINDOW=${HALF_WINDOW:-500}
OUT_WINDOW=$(( 2 * HALF_WINDOW ))
IN_WINDOW=$(( OUT_WINDOW + 2 * TRIMMING ))
COUNT_LOSS_WEIGHT=${COUNT_LOSS_WEIGHT:-10}
MAX_EPOCHS=${MAX_EPOCHS:-100}
MAX_NEGATIVES=${MAX_NEGATIVES:-50000}

GENOME="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
ATAC_BW="$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac_5p.bw"
ELEMENTS="$PROJECT_DIR/reference/K562_DNase_candidate_elements.narrowPeak"
NEGATIVES="$PROJECT_DIR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed"
FOLDS="$PROJECT_DIR/reference/hg38_five_folds.json"

OUT_DIR="$PROJ/models/multimodal5p_accs5p_ep${MAX_EPOCHS}_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"
mkdir -p "$OUT_DIR" "$PROJ/log"

echo "fold=$FOLD max_epochs=$MAX_EPOCHS (early stopping set to the budget, so it cannot fire)"
echo "out_dir=$OUT_DIR"

cd "$PROJ"
$PY "$PROJECT_DIR/scripts/train_multimodal_bpnet.py" \
    --mode multimodal --peaks "$ELEMENTS" --negatives "$NEGATIVES" \
    --genome "$GENOME" \
    --signal-plus-bw "$PROJ/data/h3k27ac_5p_plus.bw" \
    --signal-minus-bw "$PROJ/data/h3k27ac_5p_minus.bw" \
    --accessibility-bw "$ATAC_BW" \
    --fold "$FOLDS" --fold-key "$FOLD" --output-dir "$OUT_DIR" \
    --in-window "$IN_WINDOW" --out-window "$OUT_WINDOW" --n-layers "$N_LAYERS" \
    --count-loss-weight "$COUNT_LOSS_WEIGHT" --max-negatives "$MAX_NEGATIVES" \
    --negative-ratio 0.1 \
    --max-epochs "$MAX_EPOCHS" --early-stopping "$MAX_EPOCHS"

echo "Done: $OUT_DIR"
