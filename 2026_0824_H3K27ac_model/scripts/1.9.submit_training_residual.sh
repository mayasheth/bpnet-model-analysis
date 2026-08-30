#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 24:00:00
#SBATCH --mem=120G
#SBATCH -o log/trainres_%x.%j.txt
#SBATCH -e log/trainres_%x.%j.txt
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Train a model on the RESIDUAL objective: predict `observed - atac_pred`, where atac_pred
# is the held-out prediction of models/atac5p_hw500_clw10, applied to the count loss only.
#
# Fills in the residual grid. Until now exactly ONE residual model existed
# (residual5pFIXED_hw500_clw10, mode=sequence), so we could not separate "the residual is
# sequence-predictable" from "residual training helps".
#
#   multimodal  the informative cell. If sequence+ATAC trained on the residual beats the
#               jointly-trained multimodal model's residual r of 0.551, then the OBJECTIVE
#               matters even when ATAC is already an input. If it does not, joint training
#               is already optimal and residual training is purely an interpretability tool.
#   atac        negative control. Expected ~0 by construction: it would predict an ATAC
#               model's own errors from the same ATAC input. A clearly non-zero result would
#               mean the residual metric is measuring something other than what we think.
#
# ATAC INPUT: deliberately the CURRENT data/atac.bw, the same track the existing
# sequence5p / multimodal5p / residual5pFIXED models used -- NOT a ChromBPNet-style 5'-end
# rebuild. These runs exist to complete a comparison against those numbers, so the input
# must match them. The 5'-end ATAC rebuild is a separate, larger piece of work that
# invalidates all ATAC-input models at once and must be done as a set.
#
# Usage: sbatch 1.9.submit_training_residual.sh MODE FOLD
#   MODE  multimodal | atac | sequence
#   FOLD  0-4

set -euo pipefail
export PYTHONUNBUFFERED=1

MODE=${1:?Usage: sbatch 1.9.submit_training_residual.sh MODE FOLD}
FOLD=${2:?Usage: sbatch 1.9.submit_training_residual.sh MODE FOLD}
case "$MODE" in sequence|multimodal|atac) ;; *) echo "bad MODE '$MODE'" >&2; exit 1 ;; esac

PROJECT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
PROJ="$PROJECT_DIR/2026_0824_H3K27ac_model"
PY="$PROJECT_DIR/.pixi/envs/multimodal/bin/python"

N_LAYERS=8
TRIMMING=$(( 47 + 2 + 4 + 8 + 16 + 32 + 64 + 128 + 256 ))
HALF_WINDOW=${HALF_WINDOW:-500}
OUT_WINDOW=$(( 2 * HALF_WINDOW ))
IN_WINDOW=$(( OUT_WINDOW + 2 * TRIMMING ))
COUNT_LOSS_WEIGHT=${COUNT_LOSS_WEIGHT:-10}
MAX_NEGATIVES=${MAX_NEGATIVES:-50000}

GENOME="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
ATAC_BW="$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac.bw"
ELEMENTS="$PROJECT_DIR/reference/K562_DNase_candidate_elements.narrowPeak"
NEGATIVES="$PROJECT_DIR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed"
FOLDS="$PROJECT_DIR/reference/hg38_five_folds.json"
OFFSET_MODEL="$PROJ/models/atac5p_hw500_clw10"

OUT_DIR="$PROJ/models/residual5p_${MODE}_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"
mkdir -p "$OUT_DIR" "$PROJ/log"

echo "RESIDUAL mode=$MODE fold=$FOLD hw=$HALF_WINDOW clw=$COUNT_LOSS_WEIGHT"
echo "offset_model=$OFFSET_MODEL"
echo "out_dir=$OUT_DIR"

# --count-offset-model requires --accessibility-bw even in sequence mode, since the offset
# model reads accessibility even when this model does not.
GENOME_ARG=(); [[ "$MODE" != "atac" ]] && GENOME_ARG=(--genome "$GENOME")

cd "$PROJ"
$PY "$PROJECT_DIR/scripts/train_multimodal_bpnet.py" \
    --mode "$MODE" --peaks "$ELEMENTS" --negatives "$NEGATIVES" \
    ${GENOME_ARG[@]+"${GENOME_ARG[@]}"} \
    --signal-plus-bw "$PROJ/data/h3k27ac_5p_plus.bw" \
    --signal-minus-bw "$PROJ/data/h3k27ac_5p_minus.bw" \
    --accessibility-bw "$ATAC_BW" \
    --count-offset-model "$OFFSET_MODEL" \
    --fold "$FOLDS" --fold-key "$FOLD" --output-dir "$OUT_DIR" \
    --in-window "$IN_WINDOW" --out-window "$OUT_WINDOW" --n-layers "$N_LAYERS" \
    --count-loss-weight "$COUNT_LOSS_WEIGHT" --max-negatives "$MAX_NEGATIVES" \
    --negative-ratio 0.1

echo "Done: $OUT_DIR"
