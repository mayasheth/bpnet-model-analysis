#!/bin/bash
#SBATCH -p gpu
#SBATCH -t 24:00:00
#SBATCH --mem=120G
#SBATCH -o log/trainresgm_%x.%j.txt
#SBATCH -e log/trainresgm_%x.%j.txt
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Residual-objective training in GM12878: predict `observed - atac_pred`, offset from
# models/gm12878_atac5p_hw500_clw10, applied to the count loss only.
#
# WHY REPEAT THE GRID HERE. GM12878 is the harder problem in the way that matters for this
# question, and the two relevant facts pull in OPPOSITE directions, so the outcome is not a
# formality:
#   * Model-free ATAC-H3K27ac coupling is LOWER (top-quintile 0.409 vs K562's 0.510), so
#     accessibility explains less of the target and there is MORE residual left for
#     sequence to explain.
#   * Its inter-replicate ceiling is HIGHER (0.832 vs 0.760), so the target is measured more
#     cleanly and less of that residual is noise.
# In K562 the finding was that the residual objective helps only a blind input: +0.310 for
# sequence-only, -0.037 for multimodal, ~0 for the ATAC control. If that is a property of
# the objective it should reproduce here; if it is a property of K562's unusually tight
# ATAC-H3K27ac coupling it should not.
#
# The three standard GM12878 models already exist, so only the residual arms are trained.
#
# ATAC INPUT: the current $TRANS/data/atac.bw, matching the existing GM12878 models. NOT the
# new ChromBPNet-style atac_5p.bw -- these must be comparable to the existing GM12878
# numbers, and the 5' ATAC switch invalidates every ATAC-input model at once.
#
# gpu partition only: 2 of the 10 K562 residual folds were preempted on `owners`, and the
# trainer has no resume.
#
# Usage: sbatch 1.10.submit_training_residual_gm12878.sh MODE FOLD
#   MODE  sequence | multimodal | atac
#   FOLD  0-4

set -euo pipefail
export PYTHONUNBUFFERED=1

MODE=${1:?Usage: sbatch 1.10.submit_training_residual_gm12878.sh MODE FOLD}
FOLD=${2:?Usage: sbatch 1.10.submit_training_residual_gm12878.sh MODE FOLD}
case "$MODE" in sequence|multimodal|atac) ;; *) echo "bad MODE '$MODE'" >&2; exit 1 ;; esac

PROJECT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
PROJ="$PROJECT_DIR/2026_0824_H3K27ac_model"
TRANS="$PROJECT_DIR/2026_0606_GM12878_transferability"
PY="$PROJECT_DIR/.pixi/envs/multimodal/bin/python"

N_LAYERS=8
TRIMMING=$(( 47 + 2 + 4 + 8 + 16 + 32 + 64 + 128 + 256 ))
HALF_WINDOW=${HALF_WINDOW:-500}
OUT_WINDOW=$(( 2 * HALF_WINDOW ))
IN_WINDOW=$(( OUT_WINDOW + 2 * TRIMMING ))
COUNT_LOSS_WEIGHT=${COUNT_LOSS_WEIGHT:-10}
MAX_NEGATIVES=${MAX_NEGATIVES:-50000}

GENOME="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
ATAC_BW="$TRANS/data/atac.bw"
ELEMENTS="$TRANS/reference/GM12878_candidate_elements.narrowPeak"
NEGATIVES="$PROJECT_DIR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed"
FOLDS="$PROJECT_DIR/reference/hg38_five_folds.json"
OFFSET_MODEL="$PROJ/models/gm12878_atac5p_hw500_clw10"

OUT_DIR="$PROJ/models/gm12878_residual5p_${MODE}_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"
mkdir -p "$OUT_DIR" "$PROJ/log"

echo "GM12878 RESIDUAL mode=$MODE fold=$FOLD hw=$HALF_WINDOW clw=$COUNT_LOSS_WEIGHT"
echo "offset_model=$OFFSET_MODEL"
echo "out_dir=$OUT_DIR"

# --count-offset-model requires --accessibility-bw even in sequence mode, since the offset
# model reads accessibility even when this model does not.
GENOME_ARG=(); [[ "$MODE" != "atac" ]] && GENOME_ARG=(--genome "$GENOME")

cd "$PROJ"
$PY "$PROJECT_DIR/scripts/train_multimodal_bpnet.py" \
    --mode "$MODE" --peaks "$ELEMENTS" --negatives "$NEGATIVES" \
    ${GENOME_ARG[@]+"${GENOME_ARG[@]}"} \
    --signal-plus-bw "$PROJ/data/gm12878_h3k27ac_5p_plus.bw" \
    --signal-minus-bw "$PROJ/data/gm12878_h3k27ac_5p_minus.bw" \
    --accessibility-bw "$ATAC_BW" \
    --count-offset-model "$OFFSET_MODEL" \
    --fold "$FOLDS" --fold-key "$FOLD" --output-dir "$OUT_DIR" \
    --in-window "$IN_WINDOW" --out-window "$OUT_WINDOW" --n-layers "$N_LAYERS" \
    --count-loss-weight "$COUNT_LOSS_WEIGHT" --max-negatives "$MAX_NEGATIVES" \
    --negative-ratio 0.1

echo "Done: $OUT_DIR"
