#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 24:00:00
#SBATCH --mem=120G
#SBATCH -o log/gmaccvar_%x.%j.txt
#SBATCH -e log/gmaccvar_%x.%j.txt
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# H3K27ac model in GM12878 with an ARBITRARY accessibility track, supplied by env var.
#
# EXACTLY THREE changes vs 1.26: the H3K27ac target, the element set and the output prefix
# all move to GM12878. The accessibility track stays a parameter. Byte-identical otherwise,
# so a K562-vs-GM12878 difference is attributable to the cell type and nothing else.
#
# WHY IT EXISTS. F-013 scored four accessibility inputs in K562 only. The two anchors already
# exist in GM12878 (gm12878_multimodal5p_accs5p for ATAC input, gm12878_multimodal5p_dnase
# for real DNase), so adding the converted and smoothed arms here completes a 2x2: each input
# variant trained in each cell type and applied to both. That is what makes a transfer claim
# possible, and F-009 is the reason it matters -- p300 transfer was strongly asymmetric, and
# nothing says H3K27ac under a DNase-family input behaves the same way in both directions.
#
# THE TWO ARMS THIS EXISTS FOR:
#   smooth250   data/gm12878_dnase_5p_smooth250.bw -- real DNase, structure erased, magnitude
#               kept (0.36). Runnable now.
#   converted   BLOCKED until the painted track's magnitude form is settled. F-013 showed the
#               K562 converted arm is resolvably WORSE than raw ATAC, and quantile mapping
#               made the compression worse rather than better, so painting GM12878 raw parts
#               (form-agnostic) is decoupled from choosing the transform on purpose.
#
# The anchors already exist and are not retrained: multimodal5p_accs5p_hw500_clw10 (ATAC
# input, 1.11) and multimodal5p_dnase_hw500_clw10 (real DNase input, 1.22).
#
# ELEMENTS are GM12878_candidate_elements, matching 1.12 and 1.23 so this arm is comparable
# to its own anchors. Same reasoning as 1.26's use of the K562 DNase-derived set.
#
# Usage: sbatch 1.27.submit_training_h3k27ac_gm12878_accvariant.sh MODE FOLD
# Env:   ACC_BW (required), ACC_TAG (required, names the output dir),
#        COUNT_LOSS_WEIGHT (default 10), HALF_WINDOW (default 500)
set -euo pipefail
export PYTHONUNBUFFERED=1

MODE=${1:?Usage: sbatch 1.27... MODE FOLD}
FOLD=${2:?Usage: sbatch 1.27... MODE FOLD}
case "$MODE" in sequence|multimodal|atac) ;; *) echo "bad MODE '$MODE'" >&2; exit 1 ;; esac
: "${ACC_BW:?set ACC_BW to the accessibility bigwig}"
: "${ACC_TAG:?set ACC_TAG to name the output directory}"
[[ -s "$ACC_BW" ]] || { echo "ERROR: ACC_BW '$ACC_BW' missing or empty" >&2; exit 1; }

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
TRANS="$PROJECT_DIR/2026_0606_GM12878_transferability"
ELEMENTS="$TRANS/reference/GM12878_candidate_elements.narrowPeak"
NEGATIVES="$PROJECT_DIR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed"
FOLDS="$PROJECT_DIR/reference/hg38_five_folds.json"

OUT_DIR="$PROJ/models/gm12878_${MODE}5p_acc-${ACC_TAG}_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"
mkdir -p "$OUT_DIR" "$PROJ/log"

echo "mode=$MODE fold=$FOLD hw=$HALF_WINDOW clw=$COUNT_LOSS_WEIGHT"
echo "accessibility=$ACC_BW"
echo "out_dir=$OUT_DIR"

ACC_ARG=(); [[ "$MODE" != "sequence" ]] && ACC_ARG=(--accessibility-bw "$ACC_BW")
GENOME_ARG=(); [[ "$MODE" != "atac" ]] && GENOME_ARG=(--genome "$GENOME")

cd "$PROJ"
$PY "$PROJECT_DIR/scripts/train_multimodal_bpnet.py" \
    --mode "$MODE" --peaks "$ELEMENTS" --negatives "$NEGATIVES" \
    ${GENOME_ARG[@]+"${GENOME_ARG[@]}"} \
    --signal-plus-bw "$PROJ/data/gm12878_h3k27ac_5p_plus.bw" \
    --signal-minus-bw "$PROJ/data/gm12878_h3k27ac_5p_minus.bw" \
    ${ACC_ARG[@]+"${ACC_ARG[@]}"} \
    --fold "$FOLDS" --fold-key "$FOLD" --output-dir "$OUT_DIR" \
    --in-window "$IN_WINDOW" --out-window "$OUT_WINDOW" --n-layers "$N_LAYERS" \
    --count-loss-weight "$COUNT_LOSS_WEIGHT" --max-negatives "$MAX_NEGATIVES" \
    --negative-ratio 0.1

echo "Done: $OUT_DIR"
