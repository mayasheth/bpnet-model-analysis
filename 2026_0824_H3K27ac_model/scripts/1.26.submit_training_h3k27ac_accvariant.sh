#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 24:00:00
#SBATCH --mem=120G
#SBATCH -o log/accvar_%x.%j.txt
#SBATCH -e log/accvar_%x.%j.txt
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# H3K27ac model with an ARBITRARY accessibility track, supplied by env var.
#
# EXACTLY TWO changes vs 1.22: the accessibility bigwig and the output directory, both now
# parameters rather than literals. 1.22 was itself 1.11 with those same two edits, so
# copying it a third and fourth time by hand would fork the training invocation four ways
# for no reason. Everything else stays byte-identical, which is what makes the accessibility
# track the only thing that can explain a difference.
#
# THE TWO ARMS THIS EXISTS FOR:
#   converted   data/abc_predicted/convdnase_prof_k562.bw -- DNase predicted from ATAC by
#               the converter, painted at base resolution by 4.20. Tests whether the
#               converter is worth anything where its strength can actually express itself:
#               unlike ABC, this model's accessibility branch reads base resolution.
#   smooth250   data/k562_dnase_5p_smooth250.bw -- real DNase with structure erased and
#               magnitude preserved (0.36). The control. F-010 showed real DNase beats ATAC
#               here but nobody measured WHICH property wins, and the answer decides whether
#               the converter can ever work. See figures/fig18_control_schematic.
#
# The anchors already exist and are not retrained: multimodal5p_accs5p_hw500_clw10 (ATAC
# input, 1.11) and multimodal5p_dnase_hw500_clw10 (real DNase input, 1.22).
#
# ELEMENTS ARE THE DNase-DERIVED SET, matching 1.11 and 1.22. This differs from the
# CONVERTER's own training (1.24 used ATAC-derived elements, so it never trained on regions
# chosen by the signal it predicts). That asymmetry is deliberate: here the comparison is
# against 1.11 and 1.22, so the element set must match THEM or the arm is not comparable to
# its own anchors.
#
# Usage: sbatch 1.26.submit_training_h3k27ac_accvariant.sh MODE FOLD
# Env:   ACC_BW (required), ACC_TAG (required, names the output dir),
#        COUNT_LOSS_WEIGHT (default 10), HALF_WINDOW (default 500)
set -euo pipefail
export PYTHONUNBUFFERED=1

MODE=${1:?Usage: sbatch 1.26... MODE FOLD}
FOLD=${2:?Usage: sbatch 1.26... MODE FOLD}
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
ELEMENTS="$PROJECT_DIR/reference/K562_DNase_candidate_elements.narrowPeak"
NEGATIVES="$PROJECT_DIR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed"
FOLDS="$PROJECT_DIR/reference/hg38_five_folds.json"

OUT_DIR="$PROJ/models/${MODE}5p_acc-${ACC_TAG}_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"
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
    --signal-plus-bw "$PROJ/data/h3k27ac_5p_plus.bw" \
    --signal-minus-bw "$PROJ/data/h3k27ac_5p_minus.bw" \
    ${ACC_ARG[@]+"${ACC_ARG[@]}"} \
    --fold "$FOLDS" --fold-key "$FOLD" --output-dir "$OUT_DIR" \
    --in-window "$IN_WINDOW" --out-window "$OUT_WINDOW" --n-layers "$N_LAYERS" \
    --count-loss-weight "$COUNT_LOSS_WEIGHT" --max-negatives "$MAX_NEGATIVES" \
    --negative-ratio 0.1

echo "Done: $OUT_DIR"
