#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 24:00:00
#SBATCH --mem=120G
#SBATCH -o log/train_%x.%j.txt
#SBATCH -e log/train_%x.%j.txt
#SBATCH -n 1
#SBATCH --ntasks 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Train an H3K27ac model on K562 DNase candidate elements.
#
# Usage:
#   sbatch --job-name=k27_seq_hw500 1.1.submit_training.sh MODE FOLD HALF_WINDOW
#
#   MODE         sequence | multimodal | atac
#   FOLD         0-4
#   HALF_WINDOW  counting half-width in bp (500 or 1000; see results/window_tradeoff.tsv)
#
# Windows: out_window = 2*HALF_WINDOW, and in_window is forced to
# out_window + 2*trimming with trimming = 47 + sum(2^i for i in 1..n_layers) = 557
# at the default n_layers=8. The trainer asserts this, so it cannot drift.
#
# Differences from the p300 setup, all deliberate:
#   * Windows are centred on candidate ELEMENTS, not ChIP peak summits — the H3K27ac
#     summit sits on a flanking nucleosome, not on the element.
#   * Target is UNSTRANDED (no --signal-minus-bw), so n_outputs=1.
#   * The profile head is kept but DOWN-WEIGHTED via --count-loss-weight, because the
#     bimodal flanking pattern is real structure worth fitting, and reproducing the
#     central dip is a check the model learned H3K27ac biology rather than
#     "accessible implies acetylated".
#   * --max-negatives is capped: 10x of ~120k elements would be ~1.2M windows (~55 GB).

set -euo pipefail

MODE=${1:?Usage: sbatch 1.1.submit_training.sh MODE FOLD HALF_WINDOW}
FOLD=${2:?Usage: sbatch 1.1.submit_training.sh MODE FOLD HALF_WINDOW}
HALF_WINDOW=${3:?Usage: sbatch 1.1.submit_training.sh MODE FOLD HALF_WINDOW}

case "$MODE" in
  sequence|multimodal|atac) ;;
  *) echo "ERROR: MODE must be sequence, multimodal, or atac (got '$MODE')" >&2; exit 1 ;;
esac

PROJECT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
PROJ="$PROJECT_DIR/2026_0824_H3K27ac_model"
PY="$PROJECT_DIR/.pixi/envs/multimodal/bin/python"   # pixi is not on PATH on Sherlock

# Window geometry — must satisfy in_window - 2*trimming == out_window
N_LAYERS=8
TRIMMING=$(( 47 + 2 + 4 + 8 + 16 + 32 + 64 + 128 + 256 ))   # 557 at n_layers=8
OUT_WINDOW=$(( 2 * HALF_WINDOW ))
IN_WINDOW=$(( OUT_WINDOW + 2 * TRIMMING ))

# Inputs
GENOME="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
H3K27AC_BW="/oak/stanford/groups/engreitz/Users/sheth/Data/share/IGV/ENCSR000AKP_coverage.bw"
ATAC_BW="$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac.bw"
ELEMENTS="$PROJECT_DIR/reference/K562_DNase_candidate_elements.narrowPeak"
NEGATIVES="$PROJECT_DIR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed"
FOLDS="$PROJECT_DIR/reference/hg38_five_folds.json"

# Hyperparameters
COUNT_LOSS_WEIGHT=${COUNT_LOSS_WEIGHT:-10}   # >1 down-weights the profile head
MAX_NEGATIVES=${MAX_NEGATIVES:-50000}
NEGATIVE_RATIO=${NEGATIVE_RATIO:-0.1}

OUT_DIR="$PROJ/models/${MODE}_hw${HALF_WINDOW}/fold${FOLD}"
mkdir -p "$OUT_DIR" "$PROJ/log"

echo "=================================================="
echo "mode=$MODE fold=$FOLD half_window=$HALF_WINDOW"
echo "in_window=$IN_WINDOW out_window=$OUT_WINDOW (trimming=$TRIMMING)"
echo "count_loss_weight=$COUNT_LOSS_WEIGHT max_negatives=$MAX_NEGATIVES"
echo "out_dir=$OUT_DIR"
echo "=================================================="

# sequence mode never reads accessibility, so skip extracting it
ACC_ARG=()
if [[ "$MODE" != "sequence" ]]; then
  ACC_ARG=(--accessibility-bw "$ATAC_BW")
fi

# atac mode needs no genome
GENOME_ARG=()
if [[ "$MODE" != "atac" ]]; then
  GENOME_ARG=(--genome "$GENOME")
fi

cd "$PROJ"
$PY "$PROJECT_DIR/scripts/train_multimodal_bpnet.py" \
    --mode "$MODE" \
    --peaks "$ELEMENTS" \
    --negatives "$NEGATIVES" \
    "${GENOME_ARG[@]}" \
    --signal-plus-bw "$H3K27AC_BW" \
    "${ACC_ARG[@]}" \
    --fold "$FOLDS" \
    --fold-key "$FOLD" \
    --output-dir "$OUT_DIR" \
    --in-window "$IN_WINDOW" \
    --out-window "$OUT_WINDOW" \
    --n-layers "$N_LAYERS" \
    --count-loss-weight "$COUNT_LOSS_WEIGHT" \
    --max-negatives "$MAX_NEGATIVES" \
    --negative-ratio "$NEGATIVE_RATIO"

echo "Done: $OUT_DIR"
