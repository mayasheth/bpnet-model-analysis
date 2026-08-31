#!/bin/bash
#
# ===== AUTO-DERIVED: ATAC input switched to ChromBPNet-convention 5' insertion counts =====
# Generated from 1.8.submit_training_gm12878.sh by scripts/make_accs5p_scripts.py.
# EXACTLY TWO changes vs the original: the accessibility bigwig, and the output directory.
# Every other hyperparameter is byte-identical, so a difference in results is attributable
# to the input definition and nothing else.
#
# Old input: genomecov -bg over the FULL tagAlign interval -> coverage scales with read
#            length (94-95 bp K562/GM12878, 35-36 bp TeloHAEC).
# New input: genomecov -bg -5 -> single-base 5' insertion counts, read-length independent,
#            matching chrombpnet/helpers/preprocessing/reads_to_bigwig.py.
# Verified: sum(old)/sum(new) equals the mean read length in all six tracks.
#

#SBATCH -p owners,gpu
#SBATCH -t 24:00:00
#SBATCH --mem=120G
#SBATCH -o log/train_gm_%x.%j.txt
#SBATCH -e log/train_gm_%x.%j.txt
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Train an H3K27ac model on GM12878, for the RECIPROCAL transfer analysis.
#
# Everything so far trains on K562 and tests on GM12878, which cannot distinguish
# "H3K27ac has little transferable sequence signal" from "K562 is a poor cell type to
# learn from". Training the same three modalities on GM12878 and testing on K562 answers
# that: if the sequence component transfers better out of GM12878, the limitation is the
# training cell type rather than the target.
#
# Deliberately identical to 1.4.submit_training_5prime.sh except for the cell type:
# same architecture, window, count-loss weight, negative cap, folds and 5'-end target
# processing, so the two directions are comparable.
#
# Usage: sbatch 1.8.submit_training_gm12878.sh MODE FOLD
#   MODE  sequence | multimodal | atac
#   FOLD  0-4

set -euo pipefail
export PYTHONUNBUFFERED=1

MODE=${1:?Usage: sbatch 1.8.submit_training_gm12878.sh MODE FOLD}
FOLD=${2:?Usage: sbatch 1.8.submit_training_gm12878.sh MODE FOLD}
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
# GM12878 inputs throughout
ATAC_BW="$TRANS/data/atac_5p.bw"
ELEMENTS="$TRANS/reference/GM12878_candidate_elements.narrowPeak"
# GC-matched negatives are genome-wide and cell-type agnostic, so shared with K562
NEGATIVES="$PROJECT_DIR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed"
FOLDS="$PROJECT_DIR/reference/hg38_five_folds.json"

OUT_DIR="$PROJ/models/gm12878_${MODE}5p_accs5p_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"
mkdir -p "$OUT_DIR" "$PROJ/log"

echo "GM12878 training | mode=$MODE fold=$FOLD in=$IN_WINDOW out=$OUT_WINDOW clw=$COUNT_LOSS_WEIGHT"
echo "out_dir=$OUT_DIR"

ACC_ARG=(); [[ "$MODE" != "sequence" ]] && ACC_ARG=(--accessibility-bw "$ATAC_BW")
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
