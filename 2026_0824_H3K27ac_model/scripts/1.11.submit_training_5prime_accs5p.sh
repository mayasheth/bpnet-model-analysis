#!/bin/bash
#
# ===== AUTO-DERIVED: ATAC input switched to ChromBPNet-convention 5' insertion counts =====
# Generated from 1.4.submit_training_5prime.sh by scripts/make_accs5p_scripts.py.
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
#SBATCH -o log/train5p_%x.%j.txt
#SBATCH -e log/train5p_%x.%j.txt
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Train an H3K27ac model on the 5'-END target (the corrected target definition).
#
# Usage: sbatch 1.4.submit_training_5prime.sh MODE FOLD
#   MODE  sequence | multimodal | atac
#   FOLD  0-4
#
# Differences from 1.1 (fragment-extended target):
#   * Target is data/h3k27ac_5p_{plus,minus}.bw — stranded 5' ends, matching p300's
#     genomecov -5 -dz convention. bpnetlite sums both channels for the count target, so
#     counts are identical to an unstranded version; the profile head gets the better task.
#   * COUNT_LOSS_WEIGHT defaults to 10, the measured optimum on this target (sweep:
#     1 -> 0.437, 3 -> 0.484, 10 -> 0.496, 100 -> 0.467, 1000 -> 0.464). The fragment
#     target wanted 1000 because MNLL was inflated ~15-24x by smearing.
#
# Window is fixed at +/-500 (out 1000 / in 2114), the choice resolved on the fragment
# target. Worth re-checking on 5' since less smearing means less neighbour bleed.

set -euo pipefail
export PYTHONUNBUFFERED=1

MODE=${1:?Usage: sbatch 1.4.submit_training_5prime.sh MODE FOLD}
FOLD=${2:?Usage: sbatch 1.4.submit_training_5prime.sh MODE FOLD}
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
ATAC_BW="$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac_5p.bw"
ELEMENTS="$PROJECT_DIR/reference/K562_DNase_candidate_elements.narrowPeak"
NEGATIVES="$PROJECT_DIR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed"
FOLDS="$PROJECT_DIR/reference/hg38_five_folds.json"

OUT_DIR="$PROJ/models/${MODE}5p_accs5p_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"
mkdir -p "$OUT_DIR" "$PROJ/log"

echo "mode=$MODE fold=$FOLD hw=$HALF_WINDOW in=$IN_WINDOW out=$OUT_WINDOW clw=$COUNT_LOSS_WEIGHT"
echo "out_dir=$OUT_DIR"

ACC_ARG=(); [[ "$MODE" != "sequence" ]] && ACC_ARG=(--accessibility-bw "$ATAC_BW")
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
