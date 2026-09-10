#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 24:00:00
#SBATCH --mem=120G
#SBATCH -o log/convgm_%x.%j.txt
#SBATCH -e log/convgm_%x.%j.txt
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# THE ATAC->DNase CONVERTER, TRAINED IN GM12878. Sequence + ATAC in, DNase out.
#
# EXACTLY THREE changes vs 1.24: the accessibility input, the target bigwigs and the element
# set all move to GM12878, plus the output directory. Every hyperparameter is byte-identical
# to the K562 converter, so the cell type is the only thing that can explain a difference.
#
# WHY IT EXISTS. The CRISPR benchmark is K562-only, so the only way to ask whether a
# converter trained somewhere else is usable is to train one somewhere else and apply it
# here. That is the deployment case: you have ATAC in your cell type, a converter from
# someone else's, and no DNase anywhere. F-009 found the GM12878 -> K562 direction retains
# nothing for p300, so this arm is a real test rather than a formality.
#
# WHY THIS EXISTS. DNase is the better accessibility representation to be in -- it is the
# one input change that survived transfer (+0.056 over GM12878's own ATAC-only floor,
# p=0.020, where the ATAC-input model manages +0.013, p=0.25). But most cell types have
# ATAC and no DNase. A converter turns that result into something deployable.
#
# THE PROFILE IS THE DELIVERABLE, NOT THE COUNTS. This is the first target in the project
# where that is true. The downstream model reads base-resolution accessibility over a
# 2,114 bp window, so a converter that predicts per-element DNase counts and paints them
# flat would strip exactly what the accessibility branch uses. 0.33 is what unblocked this:
# DNase has a real base-resolution profile (1 bp top-quintile ceiling 0.848 K562, 0.686
# GM12878) where H3K27ac has none (0.21, 0.18).
#
# COUNT_LOSS_WEIGHT IS SETTLED. The K562 pilot found clw=1 and clw=10 differ by 0.001 in
# profile shape at every bin size, with clw=10 slightly better on counts, so the project
# standard is used and no sweep is repeated here.
#
# ELEMENT SET ASYMMETRY, STATED BECAUSE IT IS NOT FIXABLE CHEAPLY. The K562 converter trains
# on ATAC-derived elements (1.24), so it never sees regions chosen by the signal it predicts.
# No ATAC-derived GM12878 candidate set exists: 0.26 derives K562's from a specific rE2G ATAC
# run with no GM12878 counterpart, and GM12878_candidate_elements.narrowPeak is what every
# GM12878 model here has used. So this arm trains on a differently derived element set than
# its K562 counterpart. Element derivation made no difference on the H3K27ac task (p=0.83),
# which is the basis for accepting it, but if this arm underperforms the asymmetry is a live
# alternative to a transfer explanation and must be reported as one.
#
# Usage: sbatch 1.25.submit_training_converter_gm12878.sh MODE FOLD
#   MODE  multimodal | atac | sequence
#   FOLD  0-4
# Env:  COUNT_LOSS_WEIGHT (default 10), HALF_WINDOW (default 500)
set -euo pipefail
export PYTHONUNBUFFERED=1

MODE=${1:?Usage: sbatch 1.25.submit_training_converter_gm12878.sh MODE FOLD}
FOLD=${2:?Usage: sbatch 1.25.submit_training_converter_gm12878.sh MODE FOLD}
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
TRANS="$PROJECT_DIR/2026_0606_GM12878_transferability"
# Input stays ATAC. That is the whole point: the converter must run where DNase does not exist.
ATAC_BW="$TRANS/data/atac_5p.bw"

ELEMENTS="$TRANS/reference/GM12878_candidate_elements.narrowPeak"
NEGATIVES="$PROJECT_DIR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed"
FOLDS="$PROJECT_DIR/reference/hg38_five_folds.json"

# Pooled stranded DNase from 0.34 -- NOT the per-replicate tracks 0.33 wrote for the
# ceiling, and not the unstranded 0.32 track, which is an input and has no strand channels.
SIG_PLUS="$PROJ/data/gm12878_dnase_5p_plus.bw"
SIG_MINUS="$PROJ/data/gm12878_dnase_5p_minus.bw"
for f in "$SIG_PLUS" "$SIG_MINUS"; do
    [[ -s "$f" ]] || { echo "ERROR: $f missing; run 0.34 first" >&2; exit 1; }
done

TAG="conv_atac2dnase_gm12878"
OUT_DIR="$PROJ/models/gm12878_${MODE}5p_${TAG}_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"
mkdir -p "$OUT_DIR" "$PROJ/log"

echo "mode=$MODE fold=$FOLD hw=$HALF_WINDOW in=$IN_WINDOW out=$OUT_WINDOW clw=$COUNT_LOSS_WEIGHT"
echo "elements=$ELEMENTS"
echo "target=$(basename "$SIG_PLUS") / $(basename "$SIG_MINUS")"
echo "out_dir=$OUT_DIR"

ACC_ARG=(); [[ "$MODE" != "sequence" ]] && ACC_ARG=(--accessibility-bw "$ATAC_BW")
GENOME_ARG=(); [[ "$MODE" != "atac" ]] && GENOME_ARG=(--genome "$GENOME")

cd "$PROJ"
$PY "$PROJECT_DIR/scripts/train_multimodal_bpnet.py" \
    --mode "$MODE" --peaks "$ELEMENTS" --negatives "$NEGATIVES" \
    ${GENOME_ARG[@]+"${GENOME_ARG[@]}"} \
    --signal-plus-bw "$SIG_PLUS" \
    --signal-minus-bw "$SIG_MINUS" \
    ${ACC_ARG[@]+"${ACC_ARG[@]}"} \
    --fold "$FOLDS" --fold-key "$FOLD" --output-dir "$OUT_DIR" \
    --in-window "$IN_WINDOW" --out-window "$OUT_WINDOW" --n-layers "$N_LAYERS" \
    --count-loss-weight "$COUNT_LOSS_WEIGHT" --max-negatives "$MAX_NEGATIVES" \
    --negative-ratio 0.1

echo "Done: $OUT_DIR"
