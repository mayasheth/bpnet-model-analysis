#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 24:00:00
#SBATCH --mem=120G
#SBATCH -o log/conv_%x.%j.txt
#SBATCH -e log/conv_%x.%j.txt
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# THE ATAC->DNase CONVERTER. Sequence + ATAC in, DNase out, trained in K562.
#
# EXACTLY THREE changes vs 1.11 (models/multimodal5p_accs5p_hw500_clw10): the target
# bigwigs, the element set, and the output directory. Every hyperparameter is otherwise
# byte-identical, so the target swap is the only thing that can explain a difference.
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
# COUNT_LOSS_WEIGHT IS NOT SETTLED HERE. clw=10 is the measured optimum for the H3K27ac 5'
# target and is kept as the anchor arm so this model is comparable to every other model in
# the project. But that sweep scored the COUNTS head on a target with no profile to learn.
# With the profile as the deliverable the optimum may sit lower, so clw is an env var and
# 2.31 scores profile_pearson_topq across arms rather than assuming.
#
# RUN THE atac MODE TOO. It is the control that decides whether this is a converter at all:
# if ATAC alone predicts DNase as well as sequence+ATAC does, the model is a learned
# rescaling of the input and adds nothing a normalisation could not.
#
# Usage: sbatch 1.24.submit_training_converter.sh MODE FOLD
#   MODE  multimodal | atac | sequence
#   FOLD  0-4
# Env:  COUNT_LOSS_WEIGHT (default 10), HALF_WINDOW (default 500),
#       ELEMENTS_SET (default atac; atac | dnase)
set -euo pipefail
export PYTHONUNBUFFERED=1

MODE=${1:?Usage: sbatch 1.24.submit_training_converter.sh MODE FOLD}
FOLD=${2:?Usage: sbatch 1.24.submit_training_converter.sh MODE FOLD}
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
# Input stays ATAC. That is the whole point: the converter must run where DNase does not exist.
ATAC_BW="$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac_5p.bw"

# ELEMENTS ARE ATAC-DERIVED BY DEFAULT, AND THAT IS A DELIBERATE DEPARTURE FROM 1.11.
# Every other K562 model in this project trains on K562_DNase_candidate_elements. For the
# converter that would be circular in a way that matters: at deployment in a cell type with
# no DNase the element set can only come from ATAC, so training on DNase-derived elements
# means the model never sees the regions ATAC nominates and DNase does not -- exactly the
# regions where a conversion is doing work rather than copying. The swap is free: element
# derivation, ATAC- vs DNase-derived, made no difference on the H3K27ac task (p=0.83).
case "${ELEMENTS_SET:-atac}" in
  atac)  ELEMENTS="$PROJECT_DIR/reference/K562_ATAC_candidate_elements.narrowPeak" ;;
  dnase) ELEMENTS="$PROJECT_DIR/reference/K562_DNase_candidate_elements.narrowPeak" ;;
  *) echo "bad ELEMENTS_SET '${ELEMENTS_SET}' (atac|dnase)" >&2; exit 1 ;;
esac
NEGATIVES="$PROJECT_DIR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed"
FOLDS="$PROJECT_DIR/reference/hg38_five_folds.json"

# Pooled stranded DNase from 0.34 -- NOT the per-replicate tracks 0.33 wrote for the
# ceiling, and not the unstranded 0.32 track, which is an input and has no strand channels.
SIG_PLUS="$PROJ/data/k562_dnase_5p_plus.bw"
SIG_MINUS="$PROJ/data/k562_dnase_5p_minus.bw"
for f in "$SIG_PLUS" "$SIG_MINUS"; do
    [[ -s "$f" ]] || { echo "ERROR: $f missing; run 0.34 first" >&2; exit 1; }
done

TAG="conv_atac2dnase_${ELEMENTS_SET:-atac}el"
OUT_DIR="$PROJ/models/${MODE}5p_${TAG}_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"
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
