#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 24:00:00
#SBATCH --mem=120G
#SBATCH -o log/thp1_%x.%j.txt
#SBATCH -e log/thp1_%x.%j.txt
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# THP-1 H3K27ac with a DNase accessibility input. The third cell type.
#
# EXACTLY THREE changes vs 1.22 (K562 DNase input): the accessibility bigwig, the H3K27ac
# target, and the element set, all three being THP-1's. Every hyperparameter is
# byte-identical to K562's and GM12878's arms, which is what lets the three sit in one
# transfer matrix.
#
# WHY A THIRD CELL TYPE, AND WHY NOW. Two findings are stuck on the same gap. F-009 (p300)
# and F-014 (DNase inputs) both report a strongly asymmetric transfer: K562-trained models
# keep their advantage in GM12878, GM12878-trained models lose it in K562. With two cell
# types that is unattributable -- "K562 is a good TRAINING cell type" and "K562->GM12878 is
# a good PAIR" predict identical numbers. A third cell type separates them, and it is the
# only thing that does.
#
# WHY THP-1 AND NOT TeloHAEC. The surviving DNase results need REAL DNase: DNase as an input
# (F-010, F-014) and DNase as the ABC activity term (+0.0709, F-012). TeloHAEC has ATAC and no
# DNase, so it serves the ATAC-only panel rather than this one, and THP-1 is the only other
# cell type here with DNase.
#
# THIS PANEL IS NOT A DEPLOYMENT ARGUMENT, and an earlier version of this comment wrongly
# implied it was by claiming DNase exists in more cell types than ATAC. The opposite is true:
# **DNase is the SCARCER assay**, which is the whole reason the application target is
# ATAC-only and the reason a converter was attempted at all (closed, F-014). What this panel
# establishes is scientific -- whether a DNase-input model's advantage is a property of the
# training cell type, the pair, or the target (F-016, F-017) -- not that DNase-input models
# are the deployable choice. They are the less deployable one.
#
# WHAT THP-1 CANNOT DO, so nobody reads more into this panel than it holds:
#   * One DNase replicate, so there is NO DNase inter-replicate ceiling here and THP-1 cannot
#     be a converter target.
#   * No ATAC, so it cannot join an ATAC-input comparison and cannot train a converter.
#   * It IS a valid DNase-input panel member, which is exactly what is wanted.
# The H3K27ac inter-replicate ceiling DOES exist (0.974 raw, 0.993 Spearman-Brown corrected,
# both replicates paired-end, read 1 only), and that is what these arms are read against.
#
# TWO MODES, AND THE FLOOR IS NOT OPTIONAL. `multimodal` is the arm of interest.
# `atac` mode here means accessibility-ONLY (DNase, no sequence) and supplies the
# `atac_LOCAL` floor every transfer table in this project is read against; without it a
# transferred model has nothing to clear.
#
# Usage: sbatch 1.30.submit_training_dnase_thp1.sh MODE FOLD
#   MODE  sequence | multimodal | atac
#   FOLD  0-4
# Env:  ACC_BW + ACC_TAG  optional, to train on a different accessibility track. Both must
#       be set together, and ACC_TAG names the output directory so two input variants can
#       never collide in one directory. Used for the depth-matched arms of F-016.
set -euo pipefail
export PYTHONUNBUFFERED=1

MODE=${1:?Usage: sbatch 1.30.submit_training_dnase_thp1.sh MODE FOLD}
FOLD=${2:?Usage: sbatch 1.30.submit_training_dnase_thp1.sh MODE FOLD}
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
DNASE_BW="$PROJ/data/thp1_dnase_5p.bw"
TAG="dnase"
if [[ -n "${ACC_BW:-}${ACC_TAG:-}" ]]; then
    : "${ACC_BW:?set ACC_BW and ACC_TAG together}"
    : "${ACC_TAG:?set ACC_BW and ACC_TAG together}"
    DNASE_BW="$ACC_BW"
    TAG="$ACC_TAG"
fi
ELEMENTS="$PROJECT_DIR/reference/THP1_DNase_candidate_elements.narrowPeak"
NEGATIVES="$PROJECT_DIR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed"
FOLDS="$PROJECT_DIR/reference/hg38_five_folds.json"
SIG_PLUS="$PROJ/data/thp1_h3k27ac_5p_plus.bw"
SIG_MINUS="$PROJ/data/thp1_h3k27ac_5p_minus.bw"

for f in "$DNASE_BW" "$ELEMENTS" "$SIG_PLUS" "$SIG_MINUS"; do
    [[ -s "$f" ]] || { echo "ERROR: '$f' missing or empty; run 0.38 first" >&2; exit 1; }
done

OUT_DIR="$PROJ/models/thp1_${MODE}5p_${TAG}_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"
mkdir -p "$OUT_DIR" "$PROJ/log"

echo "mode=$MODE fold=$FOLD hw=$HALF_WINDOW in=$IN_WINDOW out=$OUT_WINDOW clw=$COUNT_LOSS_WEIGHT"
echo "accessibility=$(basename "$DNASE_BW")  elements=$(basename "$ELEMENTS")"
echo "out_dir=$OUT_DIR"

ACC_ARG=(); [[ "$MODE" != "sequence" ]] && ACC_ARG=(--accessibility-bw "$DNASE_BW")
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
