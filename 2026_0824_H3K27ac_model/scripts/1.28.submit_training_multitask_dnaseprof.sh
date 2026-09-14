#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 24:00:00
#SBATCH --mem=120G
#SBATCH -o log/mtprof_%x.%j.txt
#SBATCH -e log/mtprof_%x.%j.txt
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# MULTI-TASK ARM: the profile head predicts DNase shape, the counts head predicts H3K27ac.
#
# IDENTICAL TO 1.11 except for the profile head's target, its loss weight, and the output
# directory. That is what makes `multimodal5p_accs5p_hw500_clw10` the paired baseline: same
# elements, same negatives, same folds, same ATAC input, same counts target, same window,
# same count_loss_weight.
#
# WHY. The H3K27ac profile head currently trains against a target whose 1 bp
# inter-replicate ceiling is 0.21, and `profile_pearson` comes out at 0.063-0.064 across all
# four accessibility inputs (F-013), so the head learns nothing whatever it is fed. K562
# DNase's 1 bp top-quintile ceiling is 0.848. Hypothesis: a LEARNABLE base-resolution task
# makes the shared trunk better at H3K27ac COUNTS, which is the quantity of interest.
#
# K562 AND NOT GM12878, AND F-014 IS THE REASON. An auxiliary base-resolution task is only
# learnable where the target has reliable base-resolution shape. GM12878's DNase profile
# ceiling is 0.686 on a 5.6x shallower library, and in-cell GM12878 destroying ALL sub-250 bp
# DNase structure costs the model nothing (+0.0011, p=0.72). There the head would fit the
# same noise it already fits.
#
# THIS IS NOT THE CONVERTER. DNase is used as a TARGET at training time; inference needs
# sequence and ATAC only, exactly as before. The converter's failure (F-014) says nothing
# about it, and it is deployment-legal by the same argument that prefers a multi-head model
# over feeding p300 in.
#
# PROFILE_LOSS_WEIGHT IS NOT OPTIONAL. MNLL scales with the target's read depth, and K562
# DNase carries 17.8x the reads of K562 H3K27ac over these 1 kb windows (910.2 against 51.1
# mean total). Unweighted, the profile term would grow ~18x and, since
# `loss = profile + clw * count`, the effective count weight would fall ~18x: a loss on
# counts would then be uninterpretable, because the arm would have been trained to care
# about counts far less than its baseline did. 0.0561 is the measured depth ratio, so the
# profile term STARTS at the baseline's magnitude and `clw=10` keeps its meaning.
#
# Usage: sbatch 1.28.submit_training_multitask_dnaseprof.sh FOLD
# Env:   PROFILE_LOSS_WEIGHT (default 0.0561, the measured depth ratio),
#        COUNT_LOSS_WEIGHT (default 10), HALF_WINDOW (default 500)
set -euo pipefail
export PYTHONUNBUFFERED=1

FOLD=${1:?Usage: sbatch 1.28.submit_training_multitask_dnaseprof.sh FOLD}

PROJECT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
PROJ="$PROJECT_DIR/2026_0824_H3K27ac_model"
PY="$PROJECT_DIR/.pixi/envs/multimodal/bin/python"

N_LAYERS=8
TRIMMING=$(( 47 + 2 + 4 + 8 + 16 + 32 + 64 + 128 + 256 ))
HALF_WINDOW=${HALF_WINDOW:-500}
OUT_WINDOW=$(( 2 * HALF_WINDOW ))
IN_WINDOW=$(( OUT_WINDOW + 2 * TRIMMING ))
COUNT_LOSS_WEIGHT=${COUNT_LOSS_WEIGHT:-10}
PROFILE_LOSS_WEIGHT=${PROFILE_LOSS_WEIGHT:-0.0561}
MAX_NEGATIVES=${MAX_NEGATIVES:-50000}

GENOME="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
ATAC_BW="$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac_5p.bw"
ELEMENTS="$PROJECT_DIR/reference/K562_DNase_candidate_elements.narrowPeak"
NEGATIVES="$PROJECT_DIR/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed"
FOLDS="$PROJECT_DIR/reference/hg38_five_folds.json"

# Pooled stranded DNase, built by 0.34, which asserts plus+minus matches the pooled
# single-track total to within 0.5% (exactly 301,112,388 reads against 0.32's).
PROF_PLUS="$PROJ/data/k562_dnase_5p_plus.bw"
PROF_MINUS="$PROJ/data/k562_dnase_5p_minus.bw"

TAG="dnaseprof$(printf '%s' "$PROFILE_LOSS_WEIGHT" | tr -d '.')"
OUT_DIR="$PROJ/models/multimodal5p_accs5p_${TAG}_hw${HALF_WINDOW}_clw${COUNT_LOSS_WEIGHT}/fold${FOLD}"
mkdir -p "$OUT_DIR" "$PROJ/log"

for f in "$PROF_PLUS" "$PROF_MINUS" "$ATAC_BW"; do
    [[ -s "$f" ]] || { echo "ERROR: '$f' missing or empty" >&2; exit 1; }
done

echo "fold=$FOLD hw=$HALF_WINDOW clw=$COUNT_LOSS_WEIGHT plw=$PROFILE_LOSS_WEIGHT"
echo "counts head  -> h3k27ac_5p_{plus,minus}.bw"
echo "profile head -> $(basename "$PROF_PLUS") / $(basename "$PROF_MINUS")"
echo "out_dir=$OUT_DIR"

cd "$PROJ"
$PY "$PROJECT_DIR/scripts/train_multimodal_bpnet.py" \
    --mode multimodal --peaks "$ELEMENTS" --negatives "$NEGATIVES" \
    --genome "$GENOME" \
    --signal-plus-bw "$PROJ/data/h3k27ac_5p_plus.bw" \
    --signal-minus-bw "$PROJ/data/h3k27ac_5p_minus.bw" \
    --profile-target-plus-bw "$PROF_PLUS" \
    --profile-target-minus-bw "$PROF_MINUS" \
    --profile-loss-weight "$PROFILE_LOSS_WEIGHT" \
    --accessibility-bw "$ATAC_BW" \
    --fold "$FOLDS" --fold-key "$FOLD" --output-dir "$OUT_DIR" \
    --in-window "$IN_WINDOW" --out-window "$OUT_WINDOW" --n-layers "$N_LAYERS" \
    --count-loss-weight "$COUNT_LOSS_WEIGHT" --max-negatives "$MAX_NEGATIVES" \
    --negative-ratio 0.1

echo "Done: $OUT_DIR"
