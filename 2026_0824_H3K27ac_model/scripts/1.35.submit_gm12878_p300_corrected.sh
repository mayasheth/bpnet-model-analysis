#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -G 1
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH -o log/gmp300fix.%A_%a.txt
#SBATCH -e log/gmp300fix.%A_%a.txt
#SBATCH --job-name=gm_p300_fix
#SBATCH --array=0-4
#
# Retrain the GM12878 p300 model on the CORRECTED target. F-025.
#
# WHAT WAS WRONG. The 2026-06-08 model at
# 2026_0606_GM12878_transferability/GM12878_multimodal_BPNet/models/atac was trained with
# signal_plus_bw = ENCFF960OFK_plus.bw, whose ENCODE output_type is "predicted signal
# profile (plus strand)" from BPNet-model annotation ENCSR038OGP. Its minus strand,
# ENCFF941MGK, was the genuine observed profile. So one strand of the target was another
# model's output.
#
# EXACTLY ONE THING CHANGES HERE. Every hyperparameter, the peaks, the negatives, the
# accessibility track and the folds are copied from the original
# 1.1.submit_training_gm12878_multimodal.sh. Measured against what the original consumed:
#   rebuilt minus vs ENCFF941MGK   r = 1.0000, means 0.0611 / 0.0611  -> unchanged
#   rebuilt plus  vs ENCFF960OFK   r = 0.2734, means 0.0607 / 0.0528  -> this is the fix
# Both strands are taken from the 0.45 rebuild so the two are built consistently, but since
# the minus is identical either way, the only difference from the original run is the plus
# strand. Any change in the result is attributable to the target and to nothing else.
#
# WHAT IT IS FOR. F-021 called this model's transferred arm the best ATAC-input arm on the
# CRISPR benchmark, 0.4771 against a 0.4680 floor, and put its transfer penalty at +0.0458.
# Both numbers were produced by the corrupted model. Retraining and re-running the ABC and
# CRISPR arms says whether either survives.
#
# NOT the same as the multi-cell panel's GM12878 arm: this is a SINGLE-cell-type model, the
# direct replacement for the 2026-06-08 one.
set -euo pipefail
export PYTHONUNBUFFERED=1
FOLD=${SLURM_ARRAY_TASK_ID:-0}
OAK=/oak/stanford/groups/engreitz
D=$OAK/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PANEL=$P/data/p300_panel
cd "$P"; mkdir -p log

GENOME=$OAK/Users/sheth/hg38_resources/hg38.fa
SIGNAL_PLUS_BW=$PANEL/GM12878_ep300_5p_plus.bw      # CORRECTED: observed, from ENCSR000DZG
SIGNAL_MINUS_BW=$PANEL/GM12878_ep300_5p_minus.bw    # identical to the original's minus
ATAC_BW=$D/2026_0606_GM12878_transferability/data/atac.bw
PEAKS=$OAK/Users/sheth/Data/ENCODE/GM12878/EP300/ENCFF926AKK.bed.gz
NEGATIVES=$D/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed
FOLD_JSON=$D/reference/hg38_five_folds.json
OUTPUT_DIR=$P/models/gm12878_p300_corrected/fold${FOLD}
mkdir -p "$OUTPUT_DIR"

for f in "$SIGNAL_PLUS_BW" "$SIGNAL_MINUS_BW" "$ATAC_BW" "$PEAKS" "$NEGATIVES" "$FOLD_JSON"; do
    [[ -s "$f" ]] || { echo "ERROR missing $f" >&2; exit 1; }
done

echo "$(date): GM12878 p300 CORRECTED target, fold $FOLD"
echo "  plus : $SIGNAL_PLUS_BW"
echo "  minus: $SIGNAL_MINUS_BW"
echo "  out  : $OUTPUT_DIR"

"$D/.pixi/envs/multimodal/bin/python" "$D/scripts/train_multimodal_bpnet.py" \
    --peaks "$PEAKS" \
    --negatives "$NEGATIVES" \
    --genome "$GENOME" \
    --signal-plus-bw "$SIGNAL_PLUS_BW" \
    --signal-minus-bw "$SIGNAL_MINUS_BW" \
    --accessibility-bw "$ATAC_BW" \
    --fold "$FOLD_JSON" \
    --fold-key "$FOLD" \
    --output-dir "$OUTPUT_DIR" \
    --n-filters 64 \
    --n-acc-filters 8 \
    --n-layers 8 \
    --count-loss-weight 1.0 \
    --batch-size 64 \
    --max-epochs 100 \
    --early-stopping 10 \
    --lr 1e-3 \
    --negative-ratio 0.1 \
    --max-jitter 50 \
    --device cuda

echo "$(date): done fold $FOLD -> $OUTPUT_DIR"
