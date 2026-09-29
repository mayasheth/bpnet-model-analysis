#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -G 1
#SBATCH -t 16:00:00
#SBATCH --mem=96G
#SBATCH -c 8
#SBATCH -o log/mcteach.%A_%a.txt
#SBATCH -e log/mcteach.%A_%a.txt
#SBATCH --job-name=mc_teach
#SBATCH --array=0-14
#
# Teaching the Test 1 model to USE its TF-activity channel (F-035: with oracle weights the
# channel removed only a third of the GATA under-ranking and moved ABC by +0.004).
#
# ARMS (task // 5), all with oracle weights, so the ceiling is tested before the deployable
# weighting:
#   head      the channel bypasses the trunk and enters only through a 4-weight additive term
#             on the counts (train_multimodal_bpnet.py --motif-head): the same rule in every
#             cell type, no room to use the channel as a cell fingerprint.
#   pair      the channel in the trunk as in 1.37, plus a cross-cell paired loss: the same
#             windows in two training cell types, predicted vs observed log-count difference
#             (--pair-loss-weight). Identical DNA, so only ATAC and the channel can explain it.
#   headpair  both.
#
# DESIGN. Identical to 1.37 (panel, peaks, negatives, folds, windows, count_loss_weight 10,
# epochs, early stopping, no ranking term) except for the arm's change, so
# models/p300_multicell_holdK562_hw500_clw10_motifact_oracle is the paired comparator.
# Pair weight = the count loss weight (10), so the two log-count terms start on one scale.
#
# SMOKE=1 runs fold 0 of ARMS_OVERRIDE (default headpair) for 1 epoch on small sets into
# $SCRATCH, and runs the loss regression gate first.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
W=${WEIGHTING:-oracle}
PANEL=$P/config/p300_panel_motifact_${W}.json
FOLDS=$P/../reference/hg38_five_folds.json
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
cd "$P"; mkdir -p log
[[ -f "$PANEL" ]] || { echo "ERROR: $PANEL missing (make it with 4.46)" >&2; exit 1; }

ARMS=(head pair headpair)
i=${SLURM_ARRAY_TASK_ID:-0}
ARM=${ARMS[$(( i / 5 ))]}
FOLD=$(( i % 5 ))
PAIRW=10
case "$ARM" in
    head)     TEACH=(--motif-head) ;;
    pair)     TEACH=(--pair-loss-weight "$PAIRW") ;;
    headpair) TEACH=(--motif-head --pair-loss-weight "$PAIRW") ;;
esac
EXTRA=(--max-epochs 100 --early-stopping 10 --max-negatives 20000)
OUT_DIR="$P/models/p300_multicell_holdK562_hw500_clw10_motifact_${W}_${ARM}/fold${FOLD}"
if [[ "${SMOKE:-0}" == 1 ]]; then
    ARM=${ARMS_OVERRIDE:-headpair}; FOLD=0
    TEACH=(--motif-head --pair-loss-weight "$PAIRW")
    [[ "$ARM" == head ]] && TEACH=(--motif-head)
    [[ "$ARM" == pair ]] && TEACH=(--pair-loss-weight "$PAIRW")
    EXTRA=(--max-epochs 1 --early-stopping 1 --max-negatives 2000
           --pair-n-train 2000 --pair-n-val 500)
    OUT_DIR="$SCRATCH/smoke_teach_${ARM}/fold${FOLD}"
    $PY "$D/scripts/test_asymmetric_loss.py"
fi
HALF_WINDOW=500; TRIMMING=557
OUT_WINDOW=$(( 2 * HALF_WINDOW )); IN_WINDOW=$(( OUT_WINDOW + 2 * TRIMMING ))
mkdir -p "$OUT_DIR"
echo "task $i: arm $ARM, fold $FOLD, weighting $W, panel $PANEL -> $OUT_DIR"

$PY "$D/scripts/train_multimodal_bpnet.py" \
    --cell-types-json "$PANEL" --holdout-cell K562 \
    --mode multimodal --genome "$GEN" \
    --negatives "$D/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed" \
    --peaks "$D/reference/ENCSR000EGE_peaks_inliers.narrowPeak" \
    --signal-plus-bw "$D/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig" \
    --fold "$FOLDS" --fold-key "$FOLD" --output-dir "$OUT_DIR" \
    --in-window "$IN_WINDOW" --out-window "$OUT_WINDOW" \
    --count-loss-weight 10 "${TEACH[@]}" "${EXTRA[@]}"
echo "Done: $OUT_DIR"
