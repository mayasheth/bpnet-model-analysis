#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -G 1
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH -o log/mctwostage.%A_%a.txt
#SBATCH -e log/mctwostage.%A_%a.txt
#SBATCH --job-name=mc_2stage
#SBATCH --array=0-4
#
# Frozen two-stage Test 1 model (train_multimodal_bpnet.py --base-model): the MSE-only
# multi-cell model (K562 held out, no channel) is FROZEN, and only a small DNA-free network on
# the TF-activity channel + ATAC is trained, adding to its log counts. Within a cell type the
# channel carries nothing the DNA does not, so a model that sees both can learn TF effects from
# the DNA and misapply them in K562; here the part that learns from the channel cannot see DNA.
#
# WEIGHTING=oracle (ceiling) or assoc_noctcf (deployable). The frozen model is fold-matched
# (fold k corrects fold k). Same panel, peaks, negatives, windows, epochs and early stopping as
# 1.37, so models/p300_multicell_holdK562_hw500_clw10_motifact_<W> is the comparator.
#
# SMOKE=1 runs fold 0 for 1 epoch on 2,000 negatives into $SCRATCH.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
W=${WEIGHTING:-oracle}
PANEL=$P/config/p300_panel_motifact_${W}.json
BASE=${BASE_MODEL:-$P/models/p300_multicell_holdK562_hw500_clw10}
FOLDS=$P/../reference/hg38_five_folds.json
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
cd "$P"; mkdir -p log
[[ -f "$PANEL" ]] || { echo "ERROR: $PANEL missing (make it with 4.46)" >&2; exit 1; }

FOLD=${SLURM_ARRAY_TASK_ID:-0}
EXTRA=(--max-epochs 100 --early-stopping 10 --max-negatives 20000)
OUT_DIR="$P/models/p300_multicell_holdK562_hw500_clw10_motifact_${W}_twostage/fold${FOLD}"
if [[ "${SMOKE:-0}" == 1 ]]; then
    EXTRA=(--max-epochs 1 --early-stopping 1 --max-negatives 2000)
    OUT_DIR="$SCRATCH/smoke_twostage_${W}/fold${FOLD}"
fi
HALF_WINDOW=500; TRIMMING=557
OUT_WINDOW=$(( 2 * HALF_WINDOW )); IN_WINDOW=$(( OUT_WINDOW + 2 * TRIMMING ))
mkdir -p "$OUT_DIR"
echo "fold $FOLD, weighting $W, frozen $BASE -> $OUT_DIR"

$PY "$D/scripts/train_multimodal_bpnet.py" \
    --cell-types-json "$PANEL" --holdout-cell K562 \
    --mode multimodal --genome "$GEN" \
    --negatives "$D/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed" \
    --peaks "$D/reference/ENCSR000EGE_peaks_inliers.narrowPeak" \
    --signal-plus-bw "$D/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig" \
    --fold "$FOLDS" --fold-key "$FOLD" --output-dir "$OUT_DIR" \
    --in-window "$IN_WINDOW" --out-window "$OUT_WINDOW" \
    --count-loss-weight 10 --base-model "$BASE" "${EXTRA[@]}"
echo "Done: $OUT_DIR"
