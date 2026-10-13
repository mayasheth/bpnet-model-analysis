#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -G 1
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH -o log/mcpanelmean.%A_%a.txt
#SBATCH -e log/mcpanelmean.%A_%a.txt
#SBATCH --job-name=mc_panelmean
#SBATCH --array=0-4
#
# Corgi idea 1 (decisions [2026-10-13]): the multi-cell p300 model, K562 held out, trained as a
# deviation from the panel mean (train_multimodal_bpnet.py --panel-mean). Leave-one-out mean over
# the other training cells as two extra input channels and as the count offset.
#
# DESIGN. Identical to 1.37 WEIGHTING=none, the MSE-only hold-out-K562 model (config/p300_panel.json,
# peaks, negatives, folds, windows, count_loss_weight 10, 100 epochs, early stopping 10, 20,000
# negatives) except for --panel-mean. Paired comparators: models/p300_multicell_holdK562_hw500_clw10
# {,_seed1,_seed2}.
#
# SEED=n makes a replicate (--seed n) in a *_seed<n> directory. SMOKE=1 runs fold 0 for 1 epoch on
# 2,000 negatives into $SCRATCH.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
PANEL=$P/config/p300_panel.json
SEEDARG=(); SFX=""
[[ -n "${SEED:-}" ]] && { SEEDARG=(--seed "$SEED"); SFX="_seed${SEED}"; }
FOLDS=$P/../reference/hg38_five_folds.json
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
cd "$P"; mkdir -p log

FOLD=${SLURM_ARRAY_TASK_ID:-0}
EXTRA=(--max-epochs 100 --early-stopping 10 --max-negatives 20000)
OUT_DIR="$P/models/p300_multicell_panelmean_holdK562_hw500_clw10${SFX}/fold${FOLD}"
if [[ "${SMOKE:-0}" == 1 ]]; then
    EXTRA=(--max-epochs 1 --early-stopping 1 --max-negatives 2000)
    OUT_DIR="$SCRATCH/smoke_panelmean/fold${FOLD}"
fi
if [[ -f "$OUT_DIR/training_complete.json" ]]; then
    echo "already complete: $OUT_DIR"; exit 0
fi
HALF_WINDOW=500; TRIMMING=557
OUT_WINDOW=$(( 2 * HALF_WINDOW )); IN_WINDOW=$(( OUT_WINDOW + 2 * TRIMMING ))
mkdir -p "$OUT_DIR"
echo "fold $FOLD, panel $PANEL, panel mean -> $OUT_DIR"

$PY "$D/scripts/train_multimodal_bpnet.py" \
    --cell-types-json "$PANEL" --holdout-cell K562 --panel-mean \
    --mode multimodal --genome "$GEN" \
    --negatives "$D/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed" \
    --peaks "$D/reference/ENCSR000EGE_peaks_inliers.narrowPeak" \
    --signal-plus-bw "$D/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig" \
    --fold "$FOLDS" --fold-key "$FOLD" --output-dir "$OUT_DIR" \
    --in-window "$IN_WINDOW" --out-window "$OUT_WINDOW" \
    --count-loss-weight 10 "${EXTRA[@]}" ${SEEDARG[@]+"${SEEDARG[@]}"}
echo "Done: $OUT_DIR"
