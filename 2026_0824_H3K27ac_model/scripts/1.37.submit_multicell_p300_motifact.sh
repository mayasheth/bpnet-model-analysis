#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -G 1
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH -o log/mcmotif.%A_%a.txt
#SBATCH -e log/mcmotif.%A_%a.txt
#SBATCH --job-name=mc_motif
#SBATCH --array=0-4
#
# Test 1: the multi-cell p300 model, K562 held out, with a TF-ACTIVITY input channel.
#
# Each training cell type's channel is its own motif-activity track (4.45): at each base, the
# sum over FIMO motif hits of that motif's weight in the cell. WEIGHTING=oracle uses each
# cell's true p300 motif importances (4.43 truth, from its own p300), so this is the upper
# bound on what telling the model which TFs are active can buy. A deployable weighting from
# ATAC alone comes later.
#
# DESIGN. Identical to 1.34's hold-out-K562 MSE-only model (panel, peaks, negatives, folds,
# windows, count_loss_weight 10, epochs, early stopping, no ranking term) except for the
# channel, so models/p300_multicell_holdK562_hw500_clw10 is the paired comparator. The ATAC
# channels are normalised exactly as in the comparator (train_multimodal_bpnet.py splits the
# activity channel off before normalisation).
#
# SMOKE=1 runs fold 0 for 1 epoch on 2,000 negatives into a scratch directory.
#
# SEED=n makes a replicate (train_multimodal_bpnet.py --seed n) in a *_seed<n> directory.
# WEIGHTING=none trains the MSE-only comparator itself (config/p300_panel.json, no channel), so
# replicates of the comparator and of the channel models come from one script.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
W=${WEIGHTING:-oracle}
PANEL=$P/config/p300_panel_motifact_${W}.json
[[ "$W" == none ]] && PANEL=$P/config/p300_panel.json
SEEDARG=(); SFX=""
[[ -n "${SEED:-}" ]] && { SEEDARG=(--seed "$SEED"); SFX="_seed${SEED}"; }
FOLDS=$P/../reference/hg38_five_folds.json
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
cd "$P"; mkdir -p log
[[ -f "$PANEL" ]] || { echo "ERROR: $PANEL missing (make it with 4.46)" >&2; exit 1; }

FOLD=${SLURM_ARRAY_TASK_ID:-0}
EXTRA=(--max-epochs 100 --early-stopping 10 --max-negatives 20000)
OUT_DIR="$P/models/p300_multicell_holdK562_hw500_clw10_motifact_${W}${SFX}/fold${FOLD}"
[[ "$W" == none ]] && OUT_DIR="$P/models/p300_multicell_holdK562_hw500_clw10${SFX}/fold${FOLD}"
if [[ "${SMOKE:-0}" == 1 ]]; then
    EXTRA=(--max-epochs 1 --early-stopping 1 --max-negatives 2000)
    OUT_DIR="$SCRATCH/smoke_motifact_${W}/fold${FOLD}"
fi
HALF_WINDOW=500; TRIMMING=557
OUT_WINDOW=$(( 2 * HALF_WINDOW )); IN_WINDOW=$(( OUT_WINDOW + 2 * TRIMMING ))
mkdir -p "$OUT_DIR"
echo "fold $FOLD, weighting $W, panel $PANEL -> $OUT_DIR"

$PY "$D/scripts/train_multimodal_bpnet.py" \
    --cell-types-json "$PANEL" --holdout-cell K562 \
    --mode multimodal --genome "$GEN" \
    --negatives "$D/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed" \
    --peaks "$D/reference/ENCSR000EGE_peaks_inliers.narrowPeak" \
    --signal-plus-bw "$D/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig" \
    --fold "$FOLDS" --fold-key "$FOLD" --output-dir "$OUT_DIR" \
    --in-window "$IN_WINDOW" --out-window "$OUT_WINDOW" \
    --count-loss-weight 10 "${EXTRA[@]}" ${SEEDARG[@]+"${SEEDARG[@]}"}
echo "Done: $OUT_DIR"
