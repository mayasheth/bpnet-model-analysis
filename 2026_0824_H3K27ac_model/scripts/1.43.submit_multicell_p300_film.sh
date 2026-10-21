#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -G 1
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH -o log/mcfilm.%A_%a.txt
#SBATCH -e log/mcfilm.%A_%a.txt
#SBATCH --job-name=mc_film
#SBATCH --array=0-4
#
# Bolero idea 2 (2026-10-20): Corgi (1.42, --panel-mean) plus FiLM conditioning of every dilated
# layer on a per-cell vector from the cell's own ATAC (bulk chromVAR 4.63 -> PCs 4.64). K562 held
# out; its row is used only at inference (4.1 --cond-cell K562).
#
# PANEL_TAG names the panel JSON config/p300_panel_<PANEL_TAG>.json ("" = config/p300_panel.json).
# COND_TAG names the table config/cell_cond_<COND_TAG>.tsv; COND_TAG="" trains the matched
# comparator, Corgi on the same panel with no conditioning. Both tags are in the model directory
# name, so arms cannot overwrite each other. SEED=n makes a replicate. SMOKE=1: fold 0, 1 epoch,
# 2,000 negatives, into $SCRATCH.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
PANEL_TAG=${PANEL_TAG:-}
COND_TAG=${COND_TAG:-}
PANEL=$P/config/p300_panel${PANEL_TAG:+_$PANEL_TAG}.json
[[ -f "$PANEL" ]] || { echo "error: no panel $PANEL"; exit 1; }
CONDARG=(); ARM="panelmean"
if [[ -n "$COND_TAG" ]]; then
    COND=$P/config/cell_cond_${COND_TAG}.tsv
    [[ -f "$COND" ]] || { echo "error: no conditioning table $COND"; exit 1; }
    CONDARG=(--cell-cond-tsv "$COND" --cell-cond-dropout "${COND_DROPOUT:-0.1}")
    ARM="panelmean_film_${COND_TAG}_drop${COND_DROPOUT:-0.1}"
fi
SEEDARG=(); SFX=""
[[ -n "${SEED:-}" ]] && { SEEDARG=(--seed "$SEED"); SFX="_seed${SEED}"; }
FOLDS=$P/../reference/hg38_five_folds.json
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
cd "$P"; mkdir -p log

FOLD=${SLURM_ARRAY_TASK_ID:-0}
EXTRA=(--max-epochs 100 --early-stopping 10 --max-negatives 20000)
NAME="p300_multicell_${ARM}${PANEL_TAG:+_panel$PANEL_TAG}_holdK562_hw500_clw10${SFX}"
OUT_DIR="$P/models/$NAME/fold${FOLD}"
if [[ "${SMOKE:-0}" == 1 ]]; then
    EXTRA=(--max-epochs 1 --early-stopping 1 --max-negatives 2000)
    OUT_DIR="$SCRATCH/smoke_film/$NAME/fold${FOLD}"
fi
if [[ -f "$OUT_DIR/training_complete.json" ]]; then
    echo "already complete: $OUT_DIR"; exit 0
fi
HALF_WINDOW=500; TRIMMING=557
OUT_WINDOW=$(( 2 * HALF_WINDOW )); IN_WINDOW=$(( OUT_WINDOW + 2 * TRIMMING ))
mkdir -p "$OUT_DIR"
echo "fold $FOLD, panel $PANEL, conditioning ${COND_TAG:-none} -> $OUT_DIR"

$PY "$D/scripts/train_multimodal_bpnet.py" \
    --cell-types-json "$PANEL" --holdout-cell K562 --panel-mean \
    --mode multimodal --genome "$GEN" \
    --negatives "$D/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed" \
    --peaks "$D/reference/ENCSR000EGE_peaks_inliers.narrowPeak" \
    --signal-plus-bw "$D/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig" \
    --fold "$FOLDS" --fold-key "$FOLD" --output-dir "$OUT_DIR" \
    --in-window "$IN_WINDOW" --out-window "$OUT_WINDOW" \
    --count-loss-weight 10 "${EXTRA[@]}" ${CONDARG[@]+"${CONDARG[@]}"} ${SEEDARG[@]+"${SEEDARG[@]}"}
echo "Done: $OUT_DIR"
