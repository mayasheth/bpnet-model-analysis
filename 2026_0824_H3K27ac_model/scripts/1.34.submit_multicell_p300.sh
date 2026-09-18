#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -G 1
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH -o log/mcp300.%A_%a.txt
#SBATCH -e log/mcp300.%A_%a.txt
#SBATCH --job-name=mc_p300
#SBATCH --array=0-24%8
#
# Leave-one-out multi-cell-type p300 training: 5 held-out cell types x 5 chromosome folds.
#
# THE ARRAY INDEX ENCODES BOTH. index = holdout_i * 5 + fold, so task 0-4 holds out the
# first cell type across its five folds, 5-9 the second, and so on. Written out rather than
# nested because a 25-task array retries one (holdout, fold) pair on preemption, and
# `owners` preempts.
#
# CHROMOSOME FOLDS ARE APPLIED WITHIN EVERY CELL TYPE, by the trainer, so a held-out
# chromosome is held out in all of them at once. Holding a chromosome out of only one cell
# type would let the model see it in another and report a leaked validation score.
#
# WHAT EACH MODEL IS FOR. The hold-out-K562 model is the one that goes through ABC and the
# CRISPR benchmark, because that benchmark is K562-only. The other four exist for the
# upstream transfer matrix: each is scored on the cell type it never saw, which is the only
# way to get more than one transfer direction out of a K562-only benchmark.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
PANEL=${PANEL_JSON:-$P/config/p300_panel.json}
FOLDS=$P/../reference/hg38_five_folds.json
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
cd "$P"; mkdir -p log

[[ -s "$PANEL" ]] || { echo "ERROR: no panel json at $PANEL (run 1.33)" >&2; exit 1; }
mapfile -t CELLS < <($PY -c "import json,sys; print('\n'.join(json.load(open('$PANEL'))))")
N=${#CELLS[@]}
i=${SLURM_ARRAY_TASK_ID:-0}
HOLD=${CELLS[$(( i / 5 ))]}
FOLD=$(( i % 5 ))
SAFE=$(echo "$HOLD" | tr '-' '_')

HALF_WINDOW=${HALF_WINDOW:-500}
TRIMMING=${TRIMMING:-557}
OUT_WINDOW=$(( 2 * HALF_WINDOW ))
IN_WINDOW=$(( OUT_WINDOW + 2 * TRIMMING ))
CLW=${COUNT_LOSS_WEIGHT:-10}
OUT_DIR="$P/models/p300_multicell_hold${SAFE}_hw${HALF_WINDOW}_clw${CLW}/fold${FOLD}"
mkdir -p "$OUT_DIR"

echo "task $i of $(( N * 5 )): holding out $HOLD, fold $FOLD"
echo "panel: $PANEL"
echo "out:   $OUT_DIR"

$PY "$D/scripts/train_multimodal_bpnet.py" \
    --cell-types-json "$PANEL" --holdout-cell "$HOLD" \
    --mode multimodal --genome "$GEN" \
    --negatives "$D/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed" \
    --peaks "$D/reference/ENCSR000EGE_peaks_inliers.narrowPeak" \
    --signal-plus-bw "$D/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig" \
    --fold "$FOLDS" --fold-key "$FOLD" --output-dir "$OUT_DIR" \
    --in-window "$IN_WINDOW" --out-window "$OUT_WINDOW" \
    --count-loss-weight "$CLW" --max-negatives "${MAX_NEGATIVES:-20000}" \
    --max-epochs "${MAX_EPOCHS:-100}" --early-stopping "${EARLY_STOPPING:-10}"

echo "Done: $OUT_DIR"
