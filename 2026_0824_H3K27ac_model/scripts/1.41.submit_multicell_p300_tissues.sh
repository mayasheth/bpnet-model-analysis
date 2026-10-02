#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -G 1
#SBATCH -t 12:00:00
#SBATCH --mem=96G
#SBATCH -c 8
#SBATCH -o log/mctissue.%A_%a.txt
#SBATCH -e log/mctissue.%A_%a.txt
#SBATCH --job-name=mc_tissue
#SBATCH --array=0-4
#
# The multi-cell p300 model, K562 held out, trained on the five-cell panel PLUS three ENCODE
# tissues (stomach, transverse colon, ovary): seven training cell types instead of four.
#
# WHY. F-032 traced the transfer error to TF activity (GATA under-predicted in K562), and 4.49
# found that no panel cell carries GATA activity; stomach is the only candidate that does. The
# question is whether more, and more TF-diverse, training cell types close that gap.
#
# DESIGN. Identical to 1.37 WEIGHTING=none, the MSE-only hold-out-K562 model (peaks, negatives,
# folds, windows, count_loss_weight 10, 100 epochs, early stopping 10, 20,000 negatives), except
# for the panel: config/p300_panel_tissues.json (1.33 --with-tissues). Its paired comparators
# are models/p300_multicell_holdK562_hw500_clw10{,_seed1,_seed2}.
#
# Memory: the 4-training-cell run peaked at 37 GB, so 7 cells get 96 GB.
#
# SEED=n makes a replicate (--seed n) in a *_seed<n> directory. SMOKE=1 runs fold 0 for 1
# epoch on 2,000 negatives into $SCRATCH.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
PANEL=$P/config/p300_panel_tissues.json
SEEDARG=(); SFX=""
[[ -n "${SEED:-}" ]] && { SEEDARG=(--seed "$SEED"); SFX="_seed${SEED}"; }
FOLDS=$P/../reference/hg38_five_folds.json
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
cd "$P"; mkdir -p log
[[ -f "$PANEL" ]] || { echo "ERROR: $PANEL missing (1.33 --with-tissues)" >&2; exit 1; }

FOLD=${SLURM_ARRAY_TASK_ID:-0}
EXTRA=(--max-epochs 100 --early-stopping 10 --max-negatives 20000)
OUT_DIR="$P/models/p300_multicell_tissues_holdK562_hw500_clw10${SFX}/fold${FOLD}"
if [[ "${SMOKE:-0}" == 1 ]]; then
    EXTRA=(--max-epochs 1 --early-stopping 1 --max-negatives 2000)
    OUT_DIR="$SCRATCH/smoke_mctissue/fold${FOLD}"
fi
if [[ -f "$OUT_DIR/training_complete.json" ]]; then
    echo "already complete: $OUT_DIR"; exit 0
fi
HALF_WINDOW=500; TRIMMING=557
OUT_WINDOW=$(( 2 * HALF_WINDOW )); IN_WINDOW=$(( OUT_WINDOW + 2 * TRIMMING ))
mkdir -p "$OUT_DIR"
echo "fold $FOLD, panel $PANEL -> $OUT_DIR"

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
