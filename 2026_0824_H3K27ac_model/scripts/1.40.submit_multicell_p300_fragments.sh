#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -G 1
#SBATCH -t 12:00:00
#SBATCH --mem=96G
#SBATCH -c 8
#SBATCH -o log/mcfrag.%A_%a.txt
#SBATCH -e log/mcfrag.%A_%a.txt
#SBATCH --job-name=mc_frag
#SBATCH --array=0-4
#
# Fragment-size ATAC channels for the multi-cell p300 model (K562 held out). F-005 found
# fragment-size channels gain in-cell and are null transferred between two cell types for
# H3K27ac; this asks whether training across FOUR libraries lets the model use them in a
# library-independent way, scored downstream.
#
# ARM=flat5p  input = the cell's flat 5' Tn5-insertion track only (the comparator: the
#             fragment channels are 5' insertion counts, the existing p300 models read a
#             coverage track, and the two conventions should not be mixed, see 0.23)
# ARM=frag5p  input = [flat 5', sub <=139, mono 140-329, di 330-620, poly >=621 bp] 5'
#             insertions split by fragment length (0.23); the first convolution can zero the
#             four split channels, so it can reproduce flat5p exactly
# SEED=n      replicate in a *_seed<n> directory. Otherwise identical to 1.37 / 1.34.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
ARM=${ARM:?set ARM=flat5p or frag5p}
PANEL=$P/config/p300_panel_${ARM}.json
FOLDS=$P/../reference/hg38_five_folds.json
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
cd "$P"; mkdir -p log
SEEDARG=(); SFX=""
[[ -n "${SEED:-}" ]] && { SEEDARG=(--seed "$SEED"); SFX="_seed${SEED}"; }
FOLD=${SLURM_ARRAY_TASK_ID:-0}
OUT_DIR="$P/models/p300_multicell_holdK562_hw500_clw10_${ARM}${SFX}/fold${FOLD}"
HALF_WINDOW=500; TRIMMING=557
OUT_WINDOW=$(( 2 * HALF_WINDOW )); IN_WINDOW=$(( OUT_WINDOW + 2 * TRIMMING ))
mkdir -p "$OUT_DIR"
echo "fold $FOLD, arm $ARM, panel $PANEL -> $OUT_DIR"
$PY "$D/scripts/train_multimodal_bpnet.py" \
    --cell-types-json "$PANEL" --holdout-cell K562 \
    --mode multimodal --genome "$GEN" \
    --negatives "$D/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed" \
    --peaks "$D/reference/ENCSR000EGE_peaks_inliers.narrowPeak" \
    --signal-plus-bw "$D/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig" \
    --fold "$FOLDS" --fold-key "$FOLD" --output-dir "$OUT_DIR" \
    --in-window "$IN_WINDOW" --out-window "$OUT_WINDOW" \
    --count-loss-weight 10 --max-epochs 100 --early-stopping 10 --max-negatives 20000 \
    ${SEEDARG[@]+"${SEEDARG[@]}"}
echo "Done: $OUT_DIR"
