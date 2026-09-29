#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -G 1
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH -o log/mcrank.%A_%a.txt
#SBATCH -e log/mcrank.%A_%a.txt
#SBATCH --job-name=mc_rank
#SBATCH --array=0-9%8
#
# Multi-cell p300, K562 held out, with a pairwise RANKING term added to the counts loss.
#
# WHY. ABC's rank qnorm keeps only each element's within-class rank of the activity track
# (F-022), so the counts head is judged downstream on ordering alone. log1p MSE also spends
# capacity on scale. This adds a RankNet term over peak pairs within each cell type
# (multimodal_bpnet._pairwise_rank_loss), KEEPING the MSE term so negatives and scale stay
# anchored -- ABC ranks all 153,545 candidate regions, most of which are not EP300 peaks.
#
# DESIGN. Identical to the hold-out-K562 tasks of 1.34 (panel, peaks, negatives, folds,
# windows, count_loss_weight 10, epochs, early stopping) except for the ranking term, so the
# existing models/p300_multicell_holdK562_hw500_clw10 is the paired comparator.
#
# TWO WEIGHTS, a dose-response. After warm-up the weighted count term is about 2 (MSE ~0.2 x
# weight 10) and the ranking term starts at log 2 = 0.69. Weight 2 puts ranking at roughly
# half the count term; weight 10 makes it dominant. index = weight_i * 5 + fold.
#
# NOTE ON EARLY STOPPING. The validation loss includes the ranking term, so the stopping
# epoch will differ from the MSE baseline's. F-019 showed training length alone moves the
# benchmark by +0.0086, so report the epochs each arm stopped at alongside the result.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
PANEL=${PANEL_JSON:-$P/config/p300_panel.json}
FOLDS=$P/../reference/hg38_five_folds.json
GEN=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa
cd "$P"; mkdir -p log

# RANK_WEIGHTS overrides the list, e.g. RANK_WEIGHTS="20 50 100" with --array=0-14.
WEIGHTS=(${RANK_WEIGHTS:-2 10})
i=${SLURM_ARRAY_TASK_ID:-0}
RW=${WEIGHTS[$(( i / 5 ))]}
FOLD=$(( i % 5 ))
HOLD=K562

HALF_WINDOW=500
TRIMMING=557
OUT_WINDOW=$(( 2 * HALF_WINDOW ))
IN_WINDOW=$(( OUT_WINDOW + 2 * TRIMMING ))
CLW=10
# TAG names a variant of the panel (PANEL_JSON), e.g. TAG=motifact_oracle for the ranking
# term plus the Test 1 TF-activity channel; it keeps the output apart from the plain rank model.
OUT_DIR="$P/models/p300_multicell_hold${HOLD}_hw${HALF_WINDOW}_clw${CLW}_rank${RW}${TAG:+_$TAG}/fold${FOLD}"
mkdir -p "$OUT_DIR"

echo "task $i: hold out $HOLD, fold $FOLD, rank weight $RW"
echo "out: $OUT_DIR"

$PY "$D/scripts/train_multimodal_bpnet.py" \
    --cell-types-json "$PANEL" --holdout-cell "$HOLD" \
    --mode multimodal --genome "$GEN" \
    --negatives "$D/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed" \
    --peaks "$D/reference/ENCSR000EGE_peaks_inliers.narrowPeak" \
    --signal-plus-bw "$D/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig" \
    --fold "$FOLDS" --fold-key "$FOLD" --output-dir "$OUT_DIR" \
    --in-window "$IN_WINDOW" --out-window "$OUT_WINDOW" \
    --count-loss-weight "$CLW" --max-negatives 20000 \
    --max-epochs 100 --early-stopping 10 \
    --rank-loss-weight "$RW" --rank-min-gap 0.1

echo "Done: $OUT_DIR"
