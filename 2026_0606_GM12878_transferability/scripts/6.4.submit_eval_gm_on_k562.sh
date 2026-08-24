#!/bin/bash
# Evaluate all three GM12878-trained models on K562 candidate elements
# (mirror direction of the GM12878 transferability analysis).
# Run after the 6.1/6.2/6.3 prediction jobs have completed.
# Usage: bash scripts/6.4.submit_eval_gm_on_k562.sh

set -euo pipefail

OAK=/oak/stanford/groups/engreitz
PROJECT_DIR=$OAK/Users/sheth/EP300_BPNet
THIS_DIR=$PROJECT_DIR/2026_0606_GM12878_transferability
SCRIPTS_DIR=$PROJECT_DIR/scripts
MM_SCRIPTS=$PROJECT_DIR/2026_0529_multimodal_p300_model/scripts
LOG_DIR=$THIS_DIR/log
PRED_DIR=$THIS_DIR/predictions

K562_PEAKS=$PROJECT_DIR/reference/ENCSR000EGE_peaks_inliers.narrowPeak

# [1] Sequence-only GM12878 BPNet on K562 elements
sbatch --partition=owners,engreitz,normal --time=2:00:00 --mem=32G \
  --job-name=eval_gm_bpnet_k562 \
  --output=$LOG_DIR/eval_gm_bpnet_on_k562.%j.txt \
  --error=$LOG_DIR/eval_gm_bpnet_on_k562.%j.txt \
  --wrap="module load devel pixi/0.53.0 && pixi run -e ism python \
    $SCRIPTS_DIR/2.3.compute_prediction_performance.py \
      --mean-pred-dir $PRED_DIR/gm12878_bpnet_on_k562/mean \
      --peaks $K562_PEAKS \
      --overlap-col EP300_peak_overlap \
      --h5-name ENCSR000DZG_split000_predictions.h5"

# [2] ATAC-only GM12878 BPNet on K562 elements
sbatch --partition=owners,engreitz,normal --time=2:00:00 --mem=32G \
  --job-name=eval_gm_atac_only_k562 \
  --output=$LOG_DIR/eval_gm_atac_only_on_k562.%j.txt \
  --error=$LOG_DIR/eval_gm_atac_only_on_k562.%j.txt \
  --wrap="module load devel pixi/0.53.0 && pixi run -e ism python \
    $MM_SCRIPTS/2.2.plot_prediction_accuracy.py \
    --predictions-dir $PRED_DIR/gm12878_atac_only_on_k562 \
    --max-counts 10"

# [3] Multimodal (sequence+ATAC) GM12878 BPNet on K562 elements
sbatch --partition=owners,engreitz,normal --time=2:00:00 --mem=32G \
  --job-name=eval_gm_mm_k562 \
  --output=$LOG_DIR/eval_gm_multimodal_on_k562.%j.txt \
  --error=$LOG_DIR/eval_gm_multimodal_on_k562.%j.txt \
  --wrap="module load devel pixi/0.53.0 && pixi run -e ism python \
    $MM_SCRIPTS/2.2.plot_prediction_accuracy.py \
    --predictions-dir $PRED_DIR/gm12878_multimodal_on_k562 \
    --max-counts 10"
