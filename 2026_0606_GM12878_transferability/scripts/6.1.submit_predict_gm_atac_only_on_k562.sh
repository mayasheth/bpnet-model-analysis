#!/bin/bash
#SBATCH -p owners,engreitz,normal
#SBATCH -t 08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --job-name=gm_atac_only_on_k562
#SBATCH -o /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0606_GM12878_transferability/log/predict_gm_atac_only_on_k562.%j.txt
#SBATCH -e /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0606_GM12878_transferability/log/predict_gm_atac_only_on_k562.%j.txt

# Mirror direction of the GM12878 transferability analysis: predict K562 candidate
# elements with the GM12878-trained ATAC-only BPNet (cross-cell-type).
# Usage: sbatch scripts/6.1.submit_predict_gm_atac_only_on_k562.sh

set -euo pipefail

OAK=/oak/stanford/groups/engreitz
PROJECT_DIR=$OAK/Users/sheth/EP300_BPNet
THIS_DIR=$PROJECT_DIR/2026_0606_GM12878_transferability
MM_SCRIPTS=$PROJECT_DIR/2026_0529_multimodal_p300_model/scripts

ELEMENTS=$PROJECT_DIR/reference/K562_DNase_candidate_elements.narrowPeak
SIGNAL_PLUS_BW=$PROJECT_DIR/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig
SIGNAL_MINUS_BW=$PROJECT_DIR/2025_0703_retrain_p300_model/data/ENCSR000EGE_minus.bigWig
ATAC_BW=$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac.bw
MODEL_DIR=$THIS_DIR/GM12878_ATAC_only_BPNet/models/atac_only
FOLD_JSON=$PROJECT_DIR/reference/hg38_five_folds.json
PEAKS=$PROJECT_DIR/reference/ENCSR000EGE_peaks_inliers.narrowPeak
OUTPUT_DIR=$THIS_DIR/predictions/gm12878_atac_only_on_k562

module load devel pixi/0.53.0

pixi run -e multimodal python "$MM_SCRIPTS/2.1.predict_multimodal.py" \
    --mode atac \
    --elements "$ELEMENTS" \
    --signal-plus-bw "$SIGNAL_PLUS_BW" \
    --signal-minus-bw "$SIGNAL_MINUS_BW" \
    --accessibility-bw "$ATAC_BW" \
    --model-dir "$MODEL_DIR" \
    --fold-json "$FOLD_JSON" \
    --peaks "$PEAKS" \
    --output-dir "$OUTPUT_DIR" \
    --batch-size 512 \
    --device cpu
