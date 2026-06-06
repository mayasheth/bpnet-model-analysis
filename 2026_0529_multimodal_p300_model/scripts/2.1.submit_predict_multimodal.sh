#!/bin/bash
#SBATCH -p owners,engreitz,normal
#SBATCH -t 6:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH -o log/predict_multimodal.%j.txt
#SBATCH -e log/predict_multimodal.%j.txt
#SBATCH --job-name=mm_predict

set -euo pipefail

ACCESSIBILITY=${1:?Usage: sbatch 2.1.submit_predict_multimodal.sh ACCESSIBILITY [atac|dnase]}

PROJECT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
MODEL_DIR="$PROJECT_DIR/2026_0529_multimodal_p300_model/models/${ACCESSIBILITY}"
OUTPUT_DIR="$PROJECT_DIR/2026_0529_multimodal_p300_model/predictions/${ACCESSIBILITY}"

ELEMENTS="$PROJECT_DIR/reference/K562_DNase_candidate_elements.narrowPeak"
GENOME="$PROJECT_DIR/../hg38_resources/hg38.fa"
SIGNAL_PLUS_BW="$PROJECT_DIR/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig"
SIGNAL_MINUS_BW="$PROJECT_DIR/2025_0703_retrain_p300_model/data/ENCSR000EGE_minus.bigWig"
FOLD_JSON="$PROJECT_DIR/reference/hg38_five_folds.json"
PEAKS="$PROJECT_DIR/reference/ENCSR000EGE_peaks_inliers.narrowPeak"

if [[ "$ACCESSIBILITY" == "atac" ]]; then
    ACC_BW="$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac.bw"
elif [[ "$ACCESSIBILITY" == "dnase" ]]; then
    ACC_BW="$PROJECT_DIR/2026_0529_multimodal_p300_model/data/dnase.bw"
else
    echo "ERROR: ACCESSIBILITY must be 'atac' or 'dnase'" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR" "$PROJECT_DIR/2026_0529_multimodal_p300_model/log"

module load devel pixi/0.53.0

pixi run -e multimodal python \
    "$PROJECT_DIR/2026_0529_multimodal_p300_model/scripts/2.1.predict_multimodal.py" \
    --elements "$ELEMENTS" \
    --genome "$GENOME" \
    --signal-plus-bw "$SIGNAL_PLUS_BW" \
    --signal-minus-bw "$SIGNAL_MINUS_BW" \
    --accessibility-bw "$ACC_BW" \
    --model-dir "$MODEL_DIR" \
    --fold-json "$FOLD_JSON" \
    --peaks "$PEAKS" \
    --output-dir "$OUTPUT_DIR" \
    --batch-size 512 \
    --device cpu

pixi run -e multimodal python \
    "$PROJECT_DIR/2026_0529_multimodal_p300_model/scripts/2.2.plot_prediction_accuracy.py" \
    --predictions-dir "$OUTPUT_DIR"

echo "Done. Outputs in $OUTPUT_DIR"
