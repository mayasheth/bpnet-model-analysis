#!/bin/bash
#SBATCH -p owners,engreitz,normal
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH -o /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0606_GM12878_transferability/log/gm12878_mm_pred.%j.txt
#SBATCH -e /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0606_GM12878_transferability/log/gm12878_mm_pred.%j.txt
#SBATCH --job-name=gm12878_mm_pred

# Predict p300 log counts from GM12878 multimodal BPNet on 154k GM12878 candidate elements.
# Outputs both CV and mean predictions to predictions/gm12878_multimodal_atac/.
# Usage: sbatch scripts/2.1.submit_predict_gm12878_multimodal.sh

set -euo pipefail

OAK=/oak/stanford/groups/engreitz
PROJECT_DIR=$OAK/Users/sheth/EP300_BPNet
THIS_DIR=$PROJECT_DIR/2026_0606_GM12878_transferability
MM_SCRIPTS=$PROJECT_DIR/2026_0529_multimodal_p300_model/scripts

ELEMENTS=$THIS_DIR/reference/GM12878_candidate_elements.narrowPeak
GENOME=$OAK/Users/sheth/hg38_resources/hg38.fa
SIGNAL_PLUS_BW=$THIS_DIR/data/EP300_plus.bw
SIGNAL_MINUS_BW=$THIS_DIR/data/EP300_minus.bw
ATAC_BW=$THIS_DIR/data/atac.bw
MODEL_DIR=$THIS_DIR/GM12878_multimodal_BPNet/models/atac
FOLD_JSON=$PROJECT_DIR/reference/hg38_five_folds.json
PEAKS=$OAK/Users/sheth/Data/ENCODE/GM12878/EP300/ENCFF926AKK.bed.gz
OUTPUT_DIR=$THIS_DIR/predictions/gm12878_multimodal_atac

mkdir -p "$OUTPUT_DIR"

module load devel pixi/0.53.0

pixi run -e multimodal python "$MM_SCRIPTS/2.1.predict_multimodal.py" \
    --elements "$ELEMENTS" \
    --genome "$GENOME" \
    --signal-plus-bw "$SIGNAL_PLUS_BW" \
    --signal-minus-bw "$SIGNAL_MINUS_BW" \
    --accessibility-bw "$ATAC_BW" \
    --model-dir "$MODEL_DIR" \
    --fold-json "$FOLD_JSON" \
    --peaks "$PEAKS" \
    --output-dir "$OUTPUT_DIR" \
    --batch-size 512 \
    --device cpu

echo "Done. Outputs in $OUTPUT_DIR"
