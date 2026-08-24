#!/bin/bash
# Mirror direction of the GM12878 transferability analysis: predict K562 candidate
# elements with the GM12878-trained sequence-only BPNet (cross-cell-type).
# Wraps the existing generic scripts/2.2.submit_mean_predict.sh (array 0-4).
# Reuses the K562 v1 predict config as-is: signal/loci/bias are all K562, only
# the model differs.
# Usage: bash scripts/6.3.submit_predict_gm_bpnet_on_k562.sh

set -euo pipefail

OAK=/oak/stanford/groups/engreitz
PROJECT_DIR=$OAK/Users/sheth/EP300_BPNet
THIS_DIR=$PROJECT_DIR/2026_0606_GM12878_transferability
SCRIPTS_DIR=$PROJECT_DIR/scripts
LOG_DIR=$THIS_DIR/log

CHR_SIZES=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
GENOME=$OAK/Users/sheth/hg38_resources/hg38.fa

GM12878_MODEL_PATTERN="$THIS_DIR/GM12878_EP300_BPNet/models/fold_{fold}/ENCSR000DZG_split000"
K562_CONFIG=$PROJECT_DIR/2025_0517_official_EP300_K562_model/config/input_data_predict.json
OUTPUT_DIR=$THIS_DIR/predictions/gm12878_bpnet_on_k562/mean

sbatch \
  --array=0-4 \
  --output=$LOG_DIR/gm_bpnet_on_k562_f%a.%A.txt \
  --error=$LOG_DIR/gm_bpnet_on_k562_f%a.%A.txt \
  --job-name=gm_bpnet_k562 \
  $SCRIPTS_DIR/2.2.submit_mean_predict.sh \
    "$GM12878_MODEL_PATTERN" \
    $CHR_SIZES $GENOME \
    $OUTPUT_DIR \
    $K562_CONFIG
