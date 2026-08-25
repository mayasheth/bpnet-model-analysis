#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 6:00:00
#SBATCH --mem=64G
#SBATCH -o log/comparison.%j.txt
#SBATCH -e log/comparison.%j.txt
#SBATCH --job-name=k27_compare
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Apples-to-apples stratified comparison of the H3K27ac models against the existing
# p300 multimodal models, all evaluated the same way: same DNase candidate elements,
# element-centered windows, same signal-quintile stratification, log counts.
#
# Caveat to carry with the numbers: the p300 models were TRAINED on ENCSR000EGE peak
# summits but are evaluated here on candidate element centers. That matches the
# existing p300 prediction pipeline, so it is consistent with previously reported
# numbers, but it is not perfectly matched to H3K27ac, which was trained AND
# evaluated element-centered. The p300 numbers may therefore be pessimistic.
#
# Usage: sbatch 2.3.submit_comparison.sh

set -euo pipefail
export PYTHONUNBUFFERED=1

PROJECT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
PROJ="$PROJECT_DIR/2026_0824_H3K27ac_model"
PY="$PROJECT_DIR/.pixi/envs/multimodal/bin/python"

cd "$PROJ"
mkdir -p log results figures

$PY scripts/2.2.evaluate_stratified.py \
    --config-json config/compare_configs.json \
    --out-prefix p300_vs_h3k27ac_ \
    --elements  "$PROJECT_DIR/reference/K562_DNase_candidate_elements.narrowPeak" \
    --genome    "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa" \
    --signal-bw "/oak/stanford/groups/engreitz/Users/sheth/Data/share/IGV/ENCSR000AKP_coverage.bw" \
    --accessibility-bw "$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac.bw" \
    --fold-json "$PROJECT_DIR/reference/hg38_five_folds.json" \
    --outdir results --figdir figures

echo "Comparison complete."
