#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 6:00:00
#SBATCH --mem=64G
#SBATCH -o log/evaluate.%j.txt
#SBATCH -e log/evaluate.%j.txt
#SBATCH --job-name=k27_eval
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Stratified evaluation of every trained H3K27ac config.
# Must run as a batch job: window extraction over ~30k elements x 5 folds x 6 configs is
# far too heavy for a login node, and it needs a GPU for inference.
#
# Usage: sbatch 2.2.submit_evaluate.sh [CONFIG ...]

set -euo pipefail
export PYTHONUNBUFFERED=1

PROJECT_DIR="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
PROJ="$PROJECT_DIR/2026_0824_H3K27ac_model"
PY="$PROJECT_DIR/.pixi/envs/multimodal/bin/python"

cd "$PROJ"
mkdir -p log results figures

CONFIG_ARG=()
if [[ $# -gt 0 ]]; then
  CONFIG_ARG=(--configs "$@")
fi

$PY scripts/2.2.evaluate_stratified.py \
    --models-dir models \
    --elements  "$PROJECT_DIR/reference/K562_DNase_candidate_elements.narrowPeak" \
    --genome    "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa" \
    --signal-bw "/oak/stanford/groups/engreitz/Users/sheth/Data/share/IGV/ENCSR000AKP_coverage.bw" \
    --accessibility-bw "$PROJECT_DIR/2026_0529_multimodal_p300_model/data/atac.bw" \
    --fold-json "$PROJECT_DIR/reference/hg38_five_folds.json" \
    --outdir results --figdir figures \
    ${CONFIG_ARG[@]+"${CONFIG_ARG[@]}"}

echo "Evaluation complete."
