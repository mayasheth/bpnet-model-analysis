#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 12:00:00
#SBATCH --mem=60G
#SBATCH -o log/shap.%j.txt
#SBATCH -e log/shap.%j.txt
#SBATCH --job-name=mm_shap
#SBATCH -n 1
#SBATCH --ntasks 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'

set -euo pipefail

MODALITY=${1:?Usage: sbatch 2.1.submit_shap.sh {atac|dnase} FOLD CHROM}
FOLD=${2:?}
CHROM=${3:?}

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$( cd "$SCRIPT_DIR/../.." && pwd )"

GENOME="/oak/stanford/groups/engreitz/Users/sheth/genome/GRCh38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
REGIONS="$PROJECT_DIR/reference/ENCSR000EGE_peaks_inliers.narrowPeak"

MODEL_DIR="$SCRIPT_DIR/../models/${MODALITY}/fold${FOLD}"
MODEL="$MODEL_DIR/multimodal_bpnet.torch"
ACC_STATS="$MODEL_DIR/acc_normalization_stats.json"
OUTPUT_PREFIX="$SCRIPT_DIR/../shap/${MODALITY}/fold${FOLD}/${CHROM}"

if [[ "$MODALITY" == "atac" ]]; then
    ACC_BW="PLACEHOLDER_ATAC.bigWig"
elif [[ "$MODALITY" == "dnase" ]]; then
    ACC_BW="PLACEHOLDER_DNASE.bigWig"
else
    echo "Error: MODALITY must be atac or dnase"; exit 1
fi

mkdir -p "$(dirname "$OUTPUT_PREFIX")"
module load devel pixi/0.53.0

pixi run -e multimodal python "$PROJECT_DIR/scripts/shap_multimodal_bpnet.py" \
    --model "$MODEL" \
    --regions "$REGIONS" \
    --genome "$GENOME" \
    --accessibility-bw "$ACC_BW" \
    --acc-stats "$ACC_STATS" \
    --output-prefix "$OUTPUT_PREFIX" \
    --chrom "$CHROM" \
    --num-shuffles 20 \
    --batch-size 8 \
    --target 1 \
    --device cuda

echo "SHAP done: $MODALITY fold$FOLD $CHROM"
