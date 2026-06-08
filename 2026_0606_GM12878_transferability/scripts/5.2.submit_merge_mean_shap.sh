#!/bin/bash
#SBATCH -p owners,engreitz,normal
#SBATCH -t 2:00:00
#SBATCH --mem=64G
#SBATCH -o /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0606_GM12878_transferability/log/merge_mean_shap.%j.txt
#SBATCH -e /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0606_GM12878_transferability/log/merge_mean_shap.%j.txt
#SBATCH --job-name=gm_shap_merge

# Merge per-chromosome SHAP h5s per fold, then average across folds.
# Run after all 5 folds have 24/24 DONE.txt files in shap_peaks/.
# Usage: sbatch scripts/5.2.submit_merge_mean_shap.sh

set -eo pipefail

OAK=/oak/stanford/groups/engreitz
PROJECT_DIR=$OAK/Users/sheth/EP300_BPNet
THIS_DIR=$PROJECT_DIR/2026_0606_GM12878_transferability
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet
SHAP_DIR=$THIS_DIR/shap_peaks

source ~/.bashrc
conda activate bpnet_37
module load cuda/11.1.1 cudnn/8.1.1.33

echo "$(date): Merging per-chromosome SHAP h5s"

for FOLD in 0 1 2 3 4; do
    echo "  Merging fold $FOLD..."
    python $BPNET_DIR/utils/merge_shap_across_chrom.py \
        --input-dir $SHAP_DIR/fold${FOLD} \
        --h5-filename counts_scores.h5 \
        --output-file $SHAP_DIR/fold${FOLD}/shap_counts_merged.h5
done

echo "$(date): Computing mean SHAP across folds"

SHAP_H5_LIST=$(printf "$SHAP_DIR/fold%d/shap_counts_merged.h5," {0..4} | sed 's/,$//')
mkdir -p $SHAP_DIR/all_folds

python $BPNET_DIR/utils/mean_shap_plus_peaks.py \
    --counts_shaps $SHAP_H5_LIST \
    --output_dir $SHAP_DIR/all_folds

echo "$(date): Done. Output: $SHAP_DIR/all_folds/"
