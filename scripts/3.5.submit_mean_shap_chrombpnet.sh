#!/bin/bash
# Compute mean SHAP scores across all folds for ChromBPNet.
# Run this after all per-fold shap_counts_merged.h5 files exist.
#
# Usage: bash 3.5.submit_mean_shap_chrombpnet.sh <shap_dir>

PROJECT_DIR=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet

if [[ $# -lt 1 ]]; then
    echo "Usage: bash $0 <shap_dir>"
    exit 1
fi

shap_dir=$1
out_dir=$shap_dir/all_folds

# Collect per-fold merged h5 files
shap_h5s=$(IFS=,; echo "${shap_dir}/fold0/shap_counts_merged.h5,${shap_dir}/fold1/shap_counts_merged.h5,${shap_dir}/fold2/shap_counts_merged.h5,${shap_dir}/fold3/shap_counts_merged.h5,${shap_dir}/fold4/shap_counts_merged.h5")

sbatch -p owners,normal -t 4:00:00 --mem=200G \
    -o "$PROJECT_DIR/slurm_logs/shap_mean.%j.txt" \
    -e "$PROJECT_DIR/slurm_logs/shap_mean.%j.txt" \
    --job-name="shap_mean_atac" \
    --wrap="module load devel pixi/0.53.0 && \
        pixi run python $PROJECT_DIR/scripts/3.4.mean_shap_chrombpnet.py \
            --shap-h5s \"$shap_h5s\" \
            --output-dir $out_dir"

echo "Submitted mean SHAP job. Output will go to: $out_dir"
