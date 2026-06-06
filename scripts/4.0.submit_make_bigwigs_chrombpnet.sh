#!/bin/bash
# Generate BigWig files for ChromBPNet mean SHAP scores and mean profile predictions.
#
# Usage:
#   bash 4.0.submit_make_bigwigs_chrombpnet.sh \
#       <shap_dir>        # e.g. K562_ATAC_ChromBPNet/shap
#       <pred_dir>        # e.g. K562_ATAC_ChromBPNet/predictions
#       <chrom_sizes>     # e.g. $OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
#
# Inputs:
#   <shap_dir>/all_folds/counts_mean_shap_scores.h5  (from 3.5.submit_mean_shap_chrombpnet.sh)
#   <shap_dir>/all_folds/counts_peaks_valid_scores.bed
#   <pred_dir>/mean_predictions.h5                   (from 2.5.mean_predictions_chrombpnet.py)
#
# Outputs:
#   <shap_dir>/all_folds/counts_scores.bw
#   <shap_dir>/all_folds/counts_scores.stats.txt
#   <pred_dir>/mean_predictions_profile.bw

if [[ $# -lt 3 ]]; then
    echo "Usage: bash $0 <shap_dir> <pred_dir> <chrom_sizes>"
    exit 1
fi

shap_dir=$1
pred_dir=$2
chr_sizes=$3

PROJECT_DIR=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
BPNET_DIR=/oak/stanford/groups/engreitz/Users/sheth/bpnet-refactor

mkdir -p $PROJECT_DIR/slurm_logs

# --- SHAP BigWig ---
shap_h5=$shap_dir/all_folds/counts_mean_shap_scores.h5
shap_bed=$shap_dir/all_folds/counts_peaks_valid_scores.bed
shap_bw=$shap_dir/all_folds/counts_scores.bw
shap_stats=$shap_dir/all_folds/counts_scores.stats.txt

sbatch -p owners,normal -t 2:00:00 --mem=64G \
    -o "$PROJECT_DIR/slurm_logs/shap_bw.%j.txt" \
    -e "$PROJECT_DIR/slurm_logs/shap_bw.%j.txt" \
    --job-name="shap_bw" \
    --wrap="conda run -n tfmodisco python $PROJECT_DIR/scripts/make_shap_bigwig_chrombpnet.py \
        --shap-h5 $shap_h5 \
        --regions-bed $shap_bed \
        --chrom-sizes $chr_sizes \
        --output-bw $shap_bw \
        --output-stats $shap_stats"

echo "Submitted SHAP BigWig job. Output: $shap_bw"

# --- Profile predictions BigWig ---
pred_h5=$pred_dir/mean_predictions.h5
pred_bw=$pred_dir/mean_predictions_profile.bw

sbatch -p owners,normal -t 1:00:00 --mem=32G \
    -o "$PROJECT_DIR/slurm_logs/pred_bw.%j.txt" \
    -e "$PROJECT_DIR/slurm_logs/pred_bw.%j.txt" \
    --job-name="pred_bw" \
    --wrap="conda run -n tfmodisco python $PROJECT_DIR/scripts/make_predictions_bigwig_chrombpnet.py \
        --mean-h5 $pred_h5 \
        --chrom-sizes $chr_sizes \
        --output-bw $pred_bw"

echo "Submitted profile BigWig job. Output: $pred_bw"
