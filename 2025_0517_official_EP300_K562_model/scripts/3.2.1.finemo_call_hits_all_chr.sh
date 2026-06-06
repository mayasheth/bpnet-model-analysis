#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 6:00:00  
#SBATCH --mem=64G
#SBATCH -o log/finemo.%j.txt
#SBATCH -e log/finemo.%j.txt
#SBATCH --job-name=finemo
#SBATCH -n 1
#SBATCH --ntasks 1
#SBATCH -G 2
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_MEM:16GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:P100_PCIE|GPU_SKU:V100_PCIE|GPU_SKU:TITAN_V|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

base_dir=$OAK/Users/sheth/EP300_BPNet

this_dir=$base_dir/2025_0517_official_EP300_K562_model
out_dir=$this_dir/finemo_all_chr; mkdir -p $out_dir

modisco_counts_h5=$base_dir/2025_0407_MoDISCo/ENCSR000EGE_trim20_flank5_0/modisco_results_counts.h5
counts_shap_h5=$this_dir/shap_all_chr/all_folds/counts_mean_shap_scores.h5
peaks_use=$this_dir/shap_all_chr/all_folds/counts_peaks_valid_scores.bed

regions_npz=$out_dir/contribution_regions.npz

# env
source ~/.bashrc
conda activate finemo_gpu

## process input regions
finemo extract-regions-bpnet-h5 \
        --h5s $counts_shap_h5 \
        --out-path $regions_npz \
        --peaks $peaks_use \
        --region-width 500

## call hits (required gpu)
finemo call-hits \
        --regions $regions_npz \
        --modisco-h5 $modisco_counts_h5 \
        --out-dir $out_dir \
        --cwm-trim-threshold 0.2 \
        --global-lambda 0.7 \
        --batch-size 2000

## report
finemo report \
        --regions $regions_npz \
        --hits $out_dir \
        --out-dir $out_dir \
        --modisco-region-width 400
