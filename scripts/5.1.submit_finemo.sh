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

set -eo pipefail

function timestamp {
    # Function to get the current time with the new line character
    # removed
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

## input args
out_dir=$1 # .../finemo
shap_counts_h5=$2 # shap/all_folds/counts_mean_shap_scores.h5
peaks_use=$3 # shap/all_folds/counts_peaks_valid_scores.bed
peak_width=$4 # default 1000
modisco_counts_h5=$5 # modisco_results_counts.h5
modisco_width=$6 # 400 (default) or 500


echo  $( timestamp )
echo $1 $2 $3 $4
mkdir -p $out_dir

## additional vars
regions_npz=$out_dir/contribution_regions.npz

## env
source ~/.bashrc
conda activate finemo_gpu

## 1) process input regions
echo $( timestamp ): "PROCESSING INPUT REGIONS"

finemo extract-regions-bpnet-h5 \
        --h5s $shap_counts_h5 \
        --out-path $regions_npz \
        --peaks $peaks_use \
        --region-width $peak_width

## 2) call hits (requires gpu)
echo $( timestamp ): "RUNNING FINEMO CALL-HITS"

finemo call-hits \
        --regions $regions_npz \
        --modisco-h5 $modisco_counts_h5 \
        --out-dir $out_dir \
        --cwm-trim-threshold 0.2 \
        --global-lambda 0.7 \
        --batch-size 2000

## 3) report
echo $( timestamp ): "RUNNING FINEMO REPORT"

finemo report \
        --regions $regions_npz \
        --hits $out_dir \
        --out-dir $out_dir \
        --modisco-h5 $modisco_counts_h5 \
        --modisco-region-width $modisco_width
