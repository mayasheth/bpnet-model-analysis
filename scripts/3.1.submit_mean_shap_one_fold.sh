#!/bin/bash

# Usage: bash 3.1.submit_mean_shap.sh <FOLD> <model_path_pattern> <pred_data_config> <peaks_use> <shap_dir> <genome> <chr_sizes> [bpnet_dir]

# If called with SLURM variables, this is the worker script
if [[ -n "$SLURM_ARRAY_TASK_ID" ]]; then
    # Worker script execution
    CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
            chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
    CHR=${CHROMS[$SLURM_ARRAY_TASK_ID]}
    
    # Check if already done
    OUTPUT_CHECK=$5/fold${1}/${CHR}/DONE.txt
    if [[ -f "$OUTPUT_CHECK" ]]; then
        echo "Output for $CHR already exists. Skipping."
        exit 0
    fi
    
    echo "Running SHAP for fold $1 on chromosome $CHR"
    mkdir -p $5/fold${1}/$CHR
    
    # Activate environment
    source ~/.bashrc
    conda activate bpnet_37
    module load cuda/11.1.1 cudnn/8.1.1.33
    
    # Run SHAP
    model_path="${2//\{fold\}/$1}"
    bpnet_dir=${8:-"$OAK/Users/sheth/bpnet-refactor/bpnet"}
    
    python $bpnet_dir/cli/shap_split.py \
        --reference-genome $6 \
        --model $model_path \
        --bed-file $4 \
        --chrom $CHR \
        --output-directory $5/fold${1}/$CHR \
        --input-seq-len 2114 \
        --control-len 1000 \
        --task-id 0 \
        --input-data $3 \
        --chrom-sizes $7 \
        --counts-only
    
    exit 0
fi

# Submission script logic
if [ $# -lt 7 ]; then
    echo "Usage: bash $0 <FOLD> <model_path_pattern> <pred_data_config> <peaks_use> <shap_dir> <genome> <chr_sizes> [bpnet_dir]"
    exit 1
fi

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

# Find incomplete tasks
incomplete_tasks=()
for i in {0..23}; do
    if [[ ! -f "$5/fold${1}/${CHROMS[$i]}/DONE.txt" ]]; then
        incomplete_tasks+=($i)
    fi
done

if [ ${#incomplete_tasks[@]} -eq 0 ]; then
    echo "All tasks completed."
    exit 0
fi

mkdir -p $5/fold${1}

# Submit long tasks (0,1) and short tasks separately
long_tasks=(); short_tasks=()
for task in "${incomplete_tasks[@]}"; do
    if [[ $task -eq 0 || $task -eq 1 ]]; then
        long_tasks+=($task)
    else
        short_tasks+=($task)
    fi
done

if [ ${#long_tasks[@]} -gt 0 ]; then
    long_list=$(IFS=,; echo "${long_tasks[*]}")
    sbatch -p owners,gpu -t 12:00:00 --mem=64G --array="$long_list" -G 1 \
        -o "slurm_logs/shap_fold${1}.%a.%A.txt" \
        -e "slurm_logs/shap_fold${1}.%a.%A.txt" \
        -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_MEM:16GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:P100_PCIE|GPU_SKU:V100_PCIE|GPU_SKU:TITAN_V|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2' \
        "$0" "$@"
fi

if [ ${#short_tasks[@]} -gt 0 ]; then
    short_list=$(IFS=,; echo "${short_tasks[*]}")
    sbatch -p owners,gpu -t 5:00:00 --mem=64G --array="$short_list" -G 1 \
        -o "slurm_logs/shap_fold${1}.%a.%A.txt" \
        -e "slurm_logs/shap_fold${1}.%a.%A.txt" \
        -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_MEM:16GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:P100_PCIE|GPU_SKU:V100_PCIE|GPU_SKU:TITAN_V|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2' \
        "$0" "$@"
fi

echo "Submitted $(( ${#long_tasks[@]} + ${#short_tasks[@]} )) tasks for fold $1"