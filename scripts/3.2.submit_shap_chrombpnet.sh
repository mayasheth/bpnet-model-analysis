#!/bin/bash
# Compute DeepSHAP counts-head scores for a ChromBPNet model, one chromosome per job.
#
# Usage: bash 3.2.submit_shap_chrombpnet.sh <FOLD> <model_path_pattern> <regions> <shap_dir> <genome>
#
# Arguments:
#   FOLD                Fold index (0-4)
#   model_path_pattern  Path with {fold} placeholder, e.g. /path/fold_{fold}/model.h5
#   regions             narrowPeak file of regions
#   shap_dir            Output directory (per-fold subdirs will be created)
#   genome              Genome FASTA

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

PROJECT_DIR=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet


# --- Worker (runs when dispatched by SLURM) ---
if [[ -n "$SLURM_ARRAY_TASK_ID" ]]; then
    FOLD=$1
    model_path_pattern=$2
    regions=$3
    shap_dir=$4
    genome=$5

    CHR=${CHROMS[$SLURM_ARRAY_TASK_ID]}
    out_dir=$shap_dir/fold${FOLD}/${CHR}
    out_prefix=$out_dir/shap

    if [[ -f "${out_prefix}.DONE.txt" ]]; then
        echo "Already done: ${CHR}. Skipping."
        exit 0
    fi

    mkdir -p $out_dir

    model_path="${model_path_pattern//\{fold\}/$FOLD}"
    echo "$(date) Starting SHAP fold=${FOLD} chrom=${CHR}"
    echo "Model: $model_path"
    echo "Regions: $regions"
    echo "Output: $out_prefix"

    source ~/.bashrc
    conda activate bpnet_37
    # cuda/11.2.0 required: 11.1.1 lacks libcusolver.so.10 (renamed to .so.11), which TF 2.4 needs
    module load cuda/11.2.0 cudnn/8.1.1.33

    python $PROJECT_DIR/scripts/shap_chrombpnet.py \
        -m $model_path \
        -r $regions \
        -g $genome \
        -o $out_prefix \
        --chrom $CHR

    echo "$(date) Completed SHAP fold=${FOLD} chrom=${CHR}"
    exit 0
fi


# --- Submission logic ---
if [[ $# -lt 5 ]]; then
    echo "Usage: bash $0 <FOLD> <model_path_pattern> <regions> <shap_dir> <genome>"
    exit 1
fi

FOLD=$1
shap_dir=$4

# Find incomplete chromosomes
incomplete_tasks=()
for i in "${!CHROMS[@]}"; do
    done_file=$shap_dir/fold${FOLD}/${CHROMS[$i]}/shap.DONE.txt
    if [[ ! -f "$done_file" ]]; then
        incomplete_tasks+=($i)
    fi
done

if [[ ${#incomplete_tasks[@]} -eq 0 ]]; then
    echo "All chromosomes already complete for fold ${FOLD}."
    exit 0
fi

mkdir -p $shap_dir/fold${FOLD}

# chr1 and chr2 get more time; others are shorter
long_tasks=()
short_tasks=()
for task in "${incomplete_tasks[@]}"; do
    if [[ $task -eq 0 || $task -eq 1 ]]; then
        long_tasks+=($task)
    else
        short_tasks+=($task)
    fi
done

GPU_CONSTRAINT='GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_MEM:16GB'
all_job_ids=()

if [[ ${#long_tasks[@]} -gt 0 ]]; then
    long_list=$(IFS=,; echo "${long_tasks[*]}")
    jid=$(sbatch --parsable -p owners,gpu -t 12:00:00 --mem=64G -G 1 \
        --array="$long_list" \
        -C "$GPU_CONSTRAINT" \
        -o "$PROJECT_DIR/slurm_logs/shap_atac_fold${FOLD}.%a.%A.txt" \
        -e "$PROJECT_DIR/slurm_logs/shap_atac_fold${FOLD}.%a.%A.txt" \
        --job-name="shap_atac_f${FOLD}" \
        "$0" "$@")
    all_job_ids+=($jid)
fi

if [[ ${#short_tasks[@]} -gt 0 ]]; then
    short_list=$(IFS=,; echo "${short_tasks[*]}")
    jid=$(sbatch --parsable -p owners,gpu -t 5:00:00 --mem=64G -G 1 \
        --array="$short_list" \
        -C "$GPU_CONSTRAINT" \
        -o "$PROJECT_DIR/slurm_logs/shap_atac_fold${FOLD}.%a.%A.txt" \
        -e "$PROJECT_DIR/slurm_logs/shap_atac_fold${FOLD}.%a.%A.txt" \
        --job-name="shap_atac_f${FOLD}" \
        "$0" "$@")
    all_job_ids+=($jid)
fi

echo "Submitted ${#incomplete_tasks[@]} chromosome jobs for fold ${FOLD}."

# Submit merge job to run after all chromosome jobs complete
dep=$(IFS=:; echo "afterok:${all_job_ids[*]}")
merge_out=$shap_dir/fold${FOLD}/shap_counts_merged.h5
sbatch -p owners,normal -t 2:00:00 --mem=64G \
    --dependency="$dep" \
    -o "$PROJECT_DIR/slurm_logs/shap_merge_fold${FOLD}.%j.txt" \
    -e "$PROJECT_DIR/slurm_logs/shap_merge_fold${FOLD}.%j.txt" \
    --job-name="shap_merge_f${FOLD}" \
    --wrap="conda run -n bpnet_37 python $PROJECT_DIR/scripts/3.3.merge_shap_chrombpnet.py \
            --input-dir $shap_dir/fold${FOLD} \
            --output-file $merge_out"

echo "Submitted merge job for fold ${FOLD} (runs after chromosome jobs)."
