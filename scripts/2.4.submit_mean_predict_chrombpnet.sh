#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 22:00:00
#SBATCH --mem=64G
#SBATCH -o slurm_logs/mean_pred_f.%a.%j.txt
#SBATCH -e slurm_logs/mean_pred_f.%a.%j.txt
#SBATCH --job-name=mean_pred
#SBATCH -n 1
#SBATCH --ntasks 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_MEM:16GB'
#SBATCH --array=0-4

function timestamp {
    # Function to get the current time with the new line character
    # removed
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

### SET UP ###
FOLD=$SLURM_ARRAY_TASK_ID
model_path_pattern=$1
chr_sizes=$2
genome=$3
pred_dir=$4
regions=$5


PROJECT_DIR=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
cd $PROJECT_DIR

# resolve model path
model_path="${model_path_pattern//\{fold\}/$FOLD}"
out_prefix=$pred_dir/fold${FOLD}


# Configuration
mkdir -p slurm_logs
mkdir -p $pred_dir

echo "$(date) Starting mean predict array task ${SLURM_ARRAY_TASK_ID}"
echo "Fold: ${SLURM_ARRAY_TASK_ID}"
echo "Model: $model_path"
echo "Regions: $regions"
echo "Output prefix: $out_prefix"

module load devel pixi/0.53.0
pixi run -e ism python $PROJECT_DIR/scripts/predict_chrombpnet_regions.py \
    -m $model_path \
    -r $regions \
    -g $genome \
    -op $out_prefix


echo "$(date) Completed mean predict array task ${SLURM_ARRAY_TASK_ID}"
