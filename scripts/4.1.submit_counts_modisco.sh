#!/bin/bash
#SBATCH -p engreitz,owners,normal
#SBATCH -t 24:00:00
#SBATCH --mem=200G
#SBATCH -o slurm_logs/modisco.%j.txt
#SBATCH -e slurm_logs/modisco.%j.txt
#SBATCH --job-name=modisco

# adapted from: https://github.com/viramalingam/tf-atlas-pipeline/blob/main/anvil/modisco/modisco_pipeline.sh

function timestamp {
    # Function to get the current time with the new line character
    # removed
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

modisco_counts_dir=$1
shap_counts_h5=$2
max_seqlets=$3
trim_size=$4
initial_flank_to_add=$5
final_flank_to_add=$6
motif_reference=$7
modisco_width=$8

# env
source ~/.bashrc
conda activate tfmodisco
module load cuda/11.1.1 cudnn/8.1.1.33

echo  $( timestamp )
echo $1 $2 $3 $4 $5 $6 $7

mkdir -p $modisco_counts_dir

# 1) run modisco
echo $( timestamp ): "RUNNING MODISCO MOTIFS"

modisco motifs \
    --max_seqlets $max_seqlets \
    --h5py $shap_counts_h5 \
    --output $modisco_counts_dir/counts_scores.h5 \
    --window $modisco_width \
    --trim_size $trim_size \
    --initial_flank_to_add $initial_flank_to_add \
    --final_flank_to_add $final_flank_to_add \
    --verbose

# 2) report
echo $( timestamp ): "RUNNING MODISCO REPORT"
modisco report \
    -i $modisco_counts_dir/counts_scores.h5 \
    -o $modisco_counts_dir/ \
    -m $motif_reference
