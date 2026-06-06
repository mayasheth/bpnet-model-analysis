#!/bin/bash
#SBATCH -p engreitz,owners,normal
#SBATCH -t 24:00:00
#SBATCH --mem=96G
#SBATCH -c 8
#SBATCH --array=0-4            # launches jobs with SLURM_ARRAY_TASK_ID=0,1,2,3,4
#SBATCH -o log/mean_pred_f.%a.%j.txt     # %a = array index, %j = jobid
#SBATCH -e log/mean_pred_f.%a.%j.txt
#SBATCH --job-name=mp_f%a    # makes job name

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
pred_data_config=$5

# other variables
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

# resolve model path
model_path="${model_path_pattern//\{fold\}/$FOLD}"
out_dir=$pred_dir/fold${FOLD}

# record arguments...
echo  $( timestamp )
echo "FOLD $FOLD"
echo "OUT DIR: $out_dir"
echo "RESOLVED MODEL PATH: $model_path"
echo "PRED INPUT CONFIG: $pred_data_config"
echo "GENOME: $genome"
echo "CHR SIZES: $chr_sizes"

# make output directory
mkdir -p $pred_dir
mkdir -p $out_dir

# env
source ~/.bashrc
conda activate bpnet_37
module load cuda/11.1.1 cudnn/8.1.1.33

# run predict
bpnet-predict \
        --model $model_path \
        --chrom-sizes $chr_sizes \
        --chroms $CHROMS \
        --test-indices-file None \
        --reference-genome $genome \
        --output-dir $out_dir \
        --input-data $pred_data_config \
        --sequence-generator-name BPNet \
        --input-seq-len 2114 \
        --output-len 1000 \
        --output-window-size 1000 \
        --batch-size 64 \
        --reverse-complement-average \
        --threads $SLURM_CPUS_ON_NODE