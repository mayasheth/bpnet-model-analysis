#!/bin/bash
#SBATCH -p engreitz,owners,normal
#SBATCH -t 8:00:00
#SBATCH --mem=82G
#SBATCH --array=0-4            # launches jobs with SLURM_ARRAY_TASK_ID=0,1,2,3,4
#SBATCH -o log/mean_pred_f.%a.%j.txt     # %a = array index, %j = jobid
#SBATCH -e log/mean_pred_f.%a.%j.txt
#SBATCH --job-name=mp_f%a    # makes job name

FOLD=$SLURM_ARRAY_TASK_ID
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

## define some helpful things
project_dir=$OAK/Users/sheth/EP300_BPNet
results_dir=$OAK/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model

# reference files
chr_sizes=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
genome=$OAK/Users/sheth/hg38_resources/hg38.fa

# output directory
pred_dir=$results_dir/predictions_mean; mkdir -p $pred_dir
out_dir=$pred_dir/fold${FOLD}; mkdir -p $out_dir

# inputs
model_path=$results_dir/models/fold${FOLD}/model_split000
peaks_use=$base_dir/2025_0325_K562_BPNet/data/input_peaks.narrowPeak
pred_data_config=$results_dir/config/input_data_predict.json

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
        --threads 2