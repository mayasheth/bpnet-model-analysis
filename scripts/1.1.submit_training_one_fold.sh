#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 18:00:00  
#SBATCH --mem=160G
#SBATCH -o log/train.%j.txt
#SBATCH -e log/train.%j.txt
#SBATCH --job-name=train
#SBATCH -n 1
#SBATCH --ntasks 1
#SBATCH -G 4
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_MEM:16GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:P100_PCIE|GPU_SKU:V100_PCIE|GPU_SKU:TITAN_V|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'

function timestamp {
    # Function to get the current time with the new line character
    # removed
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

### SET UP ###
# load inputs 
FOLD=$1
model_out_dir=$2 # res/models
splits=$3
input_config=$4 # config/input_data.json
params_config=$5 # config/bpnet_params.json
genome=$6 # hg38.fa
chr_sizes=$7

echo  $( timestamp )
echo "FOLD $1"
echo "OUT DIR: $2"
echo "SPLITS: $3"
echo "INPUT CONFIG: $4"
echo "PARAMS CONFIG: $5"
echo "GENOME: $6"
echo "CHR SIZES: $7"

# other vars
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

# conda env
source ~/.bashrc  # ensure Conda is available in non-interactive shells
conda activate bpnet_37
module load cuda/11.1.1 cudnn/8.1.1.33

## TRAIN MODEL
mkdir -p $model_out_dir 
mkdir -p $model_out_dir/fold${FOLD}

bpnet-train \
    --input-data $input_config \
    --output-dir $model_out_dir/fold${FOLD} \
    --reference-genome $genome \
    --chroms $CHROMS \
    --chrom-sizes $chr_sizes \
    --splits $splits \
    --model-arch-name BPNet \
    --model-arch-params-json $params_config \
    --sequence-generator-name BPNet \
    --model-output-filename model \
    --input-seq-len 2114 \
    --output-len 1000 \
    --shuffle \
    --threads 10 \
    --epochs 100 \
    --batch-size 64 \
    --reverse-complement-augmentation \
    --early-stopping-patience 10 \
    --reduce-lr-on-plateau-patience 5 \
    --learning-rate 0.001