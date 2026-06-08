#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 4:00:00
#SBATCH --mem=96G
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_MEM:16GB'
#SBATCH --job-name=gata1_mp_f%a
#SBATCH --output=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/K562_GATA1_BPNet/log/mp_f%a.%A.txt
#SBATCH --error=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/K562_GATA1_BPNet/log/mp_f%a.%A.txt

# Mean predictions resubmission for GATA1 BPNet (all 5 folds as array).
# Usage: sbatch --array=0-4 resubmit_mean_predict.sh
#
# Fixes vs shared 2.2.submit_mean_predict.sh:
#   1. GPU constraint added (previous jobs ran on CPU-only nodes)
#   2. chrM excluded from CHROMS — the candidate elements contain
#      chrM:15952-16568 whose 2114bp window (chrM:15203-17317) exceeds
#      the GATA1 BigWig chromosome size (16569 bp), causing a pyBigWig
#      RuntimeError that crashes a generator worker and deadlocks the job
#   3. --threads 2 (fixed) instead of $SLURM_CPUS_ON_NODE to avoid
#      spawning excessive generator processes on multi-CPU GPU nodes

OAK=/oak/stanford/groups/engreitz
GATA1_DIR=$OAK/Users/sheth/EP300_BPNet/K562_GATA1_BPNet

FOLD=$SLURM_ARRAY_TASK_ID
MODEL_PATH=$GATA1_DIR/models/fold${FOLD}/model_split000
CHR_SIZES=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
GENOME=$OAK/Users/sheth/hg38_resources/hg38.fa
OUT_DIR=$GATA1_DIR/predictions_mean/fold${FOLD}
CONFIG=$GATA1_DIR/config/input_data_predict.json

# chrM excluded: its last element (chrM:15952-16568) produces a 2114bp window
# that extends past the end of the GATA1 BigWig, crashing the generator
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
        chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

echo "$(date +'%Y-%m-%d_%H-%M-%S')"
echo "FOLD $FOLD (mean predictions resubmission)"
echo "OUT DIR: $OUT_DIR"
echo "MODEL: $MODEL_PATH"
echo "GPU: $CUDA_VISIBLE_DEVICES"

mkdir -p $OUT_DIR

source ~/.bashrc
conda activate bpnet_37
module load cuda/11.1.1 cudnn/8.1.1.33

bpnet-predict \
    --model $MODEL_PATH \
    --chrom-sizes $CHR_SIZES \
    --chroms $CHROMS \
    --test-indices-file None \
    --reference-genome $GENOME \
    --output-dir $OUT_DIR \
    --input-data $CONFIG \
    --sequence-generator-name BPNet \
    --input-seq-len 2114 \
    --output-len 1000 \
    --output-window-size 1000 \
    --batch-size 64 \
    --reverse-complement-average \
    --threads 2
