#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 2:00:00
#SBATCH --mem=82G
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_MEM:16GB'
#SBATCH --job-name=gata1_cvp_f2
#SBATCH --output=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/K562_GATA1_BPNet/log/cvp_f2.%j.txt
#SBATCH --error=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/K562_GATA1_BPNet/log/cvp_f2.%j.txt

# CV fold 2 resubmission for GATA1 BPNet.
# Fixes vs shared 2.1.submit_cv_predict.sh:
#   1. GPU constraint added (previous jobs ran on CPU-only nodes)
#   2. --threads 1 to avoid multiprocessing deadlock caused by a BigWig read
#      error in one generator worker (process 1 crashes on fold 2 test chroms,
#      hanging the stealer thread indefinitely with --threads 2)

OAK=/oak/stanford/groups/engreitz
SCRIPTS_DIR=$OAK/Users/sheth/EP300_BPNet/scripts
GATA1_DIR=$OAK/Users/sheth/EP300_BPNet/K562_GATA1_BPNet

MODEL_PATH=$GATA1_DIR/models/fold2/model_split000
CHR_SIZES=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
GENOME=$OAK/Users/sheth/hg38_resources/hg38.fa
OUT_DIR=$GATA1_DIR/predictions_cv/fold2
CONFIG=$GATA1_DIR/config/input_data_predict.json

TEST_CHR="chr4 chr11 chr12 chr15 chrY"

echo "$(date +'%Y-%m-%d_%H-%M-%S')"
echo "FOLD 2 (CV resubmission)"
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
    --chroms $TEST_CHR \
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
    --threads 1
