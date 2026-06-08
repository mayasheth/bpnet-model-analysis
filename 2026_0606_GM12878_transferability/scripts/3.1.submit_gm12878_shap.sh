#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_MEM:16GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:P100_PCIE|GPU_SKU:V100_PCIE|GPU_SKU:TITAN_V|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#SBATCH --job-name=gm_shap
#SBATCH --output=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0606_GM12878_transferability/log/shap_fold%a.%A.txt
#SBATCH --error=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0606_GM12878_transferability/log/shap_fold%a.%A.txt

# Usage: sbatch --array=0-4 scripts/3.1.submit_gm12878_shap.sh
# One job per fold; loops through all chromosomes sequentially within the job.
# With 21k peaks, each fold runs in ~2h vs. 25 separate array tasks.

FOLD=$SLURM_ARRAY_TASK_ID

OAK=/oak/stanford/groups/engreitz
PROJECT_DIR=$OAK/Users/sheth/EP300_BPNet
THIS_DIR=$PROJECT_DIR/2026_0606_GM12878_transferability
BPNET_DIR=$OAK/Users/sheth/bpnet-refactor/bpnet

MODEL_DIR=$THIS_DIR/GM12878_EP300_BPNet/models/fold_${FOLD}/ENCSR000DZG_split000
SHAP_DIR=$THIS_DIR/shap_peaks/fold${FOLD}
CONFIG=$THIS_DIR/config/input_data_gm12878_shap.json
PEAKS=$OAK/Users/sheth/Data/ENCODE/GM12878/EP300/ENCFF926AKK.bed.gz
GENOME=$OAK/Users/sheth/hg38_resources/hg38.fa
CHR_SIZES=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
        chr21 chr22 chrX)

source ~/.bashrc
conda activate bpnet_37
module load cuda/11.1.1 cudnn/8.1.1.33

echo "Starting SHAP for fold $FOLD on $(date)"

for CHR in "${CHROMS[@]}"; do
    DONE=$SHAP_DIR/$CHR/DONE.txt
    if [[ -f "$DONE" ]]; then
        echo "  $CHR: already done, skipping"
        continue
    fi
    echo "  Running $CHR..."
    mkdir -p $SHAP_DIR/$CHR
    python $BPNET_DIR/cli/shap_split.py \
        --reference-genome $GENOME \
        --model $MODEL_DIR \
        --bed-file $PEAKS \
        --chrom $CHR \
        --output-directory $SHAP_DIR/$CHR \
        --input-seq-len 2114 \
        --control-len 1000 \
        --task-id 0 \
        --input-data $CONFIG \
        --chrom-sizes $CHR_SIZES \
        --counts-only
done

echo "Finished SHAP for fold $FOLD on $(date)"
