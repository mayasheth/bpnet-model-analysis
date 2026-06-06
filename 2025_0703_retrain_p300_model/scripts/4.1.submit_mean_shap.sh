#!/bin/bash
#SBATCH -p owners,gpu
#SBATCH -t 5:00:00  
#SBATCH --mem=64G
#SBATCH --array=0-24
#SBATCH -n 1
#SBATCH --ntasks 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_MEM:16GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:P100_PCIE|GPU_SKU:V100_PCIE|GPU_SKU:TITAN_V|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'

# Usage: sbatch run_shap_split.sh <FOLD>
# Example: sbatch run_shap_split.sh 2

if [ -z "$1" ]; then
  echo "Error: No fold number provided."
  echo "Usage: bash $0 <FOLD>"
  exit 1
fi

FOLD=$1
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
CHR=${CHROMS[$SLURM_ARRAY_TASK_ID]}

base_dir=$OAK/Users/sheth/EP300_BPNet
results_dir=$OAK/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model
bpnet_dir=$OAK/Users/sheth/bpnet-refactor/bpnet

model_path=$results_dir/models/fold${FOLD}/model_split000
out_dir=$results_dir/mean_shap/fold${FOLD}

genome=$OAK/Users/sheth/hg38_resources/hg38.fa
chr_sizes=$OAK/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
peaks_use=$base_dir/2025_0325_K562_BPNet/data/input_peaks.narrowPeak
pred_data_config=$results_dir/config/input_data_predict.json
TASK_ID=0


# File that indicates successful completion
OUTPUT_CHECK=$out_dir/${CHR}/DONE.txt

# Skip if output exists
if [[ -f "$OUTPUT_CHECK" ]]; then
  echo "Output for $CHR in fold $FOLD already exists at $OUTPUT_CHECK. Skipping."
  exit 0
fi

# Optional print for debugging
echo "Running SHAP for fold $FOLD on chromosome $CHR"

mkdir -p $out_dir/$CHR

# Activate environment
source ~/.bashrc
conda activate bpnet_37
module load cuda/11.1.1 cudnn/8.1.1.33

# Call your SHAP script
python $bpnet_dir/cli/shap_split.py \
  --reference-genome $genome \
  --model $model_path \
  --bed-file $peaks_use \
  --chrom $CHR \
  --output-directory $out_dir/$CHR \
  --input-seq-len 2114 \
  --control-len 1000 \
  --task-id $TASK_ID \
  --input-data $pred_data_config \
  --chrom-sizes $chr_sizes \
  --counts-only 