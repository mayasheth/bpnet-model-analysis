#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 3:00:00
#SBATCH --mem=32G
#SBATCH -c 1
#SBATCH -o log/footprints.%x.%A_%a.txt
#SBATCH --array=1-23
#
# 1.1 per chromosome. MODE=bias first (all 23), then MODE=score (pools the 23 bias files).
# Usage: MODE=bias sbatch -J fp_bias 1.1...; MODE=score sbatch -J fp_score --dependency=afterok:<bias> 1.1...
set -euo pipefail
export PYTHONUNBUFFERED=1
A2=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
E=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/.pixi/envs/multimodal/bin
export PATH="$E:$PATH"
cd "$A2"
CHR=($(seq -f "chr%g" 1 22) chrX)
C=${CHR[$(( SLURM_ARRAY_TASK_ID - 1 ))]}
$E/python scripts/1.1.footprints.py "${MODE:?set MODE=bias or score}" --chrom "$C"
