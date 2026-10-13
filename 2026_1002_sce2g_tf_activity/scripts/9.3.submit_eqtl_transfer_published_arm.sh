#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 1:00:00
#SBATCH --mem=16G
#SBATCH -c 2
#SBATCH -o log/eqtl_transfer_published.%A_%a.txt
#SBATCH --job-name=eqtl_published
# 9.3 for one biosample of $SCRATCH/eqtl_transfer/matches.tsv; METHOD and ARM from the environment
set -euo pipefail
export PYTHONUNBUFFERED=1
cd /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
B=$(tail -n +2 $SCRATCH/eqtl_transfer/matches.tsv | cut -f1 | awk '!s[$0]++' | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
[[ -n "$B" ]] || { echo "no biosample for task $SLURM_ARRAY_TASK_ID"; exit 1; }
.pixi/envs/tobias/bin/python scripts/9.3.eqtl_transfer_published_arm.py --biosample "$B" --method "${METHOD:-scE2G_ATAC}" --arm "${ARM:-scatac_published}"
