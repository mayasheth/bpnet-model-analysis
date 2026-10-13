#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 2:00:00
#SBATCH --mem=48G
#SBATCH -c 2
#SBATCH -o log/eqtl_transfer_cluster.%A_%a.txt
#SBATCH --job-name=eqtl_cluster
# 9.1 for one biosample of $SCRATCH/eqtl_transfer/matches.tsv (array index = row of its unique biosamples)
set -euo pipefail
export PYTHONUNBUFFERED=1
A2=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
cd $A2
B=$(tail -n +2 $SCRATCH/eqtl_transfer/matches.tsv | cut -f1 | awk '!s[$0]++' | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
[[ -n "$B" ]] || { echo "no biosample for task $SLURM_ARRAY_TASK_ID"; exit 1; }
.pixi/envs/tobias/bin/python scripts/9.1.eqtl_transfer_cluster.py --biosample "$B"
