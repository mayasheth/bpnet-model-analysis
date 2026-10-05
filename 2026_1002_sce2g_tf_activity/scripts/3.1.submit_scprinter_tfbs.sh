#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 12:00:00
#SBATCH --mem=128G
#SBATCH -c 32
#SBATCH -o log/scprinter_tfbs.%j.txt
#SBATCH --job-name=scprinter_tfbs
# 3.1 PRINT TF binding scores at the TOBIAS sites. Smoke test: sbatch ... --max-sites 2000
set -euo pipefail
export PYTHONUNBUFFERED=1
export SCPRINTER_DATA=$SCRATCH/scprinter_cache CONDA_OVERRIDE_CUDA=12.0
A2=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
cd $A2
$A2/.pixi/envs/scprinter/bin/python scripts/3.1.scprinter_tfbs.py "$@"
