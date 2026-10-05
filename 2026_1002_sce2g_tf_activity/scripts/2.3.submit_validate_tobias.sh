#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 1:00:00
#SBATCH --mem=16G
#SBATCH -c 1
#SBATCH -o log/tobias_validate.%j.txt
#SBATCH --job-name=tobias_validate
# 2.3 ChIP and motif-ranking validation of the TOBIAS footprints. Run from the A2 dir.
set -euo pipefail
export PYTHONUNBUFFERED=1
A2=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
cd $A2
$A2/.pixi/envs/tobias/bin/python scripts/2.3.validate_tobias.py
