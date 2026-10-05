#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 4:00:00
#SBATCH --mem=96G
#SBATCH -c 4
#SBATCH -o log/peak_tf_kendall.%j.txt
#SBATCH --job-name=featA_kendall
# 5.1 feature A: peak-TF Kendall and its p300 gate. Args passed to 5.1.peak_tf_kendall.py. Run from the A2 dir.
set -euo pipefail
export PYTHONUNBUFFERED=1
A2=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
cd $A2
$A2/.pixi/envs/tobias/bin/python scripts/5.1.peak_tf_kendall.py "$@"
