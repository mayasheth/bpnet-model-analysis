#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 2:00:00
#SBATCH --mem=24G
#SBATCH -c 1
#SBATCH -o log/fp_p300_gate.%j.txt
#SBATCH --job-name=fp_p300_gate
# 2.4 beyond-ATAC p300 gate for the TOBIAS footprint features. Run from the A2 dir.
set -euo pipefail
export PYTHONUNBUFFERED=1
A2=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
cd $A2
$A2/.pixi/envs/tobias/bin/python scripts/2.4.footprint_p300_gate.py "$@"
