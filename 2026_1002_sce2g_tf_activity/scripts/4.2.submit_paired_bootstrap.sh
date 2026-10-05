#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 1:00:00
#SBATCH --mem=8G
#SBATCH -c 1
#SBATCH -o log/paired_bootstrap.%j.txt
#SBATCH --job-name=paired_boot
# 4.2 paired bootstrap of CRISPR AUPRC differences. Args passed to 4.2.paired_bootstrap.py. Run from the A2 dir.
set -euo pipefail
export PYTHONUNBUFFERED=1
A2=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
cd $A2
$A2/.pixi/envs/scprinter/bin/python scripts/4.2.paired_bootstrap.py "$@"
