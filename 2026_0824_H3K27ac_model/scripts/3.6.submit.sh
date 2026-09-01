#!/bin/bash
#SBATCH -p owners,normal
#SBATCH -t 4:00:00
#SBATCH --mem=32G
#SBATCH -c 2
#SBATCH -o log/coupling_fig.%j.txt
#SBATCH -e log/coupling_fig.%j.txt
#SBATCH --job-name=k27_coupfig
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
cd $D/2026_0824_H3K27ac_model
$D/.pixi/envs/multimodal/bin/python scripts/3.6.plot_coupling_comparison.py
echo ALL_DONE
