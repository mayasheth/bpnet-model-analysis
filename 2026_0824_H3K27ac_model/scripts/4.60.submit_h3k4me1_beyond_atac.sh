#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 4:00:00
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH -o log/k4beyond.%j.txt
#SBATCH --job-name=k4beyond
#
# 4.60: H3K4me1 beyond ATAC in the five panel cells, after 0.47 has built the tracks.
set -euo pipefail
export PYTHONUNBUFFERED=1
P=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model
PY=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/.pixi/envs/multimodal/bin/python
cd "$P"
$PY scripts/4.60.h3k4me1_beyond_atac.py
