#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 2:00:00
#SBATCH --mem=32G
#SBATCH -c 2
#SBATCH -o log/cgi_panel.%j.txt
#SBATCH --job-name=cgi_panel
#
# 8.8: CGI distance vs p300 / H3K27ac / H3K4me1 beyond ATAC in the five panel cells (reading fixed in the docstring).
set -euo pipefail
export PYTHONUNBUFFERED=1
A2=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
PY=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/.pixi/envs/multimodal/bin/python
cd "$A2"
$PY scripts/8.8.cgi_distance_activity_panel.py
