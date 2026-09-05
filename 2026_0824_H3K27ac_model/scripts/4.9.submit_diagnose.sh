#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 2:00:00
#SBATCH --mem=64G
#SBATCH -c 4
#SBATCH -o log/abcdiag.%j.txt
#SBATCH -e log/abcdiag.%j.txt
#SBATCH --job-name=abc_diag
#
# Why doesn't better H3K27ac help ABC? Runs after the ABC rerun so the predicted arm's
# EnhancerList exists and carries RC-averaged predictions, consistent with everything else.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
cd "$D/2026_0824_H3K27ac_model"
"$D/.pixi/envs/multimodal/bin/python" scripts/4.9.diagnose_abc_gap.py
