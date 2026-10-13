#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 4:00:00
#SBATCH --mem=64G
#SBATCH -c 2
#SBATCH -o log/eqtl_withingene_gather.%j.txt
#SBATCH --job-name=eqtl_wg_gather
# 10.2; BLIND (space-separated arms) from the environment, empty = unblinded
set -euo pipefail
export PYTHONUNBUFFERED=1
cd /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
BL=(${BLIND:-})
.pixi/envs/tobias/bin/python scripts/10.2.eqtl_withingene_gather.py ${BL[@]+--blind "${BL[@]}"}
