#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 1:00:00
#SBATCH --mem=64G
#SBATCH -c 4
#SBATCH -o log/qinv.%j.txt
#SBATCH -e log/qinv.%j.txt
#SBATCH --job-name=qnorm_inv
#
# The script MUST live on Oak, not /tmp: /tmp is node-local, so a login-node copy is invisible
# to the compute node. Same trap 4.5 documents for its mamba shim.
set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
cd $D/2026_0824_H3K27ac_model
$D/.pixi/envs/multimodal/bin/python scripts/4.30.qnorm_invariance.py
