#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 4:00:00
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH -o log/smoothctl.%j.txt
#SBATCH -e log/smoothctl.%j.txt
#SBATCH --job-name=smooth_ctl
#
# Build the shape-destroyed DNase control for both cell types. See 0.36's docstring for what
# it is for and why the width is 250 bp.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
CHR=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
W=${W:-250}
cd "$P"

for CELL in k562 gm12878; do
    OUT=$P/data/${CELL}_dnase_5p_smooth${W}.bw
    if [[ -s "$OUT" ]]; then echo "exists: $(basename "$OUT")"; continue; fi
    echo "=== $CELL, width ${W} bp ==="
    $PY scripts/0.36.smooth_dnase_control.py \
        --in-bw "$P/data/${CELL}_dnase_5p.bw" --out-bw "$OUT" \
        --chrom-sizes "$CHR" -w "$W"
    echo
done
echo done
