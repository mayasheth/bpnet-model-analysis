#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 8:00:00
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH -o log/qmap.%j.txt
#SBATCH -e log/qmap.%j.txt
#SBATCH --job-name=qmap
#
# Quantile-map the painted converter track onto GM12878 DNase's value distribution.
# See 0.37's docstring for why the reference is another cell type's DNase and not K562's.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
CHR=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
ARM=${1:-k562}
cd "$P"
$PY scripts/0.37.quantile_map_track.py \
    --in-glob "data/paint_parts/$ARM/${ARM}_*.bw" \
    --reference-bw "$P/data/gm12878_dnase_5p.bw" \
    --out-bw "$P/data/convdnase_qmap_${ARM}.bw" \
    --chrom-sizes "$CHR"
ls -la "$P/data/convdnase_qmap_${ARM}.bw"
echo done
