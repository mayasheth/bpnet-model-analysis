#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 8:00:00
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH -o log/mergepaint.%j.txt
#SBATCH -e log/mergepaint.%j.txt
#SBATCH --job-name=merge_paint
#
# Merge the per-chromosome painted parts into one accessibility track.
#
# BATCH, NOT LOGIN NODE. The first attempt was SIGKILLed interactively: the parts total
# ~22 GB and the old merge materialised every interval. 4.22 now streams, but reading 22 GB
# still does not belong on a shared login node.
#
# MIN_VALUE=0.5 IS NOT A TUNED KNOB. A real 5-prime track holds integer read counts, so a
# painted base carrying less than one read is below what the real track can represent. The
# measured consequences on chr22: sparsity goes to 74.5% against real DNase's 71.3%, the
# 1 bp top-quintile shape correlation moves 0.7814 -> 0.7793, and the lowest-quintile
# magnitude inflation falls from 2.73x to 1.54x. Raising it to 1.0 overshoots real DNase's
# sparsity (86.7% zero) and starts costing shape.
#
# chrM is absent by design: it sits in no fold's held-out split, so no leakage-free model
# exists for it. It holds 9 of ~153k candidate elements and 17 negative windows, none of
# which can appear in a validation metric, so --allow-missing is correct here rather than
# permissive.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
CHR=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
ARM=${1:?usage: sbatch 4.24.submit_merge_paint.sh ARM}
MIN_VALUE=${MIN_VALUE:-0.5}
cd "$P"

$PY scripts/4.22.merge_chrom_bigwigs.py \
    --in-glob "data/paint_parts/$ARM/${ARM}_*.bw" \
    --out-bw "data/convdnase_prof_${ARM}.bw" \
    --chrom-sizes "$CHR" --min-value "$MIN_VALUE" --allow-missing
ls -la "data/convdnase_prof_${ARM}.bw"
echo done
