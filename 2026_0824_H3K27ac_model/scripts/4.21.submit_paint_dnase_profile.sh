#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 12:00:00
#SBATCH --mem=96G
#SBATCH -o log/paintprof.%j.txt
#SBATCH -e log/paintprof.%j.txt
#SBATCH --job-name=paint_prof
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Paint converted DNase genome-wide at base resolution (4.20).
#
# GPU, NOT CPU. Unlike 4.1 this is ~3.1M forward passes rather than 153k, because the track
# has to exist everywhere the downstream model reads -- including the flanking windows and
# the genome-wide negative pool. A chr22-only trial on a login core did not finish inside
# 25 minutes, which sets the CPU estimate for the whole genome at well over a day.
#
# Memory is 96G: one float32 array per chromosome (chr1 is ~1 GB) plus a chunk of sequence
# and accessibility windows.
#
# Usage: sbatch 4.21.submit_paint_dnase_profile.sh ARM MODEL_DIR [CHROMS]
#   ARM        output stem: data/abc_predicted/convdnase_prof_<ARM>.bw
#   MODEL_DIR  relative to the project, or absolute
#   CHROMS     optional comma-separated subset for a trial
# Env: ACC_BW (accessibility the converter reads), SIGNAL_BW (bounds check only)
set -euo pipefail
export PYTHONUNBUFFERED=1
ARM=${1:?usage: sbatch 4.21... ARM MODEL_DIR [CHROMS]}
MODEL_DIR=${2:?}
CHROMS=${3:-}
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
CHR=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
ACC_BW=${ACC_BW:-$D/2026_0529_multimodal_p300_model/data/atac_5p.bw}
SIGNAL_BW=${SIGNAL_BW:-$P/data/k562_dnase_5p_plus.bw}
OUTDIR=$P/data/abc_predicted
mkdir -p "$OUTDIR" "$P/log"
cd "$P"

case "$MODEL_DIR" in /*) MP="$MODEL_DIR" ;; *) MP="$P/$MODEL_DIR" ;; esac
[[ -d "$MP" ]] || { echo "ERROR: no model dir $MP" >&2; exit 1; }

SUFFIX=""; [[ -n "$CHROMS" ]] && SUFFIX="_$(echo "$CHROMS" | tr ',' '-')"
OUT=$OUTDIR/convdnase_prof_${ARM}${SUFFIX}.bw
CH=(); [[ -n "$CHROMS" ]] && CH=(--chroms "$CHROMS")

echo "model=$MP acc=$ACC_BW out=$OUT"
$PY scripts/4.20.paint_dnase_profile.py \
    --model-dir "$MP" --mode multimodal --accessibility-bw "$ACC_BW" \
    --signal-bw "$SIGNAL_BW" --chrom-sizes "$CHR" --out-bw "$OUT" \
    ${CH[@]+"${CH[@]}"}
echo "Done: $OUT"
