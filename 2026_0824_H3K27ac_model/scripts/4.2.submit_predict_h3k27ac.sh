#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 8:00:00
#SBATCH --mem=96G
#SBATCH -o log/abcpred.%j.txt
#SBATCH -e log/abcpred.%j.txt
#SBATCH --job-name=abc_pred
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
#
# Predict H3K27ac over ABC's candidate regions and write a bigWig ABC can consume.
# Regions come from the completed July run, so every arm scores an IDENTICAL region set
# and the comparison isolates the activity term.
#
# Usage: sbatch 4.2.submit_predict_h3k27ac.sh ARM MODEL_DIR MODE [CHROMS]
#   ARM        output name, e.g. k562_multimodal
#   MODEL_DIR  relative to the H3K27ac project, e.g. models/multimodal5p_hw500_clw10
#   MODE       sequence | atac | multimodal
#   CHROMS     optional comma-separated subset for a fast test, e.g. chr22
#
# The accessibility input is the ORIGINAL atac.bw, matching what these models trained on
# (results/TARGET_PROVENANCE.md: every result except accs5p_* uses atac.bw). Feeding the
# 5-prime track to a model trained on full-interval coverage would run silently and be wrong.

set -euo pipefail
export PYTHONUNBUFFERED=1
ARM=${1:?usage: sbatch 4.2... ARM MODEL_DIR MODE [CHROMS]}
MODEL_DIR=${2:?}
MODE=${3:?}
CHROMS=${4:-}

D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
A=/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction
PY=$D/.pixi/envs/multimodal/bin/python
REG=$A/results/2026_0721_h3k27ac_counting_comparison/K562_ATAC_only/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed
SIZES=$A/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv
OUTDIR=$P/data/abc_predicted
mkdir -p "$OUTDIR" "$P/log"
cd "$P"

ACC=()
[[ "$MODE" != "sequence" ]] && ACC=(--accessibility-bw "$D/2026_0529_multimodal_p300_model/data/atac.bw")
CH=()
[[ -n "$CHROMS" ]] && CH=(--chroms "$CHROMS")

$PY scripts/4.1.predict_h3k27ac_for_abc.py \
    --regions "$REG" --model-dir "$P/$MODEL_DIR" --mode "$MODE" \
    ${ACC[@]+"${ACC[@]}"} ${CH[@]+"${CH[@]}"} \
    --chrom-sizes "$SIZES" --out-bw "$OUTDIR/predk27ac_${ARM}.bw"
