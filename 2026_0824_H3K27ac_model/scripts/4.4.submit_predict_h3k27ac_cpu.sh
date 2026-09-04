#!/bin/bash
#SBATCH -p engreitz,normal,owners
#
# PARTITION: `engreitz` first -- the lab-owned partition, 9 nodes at 24+ cores and 192 GB+
# with a 7-day limit and no GPUs, so it is the right home for CPU work and does not compete
# with the general GPU queues. `normal` and `owners` follow as fallbacks; SLURM starts the
# job wherever a slot frees first.
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -c 16
#SBATCH -o log/abcpredcpu.%j.txt
#SBATCH -e log/abcpredcpu.%j.txt
#SBATCH --job-name=abc_pred_cpu
#
# CPU twin of 4.2. Identical arguments and identical output, no GPU requested.
#
# WHY: this is pure inference on an 8-layer, 64-filter model over a 2,114 bp window, which
# does not need a GPU. Requesting one puts the job behind a depleted GPU fairshare while
# hundreds of CPU nodes sit idle. 4.1 selects its device with
# `torch.cuda.is_available()`, so it falls back to CPU with no code change.
#
# Cost estimate from the chr22 test: ~3,500 regions took a couple of minutes on a login-node
# core, so 153,545 regions across 16 threads should land well inside the 8 h limit.
#
# Do NOT run this concurrently with 4.2 for the same ARM -- both write the same bigwig.
#
# Usage: sbatch 4.4.submit_predict_h3k27ac_cpu.sh ARM MODEL_DIR MODE [CHROMS]

set -euo pipefail
export PYTHONUNBUFFERED=1
ARM=${1:?usage: sbatch 4.4... ARM MODEL_DIR MODE [CHROMS]}
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

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

ACC=()
[[ "$MODE" != "sequence" ]] && ACC=(--accessibility-bw "$D/2026_0529_multimodal_p300_model/data/atac.bw")
CH=()
[[ -n "$CHROMS" ]] && CH=(--chroms "$CHROMS")

$PY scripts/4.1.predict_h3k27ac_for_abc.py \
    --regions "$REG" --model-dir "$P/$MODEL_DIR" --mode "$MODE" \
    ${ACC[@]+"${ACC[@]}"} ${CH[@]+"${CH[@]}"} \
    --chrom-sizes "$SIZES" --out-bw "$OUTDIR/predk27ac_${ARM}.bw"

echo "Done: $OUTDIR/predk27ac_${ARM}.bw"
