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
# Env overrides so the same script serves the p300 models, which live in another project
# directory and were trained against a different target and accessibility track. Duplicating
# this script for p300 would fork the RC-averaging and leakage-assembly logic, which is
# exactly how 2.8/2.11 drifted before 2.15 replaced them.
#   MODEL_ROOT  prefix for a relative MODEL_DIR (default: the H3K27ac project)
#   ACC_BW      accessibility track the model was trained on
#   SIGNAL_BW   any valid bigwig on this genome; used for bounds checking only
#   OUT_PREFIX  output filename stem (default predk27ac)

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

MODEL_ROOT=${MODEL_ROOT:-$P}
ACC_BW=${ACC_BW:-$D/2026_0529_multimodal_p300_model/data/atac.bw}
OUT_PREFIX=${OUT_PREFIX:-predk27ac}
case "$MODEL_DIR" in /*) MODEL_PATH="$MODEL_DIR" ;; *) MODEL_PATH="$MODEL_ROOT/$MODEL_DIR" ;; esac
[[ -d "$MODEL_PATH" ]] || { echo "ERROR: no model dir $MODEL_PATH" >&2; exit 1; }

ACC=()
[[ "$MODE" != "sequence" ]] && ACC=(--accessibility-bw "$ACC_BW")
# ACC_MEAN / ACC_STD: standardize the accessibility input on the PREDICTION regions instead
# of the model's training windows. Needed for a transferred model, whose stored statistics
# came from another cell type's library; get the pair from 4.26. Unset means the old
# behaviour, so no existing arm moves.
NRM=()
if [[ -n "${ACC_MEAN:-}" || -n "${ACC_STD:-}" ]]; then
    [[ -n "${ACC_MEAN:-}" && -n "${ACC_STD:-}" ]] || {
        echo "ERROR: set ACC_MEAN and ACC_STD together" >&2; exit 1; }
    NRM=(--acc-mean "$ACC_MEAN" --acc-std "$ACC_STD")
fi
# OFFSET_MODEL: required for a RESIDUAL model, which predicts observed minus an
# accessibility-only model's counts. Without it the painted track is a residual, not an
# activity estimate. Read it off the model's training_target.json count_offset_model field.
OFF=()
if [[ -n "${OFFSET_MODEL:-}" ]]; then
    case "$OFFSET_MODEL" in /*) OFF_PATH="$OFFSET_MODEL" ;;
                            *)  OFF_PATH="$MODEL_ROOT/models/$OFFSET_MODEL" ;; esac
    [[ -d "$OFF_PATH" ]] || { echo "ERROR: no offset model dir $OFF_PATH" >&2; exit 1; }
    OFF=(--count-offset-model "$OFF_PATH")
fi
SIG=()
[[ -n "${SIGNAL_BW:-}" ]] && SIG=(--signal-bw "$SIGNAL_BW")
CH=()
[[ -n "$CHROMS" ]] && CH=(--chroms "$CHROMS")

echo "model=$MODEL_PATH mode=$MODE acc=$ACC_BW out=${OUT_PREFIX}_${ARM}.bw"
echo "accnorm=${ACC_MEAN:-model}/${ACC_STD:-model} offset=${OFFSET_MODEL:-none}"
$PY scripts/4.1.predict_h3k27ac_for_abc.py \
    --regions "$REG" --model-dir "$MODEL_PATH" --mode "$MODE" \
    ${ACC[@]+"${ACC[@]}"} ${CH[@]+"${CH[@]}"} ${SIG[@]+"${SIG[@]}"} \
    ${NRM[@]+"${NRM[@]}"} ${OFF[@]+"${OFF[@]}"} \
    --chrom-sizes "$SIZES" --out-bw "$OUTDIR/${OUT_PREFIX}_${ARM}.bw"

echo "Done: $OUTDIR/${OUT_PREFIX}_${ARM}.bw"
