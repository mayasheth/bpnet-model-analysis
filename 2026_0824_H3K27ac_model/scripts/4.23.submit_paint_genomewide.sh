#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -t 4:00:00
#SBATCH --mem=64G
#SBATCH -o log/paintgw.%A_%a.txt
#SBATCH -e log/paintgw.%A_%a.txt
#SBATCH --job-name=paint_gw
#SBATCH -n 1
#SBATCH -G 1
#SBATCH -C 'GPU_MEM:40GB|GPU_MEM:32GB|GPU_MEM:24GB|GPU_SKU:A100_PCIE|GPU_SKU:A100_SXM4|GPU_SKU:V100_PCIE|GPU_SKU:V100S_PCIE|GPU_SKU:V100_SXM2'
# Array size MUST cover every line of the sizes file. GRCh38.main.chrom.sizes has 25
# contigs (chr1-22, X, Y, M) and an 0-23 array silently dropped one; the guard below
# exits cleanly for out-of-range tasks, so over-sizing is safe and under-sizing is not.
#SBATCH --array=0-24
#
# Genome-wide converted-DNase painting, one chromosome per array task.
#
# WHY AN ARRAY AND NOT ONE JOB. Serial is ~10 h (chr22, 1.6% of the genome, took 10 min at
# stride 250 and data loading dominates, not inference). Chromosomes are independent -- each
# has its own held-out fold and its own accumulator array -- so there is nothing to
# serialise. It is also the preemption-tolerant shape: on `owners` a single 10 h job that
# gets preempted loses everything, while one task losing its chromosome costs 20 minutes.
#
# 4.22 concatenates the per-chromosome outputs afterwards and REFUSES to write a track that
# is missing a chromosome, so a silently-dropped array task cannot become a model trained on
# zero accessibility for a whole chromosome.
#
# Usage: sbatch 4.23.submit_paint_genomewide.sh ARM MODEL_DIR
# Env:   ACC_BW, SIGNAL_BW, STRIDE (default 250)
set -euo pipefail
export PYTHONUNBUFFERED=1
ARM=${1:?usage: sbatch 4.23... ARM MODEL_DIR}
MODEL_DIR=${2:?}
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
CHR=/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GRCh38.main.chrom.sizes
ACC_BW=${ACC_BW:-$D/2026_0529_multimodal_p300_model/data/atac_5p.bw}
SIGNAL_BW=${SIGNAL_BW:-$P/data/k562_dnase_5p_plus.bw}
STRIDE=${STRIDE:-250}
OUTDIR=$P/data/paint_parts/$ARM
mkdir -p "$OUTDIR" "$P/log"
cd "$P"

case "$MODEL_DIR" in /*) MP="$MODEL_DIR" ;; *) MP="$P/$MODEL_DIR" ;; esac

# Chromosome list taken from the sizes file, so it cannot drift from what the track is built
# against. GRCh38.main.chrom.sizes INCLUDES chrM, which is deliberate elsewhere in this
# project and is kept here for the same reason.
mapfile -t CHROMS < <(cut -f1 "$CHR" | sort -u)
N=${#CHROMS[@]}
IDX=${SLURM_ARRAY_TASK_ID:-0}
if (( IDX >= N )); then echo "task $IDX beyond $N chromosomes; nothing to do"; exit 0; fi
C=${CHROMS[$IDX]}

OUT=$OUTDIR/${ARM}_${C}.bw
if [[ -s "$OUT" ]]; then echo "exists: $OUT"; exit 0; fi

echo "task $IDX/$N chrom=$C stride=$STRIDE model=$MP"
$PY scripts/4.20.paint_dnase_profile.py \
    --model-dir "$MP" --mode multimodal --accessibility-bw "$ACC_BW" \
    --signal-bw "$SIGNAL_BW" --chrom-sizes "$CHR" --out-bw "$OUT" \
    --stride "$STRIDE" --chroms "$C"
echo "Done: $OUT"
