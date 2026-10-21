#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH -o log/motifpair.%A_%a.txt
#SBATCH -e log/motifpair.%A_%a.txt
#SBATCH --job-name=motif_pair
#
# Motif-pair spacing/orientation assay on the Corgi panel-mean p300 model (6.1 python).
#
# THREE STAGES, submitted separately so each asks for what it needs:
#   STAGE=prep     CPU. Chooses the backgrounds and computes the per-cell accessibility
#                  normalisation. Submit with -p engreitz,normal,owners and no GPU.
#   STAGE=run      GPU array over chunks of backgrounds. Submit with -p gpu,owners -G 1
#                  --array=0-<N-1>. Each task writes its own partial .npz and is skipped if
#                  that file already exists, so an `owners` preemption requeue is free.
#   STAGE=combine  CPU. Pools the partials, runs the t-tests and BH-FDR, writes the tables.
#
# The partition and GPU are NOT in the #SBATCH block because they differ by stage; pass them
# on the sbatch command line. Examples:
#   sbatch -p engreitz,normal,owners -t 2:00:00 --export=ALL,STAGE=prep \
#       scripts/6.1.submit_motif_pair_spacing.sh
#   sbatch -p gpu,owners -G 1 --array=0-19 --export=ALL,STAGE=run,CHUNK=25 \
#       scripts/6.1.submit_motif_pair_spacing.sh
#   sbatch -p engreitz,normal,owners -t 0:30:00 --export=ALL,STAGE=combine \
#       scripts/6.1.submit_motif_pair_spacing.sh
#
# SMOKE=1 runs 20 backgrounds with 2 pairs (anchor GATA1 x partners EBOX_TAL1_TCF3, CTCF) into
# $SCRATCH, which is what to run before the full array.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
STAGE=${STAGE:?set STAGE=prep|run|combine}
CHUNK=${CHUNK:-25}
NBG=${NBG:-500}
cd "$P"; mkdir -p log results

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

OUT=$P/results/2026_1020_motif_pair_spacing.tsv
BG=$P/results/2026_1020_motif_pair_backgrounds.tsv
ST=$P/results/2026_1020_motif_pair_accstats.json
PARTS=$P/results/parts_2026_1020_motif_pair_spacing
SMOKE_ARGS=()
if [[ "${SMOKE:-0}" == 1 ]]; then
    NBG=20
    OUT=$SCRATCH/motif_pair_smoke/2026_1020_motif_pair_spacing.tsv
    BG=$SCRATCH/motif_pair_smoke/backgrounds.tsv
    ST=$SCRATCH/motif_pair_smoke/accstats.json
    PARTS=$SCRATCH/motif_pair_smoke/parts
    SMOKE_ARGS=(--anchors GATA1 --partners EBOX_TAL1_TCF3,CTCF --stat-sample 3000)
    mkdir -p "$SCRATCH/motif_pair_smoke"
fi
COMMON=(--bg-tsv "$BG" --stats-json "$ST" --out "$OUT" --n-bg "$NBG")

case "$STAGE" in
prep)
    $PY scripts/6.1.motif_pair_spacing_multimodal.py --stage prep \
        "${COMMON[@]}" ${SMOKE_ARGS[@]+"${SMOKE_ARGS[@]}"}
    ;;
run)
    T=${SLURM_ARRAY_TASK_ID:-0}
    PART=$PARTS/part_$(printf "%03d" "$T").npz
    if [[ -f "$PART" ]]; then
        echo "already done: $PART"; exit 0
    fi
    [[ -f "$BG" ]] || { echo "ERROR: run STAGE=prep first ($BG missing)" >&2; exit 1; }
    mkdir -p "$PARTS"
    $PY scripts/6.1.motif_pair_spacing_multimodal.py --stage run \
        "${COMMON[@]}" ${SMOKE_ARGS[@]+"${SMOKE_ARGS[@]}"} \
        --chunk-start $(( T * CHUNK )) --chunk-size "$CHUNK" --part "$PART"
    ;;
combine)
    $PY scripts/6.1.motif_pair_spacing_multimodal.py --stage combine \
        "${COMMON[@]}" --part-dir "$PARTS"
    ls -la "$OUT"
    ;;
*)
    echo "ERROR: STAGE must be prep, run or combine" >&2; exit 1
    ;;
esac
echo "Done: STAGE=$STAGE"
