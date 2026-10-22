#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH -o log/seamctl.%A_%a.txt
#SBATCH -e log/seamctl.%A_%a.txt
#SBATCH --job-name=seam_ctl
#
# Seam control for the Corgi motif-pair assay (6.2 python; F-068, decisions [2026-10-21]).
#
#   STAGE=run      GPU array, 125 backgrounds (5 gated blocks of 25) per task:
#                  sbatch -p gpu,owners -G 1 --array=0-3 --export=ALL,STAGE=run \
#                      scripts/6.2.submit_motif_pair_seam_control.sh
#   STAGE=combine  CPU:
#                  sbatch -p engreitz,normal,owners -t 0:20:00 --export=ALL,STAGE=combine \
#                      scripts/6.2.submit_motif_pair_seam_control.sh
#   STAGE=smoke    GPU, block 0 into $SCRATCH: the gate must pass against 6.1 block 0 and
#                  must FAIL against block 1; the job exits non-zero otherwise.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
S=$P/scripts/6.2.motif_pair_seam_control.py
STAGE=${STAGE:?set STAGE=run|combine|smoke}
PER=${PER:-125}
cd "$P"; mkdir -p log results
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

TAG=2026_1021_motif_pair_seam_control_corgi_k5_cs1
PARTS=$P/results/parts_$TAG
OUT=$P/results/$TAG.tsv

case "$STAGE" in
run)
    T=${SLURM_ARRAY_TASK_ID:-0}
    $PY "$S" --stage run --part-dir "$PARTS" --start $(( T * PER )) --stop $(( (T + 1) * PER ))
    ;;
combine)
    $PY "$S" --stage combine --part-dir "$PARTS" --out "$OUT"
    ;;
smoke)
    SP=$SCRATCH/seamctl_smoke
    rm -f "$SP"/pass/part_000.npz "$SP"/fail/part_000.npz
    $PY "$S" --stage run --part-dir "$SP/pass" --start 0 --stop 25
    if $PY "$S" --stage run --part-dir "$SP/fail" --start 0 --stop 25 --gate-shift 1 \
            > "$SP/fail.log" 2>&1; then
        echo "ERROR: gate did not fire against the wrong 6.1 block" >&2; exit 1
    fi
    grep "GATE FAILED" "$SP/fail.log" || { cat "$SP/fail.log"; exit 1; }
    [[ ! -f "$SP/fail/part_000.npz" ]] || { echo "ERROR: failed gate wrote a partial" >&2; exit 1; }
    $PY "$S" --stage combine --part-dir "$SP/pass" --out "$SP/smoke.tsv"
    ;;
*)
    echo "ERROR: STAGE must be run, combine or smoke" >&2; exit 1
    ;;
esac
echo "Done: STAGE=$STAGE"
