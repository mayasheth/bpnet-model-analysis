#!/bin/bash
# Submit the full H3K27ac grid: 3 modes x 2 windows x 5 folds = 30 jobs.
#
# Note --signal-minus-bw is never passed, so every model trains on an unstranded
# target (n_outputs=1).
#
# Usage:
#   bash 1.2.submit_all.sh              # submit everything
#   bash 1.2.submit_all.sh --dry-run    # print without submitting
#   MODES="sequence multimodal" bash 1.2.submit_all.sh
#   HALF_WINDOWS=500 bash 1.2.submit_all.sh

set -euo pipefail

DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODES=${MODES:-"sequence multimodal atac"}
HALF_WINDOWS=${HALF_WINDOWS:-"500 1000"}
FOLDS=${FOLDS:-"0 1 2 3 4"}

n=0
for hw in $HALF_WINDOWS; do
  for mode in $MODES; do
    for fold in $FOLDS; do
      name="k27_${mode}_hw${hw}_f${fold}"
      cmd=(sbatch --job-name="$name" "$SCRIPT_DIR/1.1.submit_training.sh" "$mode" "$fold" "$hw")
      if [[ $DRY_RUN -eq 1 ]]; then
        echo "${cmd[*]}"
      else
        "${cmd[@]}"
      fi
      n=$((n + 1))
    done
  done
done
echo "$([[ $DRY_RUN -eq 1 ]] && echo Would submit || echo Submitted) $n jobs"
