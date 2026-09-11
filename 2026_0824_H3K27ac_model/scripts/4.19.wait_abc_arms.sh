#!/bin/bash
# Gate downstream work on the CONSUMED FILES, not on the driver.
#
# The ABC snakemake driver has stalled after its last real rule three times: every child
# COMPLETED, nothing queued, driver idle indefinitely, outputs fine. Waiting on its exit
# status hangs. Waiting on the newest log file is worse and just bit me -- `ls -t
# log/abcrun.*.txt | head -1` matched a PREVIOUS run whose text already said the gate had
# passed, so the wait returned instantly on another run4s output.
#
# So: poll the specific per-arm prediction files, verify each properly (gzip -t plus a row
# count), and only then report ready. Cancel the driver afterwards.
#
# Usage: 4.19.wait_abc_arms.sh RESULTS_DIR ARM [ARM ...]
set -uo pipefail
A=/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction
RES=${1:?usage: 4.19.wait_abc_arms.sh RESULTS_DIR ARM [ARM ...]}; shift
ARMS=("$@")
F=Predictions/EnhancerPredictionsAllPutative.tsv.gz
MIN_ROWS=${MIN_ROWS:-1000000}

while true; do
    ready=0
    for arm in "${ARMS[@]}"; do
        p="$A/$RES/$arm/$F"
        [[ -s "$p" ]] || continue
        gzip -t "$p" 2>/dev/null || continue
        n=$(zcat "$p" | wc -l)
        [[ "$n" -ge "$MIN_ROWS" ]] || continue
        ready=$((ready + 1))
    done
    if [[ "$ready" -eq "${#ARMS[@]}" ]]; then break; fi
    sleep 120
done

echo "ALL_ARMS_READY"
for arm in "${ARMS[@]}"; do
    p="$A/$RES/$arm/$F"
    printf "  %-20s %s rows  %s\n" "$arm" "$(zcat "$p" | wc -l)" "$(du -h "$p" | cut -f1)"
done
