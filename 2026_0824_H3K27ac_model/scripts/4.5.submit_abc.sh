#!/bin/bash
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH --mem=128G
#SBATCH -c 16
#SBATCH -o log/abcrun.%j.txt
#SBATCH -e log/abcrun.%j.txt
#SBATCH --job-name=abc_predacts
#
# Run ABC over the nine predicted-activity arms.
#
# `normal` rather than `owners`: this is a multi-hour uncheckpointed Snakemake DAG, and a
# preemption partway through wastes every completed rule.
#
# Snakemake runs LOCALLY across the allocated cores rather than submitting its own cluster
# jobs, which keeps the whole DAG inside one allocation and avoids a second layer of
# scheduling to debug.
#
# A dry run goes first and MUST NOT list call_macs_peaks, sort_narrowpeaks or
# make_candidate_regions. Those are satisfied by the Peaks directories copied from the
# completed July run; if they appear, the mtime invariant broke and every arm would get its
# own region set, destroying the comparison. The run aborts in that case rather than
# producing nine subtly incomparable answers.

set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
A=/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction
CFG=config/mine/config_predicted_activity.yaml
mkdir -p "$D/2026_0824_H3K27ac_model/log"

source /home/groups/engreitz/Software/anaconda3/etc/profile.d/conda.sh
conda activate final-abc-env
cd "$A"

echo "=== dry run ==="
snakemake --configfile "$CFG" -n -q > /tmp/abc_dryrun.$$ 2>&1 || { cat /tmp/abc_dryrun.$$; exit 1; }
cat /tmp/abc_dryrun.$$
if grep -qE '^\s*(call_macs_peaks|sort_narrowpeaks|make_candidate_regions)\s' /tmp/abc_dryrun.$$; then
    echo "ABORT: region-calling rules are scheduled, so the copied Peaks were not honoured." >&2
    echo "Every arm must share one region set. Re-run 4.3 --copy-peaks and check mtimes." >&2
    exit 1
fi
echo "region calling is correctly skipped"
rm -f /tmp/abc_dryrun.$$

echo
echo "=== run ==="
snakemake --configfile "$CFG" -j "${SLURM_CPUS_PER_TASK:-16}" --rerun-incomplete --nolock
echo "ABC complete"
