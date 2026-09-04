#!/bin/bash
#SBATCH -p engreitz,normal,owners
#
# PARTITION: `engreitz` first -- the lab-owned partition, 9 nodes at 24+ cores and 192 GB+
# with a 7-day limit and no GPUs, so it is the right home for CPU work and does not compete
# with the general GPU queues. `normal` and `owners` follow as fallbacks; SLURM starts the
# job wherever a slot frees first.
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH -o log/abcrun.%j.txt
#SBATCH -e log/abcrun.%j.txt
#SBATCH --job-name=abc_predacts
#
# Run ABC over the nine predicted-activity arms.
#
# `owners,normal` -- the OPPOSITE call from the fragment-channel build, deliberately.
# That build is one uncheckpointed pass over three 9 GB BAMs into a temp directory under a
# cleanup trap, so preemption loses everything and it must sit on non-preemptible `normal`.
# This is a Snakemake DAG whose rule outputs persist on disk, so a preemption loses at most
# the single rule in flight and `--rerun-incomplete` cleans up any partial output. Being
# resumable, it should take the faster queue rather than the safer one.
#
# The request is 8 cores / 64 GB / 12 h rather than 16 / 128 / 24: the counting rules stream
# through bedtools rather than loading anything large, and a smaller request backfills far
# more easily against a depleted fairshare.
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
