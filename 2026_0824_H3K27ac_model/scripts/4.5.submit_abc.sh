#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 48:00:00
#SBATCH --mem=8G
#SBATCH -c 2
#SBATCH -o log/abcrun.%j.txt
#SBATCH -e log/abcrun.%j.txt
#SBATCH --job-name=abc_driver
#
# Run ABC over the nine predicted-activity arms.
#
# THIS JOB IS ONLY THE SNAKEMAKE DRIVER. With --profile, every rule instance is submitted as
# its own SLURM job (profile: 50 concurrent, partitions engreitz,owners,normal, 3 retries),
# so the driver itself needs almost nothing -- 2 cores and 8 GB -- but must outlive the whole
# DAG, hence 48 h. It sits on a non-preemptible partition: if the driver dies, running
# children finish but nothing further is submitted.
#
# --use-conda makes each rule run in the environment declared in workflow/envs/abcenv.yml.
# Snakemake resolves that against .snakemake/conda in the working directory, where a built
# env already exists, so this reuses it rather than solving a new one.
#
# A dry run goes first and MUST NOT list call_macs_peaks, sort_narrowpeaks or
# make_candidate_regions. Those are satisfied by the Peaks directories copied from the
# completed July run; if they appear, the mtime invariant broke and every arm would get its
# own region set, destroying the comparison the experiment exists to make.

set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
A=/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction
SM_ENV=/oak/stanford/groups/engreitz/Users/sheth/.conda/envs/run_snakemake
SM=$SM_ENV/bin/snakemake
PROFILE="$HOME/.config/snakemake/slurm"
# --use-conda shells out to mamba from /usr/bin/bash, which does not inherit the wrapper
# env just because snakemake was called by absolute path. Without this the run dies with
# CreateCondaEnvironmentException before submitting anything.
export PATH="$SM_ENV/bin:$PATH"
CFG=config/mine/config_predicted_activity.yaml
mkdir -p "$D/2026_0824_H3K27ac_model/log"
cd "$A"
mkdir -p .snakemake/slurm_logs

echo "snakemake: $($SM --version)"
echo "profile:   $PROFILE"
echo

echo "=== dry run ==="
$SM --configfile "$CFG" --use-conda -n -q > /tmp/abc_dryrun.$$ 2>&1 || { cat /tmp/abc_dryrun.$$; exit 1; }
cat /tmp/abc_dryrun.$$
if grep -qE '^[[:space:]]*(call_macs_peaks|sort_narrowpeaks|make_candidate_regions)[[:space:]]' /tmp/abc_dryrun.$$; then
    echo "ABORT: region-calling rules are scheduled, so the copied Peaks were not honoured." >&2
    echo "Every arm must share one region set. Re-run 4.3 --copy-peaks and check mtimes." >&2
    exit 1
fi
echo "region calling is correctly skipped"
rm -f /tmp/abc_dryrun.$$

echo
echo "=== run ==="
$SM --configfile "$CFG" --profile "$PROFILE" --use-conda
echo "ABC complete"
