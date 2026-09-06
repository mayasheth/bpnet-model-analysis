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
# --use-conda shells out to mamba from /usr/bin/bash, which does not inherit the wrapper env
# just because snakemake was called by absolute path; without mamba the run dies with
# CreateCondaEnvironmentException before submitting anything.
#
# But putting $SM_ENV/bin on PATH is NOT the fix: every rule job inherits the driver's
# environment, and the wrapper env's python (3.11, no pyranges) then shadows the rule conda
# env's python (3.10, with pyranges), so every create_neighborhoods job dies on
# `import pyranges`. Expose ONLY mamba and conda, through a shim directory appended to the
# end of PATH, so nothing else can be shadowed whatever order conda activation applies.
# The shim must live on a SHARED filesystem, not mktemp -d. With the SLURM executor each
# rule runs its own `snakemake --jobstep`, which re-checks for mamba on whatever node it
# lands on; a driver-local /tmp path is invisible there (and deleted when the driver exits),
# so the children failed with the same CreateCondaEnvironmentException the driver had.
SHIM=$D/2026_0824_H3K27ac_model/.mamba_shim
mkdir -p "$SHIM"
ln -sfn "$SM_ENV/bin/mamba" "$SHIM/mamba"
ln -sfn "$SM_ENV/bin/conda" "$SHIM/conda"
export PATH="$PATH:$SHIM"
command -v mamba >/dev/null || { echo "ERROR: mamba still not on PATH via $SHIM" >&2; exit 1; }
CFG=${ABC_CFG:-config/mine/config_predicted_activity.yaml}
mkdir -p "$D/2026_0824_H3K27ac_model/log"
# Re-stamp Peaks before every run rather than trusting whoever ran 4.3 last. If the
# predicted bigwigs have been regenerated, the existing Peaks are older than them and
# Snakemake would re-run region calling, silently giving each arm its own region set. Making
# this part of the run means the invariant cannot depend on operator order.
echo "=== re-stamping Peaks against current inputs ==="
"$D/.pixi/envs/multimodal/bin/python" \
    "$D/2026_0824_H3K27ac_model/scripts/4.3.setup_abc_arms.py" --copy-peaks --force-peaks \
    ${ABC_ARMS:+--arms $ABC_ARMS} ${ABC_RESULTS_DIR:+--results-dir $ABC_RESULTS_DIR} \
    | tail -3

cd "$A"
mkdir -p .snakemake/slurm_logs

echo "snakemake: $($SM --version)"
echo "profile:   $PROFILE"
echo

# A driver killed mid-run leaves the working directory locked, and the next attempt dies
# with LockException before doing anything. Clear it -- but only after confirming no other
# ABC driver is alive, because unlocking a live workflow lets two Snakemake instances write
# the same files.
OTHERS=$(squeue -u "$USER" -h -o '%i %j %T' | awk '$2 ~ /abc_driver/ && $3 == "RUNNING" {print $1}' | grep -v "^${SLURM_JOB_ID:-none}$" || true)
if [[ -n "$OTHERS" ]]; then
    echo "ABORT: another abc_driver is RUNNING (job(s): $OTHERS)." >&2
    echo "Cancel it before starting a new one; do not unlock a live workflow." >&2
    exit 1
fi
echo "=== clearing any stale lock ==="
$SM --configfile "$CFG" --unlock || true

echo "=== dry run ==="
# The dry run must use the SAME flags as the real run, or its verdict is about a different
# invocation. Omitting --profile meant it lacked the profile's rerun-incomplete, so it died
# with IncompleteFilesException on outputs left by a cancelled run while the real run would
# have handled them -- the gate was stricter than the thing it was gating.
$SM --configfile "$CFG" --profile "$PROFILE" --use-conda -n -q > /tmp/abc_dryrun.$$ 2>&1 || { cat /tmp/abc_dryrun.$$; exit 1; }
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
