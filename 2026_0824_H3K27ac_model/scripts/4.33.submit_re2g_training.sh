#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 48:00:00
#SBATCH --mem=8G
#SBATCH -c 2
#SBATCH -o log/re2g.%j.txt
#SBATCH -e log/re2g.%j.txt
#SBATCH --job-name=re2g_driver
#
# Train ENCODE-rE2G with p300 tracks in the activity slot (configs from 4.32).
#
# Snakemake driver only; every rule goes out as its own SLURM job through the profile. Same
# environment handling as 4.5: snakemake by absolute path, mamba and conda exposed through a
# shim APPENDED to PATH on shared storage so no rule env gets shadowed.
#
# The dry run MUST NOT schedule region calling. rE2G imports ABC's rules with an abc_ prefix,
# so the gate looks for abc_call_macs_peaks, abc_sort_narrowpeaks and
# abc_make_candidate_regions. If they appear, the staged Peaks were not honoured and each
# dataset would get its own region set.
#
# REQUIRES the SLURM 25.11 srun patch in run_snakemake's slurm_jobstep.py (see
# .living/learnings.md, 2026-09-22). Without it every rule job dies on the normal partition.

set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
E=/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G
SM_ENV=/oak/stanford/groups/engreitz/Users/sheth/.conda/envs/run_snakemake
SM=$SM_ENV/bin/snakemake
PROFILE="$HOME/.config/snakemake/slurm"
CFG=${E2G_CFG:-config/mine/config_training_p300.yaml}

grep -q 'SLURM_CPUS_PER_TASK:-1' \
    "$SM_ENV"/lib/python3.11/site-packages/snakemake/executors/slurm/slurm_jobstep.py || {
    echo "ERROR: the srun patch is missing from $SM_ENV (env rebuilt?). See learnings." >&2
    exit 1; }

SHIM=$D/2026_0824_H3K27ac_model/.mamba_shim
mkdir -p "$SHIM"
ln -sfn "$SM_ENV/bin/mamba" "$SHIM/mamba"
ln -sfn "$SM_ENV/bin/conda" "$SHIM/conda"
export PATH="$PATH:$SHIM"
command -v mamba >/dev/null || { echo "ERROR: mamba not on PATH via $SHIM" >&2; exit 1; }

cd "$E"
mkdir -p .snakemake/slurm_logs "$D/2026_0824_H3K27ac_model/log"
echo "snakemake: $($SM --version)   config: $CFG"

OTHERS=$(squeue -u "$USER" -h -o '%i %j %T' | awk '$2 ~ /re2g_driver/ && $3 == "RUNNING" {print $1}' | grep -v "^${SLURM_JOB_ID:-none}$" || true)
[[ -z "$OTHERS" ]] || { echo "ABORT: another re2g_driver is running ($OTHERS)" >&2; exit 1; }
$SM -s workflow/Snakefile_training --configfile "$CFG" --unlock || true

echo "=== dry run ==="
DRY=$D/2026_0824_H3K27ac_model/log/re2g_dryrun.${SLURM_JOB_ID:-local}.txt
$SM -s workflow/Snakefile_training --configfile "$CFG" --profile "$PROFILE" --use-conda -n -q \
    > "$DRY" 2>&1 || { cat "$DRY"; exit 1; }
cat "$DRY"
if grep -qE '^[[:space:]]*abc_(call_macs_peaks|sort_narrowpeaks|make_candidate_regions)[[:space:]]' "$DRY"; then
    echo "ABORT: region calling is scheduled; staged Peaks were not honoured. Re-run 4.32 --stage-peaks." >&2
    exit 1
fi
echo "region calling is correctly skipped"

echo "=== run ==="
$SM -s workflow/Snakefile_training --configfile "$CFG" --profile "$PROFILE" --use-conda
echo "rE2G training complete"
