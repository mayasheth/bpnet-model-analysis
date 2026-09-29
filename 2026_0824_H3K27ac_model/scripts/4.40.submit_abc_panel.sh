#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 48:00:00
#SBATCH --mem=8G
#SBATCH -c 2
#SBATCH -o log/abcpanel.%j.txt
#SBATCH --job-name=abc_driver
#
# ABC with observed p300 in the four non-K562 p300-panel cell types (config from 4.39).
#
# UNLIKE 4.5, region calling is SUPPOSED to run here: each cell type gets its own candidate
# regions from its own ATAC, the same way K562's July set was made. So there is no Peaks
# re-stamping and no region-calling gate; the dry run instead has to SHOW call_macs_peaks
# and make_candidate_regions for every arm. Environment handling (mamba shim, stale-lock
# check) is 4.5's; see that script for why each piece is there.
set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
A=/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction
SM_ENV=/oak/stanford/groups/engreitz/Users/sheth/.conda/envs/run_snakemake
SM=$SM_ENV/bin/snakemake
PROFILE="$HOME/.config/snakemake/slurm"
CFG=${ABC_CFG:-config/mine/config_p300_panel_obs.yaml}
grep -q 'SLURM_CPUS_PER_TASK:-1' \
    "$SM_ENV"/lib/python3.11/site-packages/snakemake/executors/slurm/slurm_jobstep.py || {
    echo "ERROR: the srun patch is missing from $SM_ENV (env rebuilt?). See learnings." >&2; exit 1; }
SHIM=$D/2026_0824_H3K27ac_model/.mamba_shim
mkdir -p "$SHIM"; ln -sfn "$SM_ENV/bin/mamba" "$SHIM/mamba"; ln -sfn "$SM_ENV/bin/conda" "$SHIM/conda"
export PATH="$PATH:$SHIM"
"$D/.pixi/envs/multimodal/bin/python" "$D/2026_0824_H3K27ac_model/scripts/4.39.setup_abc_panel.py"
cd "$A"; mkdir -p .snakemake/slurm_logs
OTHERS=$(squeue -u "$USER" -h -o '%i %j %T' | awk '$2 ~ /abc_driver/ && $3 == "RUNNING" {print $1}' | grep -v "^${SLURM_JOB_ID:-none}$" || true)
[[ -z "$OTHERS" ]] || { echo "ABORT: another abc_driver is RUNNING ($OTHERS); two drivers must not share the ABC workdir" >&2; exit 1; }
$SM --configfile "$CFG" --unlock || true
DRY=$D/2026_0824_H3K27ac_model/log/abcpanel_dryrun.${SLURM_JOB_ID:-local}.txt
$SM --configfile "$CFG" --profile "$PROFILE" --use-conda -n -q > "$DRY" 2>&1 || { cat "$DRY"; exit 1; }
cat "$DRY"
n=$(awk '$1=="make_candidate_regions"{print $2}' "$DRY")
[[ "${n:-0}" == 4 ]] || { echo "ABORT: expected make_candidate_regions for 4 arms, got '${n:-0}'" >&2; exit 1; }
$SM --configfile "$CFG" --profile "$PROFILE" --use-conda
echo "ABC panel complete"
