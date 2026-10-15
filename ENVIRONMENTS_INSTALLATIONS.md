# Environments & Installations

Current work runs in **pixi** (two workspaces, below). The original p300 v1 pipeline predates
the lab pixi-only standard and still uses three conda envs. There is no single env that runs
everything — pick by stage. See
`.living/conventions/engreitz-lab/environment-conventions.md` for the lab standard this
deviates from, and `.living/conventions.md` for why.

All work runs on **Sherlock**; the repo lives on Oak at
`$OAK/Users/sheth/EP300_BPNet` (`$OAK` = `/oak/stanford/groups/engreitz`).

## Pixi environments (current work)

Everything from the H3K27ac/p300 model work onward (`2026_0824_H3K27ac_model/`,
`2026_1002_sce2g_tf_activity/`) runs in pixi. Two workspaces, both with committed lockfiles:

| Workspace | Env | Used for | Lock (sha256, first 12) |
|---|---|---|---|
| `pixi.toml` (repo root) | `multimodal` | Training, prediction and scoring (`scripts/train_multimodal_bpnet.py`, `A/scripts/1.x`, `2.x`, `4.1`); Python 3.9, torch from the cu118 wheel index, bpnet-lite at `fa5bcf5` | `069960d36312` |
| `pixi.toml` (repo root) | default, `ism`, `variant-scoring` | ChromBPNet tooling (workspace name `chrombpnet-env`); the default env has no bpnetlite | same |
| `2026_1002_sce2g_tf_activity/pixi.toml` | `tobias` (Python 3.10), `scprinter` (Python 3.11) | TOBIAS and scPRINTER footprinting for report 7 | `e32b6a713472` |

Always name the env: `pixi run -e multimodal python ...`. Plain `pixi run` uses the default env
and fails on `import bpnetlite` (review F53).

```bash
cd $OAK/Users/sheth/EP300_BPNet
pixi install -e multimodal      # resolves from the committed pixi.lock
pixi run -e multimodal python <script>
```

## Snakemake wrapper envs (ABC, ENCODE-rE2G, scE2G drivers)

Not pixi. They live on Oak, not in `~/.conda`; see the CLAUDE.md section "Running Snakemake
pipelines on Sherlock" for invocation.

| Env | Version |
|---|---|
| `$OAK/Users/sheth/.conda/envs/run_snakemake` | Snakemake 7.32.4 |
| `$OAK/Users/sheth/.conda/envs/run_snakemake9` | Snakemake 9.6.0 |

**Local patch:** `run_snakemake` carries a SLURM 25.11 srun patch in its `slurm_jobstep.py`
(`.living/learnings.md`, 2026-09-22). Without it every rule job dies on the `normal` partition.
An env rebuild reverts it. `A/scripts/4.33.submit_re2g_training.sh` checks for it; 4.5, 4.7 and
the scE2G driver do not yet (review F59).

R: neither pixi workspace contains R, so the few `Rscript` steps in report 7 run on whatever R
is on `PATH` (review F57, open).

## Conda environments (legacy, p300 v1 pipeline)

The original p300 BPNet pipeline (`2025_0517_official_EP300_K562_model/`, stages 1-7) predates
pixi and still uses these.

| Env | Stages | Notes |
|---|---|---|
| `bpnet_37` | Training (1.x), prediction (2.x), SHAP (3.x) | Requires `module load cuda/11.1.1 cudnn/8.1.1.33` before activation |
| `tfmodisco` | MoDISCo (4.x), FIMO (7.x), model inference | |
| `analysis` | FiNeMo formatting (5.2+), downstream plots and stats | |

Activation pattern used throughout the workflow logs:

```bash
source ~/.bashrc && conda activate bpnet_37 && module load cuda/11.1.1 cudnn/8.1.1.33
```

The literal per-stage invocations are recorded in the workflow logs — treat those as
authoritative over this table:
`2025_0517_official_EP300_K562_model/scripts/0.0.log.sh` and
`2026_0606_GM12878_transferability/0.0.log.sh`.

## Dependencies

| Package | Used for |
|---|---|
| `torch` | Model inference; the multimodal architecture is pure PyTorch |
| `tangermeme` | SHAP attribution over the counts head |
| `memelite` | FIMO PWM scanning |
| `tfmodisco-lite` | De novo motif discovery from SHAP profiles |
| `pandas`, `numpy` | Data manipulation |
| `matplotlib`, `seaborn` | Plotting (see `scripts/STYLE_GUIDELINES.md`) |
| `scipy.stats` | Fisher's exact, Mann-Whitney, Spearman |
| `python-igraph`, `leidenalg` | Clustering (pinned via conda in the multimodal env) |

## System Dependencies

- **SLURM** — all training, prediction, and SHAP run as batch jobs; logs land in
  `slurm_logs/` and per-project `log/`.
- **CUDA 11.1.1 / cuDNN 8.1.1.33** — module-loaded for `bpnet_37`. GPU required for
  training and SHAP.
- `bedtools`, `samtools`, `bedGraphToBigWig` — BigWig generation (multimodal env).

## External codebases

`bpnet-refactor` lives **outside this repo** at `$OAK/Users/sheth/bpnet-refactor`. Local
modifications are vendored as a patch in `external/bpnet-refactor-patch/`
(`modified_files.diff`, `new_files/`, `UPSTREAM_COMMIT.txt`). Reproducing a training run
means applying that patch at the recorded upstream commit — see `README.md`
§ "External codebases".

## Known environment gotchas

- The multimodal training run hit 7 distinct SLURM failures (generator deadlocks,
  device-placement bugs, OOM from loading all genome-wide negatives before subsampling).
  Fixes are documented in `2026_0529_multimodal_p300_model/HANDOVER.md` — read it before
  re-running that pipeline.
- Mycelium's own scripts need Python 3.11+ (`datetime.UTC`). On Sherlock use
  `module load python/3.12.1`; the repo's 3.9 pixi env will not run them.
