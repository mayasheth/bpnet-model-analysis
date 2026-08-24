# Environments & Installations

This repo predates the lab pixi-only standard and uses **four environments across pipeline
stages**. There is no single env that runs everything — pick by stage. See
`.living/conventions/engreitz-lab/environment-conventions.md` for the lab standard this
deviates from, and `.living/conventions.md` for why.

All work runs on **Sherlock**; the repo lives on Oak at
`$OAK/Users/sheth/EP300_BPNet` (`$OAK` = `/oak/stanford/groups/engreitz`).

## Primary Environment

- **Manager**: pixi (`pixi.toml` + `pixi.lock`, both committed)
- **Python version**: 3.9.*
- **Created**: 2025-04-07 (`.pixi/`), lockfile last updated 2026-07-08

### Setup from scratch

```bash
cd $OAK/Users/sheth/EP300_BPNet
pixi install          # resolves from the committed pixi.lock
pixi run <command>
```

## Conda environments (per pipeline stage)

| Env | Stages | Notes |
|---|---|---|
| `bpnet_37` | Training (1.x), prediction (2.x), SHAP (3.x) | Requires `module load cuda/11.1.1 cudnn/8.1.1.33` before activation |
| `tfmodisco` | MoDISCo (4.x), FIMO (7.x), model inference | |
| `analysis` | FiNeMo formatting (5.2+), downstream plots and stats | |
| `multimodal` (pixi) | Multimodal model training/prediction | Local to `2026_0529_multimodal_p300_model/`. bedtools/samtools/bedGraphToBigWig plus leidenalg/igraph pulled via conda to avoid source compilation |

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
