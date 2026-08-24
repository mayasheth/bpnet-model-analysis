# Repo-Specific Conventions

Overrides to mycelium defaults or convention pack conventions.

## Layout: analyses live at the repo root, not under `analysis/`

Mycelium was initialized here on **2026-08-24**, years into the project. The canonical
mycelium directories (`analysis/`, `data/`, `algorithms/`, `reference_material/`) were
created but are **intentionally empty except for their manifests** — no files were moved,
because doing so would break every absolute path in the workflow logs
(`*/0.0.log.sh`), the SLURM submit scripts, and `CLAUDE.md`.

The real layout is:

| Mycelium expects | This repo uses |
|---|---|
| `analysis/<name>/` | dated dirs at root: `2025_0517_official_EP300_K562_model/`, `2026_0606_GM12878_transferability/`, … and model dirs `K562_*_BPNet/`, `K562_*_ChromBPNet/` |
| `data/raw`, `data/processed` | `reference/` for shared inputs; each analysis's own `data/` for project-specific signal |
| `algorithms/` | `scripts/` at root (shared) and each analysis's `scripts/` |
| `reference_material/` | `external/` for vendored third-party code; `reference/` for data |

**The manifests are the index into that real layout** — every entry carries a `path:` field
pointing at the actual location. Read the manifest, not the directory tree.

Naming convention for analyses: `YYYY_MMDD_short_description/`, where the date is when the
analysis started. Model directories predating that convention use `<CELLTYPE>_<TARGET>_<ARCH>/`.

## Reports

`/engreitzlab-report` (→ `/analysis-report`, Quarto → self-contained HTML) supersedes
mycelium's `report-generator` pack for anything written up here. The `report:` field in
`analysis/ANALYSIS_MANIFEST.md` is `null` for every entry because no analysis has been
routed through that path yet — manuscript figures are tracked in `README.md` instead.

## Figures

`scripts/STYLE_GUIDELINES.md` is the binding style reference for this project and takes
precedence over `.living/conventions/engreitz-lab/figure-conventions.md` where they differ:
sentence case axis titles, "p300" never "P300", no gridlines, black axes, diverging
`managua` / sequential `PuBu`, p300 status `#792374` (p300+) and `#49bcbc` (p300-).

## Environments

This repo predates the lab pixi standard and uses a mix: `pixi.toml` at root (Python 3.9)
plus three conda envs — `bpnet_37` (training/prediction/SHAP, needs
`module load cuda/11.1.1 cudnn/8.1.1.33`), `tfmodisco` (FIMO, MoDISCo, inference), and
`analysis` (FiNeMo formatting, downstream plots). The multimodal project has its own pixi
env named `multimodal`. Do not assume a single environment works across stages.
