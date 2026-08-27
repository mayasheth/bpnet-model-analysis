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

## SLURM submit scripts (mandatory)

Every new submit script in this repo must:

1. **`export PYTHONUNBUFFERED=1`.** Python block-buffers stdout to a file, so a running
   job shows an empty log and looks hung. Without this, a 35-minute job is
   undiagnosable while alive.
2. **Guard empty array expansion.** `set -u` plus `"${ARR[@]}"` on an empty array is an
   unbound-variable error in bash < 4.4. Use `${ARR[@]+"${ARR[@]}"}` for any optional
   argument array.
3. **Put the varied hyperparameter in the output path.** A sweep whose runs share an
   output directory silently overwrites itself. See
   `2026_0824_H3K27ac_model/scripts/1.1.submit_training.sh`, where the counting window
   and count-loss weight both appear in `OUT_DIR`.
4. **Bound memory explicitly when a set size scales with the data.** Anything of the
   form `n_items * k` needs a cap; see `--max-negatives`.

## Never run heavy work on a login node

Anything beyond a few seconds goes through `sbatch`. A long python process started over
`ssh` on a login node is killed when the calling shell detaches — it exits 0 and writes
no output, which reads as success. Login-node work is limited to inspecting files and
small tabulations; the one window-profiling script that does run there is explicitly
memory-bounded (264 MB) and documented as such.

## Model evaluation must be stratified by signal level

`reference/K562_DNase_candidate_elements.narrowPeak` holds 150,528 elements and most
carry little signal, so any correlation computed over all of them is dominated by the
dead-vs-active contrast. Report the top signal quintile alongside the overall number —
this repo's equivalent of the all-vs-p300+ split already used for p300. See
`2026_0824_H3K27ac_model/scripts/2.2.evaluate_stratified.py`. Concretely: ATAC-only
predicts H3K27ac at 0.746 over all elements but only 0.543 on the top quintile.

Ceilings derived from replicate agreement need both the Spearman-Brown and sqrt
corrections before they bound model performance — see `.living/learnings.md`.

## Decision-log entries use `###`, not `##`

Mycelium's `decision-log-entry.md` template says `## [YYYY-MM-DD] Title`, but
`generate_index.py` only counts `### ` headings in `decisions.md`. Entries written per
the template are silently absent from `.living/INDEX.md`. Use `### ` for decisions here.
See `.living/learnings.md` (2026-08-25).

## Never take a model target from a display directory

Tracks under `$OAK/Users/sheth/Data/share/IGV/` are built for visualization. They may be
raw-count and unnormalized — and so look usable — while carrying processing choices that
are wrong for modelling. `ENCSR000AKP_coverage.bw` extends 36 bp reads to a fixed 250 bp
fragment, which flattens nucleosome structure and breaks the multinomial assumption
behind bpnetlite's profile loss.

Derive training targets from the source BAMs inside the analysis project, so the
processing is explicit, versioned, and reviewable. For ChIP-seq in this repo that means
`bedtools genomecov -5 -dz` (single-base 5' ends, stranded), matching
`scripts/0.3.make_training_bw.sh`. If a display track is used as a stopgap, record the
exact command that produced it and treat replacing it as an open item.

## H3K27ac target is 5-prime ends, full stop

`2026_0824_H3K27ac_model/data/h3k27ac_5p_{plus,minus}.bw`. Nothing new is trained or
evaluated on the 250 bp fragment-extended track — it breaks MNLL's read-count assumption
and bleeds signal between neighbouring elements. See
`2026_0824_H3K27ac_model/results/TARGET_PROVENANCE.md` for which existing results used
which target, and which conclusions survive.

Enforcement: `scripts/train_multimodal_bpnet.py` writes `training_target.json` into every
model directory and refuses a `--count-offset-model` whose signal tracks or output window
differ from the current run. This exists because a fragment-trained offset model was
silently used against a 5-prime target on 2026-08-26 and trained to a plausible-looking
but meaningless result.

## Sample definition for multi-dataset comparisons

Two rules when assembling panels from the lab data inventory
(`2023_0701 Maya-s sequences.xlsx`, "ENCODE data" sheet):

1. **No Mint-ChIP.** Different assay chemistry from standard ChIP-seq and not comparable.
   Excludes WTC11 H3K27ac/H3K27me3/H3K4me1 (ENCSR146DPQ, ENCSR418YEV, ENCSR473GYW).
2. **A sample is one (experiment, processing) pair.** Replicates may be pooled, or
   compared as an inter-replicate ceiling, only within a single experiment AND a single
   processing. Rows in the sheet sometimes pool two experiments (all four TeloHAEC
   H3K27ac rows list GSE210489 and GSE210491) or mix run types across replicates
   (HCT116 rows are `se, pe`); both must be split before use. Pooling across either
   folds batch or assay differences into what gets reported as biological
   reproducibility, which understates the ceiling.

Recorded because the ceiling is the denominator for every performance number in the
H3K27ac work, so contaminating it silently rescales all of them.
