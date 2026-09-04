# EP300_BPNet project

> **Knowledge index (read first):** [`.living/INDEX.md`](.living/INDEX.md) is an auto-generated map of tag clusters, most-recent entries, and a tag → entry-ID inverted index. The SessionStart hook keeps it fresh — trust it. For targeted lookup: `module load python/3.12.1 && python3 $OAK/Users/sheth/EngreitzLabAgents/mycelium-upstream/skills/core/scripts/recall_lessons.py --living-dir .living/ --tag <tag>` (also `--id L-42`, `--since YYYY-MM-DD`).

## Overview

This project interprets what BPNet sequence models have learned about **TF motif SYNTAX** (not just identity) that drives p300 coactivator binding. The core approach uses deep learning models trained on ChIP-seq data, combined with interpretability methods (SHAP, MoDISCo) and in silico motif insertion experiments to understand how motif arrangements, spacing, and orientations affect p300 recruitment.

## Models

### p300 v1 (primary)
Main model trained on p300 ChIP-seq peaks (ENCSR000EGE) vs GC-matched negatives. 5-fold cross-validation.

**Directory**: `2025_0517_official_EP300_K562_model/`
```
├── models/release_run_1/fold{0-4}/ENCSR000EGE/  # TensorFlow SavedModel format
├── models_h5/fold{0-4}/                          # H5 format for inference
├── predictions_cv/                               # Cross-validation predictions
├── predictions_mean/                             # Mean predictions across folds
├── shap/                                         # SHAP values (all candidate elements)
├── shap_peaks/                                   # SHAP values (peaks only)
├── modisco/                                      # MoDISCo results
├── modisco_peaks/                                # MoDISCo on peaks only
├── finemo/                                       # FiNeMo hit calling results
├── motif_spacing/                                # Motif spacing experiment results
├── scripts/0.0.log.sh                            # Complete workflow log with commands
└── scripts/                                      # Model-specific scripts
```

### Other models
- **K562_DNase_ChromBPNet/**: DNase accessibility model for comparison
- **K562_GATA1_BPNet/**: GATA1 ChIP-seq model
- **K562_GATA2_BPNet/**: GATA2 ChIP-seq model
- **K562_ATAC_ChromBPNet/**: ATAC-seq accessibility model

## BPNet reference

For model architecture and training details, see: `/oak/stanford/groups/engreitz/Users/sheth/bpnet-refactor`

BPNet models predict:
- **Profile head**: Base-resolution ChIP-seq signal shape
- **Counts head**: Total read counts (log scale) - primary output used in analyses

## Key workflow log

The file `2025_0517_official_EP300_K562_model/scripts/0.0.log.sh` contains the complete command history for the p300 model workflow, including:
- Model download from TFAtlas
- Prediction commands (CV and mean)
- SHAP computation (all elements and peaks-only)
- MoDISCo runs with different parameters
- FiNeMo hit calling
- Environment setup (conda activate, module loads)

---

## Pipeline stages

### Stage 0: Data preparation
| Script | Purpose |
|--------|---------|
| `0.1.download_tfatlas_bpnet.sh` | Download pre-trained BPNet models from TFAtlas |
| `0.2.bed_to_narrowPeak.sh` | Convert BED to narrowPeak format for training |
| `0.3.make_training_bw.sh` | Create training BigWig files from BAM |

### Stage 1: Model training
| Script | Purpose |
|--------|---------|
| `1.1.submit_training_one_fold.sh` | Submit BPNet training jobs (5-fold CV) via SLURM |

### Stage 2: Prediction
| Script | Purpose |
|--------|---------|
| `2.1.submit_cv_predict.sh` | Cross-validation predictions |
| `2.2.submit_mean_predict.sh` | Mean predictions across CV folds |
| `2.3.compute_prediction_performance.py` | Calculate Pearson/Spearman correlations |

### Stage 3: SHAP interpretation
| Script | Purpose |
|--------|---------|
| `3.1.submit_mean_shap_one_fold.sh` | Compute SHAP values for model interpretability |

### Stage 4: Motif discovery
| Script | Purpose |
|--------|---------|
| `4.1.submit_counts_modisco.sh` | Run TF-MoDISco-lite on SHAP profiles |
| `extract_modisco_motifs.py` | Extract CWMs from MoDISCo HDF5 files |

### Stage 5: Motif hit calling
| Script | Purpose |
|--------|---------|
| `5.1.submit_finemo.sh` | Run FiNeMo motif hit calling |
| `5.2.format_finemo_hits.py` | Format FiNeMo output for analysis |
| `5.3.plot_finemo_hits.R` | Visualize FiNeMo results |
| `5.4.plot_hits_vs_chip_data.R` | Compare hits to ChIP-seq data |

### Stage 6: In silico motif experiments
| Script | Purpose |
|--------|---------|
| `6.0.motif_spacing.one_motif.py` | Test model response to varying counts of one motif |
| `6.1.motif_spacing.two_motifs.py` | Test model response to two motifs at various spacings |
| `6.2.motif_pairs.py` | Comprehensive motif pair insertion experiments |
| `6.2.1.submit_motif_pairs.sh` | SLURM submission for motif pairs |
| `plot_motif_spacing.py` | Visualization for spacing experiments |
| `plot_motif_pairs.py` | Heatmaps for motif pair synergy |

### Stage 7: FIMO analysis
| Script | Purpose |
|--------|---------|
| `run_fimo.py` | Run FIMO motif scanning (PWM-based) |
| `7.0.create_region_mapping.py` | Map FIMO regions to annotation regions |
| `7.1.fimo_motif_analysis.py` | Comprehensive FIMO analysis and plotting |

---

## Key utility module

**`motif_exp_utils.py`** - Core functions for model predictions and motif manipulations:

```python
get_model(model_path)                    # Load BPNet/ChromBPNet model
make_model_prediction(mod_encoded, ...)  # Universal prediction wrapper (returns log counts)
insert_motifs_with_orientation_general() # Insert motifs into sequences
one_hot_encode(sequences, seq_length)    # DNA to one-hot encoding
dinuc_shuffle(seq)                       # Dinucleotide-preserving shuffle
generate_motif_pairs(motif_dict)         # Generate all motif pair combinations
```

---

## Tabular file formats

### Motif pair results (`6.2.motif_pairs.py` output)

**`motif_pairs.raw_results.tsv.gz`**
| Column | Description |
|--------|-------------|
| `motif1_name` | Name of first motif |
| `motif2_name` | Name of second motif |
| `orientation` | Orientation pattern (e.g., "++", "+-", "-+", "--") |
| `spacing` | Distance between motifs (bp) |
| `pred_log_counts` | Predicted log counts from model |
| `log2_fc_vs_baseline` | Log2 fold-change vs dinucleotide shuffled background |
| `log2_synergy` | Log2 synergy score (pair effect beyond individual motifs) |

**`individual_motifs.raw_results.tsv.gz`**
| Column | Description |
|--------|-------------|
| `motif_name` | Motif name |
| `orientation` | "+" or "-" |
| `pred_log_counts` | Predicted log counts |
| `log2_fc_vs_baseline` | Log2 fold-change vs baseline |

### Motif spacing results (`6.0`, `6.1` output)

**`motif_spacing.one_motif.results.tsv`**
| Column | Description |
|--------|-------------|
| `motif_name` | Motif name |
| `n_motifs` | Number of motif copies inserted |
| `spacing` | Distance between copies (bp) |
| `mean_pred` | Mean predicted log counts |
| `std_pred` | Standard deviation |

### FIMO analysis output (`7.1.fimo_motif_analysis.py`)

**`motif_enrichment.tsv`**
| Column | Description |
|--------|-------------|
| `motif` | Motif name |
| `p300+_regions` | Count in p300+ regions |
| `p300-_regions` | Count in p300- regions |
| `odds_ratio` | Fisher's exact test odds ratio |
| `p_value` | Fisher's exact test p-value |
| `fdr` | Benjamini-Hochberg corrected p-value |

**`motif_pair_enrichment.tsv`**
| Column | Description |
|--------|-------------|
| `motif1` | First motif |
| `motif2` | Second motif |
| `p300+_cooccur` | Co-occurrence count in p300+ |
| `p300-_cooccur` | Co-occurrence count in p300- |
| `odds_ratio` | Enrichment odds ratio |
| `p_value` | Fisher's exact test p-value |

**`spacing_distributions.tsv`**
| Column | Description |
|--------|-------------|
| `motif1` | First motif |
| `motif2` | Second motif |
| `spacing` | Distance between motif centers |
| `count` | Frequency |
| `p300_status` | "p300+" or "p300-" |

### FiNeMo hits (`5.2.format_finemo_hits.py` output)

**`finemo_hits.formatted.tsv`**
| Column | Description |
|--------|-------------|
| `chrom` | Chromosome |
| `start` | Start position |
| `end` | End position |
| `motif_name` | Motif name |
| `strand` | "+" or "-" |
| `contribution_score` | SHAP-based importance score |
| `peak_id` | Associated peak identifier |

### Region mapping (`7.0.create_region_mapping.py`)

**`fimo_to_annotation_mapping.tsv`**
| Column | Description |
|--------|-------------|
| `fimo_region` | FIMO region name (chr:start-end) |
| `annot_region` | Matched annotation region |
| `overlap` | Overlap length (bp) |

---

## Key input files

| File | Description | Example path |
|------|-------------|--------------|
| Training peaks | narrowPeak format | `/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/reference/ENCSR000EGE_peaks_inliers.narrowPeak` |
| Candidate elements | FASTA sequences | `/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/reference/K562_DNase_candidate_elements.fa` |
| Motif compendium | MEME format | `/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/reference/MotifCompendium-Database-Human.meme.txt` |
| Chromatin annotations | TSV with p300 status | `/oak/stanford/groups/engreitz/Users/sheth/TF_analysis/2025_0609_K562_TF_annotations/finemo_peaks_all_chr.chromatin_annotations.tsv` |

---

## SLURM partitions — pick deliberately

**`engreitz` is the lab's own partition. Use it for everything that does not need a GPU.**
9 nodes, 24+ cores, 192 GB+, 7-day limit, no GPUs. It does not compete with the general GPU
queues, so it sidesteps a fairshare depleted by training runs.

```bash
#SBATCH -p engreitz,normal,owners   # CPU work: lab partition first, fallbacks after
#SBATCH -p gpu,owners               # training only
```

| job | partition | why |
|---|---|---|
| Training (per fold, ~1–3 h) | `gpu,owners` | needs a GPU; short enough that preemption is cheap |
| Inference / prediction | `engreitz,normal,owners` | no GPU needed. `torch.cuda.is_available()` makes the fallback automatic, so the same script runs either way |
| Track building, evaluation, plotting | `engreitz,normal,owners` | CPU only |
| One long uncheckpointed step | `engreitz,normal` — **never `owners`** | preemption loses the whole thing. The fragment-channel build (one pass over three 9 GB BAMs into a temp dir under a cleanup trap) was preempted at 40 minutes |
| Resumable Snakemake DAG | `engreitz,normal,owners` | rule outputs persist, so preemption costs one rule; `--rerun-incomplete` handles partials. Take the faster queue |

Two further habits that cost time on 2026-09-03:

- **Request less.** 8 cores / 64 GB backfills far more easily than 16 / 128, and the counting
  and inference steps here stream rather than holding anything large.
- **Submit in dependency order, and `scontrol hold` the speculative batch.** Twenty training
  jobs submitted inside an hour drained fairshare and left the work that mattered queued
  behind exploratory follow-ups.

## Ask before running a broad search

**Before any "search for anything" tool call — a repo-wide `grep -r`, a `find` over Oak, or
hunting for where something lives — ask Maya first.** She usually knows the path outright or
can narrow it to one directory, and that is faster and more reliable than guessing.

Two failures on 2026-09-03/04 make the case. A `find` over Oak for an element set timed out
and returned nothing, which would read as "the file does not exist" if it had not been
recorded as a non-result. And a `grep -rl 'executor:'` across four lab repos found no
snakemake9 profile, so this file briefly asserted that none existed; one question produced
`seq-processing-snakemake`, whose profile sits in a hidden `.snakemake_profile/` that the
search patterns never covered.

A timed-out or empty search is not evidence of absence. Targeted `ls` of a named directory
is fine and does not need asking; open-ended discovery does.

## Running Snakemake pipelines on Sherlock

Do not reach for whatever `snakemake` happens to be on `PATH`. Several stale ones exist in
`/home/groups/engreitz/Software/anaconda3/envs` (5.5.4, 5.10, 5.32), and
`~/.conda/envs/final-abc-env` has the bioinformatics tools but no snakemake and a python that
segfaults on import. Use the dedicated wrapper envs, which live on **Oak**, not in `~/.conda`:

| pipeline targets | env | version |
|---|---|---|
| Snakemake 7 | `$OAK/Users/sheth/.conda/envs/run_snakemake` | 7.32.4 |
| Snakemake 9 | `$OAK/Users/sheth/.conda/envs/run_snakemake9` | 9.6.0 |

Both carry `conda` and `mamba`, so `--use-conda` works from either.

**Always pass a `--profile`.** Without one, Snakemake runs every rule inside the submitting
allocation; with one, each rule instance becomes its own SLURM job, which is both faster and
the only way a large DAG fits sensible resource requests.

```bash
SM_ENV=$OAK/Users/sheth/.conda/envs/run_snakemake
export PATH="$SM_ENV/bin:$PATH"     # REQUIRED for --use-conda
$SM_ENV/bin/snakemake --configfile <cfg> --profile ~/.config/snakemake/slurm --use-conda
```

The `PATH` line is not optional. `--use-conda` shells out to `mamba` from `/usr/bin/bash`,
which does not inherit the wrapper env just because snakemake was invoked by absolute path;
without it the run dies with `CreateCondaEnvironmentException` before submitting anything.
`--conda-frontend conda` is the alternative if mamba is genuinely unavailable.

`~/.config/snakemake/slurm` — 50 concurrent jobs, `slurm_partition=engreitz,owners,normal`,
`slurm_account=engreitz`, 6 h default runtime, 3 retries, `rerun-incomplete`.
`~/.config/snakemake/slurm_long` is the same with 100 jobs and 48 h.

**The v7 and v9 profiles are NOT interchangeable.** The two under `~/.config/snakemake`
(`slurm`, `slurm_long`) are v7-style: submission goes through a `cluster:` command string.
Snakemake 9 uses the executor plugin (`snakemake_executor_plugin_slurm` 1.4.0, installed in
`run_snakemake9`), which needs `executor: slurm` and a nested `cluster:` mapping of
`submit-cmd` / `status-cmd` / `cancel-cmd`, and ignores a v7 `cluster:` string.

v9 pipelines keep their profile **inside the repo**, not in `~/.config/snakemake`. Working
template, copy this rather than writing one:

```
$OAK/Users/sheth/seq-processing-snakemake/.snakemake_profile/slurm/config.yaml
```

and its README gives the invocation:

```bash
conda activate run_snakemake9
snakemake --configfile config/config.yml --profile .snakemake_profile/slurm
```

That profile also bakes in `use-conda: true`, `conda-frontend: mamba`, `keep-going: true`
and `default-resources` with `slurm_partition=engreitz`, so those need not be passed.

**Submit the driver as its own small job.** With a profile the driver only orchestrates, so
2 cores and 8 GB is plenty, but give it a long wall clock (48 h) and a non-preemptible
partition: if the driver dies, running children finish and nothing further is submitted.

**Gate the run on a dry run** (`-n -q`) whenever correctness depends on rules being skipped.
The ABC arms share one candidate-region set by pre-populating each arm's `Peaks/`; if
`call_macs_peaks` ever appears in the dry run, the arms would silently get different region
sets and every cross-arm comparison would be meaningless. `scripts/4.5.submit_abc.sh` aborts
in that case rather than producing nine incomparable answers.

## Environment and dependencies

Primary conda environments:
- `bpnet_37` - For BPNet training, prediction, and SHAP (requires `module load cuda/11.1.1 cudnn/8.1.1.33`)
- `tfmodisco` - For FIMO, MoDISCo, and model inference
- `analysis` - For FiNeMo formatting and downstream analysis

Key packages:
- `torch` - Model inference
- `tangermeme` - SHAP computation
- `pandas`, `numpy` - Data manipulation
- `matplotlib`, `seaborn` - Visualization
- `scipy.stats` - Statistical tests (Fisher's exact, Mann-Whitney, Spearman)
- `memelite` - FIMO motif scanning

---

## Style guidelines

See `scripts/STYLE_GUIDELINES.md` for plot formatting:
- Sentence case for axis titles
- Use "p300" not "P300"
- No gridlines, black axes
- Color palettes:
  - Diverging: 'managua'
  - Sequential: 'PuBu'
  - p300 status: #792374 (p300+), #49bcbc (p300-)

---

## Directory structure

```
EP300_BPNet/
├── scripts/                              # Shared analysis scripts (this doc focuses on these)
├── reference/                            # Input files (peaks, FASTA, motifs)
├── FIMO/                                 # FIMO scan results
│   └── elements_v1/
│       └── analysis_v1/
│
├── 2025_0517_official_EP300_K562_model/  # Primary p300 model (see Models section)
│   ├── models/                           # TensorFlow models (5 folds)
│   ├── models_h5/                        # H5 models for inference
│   ├── predictions_cv/                   # CV predictions
│   ├── predictions_mean/                 # Mean predictions
│   ├── shap/                             # SHAP scores (all elements)
│   ├── shap_peaks/                       # SHAP scores (peaks only)
│   ├── modisco/                          # MoDISCo results
│   ├── finemo/                           # FiNeMo hit calls
│   ├── motif_spacing/                    # Spacing experiments
│   ├── config/                           # Config files
│   ├── data/                             # Model-specific data
│   └── scripts/0.0.log.sh                # Workflow log with all commands
│
├── K562_DNase_ChromBPNet/                # DNase model
├── K562_GATA1_BPNet/                     # GATA1 model
├── K562_GATA2_BPNet/                     # GATA2 model
├── K562_ATAC_ChromBPNet/                 # ATAC model
├── 2025_0703_retrain_p300_model/         # p300 model v2
└── 2025_1016_p300_model_v3/              # p300 model v3
```


## Installed Convention Packs

- **engreitz-lab** — See `.living/conventions/engreitz-lab/analysis-conventions.md`

- **bioinformatics** — See `.living/conventions/bioinformatics/analysis-conventions.md`

- **idea-generator** — See `.living/conventions/idea-generator/analysis-conventions.md`

- **report-generator** — See `.living/conventions/report-generator/analysis-conventions.md`

- **robust-analysis** — See `.living/conventions/robust-analysis/analysis-conventions.md`
