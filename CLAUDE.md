# EP300_BPNet project

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
