# FiNeMo motif analysis — HANDOVER

Last updated: 2026-06-11

This document tracks the motif logo figures and TF-correlation analysis for FiNeMo runs
on the primary EP300 K562 model.

---

## FiNeMo runs

| Run | Directory | Motifs | Notes |
|-----|-----------|--------|-------|
| v1 (original) | `finemo/pkw_500_curated_motifs` | 9 (GATA, AP1_1, REPEAT_G, GATA_TAL1, ETS_1, STAT_2, CREB_ATF_3, CTCF, CREB_ATF_1) | Original curated motifs |
| v2 (current) | `finemo/pkw_500_curated_motifs_v2` | 8 active (FOS_JUN excluded — 0 hits) | Updated motif set; all analyses below use v2 |

---

## v2 motif order (by total hits, active motifs only)

| Rank | Motif | Total hits |
|------|-------|-----------|
| 1 | GATA | 182,579 |
| 2 | REPEAT_G | 91,079 |
| 3 | NF2L_NFE | 59,157 |
| 4 | GATA_TAL1_9BP | 34,057 |
| 5 | ELF | 27,303 |
| 6 | STAT | 19,758 |
| 7 | ATF_CEBP | 8,609 |
| — | GATA_TAL1_8BP | 28,873 (excluded — redundant with GATA_TAL1_9BP) |
| — | FOS_JUN | 0 (excluded — no hits) |

---

## Inputs

| File | Description |
|------|-------------|
| `finemo/pkw_500_curated_motifs_v2/hits.tsv` | FiNeMo hit calls (motif × region) |
| `finemo/pkw_500_curated_motifs_v2/peaks_qc.tsv` | Peak region coordinates |
| `finemo/pkw_500_curated_motifs_v2/motif_report.tsv` | Per-motif hit counts |
| `finemo/pkw_500_curated_motifs_v2/motif_cwms.npy` | CWM matrices (N, 4, L) |
| `finemo/pkw_500_curated_motifs_v2/motif_data.tsv` | Restricted-region trim coords per motif |
| `finemo/pkw_500_curated_motifs_v2/contribution_regions.npz` | SHAP contributions per region |
| `/oak/stanford/groups/engreitz/Users/agschwin/.../enhancer_activity_features.tsv.gz` | ENCODE ChIP-seq signal RPMs (523 signals, ~150k elements) |
| `/oak/stanford/groups/engreitz/Users/sheth/TF_analysis/2025_0609_K562_TF_annotations/finemo_peaks_all_chr.chromatin_annotations.tsv` | peak_id → EP300_peak_overlap, true/pred logcounts |

Note: v2 and v1 peaks have identical coordinates and peak_ids (both are the same 150,519 K562 accessible elements). The chromatin annotations file uses 1000bp windows but the same peak_id indexing — join directly on `peak_id`.

---

## Scripts (all parameterized with `--finemo-dir`)

| Script | Description | Default `--finemo-dir` |
|--------|-------------|----------------------|
| `scripts/plot_finemo_motif_logos_row.py` | Horizontal row of forward CWM logos (restricted window, PDF only) | v2 |
| `scripts/finemo_tf_correlation.py` | Motif × ChIP-seq Spearman + Pearson + Pearson-log1p correlations | v2 |
| `scripts/plot_finemo_explained_fraction.py` | Fraction of total SHAP explained by motif hits | v2 |
| `scripts/plot_finemo_hit_importance_summary.py` | Hit importance violin + % hits vs % importance bar chart | v2 |
| `scripts/5.3.plot_finemo_hits.R` | Bubble plot, frequency/importance bars, percentile lines, heatmaps | args only |
| `scripts/submit_finemo_tf_correlation.sh` | SLURM wrapper for correlation analysis (48 GB, 2 h) | v2 |
| `scripts/submit_finemo_v2_summary_plots.sh` | SLURM wrapper for all summary plots (32 GB, 1 h) | v2 |

### Key script changes (2026-06-11)
- All Python scripts now accept `--finemo-dir` — reusable across FiNeMo runs
- Logos use `motif_cwms.npy` + `motif_data.tsv` restricted window (consistent with individual SVG logos)
- `EXCLUDE_MOTIFS = {"FOS_JUN", "GATA_TAL1_8BP"}` applied in all Python scripts and R input prep
- `finemo_tf_correlation.py`: added `pearson_r` and `pearson_log1p_r` columns to output alongside `spearman_r`; top10_per_motif.pdf now shows 3 bar panels per motif (one per metric)
- `5.3.plot_finemo_hits.R`: color palette now derived dynamically from hits data (sorted by frequency) instead of hardcoded motif names

### Running order for a new FiNeMo run
```bash
# logos row (fast, login node OK if pixi available)
pixi run -e ism python scripts/plot_finemo_motif_logos_row.py --finemo-dir <dir>

# TF correlations (needs ~48 GB, 2 h)
sbatch scripts/submit_finemo_tf_correlation.sh   # edit FINEMO_DIR inside first

# summary plots: explained fraction + violin + bubble plot (needs ~32 GB, 1 h)
sbatch scripts/submit_finemo_v2_summary_plots.sh  # edit FINEMO_V2 inside first
```

---

## Output files — v2

### `finemo/pkw_500_curated_motifs_v2/`
| File | Description | Status |
|------|-------------|--------|
| `motif_logos_row.pdf` | Restricted-window CWM logos, 7 motifs (FOS_JUN + GATA_TAL1_8BP excluded) | Done (2026-06-11) |
| `motif_importance_violin.pdf` | Hit importance distributions per motif | Done (2026-06-11, job 29079885) |
| `motif_importance_vs_frequency.pdf` | % total hits vs % total importance (paired bar) | Done (2026-06-11, job 29079885) |

### `finemo/pkw_500_curated_motifs_v2/annotated_motifs/`
| File | Description | Status |
|------|-------------|--------|
| `hits_renamed.tsv` | v2 hits.tsv filtered (FOS_JUN + GATA_TAL1_8BP excluded) — R script input | Done (2026-06-11, job 29079885) |
| `hits_per_peak.with_predictions.tsv.gz` | Peaks + true/pred logcounts + n_hits — R script input | Done (2026-06-11, job 29079885) |
| `motif_explained_fraction.pdf` | Bar chart: % total SHAP per motif + unexplained | Done (2026-06-11, job 29079885) |
| `motif_explained_fraction.motifs_only.pdf` | Same, motifs only (zoomed) | Done (2026-06-11, job 29079885) |

### `finemo/pkw_500_curated_motifs_v2/plot_annotated_motifs/`
| File | Description | Status |
|------|-------------|--------|
| `motif_bubble_plot.pdf` | Frequency vs avg importance scatter (bubble size = total importance) | Done (2026-06-11, job 29079885) |
| `motif_frequency_importance.pdf` | Frequency + total importance side-by-side bars | Done (2026-06-11, job 29079885) |
| `motif_frequency_by_percentile.pdf` | Motif frequency vs p300 signal percentile lines | Done (2026-06-11, job 29079885) |
| `motif_scores_by_percentile.pdf` | Avg hit scores vs p300 signal percentile | Done (2026-06-11, job 29079885) |
| `motif_importance_heatmap.pdf` | Avg importance by p300 signal bin (heatmap) | Done (2026-06-11, job 29079885) |
| `motif_similarity_heatmap.pdf` | Avg similarity by p300 signal bin (heatmap) | Done (2026-06-11, job 29079885) |
| `counts_per_peak.pdf` | Hits per peak distribution + observed/predicted logcounts | Done (2026-06-11, job 29079885) |
| `counts_per_peak.top50pct.pdf` | 3-panel: hits/peak histogram + predicted/observed logcount violins; top-50% hits by hit_coefficient_global; all 150k elements | Done (2026-06-12) |
| `counts_per_peak.top50pct.p300plus.pdf` | Same; EP300_peak_overlap==1 elements only (n=27,340; ~40% have true_logcounts=0 — window mismatch) | Done (2026-06-12) |
| `counts_per_peak.top50pct.top20pct_obs.pdf` | Same; union of top-20% by true_logcounts OR EP300_peak_overlap==1 (n=41,858) — preferred p300+ definition | Done (2026-06-12) |

**Note on p300+ definition:** `EP300_peak_overlap` is based on 1000bp windows; `true_logcounts` is measured over the 500bp BPNet window. ~40% of overlap-flagged elements have zero counts. The union (top-20%-obs OR overlap flag) is the preferred definition. See also `TODO.md` for list of figures to update if this definition is adopted globally.

Script: `scripts/plot_finemo_counts_per_peak_top50pct.py`

### `finemo/pkw_500_curated_motifs_v2/tf_correlation/`
| File | Description | Status |
|------|-------------|--------|
| `correlations.tsv.gz` | Spearman r, Pearson r, Pearson log1p r — all motif × signal pairs | Done (2026-06-11, job 29081409) |
| `correlation_density.pdf` | KDE distributions, 3-panel (one per metric) | Done (2026-06-11, job 29081409) |
| `top10_per_motif.pdf` | Logo + top-10 bar charts for all three metrics | Done (2026-06-11, job 29081409) |
| `top50pct_hit_rpm.pdf` | Violin: top-50% hit regions vs rest per motif | Done (2026-06-11, job 29081409) |

---

## Analysis design notes

**Region merge (TF correlations):**
- pybedtools intersect maps FiNeMo peaks → signal regions; largest-overlap wins for many-to-many
- Merged universe: ~132k regions (out of 150,519 FiNeMo peaks)

**Motif activity per region:**
- `hit_importance` summed across all hits (both strands) for a given motif × region; 0 if no hits

**Correlation metrics:**
- `spearman_r` — rank-based; robust to zero-inflation
- `pearson_r` — standard Pearson on raw signal values
- `pearson_log1p_r` — Pearson after log1p transform of both vectors; compare to assess sensitivity to distribution shape

**Top-50% hit analysis:**
- Per motif: regions with ≥1 hit at ≥50th percentile of `hit_coefficient_global`
- Mann-Whitney U test vs remaining regions

---

## Open questions / next steps

- [x] Check job 29081409 (TF correlations) outputs — done
- [ ] Review v2 logos — confirm motif identities are sensible (especially REPEAT_G, NF2L_NFE, ELF)
- [ ] Decide whether to exclude additional motifs from publication figures (REPEAT_G is a likely candidate)
- [ ] Compare Spearman vs Pearson vs Pearson-log1p results — note which metric best separates expected TF-motif pairs (e.g. GATA → GATA1/2, ELF → ELF1)
- [ ] Consider restricting correlation universe to p300+ elements only
- [ ] Potential additional figure: motif × TF correlation heatmap (all 8 × top-N TFs)
- [ ] Update `P300_INTERACTORS` set in `plot_finemo_composite_figure.py` after BioGRID review

---

## v1 outputs (for reference)

All v1 outputs remain under `finemo/pkw_500_curated_motifs/`. The v1 `annotated_motifs/` and `tf_correlation/` directories contain analogous figures for the original 9-motif set.
