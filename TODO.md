# TODO

Last updated: 2026-06-08 (end of day)

---

## Model improvements / ongoing

- [ ] Consider hashFrag chromosome splits to reduce train/test sequence-similarity leakage
- [ ] Run SHAP on K562 ATAC multimodal model (not yet submitted)
- [ ] Train DNase multimodal variant (BigWig not yet generated)
- [ ] ATAC ChromBPNet true performance — using manuscript Pearson r = 0.70 as placeholder in Fig 1d

---

## Publication figures — status

### Section 1: BPNet predicts p300 binding from DNA sequence

| Fig | Description | Status |
|-----|-------------|--------|
| 1a | BPNet architecture schematic (BioRender) | Not started |
| 1b | CV scatter — all 150k accessible elements | Done: `2025_0517.../predictions_mean/all_folds/mean_predictions.pdf` |
| 1c | CV scatter — p300+ elements only | Same PDF as 1b (panel) |
| 1d | Model comparison bar chart | Done: `figures/model_comparison.pdf` + subset split PDFs |
| S1a | Per-fold CV scatter (5 panels) | Done: `cv_predictions_by_fold.pdf`; check formatting |
| S1b | Training region comparison (v1/v2/v3) | Done: `figures/training_region_comparison.pdf` + subset split PDFs |

Key values (CV Pearson r):
- GATA1 BPNet: 0.597 all, 0.544 GATA1+
- p300 BPNet v1 (seq only): 0.651 all, 0.521 p300+
- ATAC ChromBPNet: ~0.70 all (manuscript placeholder; no p300+ number)
- Inter-replicate ceiling (K562): 0.876 all, 0.746 p300+

### Section 2: Multimodal model and cross-cell-type transferability

| Fig | Description | Status |
|-----|-------------|--------|
| 2a | Multimodal architecture schematic (BioRender) | Not started |
| 2b | CV scatter — K562 multimodal, all elements | Done: `2026_0529.../predictions/atac/mean_predictions.pdf` |
| 2c | CV scatter — K562 multimodal, p300+ only | Same PDF as 2b (panel) |
| 2d | Transferability bar chart (GM12878) | Done: `2026_0606_GM12878_transferability/figures/transferability_bar.pdf` + subset split PDFs |
| S2a | 3-panel scatter, all GM12878 elements | Done: `transferability_scatter_all.pdf` |
| S2b | 3-panel scatter, p300+ GM12878 elements | Done: `transferability_scatter_peaks.pdf` |
| 2e | Reverse transferability bar chart (GM12878-trained models -> K562) | Done: `2026_0606_GM12878_transferability/figures/transferability_bar_on_k562.pdf` + subset split PDFs |
| S2c | 3-panel scatter, all K562 elements (GM-trained models) | Done: `transferability_scatter_on_k562_all.pdf` |
| S2d | 3-panel scatter, p300+ K562 elements (GM-trained models) | Done: `transferability_scatter_on_k562_peaks.pdf` |

Key values (K562 CV Pearson r):
- ATAC counts correlation: 0.590 all, 0.335 p300+
- ATAC-only BPNet: 0.601 all, 0.428 p300+
- Sequence-only BPNet: 0.651 all, 0.521 p300+
- Multimodal BPNet: 0.785 all, 0.663 p300+
- Inter-replicate ceiling: 0.876 all, 0.746 p300+

Key values (GM12878 mean Pearson r):
- GM12878 BPNet (in-cell-type seq ceiling): 0.432 all, 0.328 p300+
- K562 seq-only -> GM12878: 0.277 all, 0.114 p300+
- K562 ATAC-only -> GM12878: 0.717 all, 0.467 p300+
- K562 multimodal -> GM12878: 0.793 all, 0.628 p300+
- GM12878 ATAC-only (in-cell-type): 0.683 all, 0.579 p300+
- GM12878 multimodal (in-cell-type): 0.821 all, 0.760 p300+
- Inter-replicate ceiling (GM12878): 0.881 all, 0.835 p300+

Key values, reverse direction (K562 mean Pearson r, computed 2026-07-09):
- GM12878 seq-only BPNet -> K562: 0.535 all, 0.337 p300+
- GM12878 ATAC-only BPNet -> K562: 0.597 all, 0.378 p300+
- GM12878 multimodal BPNet -> K562: 0.684 all, 0.451 p300+
- (for reference) K562 seq-only/ATAC-only/multimodal in-cell-type: 0.651/0.601/0.785 all; 0.521/0.428/0.663 p300+
- (for reference) K562 inter-replicate ceiling: 0.876 all, 0.746 p300+

---

## Completed this session (2026-06-08)

- [x] All bar chart and scatter PDFs regenerated with TrueType fonts (pdf.fonttype=42) for Illustrator compatibility
- [x] Subset bar charts split into two separate PDFs (_all_elements + _p300plus) for publication
- [x] K562 ATAC-only BPNet — trained, predicted (K562 + GM12878), all evals complete
- [x] GM12878 ATAC-only BPNet — trained, predicted, eval complete
- [x] GM12878 multimodal ATAC BPNet — trained, predicted, eval complete
- [x] GM12878 MoDISCo — 26 motifs; logos in modisco/max_seqlets_250k_30_10_0/logos/
- [x] p300 BPNet v2 + v3 CV performance computed; v1 recomputed for consistency (27,176 p300+ elements)
- [x] GATA1 BPNet performance computed
- [x] Three-figure structure: Fig 1d (k562-default), Fig 2 K562 (k562-multimodal), Fig 2d (transferability)
- [x] CV numbers corrected for multimodal + ATAC-only (were using mean predictions, now CV)
- [x] figures/figures_log.sh created — comprehensive log of all figure generation commands

---

## Completed this session (2026-06-10)

- [x] `scripts/plot_finemo_explained_fraction.py` — computes fraction of total |SHAP| explained by FiNeMo motif hits
  - Results: 8.4% overall, 11.4% p300+, 7.6% p300−; GATA=3.3%, GATA_TAL1=1.8%, AP1=1.6%
  - Outputs: `finemo/pkw_500_curated_motifs/annotated_motifs/motif_explained_fraction.pdf` (with unexplained) and `motif_explained_fraction.motifs_only.pdf` (zoomed)
- [x] `scripts/plot_activity_ep300_correlation.py` — Spearman r of every TF RPM vs EP300 RPM across all 150k elements
  - Top hits: TBL1XR1 (r=0.784), JUND (r=0.774), RCOR1 (r=0.737), SMARCA4 (r=0.695); chromatin remodelers dominate
  - Outputs: `finemo/pkw_500_curated_motifs/activity_ep300_correlation/top_corr_ep300.pdf`, `correlations.tsv.gz`
- [x] `scripts/plot_motif_ep300_correlation.py` — bar chart: each motif's hit activity correlation with EP300 RPM
  - Output: `tf_correlation/motif_ep300_correlation.pdf` (needs to be run)
- [x] `scripts/plot_finemo_composite_figure.py` — 3-motif composite panel (CTCF/AP-1/GATA) with CWM logos + top TF correlations + p300 interaction indicator (★)
  - Output: `tf_correlation/composite_motif_tf_figure.pdf` (job 28842356 running)

## Pending / next steps

- [ ] Update `P300_INTERACTORS` set in `plot_finemo_composite_figure.py` after literature/BioGRID review
- [x] Commit all changes to git
- [ ] Section 3+ figures: MoDISCo motif logos, FiNeMo hits, motif spacing/pair experiments
  - [x] Individual motif insertions violin plot (`scripts/plot_individual_motif_insertions.py`)
  - [x] Motif pair heatmaps — max log2FC + synergy (`scripts/plot_motif_pair_heatmaps.py`)

### Revisit p300+ definition across all figures
The current `EP300_peak_overlap` flag (from `finemo_peaks_all_chr.chromatin_annotations.tsv`) is based on 1000bp windows overlapping a p300 ChIP-seq peak call. ~40% of "p300+" elements by this definition have `true_logcounts = 0` — likely because the peak call overlaps only the edge of the 1000bp window, outside the 500bp BPNet input window where reads are counted. Consider replacing `EP300_peak_overlap == 1` with **top 20% of elements by observed p300 counts** (`true_logcounts >= 80th percentile`) as the p300+ definition in all figures. Figures to update if definition changes:
- `figures/model_comparison.pdf` (Fig 1d) — p300+ bars
- `figures/k562_multimodal_comparison.pdf` (Fig 2 K562) — p300+ bars
- `figures/training_region_comparison.pdf` (Fig S1b) — p300+ bars
- `2026_0606_GM12878_transferability/figures/transferability_bar.pdf` (Fig 2d) — p300+ bars
- All `*_p300plus.pdf` subset split PDFs
- Scatter plots (`mean_predictions.pdf`) — p300+ subset
