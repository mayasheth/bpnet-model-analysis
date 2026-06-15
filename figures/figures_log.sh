#!/bin/bash
# Log of commands to regenerate all publication figures.
# Run from the project root: /oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
# Environment: module load devel pixi/0.53.0 && pixi run -e ism python ...

PROJECT_DIR=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
cd $PROJECT_DIR

# ============================================================
# Fig 1d — K562 model comparison (default)
# Models: GATA1 BPNet, p300 seq-only, ATAC ChromBPNet, inter-replicate ceiling
# Script: scripts/plot_model_comparison.py --figure k562-default
# ============================================================

## [Fig 1d] Grouped by model — COMPLETE 2026-06-08
# pixi run -e ism python scripts/plot_model_comparison.py \
#     --figure k562-default \
#     --output figures/model_comparison.pdf

## [Fig 1d, alt layout] Grouped by subset — two separate PDFs — COMPLETE 2026-06-08
# pixi run -e ism python scripts/plot_model_comparison.py \
#     --figure k562-default --group-by-subset \
#     --output figures/model_comparison_by_subset.pdf
# Outputs: figures/model_comparison_by_subset_all_elements.pdf
#          figures/model_comparison_by_subset_p300plus.pdf


# ============================================================
# Fig 2 (K562) — Multimodal model comparison
# Models: ATAC counts correlation, ATAC-only BPNet, seq-only BPNet,
#         multimodal BPNet, inter-replicate ceiling
# Script: scripts/plot_model_comparison.py --figure k562-multimodal
# ============================================================

## [Fig 2 K562] Grouped by model — COMPLETE 2026-06-08
# pixi run -e ism python scripts/plot_model_comparison.py \
#     --figure k562-multimodal \
#     --output figures/k562_multimodal_comparison.pdf

## [Fig 2 K562, alt layout] Grouped by subset — two separate PDFs — COMPLETE 2026-06-08
# pixi run -e ism python scripts/plot_model_comparison.py \
#     --figure k562-multimodal --group-by-subset \
#     --output figures/k562_multimodal_comparison_by_subset.pdf
# Outputs: figures/k562_multimodal_comparison_by_subset_all_elements.pdf
#          figures/k562_multimodal_comparison_by_subset_p300plus.pdf

# Key values (K562, CV Pearson r):
#   ATAC counts correlation:  0.590 (all), 0.335 (p300+)
#   ATAC-only BPNet:          0.601 (all), 0.428 (p300+)  [CV]
#   Sequence-only BPNet:      0.651 (all), 0.521 (p300+)  [CV]
#   Multimodal BPNet:         0.785 (all), 0.663 (p300+)  [CV]
#   Inter-replicate ceiling:  0.876 (all), 0.746 (p300+)


# ============================================================
# Fig 2d — GM12878 transferability bar chart
# Models: GM12878 BPNet (ceiling), K562 seq-only, K562 multimodal,
#         GM12878 inter-replicate ceiling
# Script: scripts/plot_transferability.py
# Output: 2026_0606_GM12878_transferability/figures/
# ============================================================

## [Fig 2d] Grouped by model — COMPLETE 2026-06-08
# pixi run -e ism python scripts/plot_transferability.py \
#     --output-dir 2026_0606_GM12878_transferability/figures/

# Outputs:
#   transferability_bar.pdf              (grouped by model)
#   transferability_bar_all_elements.pdf (subset split — all elements)
#   transferability_bar_p300plus.pdf     (subset split — p300+ elements)
#   transferability_scatter_all.pdf      (Fig S2a)
#   transferability_scatter_peaks.pdf    (Fig S2b)

# Key values (GM12878, mean Pearson r):
#   GM12878 BPNet (in-cell-type ceiling):       0.432 (all), 0.328 (p300+)
#   K562 seq-only (cross-cell-type):            0.277 (all), 0.114 (p300+)
#   K562 multimodal (cross-cell-type):          0.793 (all), 0.628 (p300+)
#   GM12878 inter-replicate ceiling:            0.881 (all), 0.835 (p300+)
#   K562 ATAC-only → GM12878 (cross-cell-type): 0.717 (all), 0.467 (p300+)
#   GM12878 ATAC-only (in-cell-type):           0.683 (all), 0.579 (p300+)


# ============================================================
# Fig S1b — Training region strategy comparison (supplemental)
# Models: p300 BPNet v1, v2, v3 + inter-replicate ceiling
# Script: scripts/plot_training_region_comparison.py
# ============================================================

## [Fig S1b] Both layouts — COMPLETE 2026-06-08
# pixi run -e ism python scripts/plot_training_region_comparison.py \
#     --output-dir figures/

# Outputs:
#   training_region_comparison.pdf               (grouped by model)
#   training_region_comparison_all_elements.pdf  (subset split — all elements)
#   training_region_comparison_p300plus.pdf      (subset split — p300+ elements)

# Key values (CV Pearson r):
#   v1 (GC-matched negatives):          0.651 (all), 0.521 (p300+)
#   v2 (TF-annotation negatives):       0.621 (all), 0.464 (p300+)
#   v3 (v2 + extended training):        0.649 (all), 0.500 (p300+)
#   Inter-replicate ceiling:            0.876 (all), 0.746 (p300+)


# ============================================================
# Motif logos (SVG/PDF)
# Script: scripts/plot_motif_logos.py
# ============================================================

## K562 p300 v1 MoDISCo logos — COMPLETE 2026-06-08
# pixi run -e ism python scripts/plot_motif_logos.py \
#     --modisco-h5 2025_0517_official_EP300_K562_model/modisco/max_seqlets_250k_20_5_10/counts_scores.h5 \
#     --output-dir 2025_0517_official_EP300_K562_model/modisco/max_seqlets_250k_20_5_10/logos/ \
#     --format svg

## K562 p300 v1 FiNeMo curated motif logos — COMPLETE 2026-06-08
# pixi run -e ism python scripts/plot_motif_logos.py \
#     --finemo-dir 2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2/ \
#     --output-dir 2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2/logos/ \
#     --format svg

## GM12878 MoDISCo logos — COMPLETE 2026-06-08
# pixi run -e ism python scripts/plot_motif_logos.py \
#     --modisco-h5 2026_0606_GM12878_transferability/modisco/max_seqlets_250k_30_10_0/counts_scores.h5 \
#     --output-dir 2026_0606_GM12878_transferability/modisco/max_seqlets_250k_30_10_0/logos/ \
#     --format svg


# ============================================================
# Section 3: FiNeMo motif analysis figures — 2026-06-10
# ============================================================

# ── Explained fraction of SHAP signal by FiNeMo motifs ──────────────────────
# Script: scripts/plot_finemo_explained_fraction.py
# Submit: sbatch scripts/submit_finemo_explained_fraction.sh  (conda: analysis)
# Input:  2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs/contribution_regions.npz
#         .../annotated_motifs/hits_renamed.tsv
#         .../annotated_motifs/hits_per_peak.with_predictions.tsv.gz
# Outputs: .../annotated_motifs/motif_explained_fraction.pdf           (with unexplained bar)
#          .../annotated_motifs/motif_explained_fraction.motifs_only.pdf  (zoomed, motifs only)
#
# Key values (fraction of total |SHAP| signal):
#   Overall: 8.4% | p300+: 11.4% | p300−: 7.6%
#   GATA 3.3%, GATA_TAL1 1.8%, AP1_1 1.6%, REPEAT_G 1.0%, STAT_2 0.3%

# ── TF ChIP RPM correlation with EP300 RPM (all 522 TFs) ────────────────────
# Script: scripts/plot_activity_ep300_correlation.py
# Submit: sbatch scripts/submit_activity_ep300_correlation.sh  (conda: analysis)
# Input:  /oak/stanford/groups/engreitz/Users/agschwin/distal_regulation_paper/
#           predictors/enhancer_activity/results/bigWig/K562/enhancer_activity_features.tsv.gz
# Output: .../activity_ep300_correlation/top_corr_ep300.pdf
#         .../activity_ep300_correlation/correlations.tsv.gz
#
# Top correlates with EP300 RPM: TBL1XR1 r=0.784, JUND r=0.774, RCOR1 r=0.737,
#   DPF2 r=0.734, DNase-seq r=0.730, SMARCA4 r=0.695 — chromatin remodelers dominate

# ── Motif hit activity vs EP300 RPM (per-motif bar chart) ───────────────────
# Script: scripts/plot_motif_ep300_correlation.py
# Run: conda activate analysis && python3 scripts/plot_motif_ep300_correlation.py
# Input:  .../tf_correlation/correlations.tsv.gz
# Output: .../tf_correlation/motif_ep300_correlation.pdf
# NOTE: not yet run — run from project root

# ── Composite figure: 3-motif TF identity panel ─────────────────────────────
# Script: scripts/plot_finemo_composite_figure.py
# Submit: sbatch scripts/submit_composite_figure.sh  (conda: analysis)
# Input:  .../tf_correlation/correlations.tsv.gz
#         .../CWMs/<motif>/modisco_fc.txt
# Output: .../tf_correlation/composite_motif_tf_figure.pdf
#
# Story: CTCF (single TF, r=0.185) | AP-1 (large family, FOSL1/NFE2/MAFG) | GATA (TAL1 complex)
# P300 interactors marked with ★; update P300_INTERACTORS set in script after lit review
