# Publication figures: status and key values

Moved out of the root `TODO.md` on 2026-09-12 when that file was folded into
`todo/TODOLIST.md`. This is a record rather than open work, which is why it does not
live in the TODO list. Figure status was last edited 2026-06-08; the reverse-direction
transfer numbers were added 2026-07-09. Treat the statuses as of those dates and the
Pearson values as the p300-model results they were computed from, NOT as the H3K27ac
numbers that dominate the later work.

Open items extracted from here and now tracked in `todo/TODOLIST.md`: the two BioRender
architecture schematics (Figs 1a, 2a), the ATAC ChromBPNet placeholder in Fig 1d, and
the p300+ definition question.

## Publication figures, status

### Section 1: BPNet predicts p300 binding from DNA sequence

| Fig | Description | Status |
|-----|-------------|--------|
| 1a | BPNet architecture schematic (BioRender) | Not started |
| 1b | CV scatter, all 150k accessible elements | Done: `2025_0517.../predictions_mean/all_folds/mean_predictions.pdf` |
| 1c | CV scatter, p300+ elements only | Same PDF as 1b (panel) |
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
| 2b | CV scatter, K562 multimodal, all elements | Done: `2026_0529.../predictions/atac/mean_predictions.pdf` |
| 2c | CV scatter, K562 multimodal, p300+ only | Same PDF as 2b (panel) |
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
