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

---

## H3K27ac model — next steps (added 2026-08-25)

Detail and results live in `2026_0824_H3K27ac_model/`; see
`analysis/ANALYSIS_MANIFEST.md` entry `h3k27ac-k562-models` and
`.living/learnings.md`. Baseline for all compute estimates below: `in_window=2114`,
~105 s/epoch, 25–33 epochs, peak RSS 13.8 GB, 25–65 min per job on one GPU.

**Where we are.** Top-quintile Pearson at ±500 bp, 5 folds: sequence 0.357,
ATAC-only 0.543, sequence+ATAC 0.668, against a ceiling of 0.930. Sequence adds
+0.125 over accessibility — versus +0.301 for p300 evaluated identically.

**What the profile comparison showed** (`figures/profile_comparison.png`): ATAC is a
sharp single peak at the element center (offset −12 bp, no dip). H3K27ac is bimodal with
shoulders at ±250 bp and a shallow dip (center/peak 0.96), and is *much broader* — at
±2000 bp it is still at 0.28 of its maximum where ATAC has fallen to 0.17. The two
H3K27ac tracks (250 bp fragment vs 5′ end) overlay almost exactly, so the fragment
extension was **not** distorting the shape; it only shallowed the dip slightly
(0.959 → 0.944). **The dominant ATAC/H3K27ac difference is breadth, not displacement** —
which reprioritizes the items below toward receptive field and away from special
handling of spatial offset.

### 1. Residual / difference-from-ATAC framing — DO FIRST
Cheapest, and closest to the actual goal: a model that takes ATAC + sequence and
predicts activity, where the *departure from ATAC expectation* is the quantity of
interest.

- Stage 1 is already trained (`models/atac_hw500_clw1000`, 5 folds).
- Stage 2: train the sequence model on `observed − atac_pred` (held-out per fold).
- **Compute: ~5 jobs at baseline ≈ 4 GPU-hours.** No new data, no code change to the
  architecture.

New evaluation to build alongside (`2.4.evaluate_residual.py`), because overall
correlation is the wrong headline for this question:
- **Residual correlation** — `r(observed − atac_pred, model_pred − atac_pred)`. Does the
  model's departure from the ATAC expectation track the true departure? This is the
  direct measure and should become the headline number.
- **Incremental R²** — `R²(multimodal) − R²(atac_only)`: variance explained beyond ATAC.
- **Stratify by |true residual|** — elements where H3K27ac most deviates from its
  ATAC expectation (accessible-but-unacetylated, and acetylated-beyond-accessibility) are
  the biologically interesting set. Report accuracy specifically there, and export the
  extremes as a candidate list.
- Use the ATAC-only **model prediction** as the baseline, not the raw ATAC track, so the
  residual is what accessibility genuinely cannot explain.

### 2. Cross-cell-type transfer to GM12878 — DO SECOND
Separates "model capacity" from "H3K27ac is not sequence-determinable". The ceiling we
computed is technical reproducibility; it says nothing about how much is predictable from
sequence in principle. A sequence model that transfers is learning generalizable rules.

- GM12878 H3K27ac is available: `ENCFF645BAL`, `ENCFF865OOP` (single-end) in
  `$OAK/Users/sheth/Data/ENCODE/GM12878/`, already `.filtered.sorted.bam`.
- GM12878 ATAC (`2026_0606_GM12878_transferability/data/atac.bw`) and elements
  (`reference/GM12878_candidate_elements.narrowPeak`) already exist.
- **Compute: ~15 min CPU to build 5′ BigWigs + ~10 min GPU for inference ≈ 0.5 GPU-hours.**
  Evaluating existing K562 models on GM12878 needs no training at all.
- Optional in-cell-type ceiling (train GM12878 models): +5–15 jobs, ~4–12 GPU-hours.

### 3. Wider receptive field (`--n-layers 10`) — best single architecture bet
Currently `n_layers=8` → `trimming=557` → the model sees ~1.1 kb. But H3K27ac around top
elements only reaches its plateau at ±4000 bp, and that plateau is ~8× the low-decile
level, meaning strong elements sit inside acetylation domains kilobases wide. The
breadth finding above is direct evidence the model cannot see the context that sets its
target.

- `n_layers=10` → `trimming=2093`, receptive field ~4.2 kb, so `in_window` becomes
  **5186** for a ±500 count window. Already a plumbed flag; only the window changes.
- **Compute: extraction RSS scales with `in_window` (2.45×) → ~30–35 GB, still inside the
  120 GB request. Trunk FLOPs scale ~L × n_layers ≈ 3.1×, so ~325 s/epoch → 2–4 h per
  job. 3 modes × 5 folds at ±500 only ≈ 30–60 GPU-hours.** The most expensive item here.

### 4. Predict at nucleosome resolution, not base pair
H3K27ac has no meaningful 1 bp structure; it is nucleosomal. Fitting a 1 bp profile is
fitting noise, and is why the profile term dominated the loss and needed
`count_loss_weight=1000`. Bin the profile target to ~50–150 bp so MNLL's multinomial
assumption applies at the scale where the signal is real. Pairs naturally with the 5′
switch: 5′ gives honest read counts, binning puts them at the right scale.

- **Compute: negligible — trunk unchanged, only the output head and loss shrink (very
  slightly faster). ~5–10 validation jobs ≈ 4–8 GPU-hours.** The cost is development:
  changes to both `multimodal_bpnet.py` and `train_multimodal_bpnet.py`, in code shared
  with the p300 models, so it needs the same backward-compatibility care as the
  unstranded change.

### 5. Fragment-size-stratified ATAC channels — BLOCKED ON DATA
H3K27ac requires a nucleosome to acetylate, and the model currently gets flat ATAC
coverage with no nucleosome-positioning information. Sub-nucleosomal (<100 bp) and
mono-nucleosomal (180–247 bp) ATAC channels would supply it.

- **Blocker:** the K562 ATAC source is Tn5-shifted **tagAlign** (`ENCFF077FBI`,
  `ENCFF128WZG`, `ENCFF534DCE`) — per-read entries (~94 bp), so fragment lengths are not
  recoverable. Needs the original paired-end ATAC BAM from ENCODE.
- **Compute once data exists: trivial.** `MultiModalBPNet` already accepts arbitrary
  accessibility channels via `n_acc_filters`, so this is a data-prep change, not an
  architecture change. Extraction memory grows ~20% per added channel; training cost
  rises only in the first conv layer.

### Also still open
- Recompute the ceiling and re-sweep `count_loss_weight` on the 5′ target — expect the
  weight to fall a long way from 1000 now that MNLL sees real read counts. Note the 5′
  track is sparse (max 27, ~6 reads per kb window at background), so the ceiling could go
  either way; if it drops materially, test unextended 36 bp read coverage as a middle
  option.
- Re-check the window choice on the 5′ target: less smearing means less neighbour
  bleed-through, so the contamination penalty at ±1000 may look different.
- Motif syntax (SHAP/MoDISCo/FiNeMo) not started. Per F-001, expect less sequence signal
  to attribute here than for p300.
