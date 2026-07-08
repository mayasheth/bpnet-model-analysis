# EP300_BPNet

Sequence- and chromatin-informed deep learning models (BPNet / ChromBPNet) trained on p300 ChIP-seq in K562 (and GM12878) to interpret what drives p300 coactivator recruitment. We use SHAP, TF-MoDISco, and FiNeMo to identify important sequence motifs, correlate motif activity with 523 ENCODE ChIP-seq experiments to infer candidate binding TFs, and run in silico motif insertion experiments to test how motif count, spacing, and orientation combinatorially affect predicted p300 binding. A multimodal extension (DNA + ATAC-seq) tests whether chromatin accessibility improves prediction within and across cell types.

See `CLAUDE.md` for full pipeline documentation (directory layout, numbered scripts, environments, file formats) and `TODO.md` for current work status. This document maps the manuscript figures to the scripts and output files that produced them.

---

## Figure 1 — Performance of sequence-based BPNet model at predicting p300 binding

Model: `2025_0517_official_EP300_K562_model/` (p300 BPNet v1, sequence-only, 5-fold CV)

| Panel | Description | Script | Output |
|---|---|---|---|
| a | BPNet architecture schematic | — (BioRender, adapted from Pampari 2025 / Avsec 2021) | — |
| b | Predicted vs. observed p300 log counts, all accessible elements | `scripts/2.3.compute_prediction_performance.py` | `2025_0517_official_EP300_K562_model/predictions_mean/all_folds/mean_predictions.pdf` |
| c | Predicted vs. observed p300 log counts, p300+ elements only | `scripts/2.3.compute_prediction_performance.py` | same file as (b), p300+ panel |
| d | Model comparison bar chart (Pearson's r) + inter-replicate correlation | `scripts/plot_model_comparison.py --figure k562-default` | `figures/model_comparison.pdf` (+ `model_comparison_by_subset_{all_elements,p300plus}.pdf`) |

---

## Figure 2 — Sequence motifs important for p300 prediction and associated TF binding

Model: `2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2/` (FiNeMo hit calls on curated MoDISco motifs)

| Panel | Description | Script | Output |
|---|---|---|---|
| a (top) | CWM logos, ordered by descending hit frequency | `scripts/plot_finemo_motif_logos_row.py` | `finemo/pkw_500_curated_motifs_v2/motif_logos_row.pdf` |
| a (bottom) | Top-10 ChIP-seq correlates per motif (Pearson's r, log1p) out of 523 experiments, colored by p300/CBP interaction evidence | `scripts/finemo_tf_correlation.py` → `scripts/plot_finemo_composite_figure.py` | `finemo/pkw_500_curated_motifs/tf_correlation/top10_per_motif.pdf` and `composite_motif_tf_figure.pdf` |
| b | Violin plots of per-instance hit importance by motif | `scripts/plot_finemo_hit_importance_summary.py` | `finemo/pkw_500_curated_motifs_v2/motif_importance_violin.pdf` |

**Note on 2a (bottom):** `finemo_tf_correlation.py` computes the correlation and generates a bar chart per motif (`top10_per_motif.pdf`) but does not yet apply the 3-way color code described in the legend (p300/CBP interactor / TF without evidence / other nuclear protein). That 3-category coloring currently exists only in `plot_finemo_composite_figure.py`, which is scoped to 3 example motifs (CTCF, AP-1, GATA) and uses a manually curated `P300_INTERACTORS` set that still needs a literature/BioGRID pass (see `TODO.md`) and a 3rd "other nuclear protein" category. Generalizing that coloring logic to all p300-important motifs is the remaining step to produce the final panel.

---

## Figure 3 — In silico motif insertion reveals combinatorial and spacing-dependent effects

Model: `2025_0517_official_EP300_K562_model/` (motif hits from FiNeMo; insertions via `scripts/motif_exp_utils.py`)

| Panel | Description | Script | Output |
|---|---|---|---|
| a | Distribution of motif hit counts per element (top 50% by global hit coefficient), p300-bound vs. not | `scripts/plot_finemo_counts_per_peak_top50pct.py` | `finemo/pkw_500_curated_motifs_v2/plot_annotated_motifs/counts_per_peak.top50pct.top20pct_obs.pdf` |
| b | Predicted p300 counts vs. number of motif hits | same script | same file, panel 2 |
| c | Observed p300 counts vs. number of motif hits | same script | same file, panel 3 |
| d | Log₂ FC distributions from single-motif insertion into dinucleotide-shuffled sequences | `scripts/plot_individual_motif_insertions.py` (data from `scripts/6.2.motif_pairs.py`) | `motif_spacing/motif_pairs_v1/plots/individual_motif_insertions_v2.pdf` |
| e | Heatmaps: max predicted effect (left) and max log₂ synergy (right) across all pairwise motif combinations (4 orientations × 6 spacings) | `scripts/plot_motif_pair_heatmaps.py` (data from `scripts/6.2.motif_pairs.py` / `6.2.1.submit_motif_pairs.sh`) | `motif_spacing/motif_pairs_v1/plots/pairs_max_log2fc_vs_baseline_v2.pdf`, `pairs_max_synergy_v2.pdf` |
| f | GATA × E-box spacing curve (0–50 bp, 2 orientations), per-fold lines | `scripts/plot_motif_spacing_focused.py` (data from `scripts/6.1.motif_spacing.two_motifs.py`) | `motif_spacing/GATA_Ebox_50bp/spacing_focused_v1.pdf` |
| g | Two GATA motifs, 3 orientations, p300 model only, with additive-expectation reference lines | `scripts/plot_gata_ncopy_spacing.py` (n=2, p300 panel) | `motif_spacing/GATA_n2_50bp/plots/gata_n2_spacing_by_model_split_v2.pdf` (p300 panel; one panel per model, orientation lines) |
| h | Two GATA motifs, ++/+−/−+ orientations, compared across p300, GATA1, GATA2, DNase models | `scripts/plot_gata_ncopy_spacing.py` (n=2, all models) | `motif_spacing/GATA_n2_50bp/plots/gata_n2_spacing_by_orientation_split_v2.pdf` (one panel per orientation, model lines; `_split` keeps GATA1/GATA2 as separate lines rather than merging them) |

**Note on 3g/h:** `plot_gata_ncopy_spacing.py` generates both groupings (by-model, by-orientation) x (split vs. merged GATA1/2) into `GATA_n2_50bp/plots/`. Panel g pulls the p300 panel out of the by-model grid; panel h is the by-orientation, split-GATA1/2 version. Confirm against the `_v2` files (most recently regenerated, 2026-06-14).

---

## Figure 4 — Chromatin accessibility improves within- and cross-cell-type prediction

Models: `2026_0529_multimodal_p300_model/` (K562 multimodal), `2026_0606_GM12878_transferability/` (cross-cell-type)

| Panel | Description | Script | Output |
|---|---|---|---|
| a | Multimodal architecture schematic | — (BioRender); model defined in `scripts/multimodal_bpnet.py` | — |
| b | K562: seq-only vs. ATAC-only vs. multimodal accuracy, + ATAC/replicate correlation reference | `scripts/plot_model_comparison.py --figure k562-multimodal` (reference bars from `scripts/compute_atac_p300_correlation.py`, `scripts/compute_replicate_correlations.py`) | `figures/k562_multimodal_comparison.pdf` (+ `_by_subset_{all_elements,p300plus}.pdf`) |
| c | Cross-cell-type transfer to GM12878, all elements | `scripts/plot_transferability.py` | `2026_0606_GM12878_transferability/figures/transferability_bar.pdf` (+ subset split PDFs) |
| d | Cross-cell-type transfer to GM12878, p300+ elements | same script | same output, p300+ panel/PDF |

---

## Figure S1 — Performance of p300 BPNet models trained on different element sets

Models: v1/v2/v3 in `2025_0517_official_EP300_K562_model/`, `2025_0703_retrain_p300_model/`, `2025_1016_p300_model_v3/`

| Panel | Description | Script | Output |
|---|---|---|---|
| — | Bar chart: 3 training element sets (peaks+GC-neg / +all accessible / peaks+all accessible) vs. inter-replicate correlation | `scripts/plot_training_region_comparison.py` | `figures/training_region_comparison.pdf` (+ `_by_subset_{all_elements,p300plus}.pdf`) |

---

## External codebases and dependencies

Model architecture and training:
- **BPNet reference implementation** — [github.com/kundajelab/bpnet-refactor](https://github.com/kundajelab/bpnet-refactor), checked out at `/oak/stanford/groups/engreitz/Users/sheth/bpnet-refactor` with local modifications (bug fixes to mean-predictions/SHAP/generators, plus new SHAP utility scripts used throughout this project's pipeline). Those local changes are vendored in `external/bpnet-refactor-patch/` (diff + new files + upstream commit pin) so they're reconstructable without depending on that external path.
- **bpnet-lite** (multimodal BPNet base classes) — [github.com/jmschrei/bpnet-lite](https://github.com/jmschrei/bpnet-lite), used unmodified, pinned as a git dependency directly in `pixi.toml` (`[feature.multimodal.pypi-dependencies]`) at commit `fa5bcf5`, resolved and locked in `pixi.lock`. No local filesystem path dependency.

Installed as package dependencies:
- **SHAP** — [github.com/kundajelab/shap](https://github.com/kundajelab/shap) (git dependency pinned in `pixi.toml`/`pixi.lock`)
- **tangermeme** — [github.com/jmschrei/tangermeme](https://github.com/jmschrei/tangermeme), used for SHAP computation (`bpnet_37` / `tfmodisco` envs)
- **TF-MoDISco** — [github.com/jmschrei/tfmodisco-lite](https://github.com/jmschrei/tfmodisco-lite), the version used for the MoDISco runs in this project; tfmodisco-lite has since been merged into the official [github.com/kundajelab/tfmodisco](https://github.com/kundajelab/tfmodisco)
- **FiNeMo** — [github.com/austintwang/finemo_gpu](https://github.com/austintwang/finemo_gpu), the version used for hit calling in this project; the actively maintained repo is now [github.com/kundajelab/Fi-NeMo](https://github.com/kundajelab/Fi-NeMo)
- **memelite** — used for FIMO motif scanning (`tfmodisco` env); no confirmed upstream URL

---

## Reproducing figure numbers

CV Pearson's r values referenced above (e.g., 0.651 seq-only, 0.785 multimodal, 0.876 inter-replicate ceiling) are computed by `scripts/2.3.compute_prediction_performance.py` and logged per-model in each model's `predictions_mean/` or `predictions_cv/` directory; see `2025_0517_official_EP300_K562_model/scripts/0.0.log.sh` and `figures/figures_log.sh` for the exact commands used to generate the figures above.
