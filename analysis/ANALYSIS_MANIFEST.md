# Analysis Manifest

Every analysis in this repo predates its mycelium init (2026-08-24) and lives in a
**dated or model-named directory at the repo root**, not under `analysis/`. Nothing was
moved during init, the `path:` field below is authoritative. See
`.living/conventions.md` for the layout override this records.

Cross-references: pipeline/script reference in `CLAUDE.md`, figure->script mapping in
`README.md`, forward-looking priorities in `todo/TODOLIST.md`, status snapshot in `HANDOVER.md`.

---

### p300-k562-bpnet-v1

```yaml
name: p300-k562-bpnet-v1
path: 2025_0517_official_EP300_K562_model/
status: active
created: 2025-05-17
last_updated: 2026-06-14
datasets: [encsr000ege-p300-chipseq, k562-dnase-candidate-elements, motif-compendium-human, k562-tf-chromatin-annotations]
algorithms: [bpnet-refactor, motif-exp-utils, tfmodisco-lite, finemo]
parent_analysis: null
key_findings:
  - Primary p300 model for the manuscript; 5-fold CV over ENCSR000EGE peaks vs GC-matched negatives.
  - Motif set v2 has 8 active motifs after excluding FOS_JUN and the redundant GATA_TAL1_8BP.
  - Motif pair and spacing experiments (pairwise heatmaps, GATA x E-box, GATA n-copy) largely complete.
report: null
tags: [p300, bpnet, k562, motif-syntax, shap, modisco, finemo, manuscript-fig1]
```

The primary sequence-only BPNet trained on p300 ChIP-seq (ENCSR000EGE) against GC-matched
negatives, and the source of manuscript Figs 1 and S1. Full command history is in
`scripts/0.0.log.sh`. Two subproject handovers carry the detail:
`HANDOVER_finemo_tf_analysis.md` (updated 2026-06-11) and `HANDOVER_motif_spacing.md`
(updated 2026-06-14, the most recently touched handover in the repo).

**Open methodological question:** `EP300_peak_overlap` uses a 1000 bp window while
`true_logcounts` uses BPNet's 500 bp window, so ~40% of elements called "p300+" by peak
overlap have zero observed counts. The preferred fix (top-20% by observed counts, unioned
with the overlap flag) is **not yet adopted globally**; `todo/TODOLIST.md` lists the figures that
would need regenerating.

---

### modisco-early-run

```yaml
name: modisco-early-run
path: 2025_0407_MoDISCo/
status: archived
created: 2025-04-07
last_updated: 2025-04-08
datasets: [encsr000ege-p300-chipseq]
algorithms: [tfmodisco-lite]
parent_analysis: null
key_findings:
  - Earliest MoDISCo run in the repo, superseded by the modisco/ and modisco_peaks/ outputs inside p300-k562-bpnet-v1.
report: null
tags: [modisco, p300, exploratory, superseded]
```

First-pass TF-MoDISco run (`ENCSR000EGE_trim20_flank5_0`) with separate `counts_report/`
and `profile_report/` outputs and a bare `log.txt`. Kept for provenance; current motif
discovery lives in `p300-k562-bpnet-v1`.

---

### p300-k562-bpnet-v2

```yaml
name: p300-k562-bpnet-v2
path: 2025_0703_retrain_p300_model/
status: complete
created: 2025-07-03
last_updated: 2025-10-31
datasets: [encsr000ege-p300-chipseq, k562-dnase-candidate-elements]
algorithms: [bpnet-refactor, motif-exp-utils, tfmodisco-lite, finemo]
parent_analysis: p300-k562-bpnet-v1
key_findings:
  - Retrained p300 model; full downstream stack (CV/mean predictions, SHAP, MoDISCo, FiNeMo, motif spacing) reproduced.
report: null
tags: [p300, bpnet, k562, retrain, model-comparison]
```

Retrain of the primary model with the complete downstream pipeline rerun
(`predictions_cv/`, `predictions_mean/`, `shap/`, `shap_peaks/`, `modisco/`,
`modisco_peaks/`, `finemo/`, `finemo_peaks/`, `motif_spacing/`). Canonical CV/mean Pearson
r for v1/v2/v3 is tabulated in `reference/PUBLICATION_FIGURES.md` rather than duplicated here.

---

### p300-k562-bpnet-v3

```yaml
name: p300-k562-bpnet-v3
path: 2025_1016_p300_model_v3/
status: active
created: 2025-10-16
last_updated: 2025-11-03
datasets: [encsr000ege-p300-chipseq, k562-dnase-candidate-elements]
algorithms: [bpnet-refactor, motif-exp-utils]
parent_analysis: p300-k562-bpnet-v2
key_findings:
  - Third p300 model iteration; predictions and SHAP done, no MoDISCo/FiNeMo stage yet.
report: null
tags: [p300, bpnet, k562, model-comparison]
```

Third iteration. Has `models/`, `predictions_cv/`, `predictions_mean/`, and `shap/` but no
`modisco/` or `finemo/`, the interpretability stage has not been run for v3. Performance
numbers live in `reference/PUBLICATION_FIGURES.md`.

---

### multimodal-p300-model

```yaml
name: multimodal-p300-model
path: 2026_0529_multimodal_p300_model/
status: active
created: 2026-05-29
last_updated: 2026-06-08
datasets: [encsr000ege-p300-chipseq, k562-atac-bigwig]
algorithms: [multimodal-bpnet, motif-exp-utils]
parent_analysis: p300-k562-bpnet-v1
key_findings:
  - ATAC variant trained on all 5 folds; CV Pearson 0.785 (all elements), 0.663 (p300+).
  - Adding base-resolution accessibility as an explicit input substantially outperforms sequence alone.
  - DNase variant blocked on generating data/dnase.bw.
  - SHAP on the ATAC multimodal model not yet submitted.
report: null
tags: [p300, multimodal, atac, dnase, chromatin-accessibility, pytorch]
```

A BPNet variant taking **DNA sequence + base-resolution chromatin accessibility** and
predicting stranded p300 ChIP-seq (5-channel input, middle fusion, `n_outputs=2`, PyTorch).
Architecture in `scripts/multimodal_bpnet.py`, with a full interactive write-up of every
design decision in `multimodal_bpnet_architecture.html`. `HANDOVER.md` documents the
training-bug history (7 distinct SLURM failures: generator deadlocks, device-placement
bugs, OOM from loading all genome-wide negatives before subsampling).

Note: that handover's own "next steps" still lists GM12878 multimodal predictions as
"job 28387044 submitted", that finished the same day; see `gm12878-transferability`.

---

### gm12878-transferability

```yaml
name: gm12878-transferability
path: 2026_0606_GM12878_transferability/
status: active
created: 2026-06-06
last_updated: 2026-07-09
datasets: [encsr000ege-p300-chipseq, gm12878-p300-chipseq, k562-atac-bigwig, gm12878-atac-bigwig]
algorithms: [bpnet-refactor, multimodal-bpnet]
parent_analysis: multimodal-p300-model
key_findings:
  - K562 multimodal ATAC transfers to GM12878 far better than sequence-only (Pearson 0.793 vs 0.277 on all elements; 0.628 vs 0.114 on p300+).
  - Cross-cell-type multimodal transfer (0.793) nearly matches the GM12878 in-cell-type multimodal ceiling (0.821).
  - The GM12878 sequence-only in-cell-type ceiling is low (0.432 all / 0.328 p300+), so accessibility carries most of the transferable signal.
  - Reverse direction (GM12878-trained models evaluated on K562 elements, steps 6.1-6.4) in progress as of 2026-07-09.
report: null
tags: [transferability, cross-cell-type, gm12878, k562, multimodal, atac]
```

Tests whether K562-trained p300 models generalise to GM12878, against a GM12878-trained
model as in-cell-type ceiling, and (since 2026-07-09) the reverse direction. Compares
sequence-only, ATAC-only, and multimodal variants in both cell types. Commands in
`0.0.log.sh`; results table and per-model detail in `HANDOVER.md`. Figures in
`figures/` (transferability bar and scatter panels).

---

### comparison-models-k562

```yaml
name: comparison-models-k562
path: [K562_DNase_ChromBPNet/, K562_ATAC_ChromBPNet/, K562_GATA1_BPNet/, K562_GATA2_BPNet/]
status: complete
created: null   # not encoded in directory names; earliest mtime is 2025-08
last_updated: 2026-06-08
datasets: [gata-chipseq]   # ChromBPNet accessibility inputs not catalogued; see each project log.sh
algorithms: [bpnet-refactor, motif-exp-utils]
parent_analysis: p300-k562-bpnet-v1
key_findings:
  - Accessibility (DNase, ATAC ChromBPNet) and sequence-specific TF (GATA1, GATA2) models used as reference points for what the p300 model has learned.
  - GATA1 BPNet complete as of 2026-06-08.
report: null
tags: [chrombpnet, dnase, atac, gata1, gata2, model-comparison]
```

Four comparison models grouped as one manifest entry because they serve a single purpose:
contrasting p300 coactivator binding against chromatin accessibility (DNase/ATAC
ChromBPNet) and against sequence-specific TF binding (GATA1/GATA2 BPNet). Each keeps its
own directory and scripts.

---

### fimo-motif-scan

```yaml
name: fimo-motif-scan
path: FIMO/
status: complete
created: null   # not encoded in directory name; mtime 2025-11-20
last_updated: 2025-11-20
datasets: [motif-compendium-human, k562-dnase-candidate-elements, k562-tf-chromatin-annotations]
algorithms: [fimo-memelite]
parent_analysis: p300-k562-bpnet-v1
key_findings:
  - PWM-based motif scan providing a non-deep-learning baseline for the SHAP/MoDISCo/FiNeMo motif calls.
  - Fisher's exact enrichment of motifs and motif pairs in p300+ vs p300- regions, plus spacing distributions.
report: null
tags: [fimo, pwm, motif-enrichment, baseline]
```

FIMO (via `memelite`) scan results under `elements_v1/analysis_v1/`, driven by
`scripts/run_fimo.py`, `7.0.create_region_mapping.py`, and `7.1.fimo_motif_analysis.py`.
Outputs `motif_enrichment.tsv`, `motif_pair_enrichment.tsv`, and
`spacing_distributions.tsv`, column definitions in `CLAUDE.md`.

---

### h3k27ac-k562-models

```yaml
name: h3k27ac-k562-models
path: 2026_0824_H3K27ac_model/
status: active
created: 2026-08-24
last_updated: 2026-10-14
datasets: [k562-h3k27ac-chipseq, k562-dnase-candidate-elements, k562-atac-bigwig]
algorithms: [multimodal-bpnet, motif-exp-utils]
parent_analysis: multimodal-p300-model
key_findings:
  - Accessibility dominates H3K27ac prediction. Top-quintile Pearson, 5 folds, +/-500 bp, 5-prime target - sequence 0.380 [0.323, 0.437], ATAC-only 0.548 [0.497, 0.599], sequence+ATAC 0.685 [0.667, 0.704]. These are validation-chromosome numbers (review F13). See F-001.
  - Sequence adds less for H3K27ac than for p300, so p300 is the better substrate for sequence interpretation. See F-001.
  - The residual objective helps only a sequence-blind input; joint training remains the best predictor. See F-002.
  - No H3K27ac model transfers across cell types once the floor uses the same assay as the model (six directions, three cell types). This withdraws the earlier DNase-input transfer claim and the "transfer the multimodal model" deployment recommendation. See F-010, F-016.
  - Observed p300 beats observed H3K27ac as the ABC activity term, +0.038 [+0.019, +0.058], and p300 became the project's target (libraries not depth-matched, review F19). See F-024.
  - Transferred multi-cell p300 (K562 held out) keeps about a quarter of the in-cell gain; changes to the transfer model mostly do not resolve on K562 CRISPR. See reports 5-6.
  - Corgi (deviation from a leave-one-out panel mean) beats the MSE-only multi-cell model in ABC, +0.0104 [+0.0010, +0.0204] (pre-registered PASS, holds gene-level), and is now the multi-cell base. See F-067.
report: 2026_0824_H3K27ac_model/report3_design_decisions.html
reports:   # numbered reports 1-7; h3k27ac_model_report.html is the superseded single report
  - 2026_0824_H3K27ac_model/report1_data_characterisation.html
  - 2026_0824_H3K27ac_model/report2_evaluation_methodology.html
  - 2026_0824_H3K27ac_model/report3_design_decisions.html
  - 2026_0824_H3K27ac_model/report4_atac_to_dnase_converter.html
  - 2026_0824_H3K27ac_model/report5_multicell_p300_transfer.html
  - 2026_0824_H3K27ac_model/report6_improving_transfer.html
  - 2026_0824_H3K27ac_model/report7_tf_evidence_multiome.html
tags: [h3k27ac, k562, multimodal, atac, counts-only, window-selection, chromatin]
```

Two models predicting K562 H3K27ac from DNase candidate elements: sequence-only
(`mode='sequence'`) and sequence+ATAC (`mode='multimodal'`), with ATAC-only
(`mode='atac'`) as the control that separates learned sequence syntax from a model
simply tracking the accessibility track. All three modes already exist in
`scripts/multimodal_bpnet.py`; this is a retarget, not a new architecture.

Differs from the p300 models in three ways: windows are centered on **candidate
elements** rather than ChIP peak summits (the H3K27ac summit sits on a flanking
nucleosome, not on the element); the counts window is **wider** than the element; and
the profile head is **kept but down-weighted** via `--count-loss-weight`, since the
bimodal flanking pattern is real structure worth fitting and reproducing the central
dip is a check that the model learned H3K27ac biology rather than "accessible implies
acetylated".

Two questions in order: first how well H3K27ac is predictable from sequence
(generalizably), then what motif syntax drives it.

Unstranded support landed in the shared `scripts/train_multimodal_bpnet.py` on
2026-08-24: `--signal-minus-bw` and `--accessibility-bw` are now optional and
`n_outputs` is derived (2 stranded, 1 not). The stranded p300 path was verified
bit-identical to before the change. `--max-negatives` was added in the same pass because
the old `10 x n_peaks` rule implies ~55 GB at 120k elements.

**Reports:** the single report was split into numbered reports starting 2026-09-03 (three at first, seven now, listed in `reports:` above); `h3k27ac_model_report.html` is superseded and carries a banner. **Original report (draft, 2026-08-27):** `2026_0824_H3K27ac_model/h3k27ac_model_report.html`, source `h3k27ac_model_report.qmd`. Covers window selection, ATAC-vs-H3K27ac spatial
comparison, inter-replicate ceilings in both cell types, the three-mode model
comparison, the residual-beyond-accessibility metric, the p300 comparison, and
GM12878 transfer. Scope stops before residual-model results. Rendered with the
lab `/analysis-report` skill (Quarto source, pandoc fallback on Sherlock).

**Window resolved: use +/-500** (`--out-window 1000`, `--in-window 2114`). All 30 jobs
COMPLETED; 5-fold stratified top-quintile Pearson:

| Mode | +/-500 | +/-1000 |
|---|---|---|
| sequence | 0.357 | 0.373 |
| atac | 0.543 | 0.582 |
| multimodal | **0.668** | 0.659 |

The wider window raises accessibility-only more than sequence (+0.039 vs +0.016), so the
*marginal* sequence contribution shrinks: +0.125 at +/-500 vs +0.077 at +/-1000. And
multimodal is flat-to-worse on the top quintile at +/-1000, despite looking better on the
single fold seen mid-grid. So +/-500 wins on accuracy, on sequence contribution, and on
neighbour contamination (0% vs 20%) - there is no longer a trade-off between the
prediction and syntax goals.

**Open:** prioritized next steps are in `todo/TODOLIST.md`, which is the single canonical
list for this project; the per-topic detail file it used to sit beside was folded into it on
2026-09-12. Motif-syntax analysis (SHAP/MoDISCo/FiNeMo) not started. Given F-001, expect
less sequence signal to attribute here than for p300. Also unresolved: the profile head
fits fragment-extended coverage with MNLL, which expects multinomial read counts - a
possible fundamental mismatch, not just a weighting problem. And whether counting only
the FLANKING windows (excluding the nucleosome-free center) improves the target.

---

### sce2g-tf-activity

```yaml
name: sce2g-tf-activity
path: 2026_1002_sce2g_tf_activity/
status: active
created: 2026-10-02
last_updated: 2026-10-09
datasets: [xu-k562-multiome, k562-crispr-benchmark, encsr000ege-p300-chipseq, five-cell-panel-atac-p300, ucsc-cpgislandext]
algorithms: [sce2g, tobias, scprinter-print, sce2g-replica-8.0]
parent_analysis: h3k27ac-k562-models
key_findings:
  - Element TF evidence (TOBIAS / PRINT footprints, peak-TF Kendall) tracks p300 beyond ATAC but does not move scE2G CRISPR AUPRC (F-047 to F-053).
  - Observed K562 p300 would add +0.056 to +0.060; predicted p300 (+0.001) and features chosen against p300 (+0.005) recover none of it (F-056 to F-058).
  - Features selected against CRISPR inside nested chromosome CV (CGI distance, GATA cluster, IKZF2) add about +0.021 out of sample and +0.036 [+0.022, +0.049] in the pre-registered scE2G run (F-059).
  - CGI distance tracks p300 beyond ATAC in 4 of 5 cells; its K562 gain is the far range (cell-dependent) and about half is beyond gene density (F-060 to F-062).
report: 2026_0824_H3K27ac_model/report7_tf_evidence_multiome.html
tags: [sce2g, crispr-benchmark, k562, multiome, footprints, feature-discovery, cpg-island]
```

scE2G (single-cell enhancer-gene model) evaluations of activity features on the K562 CRISPR
benchmark, using the Xu K562 multiome run. Plans in `plans/` (7: p300 feature screen; 8: CRISPR
feature discovery). Report 7 (Parts B and C) lives in the H3K27ac analysis directory with the other
numbered reports; its Part C figures come from `2026_0824_H3K27ac_model/scripts/5.6`.
