# HANDOVER — EP300_BPNet

This is the top-level handover doc: current status and context for picking up
this project, consolidated from prior dated handover notes
(`HANDOVER_2025-05-29.md`, `HANDOVER_2026-06-08.md`, `HANDOVER_2026-06-12.md`,
now removed — their content lives here). For forward-looking priorities see
`TODO.md`; for pipeline/script reference see `CLAUDE.md`; for figure-to-script
mapping see `README.md`.

Five subproject-specific handover docs remain in place with full debugging
history and are referenced, not duplicated, below.

---

## Status snapshot

| Area | Status |
|---|---|
| Git repo | Set up 2025-05-29; pushed to `git@github.com:mayasheth/bpnet-model-analysis.git` (branch `main`) |
| Manuscript figures, sections 1–2 | Done (2026-06-08), Illustrator-compatible (TrueType fonts) — see `README.md` |
| Manuscript figures, section 3+ (motif syntax) | In progress — FiNeMo/TF-correlation and motif pair/spacing analyses mostly done; open items below |
| K562 multimodal model (ATAC) | Trained, predicted, evaluated (2026-06-08); SHAP not yet run |
| K562 multimodal model (DNase) | Blocked — accessibility BigWig not yet generated |
| GM12878 cross-cell-type transferability | Complete (2026-06-08) |
| GATA1 BPNet | Complete (2026-06-08) |
| `bpnet-refactor` local edits | Vendored in `external/bpnet-refactor-patch/` (2026-07-08) — see `README.md` "External codebases" |

Canonical performance numbers (CV/mean Pearson r for all models, K562 and
GM12878) live in `TODO.md` — not repeated here to avoid drift between copies.

---

## By subproject

### 1. Primary K562 p300 BPNet v1 — `2025_0517_official_EP300_K562_model/`

- CV/mean prediction pipeline, figures for manuscript Figs 1 and S1: done.
- FiNeMo hit calling + TF-ChIP correlation analysis: see
  `2025_0517_official_EP300_K562_model/HANDOVER_finemo_tf_analysis.md`
  (updated 2026-06-11) — motif set v2 (8 active motifs after excluding
  `FOS_JUN`/redundant `GATA_TAL1_8BP`), correlation scripts, open questions
  (confirm identity of `REPEAT_G`/`NF2L_NFE`/`ELF`, review `P300_INTERACTORS`
  set in `plot_finemo_composite_figure.py` against BioGRID).
- Motif pair and spacing experiments (individual insertions, pairwise
  heatmaps, GATA×E-box spacing, GATA n-copy spacing): see
  `2025_0517_official_EP300_K562_model/HANDOVER_motif_spacing.md` (updated
  2026-06-14 — most recently updated handover in the repo). Still open:
  decide v1 vs. v2 (additive-expectation-referenced) plots for publication,
  and which n=3/n=4 orientations to highlight.
- **Open methodological question:** `EP300_peak_overlap` (1000bp window) vs.
  `true_logcounts` (500bp BPNet window) mismatch means ~40% of "p300+"
  elements by peak-overlap have zero observed counts. Preferred fix (top-20%
  by observed counts, unioned with overlap flag) is not yet adopted globally
  — see `TODO.md` for the full list of figures that would need regenerating.

### 2. Multimodal model — `2026_0529_multimodal_p300_model/`

See `2026_0529_multimodal_p300_model/HANDOVER.md` (updated 2026-06-08) for
architecture details, the full training-bug history (7 distinct SLURM
failures and fixes — generator deadlocks, device-placement bugs, OOM from
loading all genome-wide negatives before subsampling, etc.), and key file
paths.

- ATAC variant: fully trained (5 folds) and predicted; CV Pearson = 0.785
  (all elements), 0.663 (p300+) — see `TODO.md` for the up-to-date table.
- DNase variant: blocked on generating `data/dnase.bw`.
- SHAP on the ATAC multimodal model: not yet submitted.
- Note: that doc's own "next steps" list still shows GM12878 multimodal
  predictions as "job 28387044 submitted" — this was completed the same day;
  see the GM12878 transferability doc below for the finished result
  (Pearson 0.821 all / 0.760 p300+).

### 3. GM12878 cross-cell-type transferability — `2026_0606_GM12878_transferability/`

See `2026_0606_GM12878_transferability/HANDOVER.md` (updated 2026-06-08) for
the full job/failure history (conda-in-SLURM-wrap issues, SHAP job sizing,
output-directory defaults bug).

- Complete. Headline result: the K562-trained multimodal (ATAC) model
  transfers to GM12878 far better than the GM12878-trained sequence-only
  ceiling (0.793/0.628 vs. 0.432/0.328 Pearson, all/p300+), i.e. chromatin
  accessibility generalizes across cell types and more than compensates for
  the cross-cell-type training gap.
- GM12878 MoDISco run complete (26 motifs); next step noted in that doc is
  comparing GM12878 motif hits to the K562 MoDISco results — not yet done.

### 4. GATA1 BPNet — `K562_GATA1_BPNet/`

See `K562_GATA1_BPNet/HANDOVER.md` (updated 2026-06-08) for the full
debugging history, including the root-caused multiprocessing generator
deadlock in `bpnet-refactor/bpnet/generators/generators.py` (a worker crash
left `_stealer`'s unbounded `mpq.get()` blocked indefinitely, making jobs
appear to hang for hours). The sentinel-based fix and the element-level
BigWig error tolerance built for this are the origin of the
`generators.py`/`mean_predictions.py` edits now vendored in
`external/bpnet-refactor-patch/`.

- Complete: CV Pearson = 0.597 (all), 0.544 (GATA1+).

---

## Cross-cutting infrastructure notes

- **`bpnet-refactor` local edits** (bug fixes to mean-predictions, SHAP
  scoring, generators, GC-negative generation, plus new SHAP utility
  scripts used throughout every subproject's pipeline) are vendored in
  `external/bpnet-refactor-patch/` as of 2026-07-08 — see `README.md`
  "External codebases and dependencies" for how to reconstruct them.
- **pixi environments**: `default`, `ism`, `variant-scoring`, `multimodal`
  (see `pixi.toml`). `bpnet-lite` is pinned as a git dependency (not a local
  path) as of 2026-07-08.
- **Sherlock/SLURM gotchas hit during this project** (kept here since they
  aren't obvious from the code):
  - SLURM copies `sbatch` scripts to `/var/spool/slurmd/` before running
    them, so `${BASH_SOURCE[0]}`-based self-path resolution breaks — hardcode
    absolute paths in submission scripts instead.
  - Loading all genome-wide negative candidates into RAM before subsampling
    caused an OOM in multimodal training — subsample first, then extract
    windows.
  - `deeptools` cannot be installed in the `multimodal` pixi environment due
    to channel-priority conflicts with the pinned old glibc/kernel
    requirements — use `bedtools genomecov` + `bedGraphToBigWig` instead.
  - Sort steps in BigWig generation should use `-T /tmp` (node-local), not a
    temp dir on Oak — much faster and avoids Lustre metadata pressure.
  - `conda activate` inside a SLURM `--wrap` string can silently fail to take
    effect — prefer `pixi run -e <env>` for job submission where possible.

---

## Open items

See `TODO.md` for the current prioritized list. Items specifically raised in
the consolidated handover history that aren't yet in `TODO.md`:
- Compare GM12878 MoDISco motifs to K562 MoDISco motifs (noted in the
  GM12878 handover, not yet done).
- Confirm identity of `REPEAT_G`, `NF2L_NFE`, `ELF` motifs in the v2 FiNeMo
  set (noted in the FiNeMo TF-analysis handover).
