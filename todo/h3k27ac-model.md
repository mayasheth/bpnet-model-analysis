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
