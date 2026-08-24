# Algorithm Manifest

Shared code lives in **`scripts/`** at the repo root and in each analysis's own `scripts/`
subdirectory. `algorithms/` was created empty by mycelium init; nothing was moved.

---

### motif-exp-utils

```yaml
name: motif-exp-utils
implementation: scripts/motif_exp_utils.py
status: active
tags: [utilities, in-silico-motifs, prediction-wrapper]
```

The core shared module for model prediction and sequence manipulation. Key entry points:

| Function | Purpose |
|---|---|
| `get_model(model_path)` | Load a BPNet or ChromBPNet model |
| `make_model_prediction(mod_encoded, ...)` | Universal prediction wrapper, returns log counts |
| `insert_motifs_with_orientation_general()` | Insert motifs into sequences at set positions/orientations |
| `one_hot_encode(sequences, seq_length)` | DNA → one-hot |
| `dinuc_shuffle(seq)` | Dinucleotide-preserving shuffle (the baseline for log2 fold-change) |
| `generate_motif_pairs(motif_dict)` | Enumerate motif pair combinations |

Anything reporting `log2_fc_vs_baseline` or `log2_synergy` routes through here, so changes
to `dinuc_shuffle` or `make_model_prediction` invalidate those columns everywhere.

---

### multimodal-bpnet

```yaml
name: multimodal-bpnet
implementation: 2026_0529_multimodal_p300_model/scripts/multimodal_bpnet.py   # a copy also sits in scripts/ at repo root — confirm which is current before editing
trainer: 2026_0529_multimodal_p300_model/scripts/train_multimodal_bpnet.py
status: active
tags: [pytorch, multimodal, accessibility, architecture]
```

PyTorch BPNet variant: 5-channel input (4 one-hot DNA + 1 accessibility), middle fusion,
`n_outputs=2` for stranded output. Training applies log1p + z-score accessibility
normalization, reverse-complement augmentation, and 5-fold CV. Architecture rationale in
`multimodal_bpnet_architecture.html`; the 7-failure debugging history is in that project's
`HANDOVER.md`.

---

### shap-interpretation

```yaml
name: shap-interpretation
implementation: 3.1.submit_mean_shap_one_fold.sh (per-analysis scripts/)
library: tangermeme
status: active
tags: [shap, interpretability, attribution]
```

SHAP attribution over the counts head, computed separately for all candidate elements
(`shap/`) and peaks only (`shap_peaks/`). Input to MoDISCo and FiNeMo. Requires the
`bpnet_37` conda env plus `module load cuda/11.1.1 cudnn/8.1.1.33`.

---

### tfmodisco-lite

```yaml
name: tfmodisco-lite
implementation: 4.1.submit_counts_modisco.sh; extract_modisco_motifs.py
status: active
tags: [modisco, motif-discovery, cwm]
```

TF-MoDISco-lite over SHAP profiles to discover de novo motifs; `extract_modisco_motifs.py`
pulls CWMs out of the resulting HDF5. Run on both all-element and peaks-only SHAP, giving
the `modisco/` and `modisco_peaks/` output pairs. Runs in the `tfmodisco` conda env.

---

### finemo

```yaml
name: finemo
implementation: 5.1.submit_finemo.sh; 5.2.format_finemo_hits.py; 5.3.plot_finemo_hits.R; 5.4.plot_hits_vs_chip_data.R
status: active
tags: [finemo, hit-calling, motif-instances]
```

Calls motif instances against the curated motif set (`curated-motifs-finemo`, use v2).
Output columns for `finemo_hits.formatted.tsv` are documented in `CLAUDE.md`. Formatting
and downstream plots run in the `analysis` conda env.

---

### in-silico-motif-experiments

```yaml
name: in-silico-motif-experiments
implementation: 6.0.motif_spacing.one_motif.py; 6.1.motif_spacing.two_motifs.py; 6.2.motif_pairs.py; 6.2.1.submit_motif_pairs.sh
plotting: plot_motif_spacing.py; plot_motif_pairs.py
status: active
tags: [motif-syntax, spacing, orientation, synergy, causal-probe]
```

The syntax probe at the heart of the project: insert motifs into shuffled backgrounds and
read out how **count, spacing, and orientation** change predicted p300 binding. `6.0` varies
copy number of one motif, `6.1` varies spacing between two, `6.2` sweeps pairs
comprehensively. Reports `log2_fc_vs_baseline` against a dinucleotide-shuffled background
and `log2_synergy` for effects beyond the individual motifs.

**Open publication decision:** v1 vs v2 (additive-expectation-referenced) plots, and which
n=3/n=4 orientations to highlight — see `HANDOVER_motif_spacing.md`.

---

### fimo-memelite

```yaml
name: fimo-memelite
implementation: scripts/run_fimo.py; 7.0.create_region_mapping.py; 7.1.fimo_motif_analysis.py
library: memelite
status: active
tags: [fimo, pwm, enrichment, fishers-exact, baseline]
```

PWM scanning as a non-deep-learning baseline against the SHAP-derived motif calls.
`7.1` computes Fisher's exact enrichment of motifs and motif pairs in p300+ vs p300-
regions plus spacing distributions, with Benjamini-Hochberg FDR correction. Runs in the
`tfmodisco` conda env.

---

### prediction-performance

```yaml
name: prediction-performance
implementation: 2.3.compute_prediction_performance.py
status: active
tags: [evaluation, pearson, spearman]
```

Computes Pearson and Spearman correlation of predicted vs observed log counts, reported
both over all candidate elements and restricted to p300+ elements. Canonical numbers for
every model live in `TODO.md`.

**Caveat:** the p300+ subset depends on `EP300_peak_overlap` (1000 bp window) while
`true_logcounts` uses BPNet's 500 bp window — so ~40% of "p300+" elements have zero
observed counts, deflating p300+ correlations. Unresolved; see
`analysis/ANALYSIS_MANIFEST.md`.
