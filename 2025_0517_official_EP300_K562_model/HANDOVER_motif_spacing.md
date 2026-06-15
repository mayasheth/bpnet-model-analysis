# GATA n-copy motif spacing — HANDOVER

Last updated: 2026-06-14

In silico experiment: insert N copies of the GATA motif (GATAA) into dinucleotide-shuffled
background sequences at varying inter-motif spacings (0–50 bp) and orientations. Compare
model predictions across p300, GATA1, GATA2, and DNase models.

---

## Data directories

| n | p300 data | GATA1/2 data | DNase data |
|---|-----------|--------------|------------|
| 2 | `motif_spacing/GATA_n2_50bp/` | `K562_GATA1_BPNet/motif_spacing/GATA_50bp_n234/` | `K562_DNase_ChromBPNet/motif_spacing/GATA_50bp_n2/` |
| 3 | `motif_spacing/GATA_n3_50bp/` | same GATA1/2 dir | `K562_DNase_ChromBPNet/motif_spacing/GATA_50bp_n3/` |
| 4 | `motif_spacing/GATA_n4_50bp/` | same GATA1/2 dir | `K562_DNase_ChromBPNet/motif_spacing/GATA_50bp_n4/` |

Each directory contains `raw_results.tsv` with columns:
`motif_seq, motif_counts, orientation_pattern, spacing, model_fold, log2_fc_vs_baseline, log2_fc_vs_single, mean_prediction, n_sequences`

- `log2_fc_vs_baseline`: log2 FC of N-motif insertion vs dinucleotide-shuffled background
- `log2_fc_vs_single`: log2 FC of N-motif insertion vs single-motif insertion

---

## Plotting script

**`scripts/plot_gata_ncopy_spacing.py`** — run directly on the login node:

```bash
module load devel pixi/0.53.0
pixi run -e ism python scripts/plot_gata_ncopy_spacing.py --n 2 3 4
```

Generates plots under `motif_spacing/GATA_n{N}_50bp/plots/` for each n.

---

## Output plots — per n-copy directory

### By-model plots (one panel per model, orientation-colored lines)

| File | y-axis | Reference lines | Status |
|------|--------|-----------------|--------|
| `gata_n{N}_spacing_by_model_split.pdf` | log2 FC vs baseline | `--` single GATA | Done |
| `gata_n{N}_spacing_by_model_split_v2.pdf` | log2 FC vs baseline | `--` single GATA + `-.` additive (2×) | Done (2026-06-14) |
| `gata_n{N}_spacing_by_model_merged.pdf` | log2 FC vs baseline | `--` single GATA | Done |
| `gata_n{N}_spacing_by_model_merged_v2.pdf` | log2 FC vs baseline | `--` single GATA + `-.` additive (2×) | Done (2026-06-14) |

Split = p300, GATA1, GATA2, DNase as separate panels. Merged = GATA1+2 averaged into one "GATA TF" panel.

### By-orientation plots (one panel per orientation, model-colored lines)

| File | y-axis | Reference lines | Status |
|------|--------|-----------------|--------|
| `gata_n{N}_spacing_by_orientation_split.pdf` | log2 FC vs single GATA | `--` at y=0 | Done |
| `gata_n{N}_spacing_by_orientation_split_v2.pdf` | log2 FC above additive expectation | `--` at y=0 = additive | Done (2026-06-14) |
| `gata_n{N}_spacing_by_orientation_merged.pdf` | log2 FC vs single GATA | `--` at y=0 | Done |
| `gata_n{N}_spacing_by_orientation_merged_v2.pdf` | log2 FC above additive expectation | `--` at y=0 = additive | Done (2026-06-14) |

---

## Key design decisions

### v2 by-model: additive expectation line
On the `log2_fc_vs_baseline` axis, additive expectation for N identical motifs = N × `single_ref`
(assuming effects sum on the log count scale). The v2 plots add this as a dash-dot line alongside
the single-motif dashed line so both reference points are visible.

### v2 by-orientation: shifted y-axis
The original `log2_fc_vs_single` axis had the additive expectation at a different y-value for each
model (`single_ref`, which varies by model), making cross-model comparison of superadditivity
ambiguous.

The v2 axis is: `log2_fc_vs_single − single_ref` = `log2_fc_vs_baseline − 2 × single_ref`

This puts y=0 at the additive expectation for all models simultaneously:
- y > 0: superadditive
- y = 0: additive (single dashed reference line, universal)
- y < 0: subadditive / interference

Y-axis label: "Log₂ fold change above additive expectation"

---

## Open questions / next steps

- [ ] Review v2 plots — assess which models show superadditivity at short spacings and how that
      varies by orientation
- [ ] Decide which version (v1 or v2) to use for the publication figure
- [ ] Consider extending to GATA + E-box or other motif pairs (see `plot_motif_spacing_focused.py`
      for a two-motif spacing script template)
- [ ] n=3 and n=4 have more orientation combinations — decide which to highlight in the paper
