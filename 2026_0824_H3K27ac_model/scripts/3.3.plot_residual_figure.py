#!/usr/bin/env python
"""Figure 5 for the report: what each model adds beyond accessibility, on the 5' target.

Replaces the earlier fragment-target figure and adds the offset-trained residual model,
which had been trained but never evaluated.

Colour encodes INPUT modality, per the project convention: blue = ATAC only,
red = sequence only, purple = sequence + ATAC. Hatching is an orthogonal channel encoding
the TRAINING OBJECTIVE: hatched = trained on `observed - atac_pred` rather than on the
total signal. So the two red bars have identical inputs and differ only in objective,
which is the comparison the figure exists to make.

Panel c reports incremental R^2 against the OBSERVED signal, which the residual metric's
built-in artifact cannot touch: because `true_resid = obs - atac_pred`, anything merely
anti-correlated with `atac_pred` scores positive residual r for free, but nothing can fake
predicting the observed signal better out of sample.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import t as tdist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nature_style import apply_rcparams, save_fig, add_panel_label, figsize

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
COLOR = {"atac": "#2166AC", "sequence": "#B2182B", "multimodal": "#762A83"}
TCRIT = tdist.ppf(0.975, df=4)

# (config key, display label, modality colour, hatched = residual objective)
BARS = [
    ("sequence5p",          "Sequence\n(total signal)",   "sequence",   False, "Seq\n(total)"),
    ("residual_sequence",   "Sequence\n(residual)",       "sequence",   True,  "Seq\n(resid.)"),
    ("multimodal5p",        "Seq + ATAC\n(total signal)", "multimodal", False, "Seq+ATAC\n(total)"),
    ("residual_multimodal", "Seq + ATAC\n(residual)",     "multimodal", True,  "Seq+ATAC\n(resid.)"),
    ("residual_atac",       "ATAC\n(residual, control)",  "atac",       True,  "ATAC\n(ctrl)"),
]


def mean_ci(v):
    v = np.asarray(v, dtype=float)
    m = v.mean()
    half = TCRIT * v.std(ddof=1) / np.sqrt(len(v))
    return m, half


def main():
    apply_rcparams()
    pf = pd.read_csv(f"{P}/results/residual_grid_per_fold.tsv", sep="\t")
    pooled = pd.read_csv(f"{P}/results/residual_grid_residual_evaluation.tsv", sep="\t")

    fig, axes = plt.subplots(1, 3, figsize=figsize(columns=2, aspect=0.42))

    # --- panel a: residual r, mean +/- 95% CI, fold points -------------------
    ax = axes[0]
    rng = np.random.default_rng(0)
    for i, (key, lab, mode, hatched, short) in enumerate(BARS):
        vals = pf.loc[pf["config"] == key, "residual_pearson"].to_numpy()
        m, half = mean_ci(vals)
        ax.bar(i, m, width=0.62, color=COLOR[mode], zorder=2, edgecolor="white",
               linewidth=0.4, hatch="///" if hatched else None)
        ax.errorbar(i, m, yerr=half, fmt="none", ecolor="black", elinewidth=0.75,
                    capsize=2, capthick=0.75, zorder=4)
        jit = (rng.random(len(vals)) - 0.5) * 0.22
        ax.scatter(np.full(len(vals), i) + jit, vals, s=4.5, color="black",
                   alpha=0.75, linewidths=0, zorder=5)
    ax.set_xticks(range(len(BARS)))
    ax.set_xticklabels([b[4].replace("\n", " ") for b in BARS], fontsize=5,
                       rotation=30, ha="right")
    ax.set_ylabel("Residual Pearson $r$")
    ax.set_title("What the model adds beyond ATAC", fontsize=8)
    ax.set_ylim(-0.12, 0.72)
    ax.axhline(0, color="black", lw=0.6, zorder=3)
    add_panel_label(ax, "a")

    # --- panel b: residual r by |true residual| quintile ---------------------
    ax = axes[1]
    for key, lab, mode, hatched, short in BARS:
        s = pooled[(pooled["config"] == key)
                   & pooled["stratum"].str.startswith("abs_resid_q")]
        if not len(s):
            continue
        ax.plot(range(1, len(s) + 1), s["residual_pearson"], marker="o", ms=3,
                lw=1.0, color=COLOR[mode], ls=(0, (3, 1.5)) if hatched else "-",
                label=lab.replace("\n", " "))
    ax.set_xlabel("|true residual| quintile\n(5 = ATAC most wrong)")
    ax.set_ylabel("Residual Pearson $r$")
    ax.set_xticks(range(1, 6))
    ax.set_title("Sequence helps where ATAC fails", fontsize=8)
    ax.legend(frameon=False, fontsize=5, loc="upper left")
    add_panel_label(ax, "b")

    # --- panel c: incremental R^2 vs observed signal (artifact-immune) -------
    ax = axes[2]
    width = 0.36
    for i, (key, lab, mode, hatched, short) in enumerate(BARS):
        for j, (col, tag) in enumerate([("incremental_r2", "all"),
                                        ("incremental_r2_topq", "top quintile")]):
            vals = pf.loc[pf["config"] == key, col].to_numpy()
            m, half = mean_ci(vals)
            x = i + (j - 0.5) * width
            ax.bar(x, m, width=width, color=COLOR[mode], zorder=2,
                   alpha=1.0 if j == 0 else 0.55, edgecolor="white", linewidth=0.4,
                   hatch="///" if hatched else None)
            ax.errorbar(x, m, yerr=half, fmt="none", ecolor="black", elinewidth=0.75,
                        capsize=1.5, capthick=0.75, zorder=4)
    ax.axhline(0, color="black", lw=0.75, zorder=3)
    ax.set_xticks(range(len(BARS)))
    ax.set_xticklabels([b[4].replace("\n", " ") for b in BARS], fontsize=5,
                       rotation=30, ha="right")
    ax.set_ylabel("Incremental $R^2$ over ATAC-only")
    ax.set_title("Gain against the observed signal", fontsize=8)
    ax.set_ylim(-0.55, 0.30)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor="#666666", label="all elements"),
                       Patch(facecolor="#666666", alpha=0.55, label="top signal quintile")],
              frameon=False, fontsize=5.5, loc="upper left", handlelength=1.2)
    add_panel_label(ax, "c")

    fig.tight_layout()
    out = save_fig(fig, f"{P}/figures/fig5_residual_grid.png")
    print("Wrote", out, "and .pdf")

    # numbers for the legend, so the caption is never hand-typed from memory
    print("\n--- values for the figure legend ---")
    for key, lab, _, _, _ in BARS:
        g = pf[pf["config"] == key]
        bits = []
        for col in ["residual_pearson", "partial_r", "r_out_vs_atac",
                    "overall_pearson", "incremental_r2", "incremental_r2_topq"]:
            m, half = mean_ci(g[col])
            bits.append(f"{col} {m:.3f} [{m-half:.3f}, {m+half:.3f}]")
        print(f"{key}: " + "; ".join(bits))


if __name__ == "__main__":
    main()
