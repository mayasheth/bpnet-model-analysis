#!/usr/bin/env python3
"""
GATA n-copy spacing analysis — n=2, 3, 4, comparing p300, GATA1, GATA2, and DNase models.

For each n, produces four PDFs (split and merged GATA1/2 variants):
  gata_n{N}_spacing_by_model_split.pdf       — one panel per model, orientation lines
  gata_n{N}_spacing_by_orientation_split.pdf — one panel per orientation, model lines
  gata_n{N}_spacing_by_model_merged.pdf      — GATA1+2 merged, one panel per model
  gata_n{N}_spacing_by_orientation_merged.pdf

Usage:
  python3 scripts/plot_gata_ncopy_spacing.py [--n 2 3 4] [--out-dir-base ...]
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
import matplotlib.pyplot as plt

BASE = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"

DATA_PATHS = {
    2: {
        "p300":  f"{BASE}/2025_0517_official_EP300_K562_model/motif_spacing/GATA_n2_50bp/raw_results.tsv",
        "GATA1": f"{BASE}/K562_GATA1_BPNet/motif_spacing/GATA_50bp_n234/raw_results.tsv",
        "GATA2": f"{BASE}/K562_GATA2_BPNet/motif_spacing/GATA_50bp_n234/raw_results.tsv",
        "DNase": f"{BASE}/K562_DNase_ChromBPNet/motif_spacing/GATA_50bp_n2/raw_results.tsv",
    },
    3: {
        "p300":  f"{BASE}/2025_0517_official_EP300_K562_model/motif_spacing/GATA_n3_50bp/raw_results.tsv",
        "GATA1": f"{BASE}/K562_GATA1_BPNet/motif_spacing/GATA_50bp_n234/raw_results.tsv",
        "GATA2": f"{BASE}/K562_GATA2_BPNet/motif_spacing/GATA_50bp_n234/raw_results.tsv",
        "DNase": f"{BASE}/K562_DNase_ChromBPNet/motif_spacing/GATA_50bp_n3/raw_results.tsv",
    },
    4: {
        "p300":  f"{BASE}/2025_0517_official_EP300_K562_model/motif_spacing/GATA_n4_50bp/raw_results.tsv",
        "GATA1": f"{BASE}/K562_GATA1_BPNet/motif_spacing/GATA_50bp_n234/raw_results.tsv",
        "GATA2": f"{BASE}/K562_GATA2_BPNet/motif_spacing/GATA_50bp_n234/raw_results.tsv",
        "DNase": f"{BASE}/K562_DNase_ChromBPNet/motif_spacing/GATA_50bp_n4/raw_results.tsv",
    },
}

OUT_DIRS = {
    2: f"{BASE}/2025_0517_official_EP300_K562_model/motif_spacing/GATA_n2_50bp/plots",
    3: f"{BASE}/2025_0517_official_EP300_K562_model/motif_spacing/GATA_n3_50bp/plots",
    4: f"{BASE}/2025_0517_official_EP300_K562_model/motif_spacing/GATA_n4_50bp/plots",
}

# Color cycle for orientations — assigned by sorted position
ORIENT_COLOR_CYCLE = ["#49bcbc", "#e96a00", "#5496ce", "#5eb342"]  # teal, orange, blue, green

MODEL_COLORS_SPLIT  = {"p300": "#b778b3", "GATA1": "#429130", "GATA2": "#5496ce", "DNase": "#e96a00"}
MODEL_LABELS_SPLIT  = {"p300": "p300", "GATA1": "GATA1", "GATA2": "GATA2", "DNase": "DNase"}
MODEL_COLORS_MERGED = {"p300": "#b778b3", "GATA TF": "#429130", "DNase": "#e96a00"}
MODEL_LABELS_MERGED = {"p300": "p300", "GATA TF": "GATA1/2", "DNase": "DNase"}


def load_data(path, n):
    df = pd.read_csv(path, sep="\t")
    if "motif_counts" in df.columns:
        df = df[df["motif_counts"] == n].copy()
    return df


def single_gata_ref(df):
    return (df["log2_fc_vs_baseline"] - df["log2_fc_vs_single"]).mean()


def orient_colors(orientations):
    return {o: ORIENT_COLOR_CYCLE[i] for i, o in enumerate(sorted(orientations))}


def plot_lines(ax, df, orient, y_col, color, label=None):
    grp = df[df["orientation_pattern"] == orient]
    if grp.empty:
        return
    spacings = sorted(grp["spacing"].unique())
    means = grp.groupby("spacing")[y_col].mean().reindex(spacings).values

    for fold in sorted(grp["model_fold"].unique()):
        fold_vals = (
            grp[grp["model_fold"] == fold]
            .set_index("spacing")[y_col]
            .reindex(spacings).values
        )
        ax.plot(spacings, fold_vals, color=color, linewidth=0.5, alpha=0.25, zorder=2)

    ax.plot(spacings, means, color=color, linewidth=1.8, label=label, zorder=3)


def style_ax(ax):
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["bottom", "left"]:
        ax.spines[spine].set_color("black")
        ax.spines[spine].set_linewidth(0.8)
    ax.tick_params(colors="black", length=3, width=0.8, labelsize=8)
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    ax.grid(False)


def option_a(data, models, model_labels, model_colors, out_path):
    """One panel per model, orientation-colored lines, y = log2_fc_vs_baseline."""
    n_panels = len(models)
    fig, axes = plt.subplots(1, n_panels, figsize=(3.5 * n_panels, 2.8), sharey=False)
    if n_panels == 1:
        axes = [axes]

    for ax, model in zip(axes, models):
        df = data[model]
        ref = single_gata_ref(df)
        orients = sorted(df["orientation_pattern"].unique())
        oc = orient_colors(orients)

        for orient in orients:
            plot_lines(ax, df, orient, "log2_fc_vs_baseline",
                       color=oc[orient], label=orient)

        ax.axhline(ref, color="black", linewidth=1.0, linestyle="--", alpha=0.5,
                   label="Single GATA")
        ax.axhline(0, color="black", linewidth=0.4, linestyle=":", alpha=0.3)

        ax.set_xticks(range(0, 51, 10))
        ax.set_xticks(range(0, 51, 5), minor=True)
        ax.set_xlim(-1, 51)
        ax.set_title(model_labels[model], fontsize=9, pad=4)
        ax.set_xlabel("Inter-motif spacing (bp)", fontsize=8)
        if ax is axes[0]:
            ax.set_ylabel("Log$_2$ fold change vs baseline", fontsize=8)
        style_ax(ax)

    axes[-1].legend(fontsize=7, frameon=False, loc="upper right")
    plt.tight_layout(pad=1.0)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def option_a_v2(data, models, model_labels, model_colors, out_path):
    """Like option_a but adds an additive expectation line at 2× single motif.

    On the log2_fc_vs_baseline axis, additive expectation for two identical
    motifs is 2 × single_ref (assuming effects sum on the log count scale).
    """
    n_panels = len(models)
    fig, axes = plt.subplots(1, n_panels, figsize=(3.5 * n_panels, 2.8), sharey=False)
    if n_panels == 1:
        axes = [axes]

    for ax, model in zip(axes, models):
        df = data[model]
        ref = single_gata_ref(df)
        orients = sorted(df["orientation_pattern"].unique())
        oc = orient_colors(orients)

        for orient in orients:
            plot_lines(ax, df, orient, "log2_fc_vs_baseline",
                       color=oc[orient], label=orient)

        ax.axhline(ref, color="black", linewidth=1.0, linestyle="--", alpha=0.55,
                   label="Single GATA")
        ax.axhline(2 * ref, color="black", linewidth=1.0, linestyle="-.", alpha=0.55,
                   label="Additive (2×)")
        ax.axhline(0, color="black", linewidth=0.4, linestyle=":", alpha=0.3)

        ax.set_xticks(range(0, 51, 10))
        ax.set_xticks(range(0, 51, 5), minor=True)
        ax.set_xlim(-1, 51)
        ax.set_title(model_labels[model], fontsize=9, pad=4)
        ax.set_xlabel("Inter-motif spacing (bp)", fontsize=8)
        if ax is axes[0]:
            ax.set_ylabel("Log$_2$ fold change vs baseline", fontsize=8)
        style_ax(ax)

    axes[-1].legend(fontsize=7, frameon=False, loc="upper right")
    plt.tight_layout(pad=1.0)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def option_bc(data, models, model_labels, model_colors, common_orients, out_path):
    """One panel per orientation (common across models), model-colored lines,
    y = log2_fc_vs_single."""
    n_panels = len(common_orients)
    fig, axes = plt.subplots(1, n_panels, figsize=(3.5 * n_panels, 2.8), sharey=True)
    if n_panels == 1:
        axes = [axes]

    for ax, orient in zip(axes, common_orients):
        for model in models:
            plot_lines(ax, data[model], orient, "log2_fc_vs_single",
                       color=model_colors[model], label=model_labels[model])

        ax.axhline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
        ax.axhline(0, color="black", linewidth=0.4, linestyle=":", alpha=0.3)

        ax.set_xticks(range(0, 51, 10))
        ax.set_xticks(range(0, 51, 5), minor=True)
        ax.set_xlim(-1, 51)
        ax.set_title(f"Orientation: {orient}", fontsize=9, pad=4)
        ax.set_xlabel("Inter-motif spacing (bp)", fontsize=8)
        if ax is axes[0]:
            ax.set_ylabel("Log$_2$ FC vs single GATA", fontsize=8)
        style_ax(ax)

    axes[-1].legend(fontsize=7, frameon=False, loc="upper right")
    plt.tight_layout(pad=1.0)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def option_bc_v2(data, models, model_labels, model_colors, common_orients, out_path):
    """One panel per orientation, model-colored lines, y = log2 synergy above additive.

    y = log2_fc_vs_single - single_ref = log2_fc_vs_baseline - 2*single_ref

    y = 0 is the additive expectation for all models, enabling direct comparison
    of superadditivity across models on a shared scale.
      y > 0 → superadditive
      y = 0 → additive (universal reference line)
      y < 0 → subadditive / interference
    """
    # Add synergy column per model (subtract each model's own single_ref)
    synergy_data = {}
    for m in models:
        df = data[m].copy()
        ref = single_gata_ref(df)
        df["log2_synergy"] = df["log2_fc_vs_single"] - ref
        synergy_data[m] = df

    n_panels = len(common_orients)
    fig, axes = plt.subplots(1, n_panels, figsize=(3.5 * n_panels, 2.8), sharey=True)
    if n_panels == 1:
        axes = [axes]

    for ax, orient in zip(axes, common_orients):
        for model in models:
            plot_lines(ax, synergy_data[model], orient, "log2_synergy",
                       color=model_colors[model], label=model_labels[model])

        ax.axhline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5,
                   zorder=0, label="Additive")

        ax.set_xticks(range(0, 51, 10))
        ax.set_xticks(range(0, 51, 5), minor=True)
        ax.set_xlim(-1, 51)
        ax.set_title(f"Orientation: {orient}", fontsize=9, pad=4)
        ax.set_xlabel("Inter-motif spacing (bp)", fontsize=8)
        if ax is axes[0]:
            ax.set_ylabel("Log$_2$ fold change above additive expectation", fontsize=8)
        style_ax(ax)

    axes[-1].legend(fontsize=7, frameon=False, loc="upper right")
    plt.tight_layout(pad=1.0)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def run_n(n, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n=== n={n} ===")

    raw = {model: load_data(path, n) for model, path in DATA_PATHS[n].items()}
    for model, df in raw.items():
        print(f"  {model}: {len(df)} rows, orientations: {sorted(df['orientation_pattern'].unique())}")

    # Common orientations across all four base models (for Option B+C)
    common_orients = sorted(
        set.intersection(*[set(df["orientation_pattern"].unique()) for df in raw.values()])
    )
    print(f"  common orientations: {common_orients}")

    # GATA1+2 merged (unique fold IDs)
    g1 = raw["GATA1"].copy()
    g2 = raw["GATA2"].copy()
    g2["model_fold"] = g2["model_fold"] + g1["model_fold"].max() + 1
    gata_merged = pd.concat([g1, g2], ignore_index=True)

    configs = [
        ("split",
         {"p300": raw["p300"], "GATA1": raw["GATA1"], "GATA2": raw["GATA2"], "DNase": raw["DNase"]},
         ["p300", "GATA1", "GATA2", "DNase"],
         MODEL_LABELS_SPLIT, MODEL_COLORS_SPLIT),
        ("merged",
         {"p300": raw["p300"], "GATA TF": gata_merged, "DNase": raw["DNase"]},
         ["p300", "GATA TF", "DNase"],
         MODEL_LABELS_MERGED, MODEL_COLORS_MERGED),
    ]

    for suffix, data, models, labels, colors in configs:
        option_a(
            data, models, labels, colors,
            os.path.join(out_dir, f"gata_n{n}_spacing_by_model_{suffix}.pdf"),
        )
        option_a_v2(
            data, models, labels, colors,
            os.path.join(out_dir, f"gata_n{n}_spacing_by_model_{suffix}_v2.pdf"),
        )
        option_bc(
            data, models, labels, colors, common_orients,
            os.path.join(out_dir, f"gata_n{n}_spacing_by_orientation_{suffix}.pdf"),
        )
        option_bc_v2(
            data, models, labels, colors, common_orients,
            os.path.join(out_dir, f"gata_n{n}_spacing_by_orientation_{suffix}_v2.pdf"),
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", nargs="+", type=int, default=[2, 3, 4],
                        help="Which n values to plot (default: 2 3 4)")
    args = parser.parse_args()

    for n in args.n:
        run_n(n, OUT_DIRS[n])


if __name__ == "__main__":
    main()
