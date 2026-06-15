#!/usr/bin/env python3
"""
For each FiNeMo motif, correlate per-region sum of hit_importance with
ENCODE ChIP-seq signal RPMs across the merged element universe.

Three correlation metrics are computed:
  spearman_r       - Spearman rank correlation
  pearson_r        - Pearson correlation (raw signal values)
  pearson_log1p_r  - Pearson correlation after log1p transform of both vectors

Outputs (under <finemo-dir>/tf_correlation/):
  correlations.tsv.gz          - all three r values for every (motif, signal) pair
  correlation_density.pdf      - KDE distributions, one panel per correlation type
  top10_per_motif.pdf          - logo + top-10 bar charts for all three metrics
  top50pct_hit_rpm.pdf         - RPM comparison: top-50% hit vs no-top-hit regions

Usage:
    pixi run -e ism python scripts/finemo_tf_correlation.py \\
        --finemo-dir 2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2
"""

import argparse
import os
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd
from scipy.stats import rankdata, mannwhitneyu
import pybedtools
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import logomaker

warnings.filterwarnings("ignore")

SIGNAL_FILE = (
    "/oak/stanford/groups/engreitz/Users/agschwin/distal_regulation_paper/"
    "predictors/enhancer_activity/results/bigWig/K562/enhancer_activity_features.tsv.gz"
)
EXCLUDE_MOTIFS = {"FOS_JUN", "GATA_TAL1_8BP"}
BASES = ["A", "C", "G", "T"]

NATURE_COLORS = [
    "#dc6464", "#5496ce", "#e9c54e", "#c5c500",
    "#5eb342", "#49bcbc", "#b778b3", "#f29742", "#bc9678",
]


# ── helpers ────────────────────────────────────────────────────────────────────

def chrom_key(s):
    s = str(s).replace("chr", "")
    if s == "X":       return 23
    if s == "Y":       return 24
    if s in ("M", "MT"): return 25
    try: return int(s)
    except: return 99


def sort_by_chrom(df, chrom_col="chr", start_col="start"):
    return (
        df.assign(_ck=df[chrom_col].map(chrom_key))
        .sort_values(["_ck", start_col])
        .drop(columns="_ck")
        .reset_index(drop=True)
    )


def load_cwm_restricted(finemo_dir, motif_name, cwms_arr, fwd_meta):
    """Return (L_trimmed, 4) CWM slice using pre-computed restricted window."""
    row = fwd_meta.loc[motif_name]
    arr = cwms_arr[int(row["motif_id"])].T  # (L, 4)
    return pd.DataFrame(
        arr[int(row["motif_start"]):int(row["motif_end"])],
        columns=BASES,
    )


def spearman_vec_matrix(x, Y):
    """Spearman r between vector x (N,) and each column of Y (N, K). Returns (K,)."""
    n = len(x)
    xr = rankdata(x).astype(float)
    Yr = np.empty_like(Y, dtype=float)
    for j in range(Y.shape[1]):
        Yr[:, j] = rankdata(Y[:, j])
    xr -= xr.mean()
    xs = xr.std(ddof=1)
    Yr -= Yr.mean(axis=0)
    Ys = Yr.std(axis=0, ddof=1)
    if xs == 0:
        return np.zeros(Y.shape[1])
    xr /= xs
    mask = Ys > 0
    Yr[:, mask] /= Ys[mask]
    return (xr @ Yr) / (n - 1)


def pearson_vec_matrix(x, Y):
    """Pearson r between vector x (N,) and each column of Y (N, K). Returns (K,)."""
    n = len(x)
    xc = (x - x.mean()).astype(float)
    xs = xc.std(ddof=1)
    Yc = (Y - Y.mean(axis=0)).astype(float)
    Ys = Yc.std(axis=0, ddof=1)
    if xs == 0:
        return np.zeros(Y.shape[1])
    xc /= xs
    mask = Ys > 0
    Yc[:, mask] /= Ys[mask]
    return (xc @ Yc) / (n - 1)


def style_ax(ax):
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["bottom", "left"]:
        ax.spines[spine].set_color("black")
    ax.tick_params(colors="black")
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    ax.grid(False)


def tf_display_name(col):
    return col.split(".")[0]


# ── main ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--finemo-dir",
        default="/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/"
                "2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2",
    )
    return p.parse_args()


def main():
    args = parse_args()
    finemo_dir = args.finemo_dir
    out_dir = os.path.join(finemo_dir, "tf_correlation")
    os.makedirs(out_dir, exist_ok=True)

    # ── load motif metadata ────────────────────────────────────────────────────
    cwms_arr = np.load(os.path.join(finemo_dir, "motif_cwms.npy"))   # (N, 4, L)
    meta = pd.read_table(os.path.join(finemo_dir, "motif_data.tsv"))
    fwd_meta = meta[meta["strand"] == "+"].set_index("motif_name")

    report = pd.read_table(os.path.join(finemo_dir, "motif_report.tsv"))
    report = report[~report["motif_name"].isin(EXCLUDE_MOTIFS)]
    report = report.sort_values("num_hits_total", ascending=False)
    motifs = report["motif_name"].tolist()
    motif_color = dict(zip(motifs, NATURE_COLORS))

    # ── load peaks and hits ────────────────────────────────────────────────────
    print("Loading peaks_qc.tsv ...")
    peaks = pd.read_csv(os.path.join(finemo_dir, "peaks_qc.tsv"), sep="\t")
    peaks["end"] = peaks["peak_region_start"] + 500
    peaks = sort_by_chrom(peaks, chrom_col="chr", start_col="peak_region_start")

    print("Loading hits.tsv ...")
    hits = pd.read_csv(os.path.join(finemo_dir, "hits.tsv"), sep="\t")
    hits = hits[~hits["motif_name"].isin(EXCLUDE_MOTIFS)]

    # ── aggregate hit importance per region × motif ────────────────────────────
    print("Aggregating hit importance per region ...")
    importance = (
        hits.groupby(["peak_id", "motif_name"])["hit_importance"]
        .sum()
        .unstack(fill_value=0.0)
    )
    for m in motifs:
        if m not in importance.columns:
            importance[m] = 0.0
    importance = importance[motifs]

    # ── top-50% flag per motif ─────────────────────────────────────────────────
    print("Computing top-50% hit flags ...")
    top50_flags = pd.DataFrame(False, index=importance.index, columns=motifs)
    for m in motifs:
        m_hits = hits[hits["motif_name"] == m]
        if m_hits.empty:
            continue
        threshold = m_hits["hit_coefficient_global"].quantile(0.5)
        top_ids = m_hits.loc[m_hits["hit_coefficient_global"] >= threshold, "peak_id"].unique()
        top50_flags.loc[top50_flags.index.isin(top_ids), m] = True

    # ── load signal file ───────────────────────────────────────────────────────
    print("Loading signal file ...")
    signal = pd.read_csv(SIGNAL_FILE, sep="\t", compression="gzip")
    signal = sort_by_chrom(signal, chrom_col="chr", start_col="start")
    signal_cols = [c for c in signal.columns if c not in ("chr", "start", "end")]

    # ── bedtools overlap merge ─────────────────────────────────────────────────
    print("Merging via bedtools overlap ...")
    importance_reset = importance.reset_index()
    top50_reset = top50_flags.reset_index()

    finemo_bed_df = peaks[["chr", "peak_region_start", "end", "peak_id"]].copy()
    finemo_bed_df = finemo_bed_df.rename(columns={"peak_region_start": "start"})
    finemo_bed_df = sort_by_chrom(finemo_bed_df)

    signal["_sig_idx"] = range(len(signal))
    signal_bed_df = sort_by_chrom(signal[["chr", "start", "end", "_sig_idx"]].copy())

    tmp_dir = tempfile.mkdtemp()
    pybedtools.set_tempdir(tmp_dir)

    a = pybedtools.BedTool.from_dataframe(finemo_bed_df[["chr", "start", "end", "peak_id"]])
    b = pybedtools.BedTool.from_dataframe(signal_bed_df[["chr", "start", "end", "_sig_idx"]])

    isect = a.intersect(b, wa=True, wb=True)
    mapping = isect.to_dataframe(
        names=["chr_a", "start_a", "end_a", "peak_id",
               "chr_b", "start_b", "end_b", "_sig_idx"]
    )
    mapping["_overlap"] = (
        mapping[["end_a", "end_b"]].min(axis=1)
        - mapping[["start_a", "start_b"]].max(axis=1)
    )
    mapping = (
        mapping.sort_values("_overlap", ascending=False)
        .drop_duplicates("peak_id", keep="first")
        [["peak_id", "_sig_idx"]]
    )

    merged = (
        importance_reset
        .merge(top50_reset.rename(columns={m: f"{m}_top50" for m in motifs}), on="peak_id")
        .merge(mapping, on="peak_id", how="inner")
        .merge(signal[["_sig_idx"] + signal_cols], on="_sig_idx", how="inner")
    )
    pybedtools.cleanup()

    print(f"Merged universe: {len(merged):,} regions")

    motif_vals = merged[motifs].values      # (N, M)
    signal_vals = merged[signal_cols].values   # (N, K)
    log1p_signal_vals = np.log1p(signal_vals)

    # ── correlations ──────────────────────────────────────────────────────────
    print("Computing Spearman correlations ...")
    # Pre-rank signal once for Spearman
    signal_ranks = np.empty_like(signal_vals, dtype=float)
    for j in range(signal_vals.shape[1]):
        signal_ranks[:, j] = rankdata(signal_vals[:, j])
    sig_rank_mean = signal_ranks.mean(axis=0)
    sig_rank_std = signal_ranks.std(axis=0, ddof=1)
    signal_ranks_norm = np.where(sig_rank_std > 0,
                                 (signal_ranks - sig_rank_mean) / sig_rank_std, 0.0)

    n = len(merged)
    records = []
    spearman_matrix   = np.empty((len(motifs), len(signal_cols)))
    pearson_matrix    = np.empty((len(motifs), len(signal_cols)))
    pearson_l1p_matrix = np.empty((len(motifs), len(signal_cols)))

    for i, m in enumerate(motifs):
        print(f"  {m} ...")
        x = motif_vals[:, i]

        # Spearman
        xr = rankdata(x).astype(float)
        xs = xr.std(ddof=1)
        if xs > 0:
            xr_norm = (xr - xr.mean()) / xs
            sp_vec = (xr_norm @ signal_ranks_norm) / (n - 1)
        else:
            sp_vec = np.zeros(len(signal_cols))
        spearman_matrix[i] = sp_vec

        # Pearson (raw)
        pe_vec = pearson_vec_matrix(x, signal_vals)
        pearson_matrix[i] = pe_vec

        # Pearson log1p
        pe_l1p_vec = pearson_vec_matrix(np.log1p(x), log1p_signal_vals)
        pearson_l1p_matrix[i] = pe_l1p_vec

        for j, sig in enumerate(signal_cols):
            records.append({
                "motif": m,
                "signal": sig,
                "tf": tf_display_name(sig),
                "spearman_r": sp_vec[j],
                "pearson_r": pe_vec[j],
                "pearson_log1p_r": pe_l1p_vec[j],
            })

    corr_df = pd.DataFrame(records)
    out_corr = os.path.join(out_dir, "correlations.tsv.gz")
    corr_df.to_csv(out_corr, sep="\t", index=False, compression="gzip")
    print(f"Saved {out_corr}")

    # ── correlation density plot ───────────────────────────────────────────────
    print("Plotting correlation density ...")
    metrics = [
        ("spearman_r",      spearman_matrix,    "Spearman r"),
        ("pearson_r",       pearson_matrix,     "Pearson r"),
        ("pearson_log1p_r", pearson_l1p_matrix, "Pearson r (log1p)"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=False)
    for ax, (key, mat, label) in zip(axes, metrics):
        for i, m in enumerate(motifs):
            sns.kdeplot(mat[i], ax=ax, label=m, color=motif_color[m], linewidth=1.5)
        ax.set_xlabel(label, fontsize=10)
        ax.set_ylabel("Density", fontsize=10)
        ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
        ax.legend(fontsize=6, frameon=False)
        style_ax(ax)
    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "correlation_density.pdf"), bbox_inches="tight")
    plt.close()
    print("Saved correlation_density.pdf")

    # ── top-10 per motif: logo + 3 bar charts ────────────────────────────────
    print("Plotting top-10 TFs per motif ...")
    with PdfPages(os.path.join(out_dir, "top10_per_motif.pdf")) as pdf:
        for i, m in enumerate(motifs):
            fig = plt.figure(figsize=(16, 3))
            gs = gridspec.GridSpec(1, 4, width_ratios=[1, 2, 2, 2], figure=fig, wspace=0.4)

            # logo panel
            ax_logo = fig.add_subplot(gs[0])
            cwm = load_cwm_restricted(finemo_dir, m, cwms_arr, fwd_meta)
            logomaker.Logo(cwm, ax=ax_logo, color_scheme="classic", font_name="DejaVu Sans")
            ax_logo.set_title(m, fontsize=10)
            ax_logo.axhline(0, color="black", linewidth=0.5)
            ax_logo.set_xticks([])
            ax_logo.set_ylabel("CWM weight", fontsize=8)
            style_ax(ax_logo)

            # bar chart panels
            for k, (key, mat, label) in enumerate(metrics):
                ax = fig.add_subplot(gs[k + 1])
                r_vals = mat[i]
                top_idx = np.argsort(r_vals)[::-1][:10]
                top_r = r_vals[top_idx]
                top_labels = [tf_display_name(signal_cols[j]) for j in top_idx]

                ax.barh(range(len(top_r))[::-1], top_r,
                        color=motif_color[m], edgecolor="none")
                ax.set_yticks(range(len(top_r))[::-1])
                ax.set_yticklabels(top_labels, fontsize=7)
                ax.set_xlabel(label, fontsize=8)
                ax.axvline(0, color="black", linewidth=0.8)
                ax.set_title(f"Top 10 — {label}", fontsize=8)
                style_ax(ax)

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close()

    print("Saved top10_per_motif.pdf")

    # ── top-50% hit vs no-top-hit RPM comparison ──────────────────────────────
    print("Plotting top-50% hit RPM comparison ...")
    with PdfPages(os.path.join(out_dir, "top50pct_hit_rpm.pdf")) as pdf:
        for i, m in enumerate(motifs):
            flag_col = f"{m}_top50"
            r_vals = spearman_matrix[i]
            top_idx = np.argsort(r_vals)[::-1][:8]
            top_sigs = [signal_cols[j] for j in top_idx]
            top_labels = [tf_display_name(s) for s in top_sigs]

            has_top = merged[flag_col]
            n_top = has_top.sum()
            n_rest = (~has_top).sum()

            fig, axes = plt.subplots(1, len(top_sigs),
                                     figsize=(len(top_sigs) * 1.4, 3.5), sharey=False)
            if len(top_sigs) == 1:
                axes = [axes]

            for k, (sig, label, ax) in enumerate(zip(top_sigs, top_labels, axes)):
                v_top  = np.log1p(merged.loc[has_top,  sig].values)
                v_rest = np.log1p(merged.loc[~has_top, sig].values)

                vp = ax.violinplot([v_rest, v_top], positions=[0, 1],
                                   showmedians=True, showextrema=False)
                vp["bodies"][0].set_facecolor("#c5cad7")
                vp["bodies"][1].set_facecolor(motif_color[m])
                for body in vp["bodies"]:
                    body.set_alpha(0.7)
                vp["cmedians"].set_color("black")

                _, pval = mannwhitneyu(v_top, v_rest, alternative="two-sided")
                sig_str = ("***" if pval < 0.001 else
                           "**"  if pval < 0.01  else
                           "*"   if pval < 0.05  else "ns")
                ax.set_xticks([0, 1])
                ax.set_xticklabels(["rest", "top 50%"], fontsize=7, rotation=30, ha="right")
                ax.set_title(f"{label}\n{sig_str}", fontsize=7)
                if k == 0:
                    ax.set_ylabel("log1p(RPM)", fontsize=8)
                style_ax(ax)

            fig.suptitle(
                f"{m}  |  top-50% hit regions (n={n_top:,}) vs rest (n={n_rest:,})",
                fontsize=9, y=1.02,
            )
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close()

    print("Saved top50pct_hit_rpm.pdf")
    print("Done.")


if __name__ == "__main__":
    main()
