#!/usr/bin/env python
"""
Compute prediction performance for BPNet (stranded) models on candidate elements.

Reads per-fold h5 prediction files (output of bpnet-predict), correctly aggregates
stranded log counts to total log counts as:
    log1p(expm1(lc_plus) + expm1(lc_minus))
rather than the incorrect sum log1p(lc_plus) + log1p(lc_minus).

Outputs (matching the multimodal model's output format):
  <cv-output-dir>/cv_predictions.tsv.gz
      chrom, start, end, true_logcounts, pred_logcounts, <overlap_col>, fold
  <cv-output-dir>/cv_predictions.pdf
  <cv-output-dir>/cv_predictions_by_fold.pdf
  <cv-output-dir>/prediction_accuracy.tsv
  <mean-output-dir>/mean_predictions.tsv.gz  (if --mean-pred-dir given)
      chrom, start, end, true_logcounts, mean_pred_logcounts, <overlap_col>, pred_fold0..4
  <mean-output-dir>/mean_predictions.pdf

Usage:
  conda activate tfmodisco
  python scripts/2.3.compute_prediction_performance.py \\
    --cv-pred-dir 2025_0517.../predictions_cv \\
    --mean-pred-dir 2025_0517.../predictions_mean \\
    --chromatin-annot /path/to/chromatin_annotations.tsv \\
    --overlap-col EP300_peak_overlap
"""

import argparse
import os

import h5py
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cv-pred-dir", required=True,
                   help="Dir with fold{0-4}/<h5-name> (one held-out chr set per fold)")
    p.add_argument("--mean-pred-dir",
                   help="Dir with fold{0-4}/<h5-name> (all 150k elements per fold; optional)")
    p.add_argument("--h5-name", default="ENCSR000EGE_split000_predictions.h5",
                   help="H5 filename inside each fold dir (default: ENCSR000EGE_split000_predictions.h5)")
    p.add_argument("--peaks",
                   help="narrowPeak file to compute peak overlap against. "
                        "Not required when --from-tsv is set.")
    p.add_argument("--overlap-col", default="EP300_peak_overlap",
                   help="Name for the peak overlap indicator column in output (default: EP300_peak_overlap)")
    p.add_argument("--cv-output-dir",
                   help="Output dir for CV tables/plots "
                        "(default: <cv-pred-dir>/all_folds)")
    p.add_argument("--mean-output-dir",
                   help="Output dir for mean tables/plots "
                        "(default: <mean-pred-dir>/all_folds)")
    p.add_argument("--from-tsv", action="store_true",
                   help="Skip h5 loading; read existing cv_predictions.tsv.gz and "
                        "mean_predictions.tsv.gz from the output dirs and regenerate plots only.")
    p.add_argument("--max-counts", type=float, default=10,
                   help="Upper axis limit for scatter plots (default: 10)")
    return p.parse_args()


def stranded_to_total(lc):
    """(N, 2) per-strand log1p counts -> (N,) total log1p counts.

    Correct formula: log1p(expm1(lc_plus) + expm1(lc_minus))
    This inverts the log1p on each strand, sums raw counts, then re-applies log1p.
    """
    return np.log1p(np.expm1(np.clip(lc[:, 0], 0, None)) +
                    np.expm1(np.clip(lc[:, 1], 0, None)))


def load_h5(path):
    with h5py.File(path, "r") as f:
        chroms = f["coords/coords_chrom"][()].astype(str)
        starts = f["coords/coords_start"][()]
        ends   = f["coords/coords_end"][()]
        true_lc = f["predictions/true_logcounts"][()]   # (N, 2)
        pred_lc = f["predictions/pred_logcounts"][()]   # (N, 2)
    return chroms, starts, ends, true_lc, pred_lc


def load_peaks(path):
    return pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                       names=["chrom", "start", "end"])


def compute_peak_overlap(coords_df, peaks_df):
    """Return int8 array: 1 if output window overlaps any peak, 0 otherwise."""
    overlap = np.zeros(len(coords_df), dtype=np.int8)
    for chrom, chrom_peaks in peaks_df.groupby("chrom"):
        r_mask = coords_df["chrom"].values == chrom
        if r_mask.sum() == 0:
            continue
        r_idx = np.where(r_mask)[0]
        r_starts = coords_df["start"].values[r_idx]
        r_ends   = coords_df["end"].values[r_idx]
        p_starts = chrom_peaks["start"].values
        p_ends   = chrom_peaks["end"].values
        has_overlap = ((r_starts[:, None] < p_ends[None, :]) &
                       (r_ends[:, None]   > p_starts[None, :])).any(axis=1)
        overlap[r_idx] = has_overlap.astype(np.int8)
    return overlap


def plot_scatter(df, x_col, y_col, ax, color="#792374",
                 xlabel="Observed log counts", ylabel="Predicted log counts", title="",
                 max_counts=10):
    x, y = df[x_col].values, df[y_col].values
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    cmap = mcolors.LinearSegmentedColormap.from_list("d", ["#ffffff", color])
    sns.kdeplot(x=x, y=y, ax=ax, cmap=cmap, fill=True)
    lo, hi = 0, max_counts
    ax.plot([lo, hi], [lo, hi], color="black", ls=(0, (5, 5)), lw=1.2, label="y = x")
    m, b = np.polyfit(x, y, 1)
    ax.plot([lo, hi], [m * lo + b, m * hi + b], color=color, lw=1.5, label="fit")
    r, _   = pearsonr(x, y)
    rho, _ = spearmanr(x, y)
    ax.text(0.04, 0.96, f"Pearson r = {r:.3f}\nSpearman ρ = {rho:.3f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="white", alpha=0.7))
    ax.set(xlim=(lo, hi), ylim=(lo, hi), aspect="equal",
           xlabel=xlabel, ylabel=ylabel)
    ax.set_title(title, fontsize=11)
    ax.tick_params(colors="black"); ax.grid(False)
    for sp in ax.spines.values(): sp.set_edgecolor("black")
    ax.legend(loc="lower right", frameon=False, fontsize=9)


def acc_row(df, x_col, y_col, label):
    x, y = df[x_col].values, df[y_col].values
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    r, _   = pearsonr(x, y)
    rho, _ = spearmanr(x, y)
    return {"subset": label, "n": int(ok.sum()),
            "pearson_r": round(r, 4), "spearman_rho": round(rho, 4),
            "mse": round(float(np.mean((y - x) ** 2)), 4)}


def make_cv_plots(cv_df, oc, cv_out, mc):
    peaks_cv    = cv_df[cv_df[oc] == 1]
    nonpeaks_cv = cv_df[cv_df[oc] == 0]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    plot_scatter(cv_df,       "true_logcounts", "pred_logcounts", axes[0], color="#006eae",
                 title=f"CV — all elements (n={len(cv_df):,})", max_counts=mc)
    plot_scatter(peaks_cv,    "true_logcounts", "pred_logcounts", axes[1], color="#006eae",
                 title=f"CV — p300+ (n={len(peaks_cv):,})", max_counts=mc)
    plot_scatter(nonpeaks_cv, "true_logcounts", "pred_logcounts", axes[2], color="#006eae",
                 title=f"CV — p300- (n={len(nonpeaks_cv):,})", max_counts=mc)
    sns.despine(trim=True); plt.tight_layout()
    plt.savefig(os.path.join(cv_out, "cv_predictions.pdf")); plt.close()

    present_folds = sorted(cv_df["fold"].unique())
    fig, axes = plt.subplots(1, len(present_folds), figsize=(5 * len(present_folds), 5))
    if len(present_folds) == 1:
        axes = [axes]
    for ax, i in zip(axes, present_folds):
        sub = cv_df[cv_df["fold"] == i]
        plot_scatter(sub, "true_logcounts", "pred_logcounts", ax, color="#006eae",
                     title=f"CV fold {i} (n={len(sub):,})", max_counts=mc)
    sns.despine(trim=True); plt.tight_layout()
    plt.savefig(os.path.join(cv_out, "cv_predictions_by_fold.pdf")); plt.close()

    rows = [
        acc_row(cv_df,       "true_logcounts", "pred_logcounts", "CV — all elements"),
        acc_row(peaks_cv,    "true_logcounts", "pred_logcounts", "CV — p300+"),
        acc_row(nonpeaks_cv, "true_logcounts", "pred_logcounts", "CV — p300-"),
    ]
    for i in present_folds:
        sub = cv_df[cv_df["fold"] == i]
        rows.append(acc_row(sub, "true_logcounts", "pred_logcounts", f"CV — fold {i}"))
    return rows


def make_mean_plots(mean_df, oc, mean_out, mc):
    peaks_m    = mean_df[mean_df[oc] == 1]
    nonpeaks_m = mean_df[mean_df[oc] == 0]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    plot_scatter(mean_df,    "true_logcounts", "mean_pred_logcounts", axes[0],
                 title=f"Mean predictions — all elements (n={len(mean_df):,})", max_counts=mc)
    plot_scatter(peaks_m,    "true_logcounts", "mean_pred_logcounts", axes[1],
                 title=f"Mean predictions — p300+ (n={len(peaks_m):,})", max_counts=mc)
    plot_scatter(nonpeaks_m, "true_logcounts", "mean_pred_logcounts", axes[2],
                 title=f"Mean predictions — p300- (n={len(nonpeaks_m):,})", max_counts=mc)
    sns.despine(trim=True); plt.tight_layout()
    plt.savefig(os.path.join(mean_out, "mean_predictions.pdf")); plt.close()

    return [
        acc_row(mean_df,    "true_logcounts", "mean_pred_logcounts", "mean — all elements"),
        acc_row(peaks_m,    "true_logcounts", "mean_pred_logcounts", "mean — p300+"),
        acc_row(nonpeaks_m, "true_logcounts", "mean_pred_logcounts", "mean — p300-"),
    ]


def main():
    args = parse_args()
    cv_out = args.cv_output_dir or os.path.join(args.cv_pred_dir, "all_folds")
    os.makedirs(cv_out, exist_ok=True)
    oc = args.overlap_col
    mc = args.max_counts
    acc_rows = []

    if args.from_tsv:
        # ── Fast replot mode: read existing TSVs, skip h5 loading ─────────────
        cv_tsv = os.path.join(cv_out, "cv_predictions.tsv.gz")
        if not os.path.exists(cv_tsv):
            raise FileNotFoundError(f"--from-tsv set but {cv_tsv} not found")
        print(f"Reading {cv_tsv}")
        cv_df = pd.read_csv(cv_tsv, sep="\t")
        acc_rows += make_cv_plots(cv_df, oc, cv_out, mc)

        if args.mean_pred_dir:
            mean_out = args.mean_output_dir or os.path.join(args.mean_pred_dir, "all_folds")
            os.makedirs(mean_out, exist_ok=True)
            mean_tsv = os.path.join(mean_out, "mean_predictions.tsv.gz")
            if not os.path.exists(mean_tsv):
                raise FileNotFoundError(f"--from-tsv set but {mean_tsv} not found")
            print(f"Reading {mean_tsv}")
            mean_df = pd.read_csv(mean_tsv, sep="\t")
            acc_rows += make_mean_plots(mean_df, oc, mean_out, mc)
    else:
        # ── Full mode: load from h5 files ─────────────────────────────────────
        if not args.peaks:
            raise ValueError("--peaks is required unless --from-tsv is set")
        peaks_df = load_peaks(args.peaks)

        print("Building CV predictions table...")
        fold_dfs = []
        for fold in range(5):
            path = os.path.join(args.cv_pred_dir, f"fold{fold}", args.h5_name)
            if not os.path.exists(path):
                print(f"  WARNING: {path} not found, skipping"); continue
            try:
                chroms, starts, ends, true_lc, pred_lc = load_h5(path)
            except OSError as e:
                print(f"  WARNING: fold {fold} h5 unreadable ({e}), skipping"); continue
            df = pd.DataFrame({
                "chrom": chroms, "start": starts, "end": ends,
                "true_logcounts": stranded_to_total(true_lc),
                "pred_logcounts": stranded_to_total(pred_lc),
                "fold": fold,
            })
            df[oc] = compute_peak_overlap(df, peaks_df)
            fold_dfs.append(df)
            print(f"  Fold {fold}: {len(df)} regions  (p300+ = {df[oc].sum()})")

        cv_df = pd.concat(fold_dfs, ignore_index=True)
        cv_df.to_csv(os.path.join(cv_out, "cv_predictions.tsv.gz"),
                     sep="\t", index=False, compression="gzip")
        print(f"Saved cv_predictions.tsv.gz ({len(cv_df)} total regions)")
        acc_rows += make_cv_plots(cv_df, oc, cv_out, mc)

        if args.mean_pred_dir:
            mean_out = args.mean_output_dir or os.path.join(args.mean_pred_dir, "all_folds")
            os.makedirs(mean_out, exist_ok=True)
            print("Building mean predictions table...")

            # Load each fold, sort by (chrom, start) to align across folds before averaging.
            # The 5 mean h5 files may have elements in different orders depending on how
            # bpnet-predict processed the loci file.
            fold_dfs_mean = {}
            for fold in range(5):
                path = os.path.join(args.mean_pred_dir, f"fold{fold}", args.h5_name)
                if not os.path.exists(path):
                    print(f"  WARNING: {path} not found"); continue
                try:
                    chroms, starts, ends, true_lc, pred_lc = load_h5(path)
                except OSError as e:
                    print(f"  WARNING: fold {fold} h5 unreadable ({e}), skipping"); continue
                df_f = pd.DataFrame({
                    "chrom": chroms, "start": starts, "end": ends,
                    "true_logcounts": stranded_to_total(true_lc),
                    "pred": stranded_to_total(pred_lc),
                }).sort_values(["chrom", "start"]).reset_index(drop=True)
                fold_dfs_mean[fold] = df_f
                print(f"  Fold {fold}: {len(df_f)} regions")

            present_mean_folds = sorted(fold_dfs_mean)
            ref_df = fold_dfs_mean[present_mean_folds[0]]
            mean_df = ref_df[["chrom", "start", "end", "true_logcounts"]].copy()
            for fold in present_mean_folds:
                mean_df[f"pred_fold{fold}"] = fold_dfs_mean[fold]["pred"].values
            pred_cols = [f"pred_fold{f}" for f in present_mean_folds]
            mean_df["mean_pred_logcounts"] = mean_df[pred_cols].mean(axis=1)
            mean_df[oc] = compute_peak_overlap(mean_df, peaks_df)

            mean_df.to_csv(os.path.join(mean_out, "mean_predictions.tsv.gz"),
                           sep="\t", index=False, compression="gzip")
            print(f"Saved mean_predictions.tsv.gz ({len(mean_df)} regions)")
            acc_rows += make_mean_plots(mean_df, oc, mean_out, mc)

    # ── Accuracy table ────────────────────────────────────────────────────────
    acc_df = pd.DataFrame(acc_rows)
    print("\nPrediction accuracy:")
    print(acc_df.to_string(index=False))
    acc_df.to_csv(os.path.join(cv_out, "prediction_accuracy.tsv"), sep="\t", index=False)
    print(f"\nSaved prediction_accuracy.tsv to {cv_out}")


if __name__ == "__main__":
    main()
