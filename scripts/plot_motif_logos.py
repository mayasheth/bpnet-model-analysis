#!/usr/bin/env python
"""
Plot vectorized sequence logo PDFs/SVGs for MoDISCo CWMs and/or MEME PWMs.

MoDISCo mode (--modisco-h5):
  Reads contrib_scores CWMs from a MoDISCo h5 file, trims low-importance
  flanks, and writes one PDF/SVG per pattern (fwd + rc side by side).

MEME mode (--meme):
  Reads PWMs from a MEME-format file and writes one PDF/SVG per motif.

Usage:
  # MoDISCo logos
  pixi run -e ism python scripts/plot_motif_logos.py \
      --modisco-h5 path/to/counts_scores.h5 \
      --output-dir path/to/logos/ \
      --format svg

  # MEME PWM logos (subset by name)
  pixi run -e ism python scripts/plot_motif_logos.py \
      --meme reference/MotifCompendium-Database-Human.meme.txt \
      --motif-names GATA_GATA1_1 IRF_IRF4_1 \
      --output-dir path/to/logos/ \
      --format pdf

  # FiNeMo curated motifs (motif_cwms.npy + motif_data.tsv)
  pixi run -e ism python scripts/plot_motif_logos.py \
      --finemo-dir 2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs_v2/ \
      --output-dir path/to/logos/ \
      --format svg
"""

import argparse
import os
import re

import h5py
import logomaker
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASES = ["A", "C", "G", "T"]
RC_MAP = {"A": "T", "C": "G", "G": "C", "T": "A"}


# ── Helpers ───────────────────────────────────────────────────────────────────

def reverse_complement_cwm(cwm):
    """Reverse-complement a (L, 4) CWM (columns: A C G T)."""
    return cwm[::-1, [3, 2, 1, 0]]  # reverse rows, swap A↔T and C↔G


def trim_cwm(cwm, threshold=0.2):
    """Trim low-importance flanks. Returns (start, end) indices."""
    score = np.abs(cwm).sum(axis=1)
    cutoff = score.max() * threshold
    passing = np.where(score >= cutoff)[0]
    return int(passing[0]), int(passing[-1]) + 1


def cwm_to_df(cwm):
    return pd.DataFrame(cwm, columns=BASES)


def plot_logo(ax, matrix_df, title="", color_scheme="classic"):
    """Plot a CWM or PWM as a sequence logo on ax."""
    logo = logomaker.Logo(
        matrix_df,
        ax=ax,
        color_scheme=color_scheme,
        font_name="DejaVu Sans",
    )
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)
    logo.ax.set_xticks([])
    logo.ax.set_ylabel("Contribution", fontsize=8)
    logo.ax.set_title(title, fontsize=9)
    logo.ax.tick_params(colors="black", labelsize=7)
    logo.ax.grid(False)


def save_fig(fig, output_dir, name, fmt):
    os.makedirs(output_dir, exist_ok=True)
    safe = re.sub(r"[^\w\-.]", "_", name)
    path = os.path.join(output_dir, f"{safe}.{fmt}")
    fig.savefig(path, bbox_inches="tight", format=fmt)
    plt.close(fig)
    return path


# ── MoDISCo ──────────────────────────────────────────────────────────────────

def load_modisco_patterns(h5_path, pattern_groups, include_subpatterns):
    patterns = {}
    with h5py.File(h5_path, "r") as f:
        for grp in pattern_groups:
            if grp not in f:
                continue
            for key in sorted(f[grp].keys()):
                is_sub = "subpattern" in key
                if is_sub and not include_subpatterns:
                    continue
                n_seqlets = int(f[grp][key]["seqlets"]["n_seqlets"][0])
                cwm = f[grp][key]["contrib_scores"][()]   # (L, 4)
                name = f"{grp}.{key}"
                patterns[name] = {"cwm": cwm, "n_seqlets": n_seqlets}
    return patterns


def plot_modisco_logos(h5_path, output_dir, fmt, trim_threshold,
                       pattern_groups, include_subpatterns):
    patterns = load_modisco_patterns(h5_path, pattern_groups, include_subpatterns)
    print(f"Found {len(patterns)} patterns in {os.path.basename(h5_path)}")

    for name, data in patterns.items():
        cwm = data["cwm"]
        n = data["n_seqlets"]

        start, end = trim_cwm(cwm, trim_threshold)
        cwm_fwd = cwm[start:end]
        cwm_rc  = reverse_complement_cwm(cwm)[start:end]

        fig, axes = plt.subplots(1, 2, figsize=(max(4, (end - start) * 0.35 + 1), 2.5))
        plot_logo(axes[0], cwm_to_df(cwm_fwd), title=f"{name} (+)\nn={n:,}")
        plot_logo(axes[1], cwm_to_df(cwm_rc),  title=f"{name} (rc)\nn={n:,}")
        plt.tight_layout()

        path = save_fig(fig, output_dir, name, fmt)
        print(f"  Saved: {path}")


# ── FiNeMo curated motifs ─────────────────────────────────────────────────────

def plot_finemo_logos(finemo_dir, output_dir, fmt, trim_threshold):
    cwms_path = os.path.join(finemo_dir, "motif_cwms.npy")
    meta_path = os.path.join(finemo_dir, "motif_data.tsv")

    cwms = np.load(cwms_path)          # (N, 4, L)
    meta = pd.read_table(meta_path)    # motif_id, motif_name, strand, motif_start, motif_end

    # Use pre-computed trim coordinates from motif_data.tsv
    fwd_rows = meta[meta["strand"] == "+"].copy()
    rc_rows  = meta[meta["strand"] == "-"].set_index("motif_name")

    print(f"Found {len(fwd_rows)} motifs in {os.path.basename(finemo_dir)}")

    for _, row in fwd_rows.iterrows():
        name = row["motif_name"]
        idx  = int(row["motif_id"])
        s, e = int(row["motif_start"]), int(row["motif_end"])

        cwm_fwd = cwms[idx].T[s:e]    # (4, L).T → (L, 4), then trim

        # rc: use rc row's trim coords if available, else mirror
        if name in rc_rows.index:
            rc_row = rc_rows.loc[name]
            rc_idx = int(rc_row["motif_id"])
            rs, re = int(rc_row["motif_start"]), int(rc_row["motif_end"])
            cwm_rc = cwms[rc_idx].T[rs:re]
        else:
            cwm_rc = reverse_complement_cwm(cwm_fwd)

        L = cwm_fwd.shape[0]
        fig, axes = plt.subplots(1, 2, figsize=(max(4, L * 0.35 + 1), 2.5))
        plot_logo(axes[0], cwm_to_df(cwm_fwd), title=f"{name} (+)")
        plot_logo(axes[1], cwm_to_df(cwm_rc),  title=f"{name} (rc)")
        plt.tight_layout()

        path = save_fig(fig, output_dir, name, fmt)
        print(f"  Saved: {path}")


# ── MEME PWMs ─────────────────────────────────────────────────────────────────

def parse_meme(meme_path, motif_names=None):
    """Parse MEME format; return dict of motif_name → (L, 4) probability matrix."""
    motifs = {}
    with open(meme_path) as fh:
        current, rows, collecting = None, [], False
        for line in fh:
            line = line.strip()
            if line.startswith("MOTIF"):
                parts = line.split()
                current = parts[1]
                rows, collecting = [], False
            elif line.startswith("letter-probability matrix"):
                collecting = True
            elif collecting and line and not line.startswith("URL"):
                try:
                    vals = list(map(float, line.split()))
                    if len(vals) == 4:
                        rows.append(vals)
                    else:
                        collecting = False
                        if current and rows:
                            if motif_names is None or current in motif_names:
                                motifs[current] = np.array(rows)
                        current, rows = None, []
                except ValueError:
                    collecting = False
                    if current and rows:
                        if motif_names is None or current in motif_names:
                            motifs[current] = np.array(rows)
                    current, rows = None, []
            elif not line and collecting and current and rows:
                collecting = False
                if motif_names is None or current in motif_names:
                    motifs[current] = np.array(rows)
                current, rows = None, []
        if current and rows:
            if motif_names is None or current in motif_names:
                motifs[current] = np.array(rows)
    return motifs


def plot_meme_logos(meme_path, motif_names, output_dir, fmt):
    motifs = parse_meme(meme_path, motif_names if motif_names else None)
    print(f"Found {len(motifs)} motifs in {os.path.basename(meme_path)}")

    for name, pwm in motifs.items():
        # Convert probability matrix to information content for cleaner logo
        # IC = sum_i { sum_b p(b) * log2[p(b) / 0.25] }
        with np.errstate(divide="ignore", invalid="ignore"):
            ic_per_pos = np.where(
                pwm > 0,
                pwm * np.log2(pwm / 0.25),
                0
            ).sum(axis=1, keepdims=True)
        ic_matrix = pwm * ic_per_pos

        rc_pwm = pwm[::-1, [3, 2, 1, 0]]
        rc_ic  = rc_pwm * (rc_pwm * np.where(rc_pwm > 0, np.log2(rc_pwm / 0.25), 0)
                           ).sum(axis=1, keepdims=True)

        fig, axes = plt.subplots(1, 2, figsize=(max(4, len(pwm) * 0.35 + 1), 2.5))
        plot_logo(axes[0], pd.DataFrame(ic_matrix, columns=BASES), title=f"{name} (+)")
        axes[0].set_ylabel("Information content (bits)", fontsize=8)
        plot_logo(axes[1], pd.DataFrame(rc_ic, columns=BASES), title=f"{name} (rc)")
        axes[1].set_ylabel("Information content (bits)", fontsize=8)
        plt.tight_layout()

        path = save_fig(fig, output_dir, name, fmt)
        print(f"  Saved: {path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--modisco-h5", help="MoDISCo counts_scores.h5 file")
    p.add_argument("--finemo-dir",
                   help="FiNeMo output directory containing motif_cwms.npy + motif_data.tsv")
    p.add_argument("--meme", help="MEME-format motif file (for PWM logos)")
    p.add_argument("--motif-names", nargs="+",
                   help="Subset of MEME motif names to plot (default: all)")
    p.add_argument("--output-dir", required=True, help="Directory to write logo files")
    p.add_argument("--format", default="svg", choices=["svg", "pdf"],
                   help="Output format (default: svg)")
    p.add_argument("--trim-threshold", type=float, default=0.2,
                   help="CWM trim threshold as fraction of max importance (default: 0.2)")
    p.add_argument("--pattern-groups", nargs="+",
                   default=["pos_patterns"],
                   help="MoDISCo pattern groups to include (default: pos_patterns)")
    p.add_argument("--include-subpatterns", action="store_true",
                   help="Also plot subpatterns (default: top-level patterns only)")
    return p.parse_args()


def main():
    args = parse_args()

    if not args.modisco_h5 and not args.finemo_dir and not args.meme:
        raise ValueError("Provide at least one of --modisco-h5, --finemo-dir, or --meme")

    if args.modisco_h5:
        plot_modisco_logos(
            args.modisco_h5, args.output_dir, args.format,
            args.trim_threshold, args.pattern_groups, args.include_subpatterns,
        )

    if args.finemo_dir:
        plot_finemo_logos(args.finemo_dir, args.output_dir, args.format, args.trim_threshold)

    if args.meme:
        plot_meme_logos(args.meme, args.motif_names, args.output_dir, args.format)


if __name__ == "__main__":
    main()
