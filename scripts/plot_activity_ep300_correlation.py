#!/usr/bin/env python3
"""
Correlate every TF ChIP-seq signal (RPM) with EP300 RPM across all K562 elements.

Input:  enhancer_activity_features.tsv.gz  (regions × TF signals)
Output: activity_ep300_correlation/
          correlations.tsv.gz       - Spearman r for all TFs vs EP300
          top_corr_ep300.pdf        - horizontal bar chart, top/bottom N TFs

Usage:
  python3 scripts/plot_activity_ep300_correlation.py
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import rankdata
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt

SIGNAL_FILE = (
    "/oak/stanford/groups/engreitz/Users/agschwin/distal_regulation_paper"
    "/predictors/enhancer_activity/results/bigWig/K562/enhancer_activity_features.tsv.gz"
)
EP300_COL = "EP300.ENCFF513GVS"
OUT_DIR = "2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs/activity_ep300_correlation"
os.makedirs(OUT_DIR, exist_ok=True)

TOP_N = 30   # TFs to show in bar chart (top and bottom)


def tf_display_name(col):
    return col.split(".")[0]


# ── load ───────────────────────────────────────────────────────────────────────

print("Loading signal file ...", flush=True)
signal = pd.read_csv(SIGNAL_FILE, sep="\t", compression="gzip")
signal_cols = [c for c in signal.columns if c not in ("chr", "start", "end")]

if EP300_COL not in signal_cols:
    raise ValueError(f"{EP300_COL} not found in signal file")

other_cols = [c for c in signal_cols if c != EP300_COL]
print(f"  {len(signal):,} regions × {len(signal_cols)} signals", flush=True)

# ── spearman correlations ─────────────────────────────────────────────────────

print("Computing Spearman correlations ...", flush=True)

ep300_vals = signal[EP300_COL].values
ep300_rank = rankdata(ep300_vals).astype(float)
ep300_rank -= ep300_rank.mean()
ep300_std = ep300_rank.std(ddof=1)
ep300_rank /= ep300_std

n = len(ep300_vals)
records = []

# Process in chunks to avoid OOM
chunk_size = 50
for start in range(0, len(other_cols), chunk_size):
    chunk = other_cols[start : start + chunk_size]
    mat = signal[chunk].values.astype(float)
    mat_rank = np.empty_like(mat)
    for j in range(mat.shape[1]):
        r = rankdata(mat[:, j]).astype(float)
        r -= r.mean()
        s = r.std(ddof=1)
        mat_rank[:, j] = r / s if s > 0 else 0.0
    r_vec = (ep300_rank @ mat_rank) / (n - 1)
    for col, r in zip(chunk, r_vec):
        records.append({"signal": col, "tf": tf_display_name(col), "spearman_r": float(r)})
    if (start // chunk_size) % 2 == 0:
        print(f"  {start + len(chunk)}/{len(other_cols)}", flush=True)

corr_df = pd.DataFrame(records).sort_values("spearman_r", ascending=False).reset_index(drop=True)
corr_df.to_csv(os.path.join(OUT_DIR, "correlations.tsv.gz"), sep="\t", index=False, compression="gzip")
print(f"Saved correlations.tsv.gz  ({len(corr_df)} TFs)", flush=True)

# ── plot: top + bottom N ───────────────────────────────────────────────────────

print("Plotting ...", flush=True)

top    = corr_df.head(TOP_N)
bottom = corr_df.tail(TOP_N).iloc[::-1]
plot_df = pd.concat([top, bottom], ignore_index=True)

# separator between top and bottom
sep_pos = TOP_N  # gap between groups

fig, ax = plt.subplots(figsize=(5, 9))

colors = ["#5496ce" if r >= 0 else "#c5373d" for r in plot_df["spearman_r"]]
y = np.arange(len(plot_df))

# shift bottom group down to add visual gap
y_plot = np.where(y >= sep_pos, y + 1.5, y)
ax.barh(y_plot, plot_df["spearman_r"], color=colors, height=0.7)

ax.set_yticks(y_plot)
ax.set_yticklabels(plot_df["tf"], fontsize=7)
ax.set_xlabel("Spearman r with p300 ChIP RPM", fontsize=9)
ax.axvline(0, color="black", linewidth=0.8)
ax.invert_yaxis()

# divider line
mid_y = (y_plot[sep_pos - 1] + y_plot[sep_pos]) / 2
ax.axhline(mid_y, color="#888888", linewidth=0.8, linestyle="--")

for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
for spine in ["bottom", "left"]:
    ax.spines[spine].set_color("black")
ax.tick_params(colors="black")
ax.grid(False)

plt.tight_layout()
out = os.path.join(OUT_DIR, "top_corr_ep300.pdf")
fig.savefig(out, bbox_inches="tight")
plt.close()
print(f"Saved {out}")

# also print top 20 to stdout
print(f"\nTop 20 TFs correlated with EP300 RPM:")
print(corr_df.head(20)[["tf", "spearman_r"]].to_string(index=False))
