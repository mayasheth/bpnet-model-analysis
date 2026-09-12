#!/usr/bin/env python
"""Schematic: what the shape-destroyed DNase control is for.

Drawn rather than written because the logic is a two-way branch on a counterfactual, and
the prose version takes a paragraph that readers skip. Real tracks are drawn from the actual
K562 bigwigs at one element so the shapes are the data and not an artist's impression of it.
"""
import os, sys
import numpy as np
import pyBigWig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nature_style import (apply_rcparams, save_fig, add_panel_label, figsize,
                          SCHEMATIC_EDGE, SCHEMATIC_ACCENT)

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
P = f"{D}/2026_0824_H3K27ac_model"
ATAC = f"{D}/2026_0529_multimodal_p300_model/data/atac_5p.bw"
DNASE = f"{P}/data/k562_dnase_5p.bw"
W = 250

HALF = 500
# The element is CHOSEN BY SIGNAL, not typed in. An arbitrary coordinate landed on a
# near-empty element (88 reads across 1 kb) whose tracks read as noise and taught the reader
# nothing. Picking the strongest element on a mid-size chromosome makes the shape difference
# the figure is about actually visible.
ELEMENTS = f"{D}/reference/K562_ATAC_candidate_elements.narrowPeak"
PICK_CHROM = "chr8"


def box(x, w):
    """Mean filter of width w, the same operation 0.36 applies genome-wide."""
    c = np.concatenate([[0.0], np.cumsum(x)])
    n = len(x)
    half = w // 2
    lo = np.clip(np.arange(n) - half, 0, n)
    hi = np.clip(np.arange(n) - half + w, 0, n)
    return (c[hi] - c[lo]) / np.maximum(hi - lo, 1)


def pick_element():
    """Strongest DNase element on PICK_CHROM, so the illustration is of real structure."""
    import pandas as pd
    els = pd.read_csv(ELEMENTS, sep="\t", header=None, usecols=[0, 1, 2],
                      names=["chr", "start", "end"])
    els = els[els["chr"] == PICK_CHROM].reset_index(drop=True)
    b = pyBigWig.open(DNASE)
    best, best_v = None, -1.0
    for _, r in els.iterrows():
        c = (int(r["start"]) + int(r["end"])) // 2
        if c - HALF < 0:
            continue
        # exact=True is REQUIRED. Without it pyBigWig answers from zoom-level summaries,
        # which returned 79 as the maximum over all of chr8 -- lower than an arbitrary
        # element picked by hand, which is what gave the lie away.
        v = b.stats(PICK_CHROM, c - HALF, c + HALF, type="sum", exact=True)[0] or 0.0
        if v > best_v:
            best_v, best = v, c
    b.close()
    print(f"illustrating {PICK_CHROM}:{best:,} (DNase sum {best_v:,.0f} over 1 kb)")
    return best


CTR = pick_element()


PAD = W  # smoothing margin, cropped off before display


def grab(path, pad=0):
    b = pyBigWig.open(path)
    v = np.nan_to_num(b.values(PICK_CHROM, CTR - HALF - pad, CTR + HALF + pad, numpy=True),
                      nan=0.0)
    b.close()
    return v


atac = grab(ATAC)
dnase = grab(DNASE)
# Smooth on the padded window and crop, so the box filter sees real flanking data at the
# display edges. Smoothing the 1 kb slice alone left the displayed totals 0.9% apart, which
# is an artefact of the illustration and would undercut the figure's own claim that the
# control preserves magnitude.
smooth = box(grab(DNASE, PAD), W)[PAD:PAD + 2 * HALF]

apply_rcparams()
fig = plt.figure(figsize=figsize(columns=2, aspect=0.52))
gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1], hspace=0.5, wspace=0.28)

# --- a: the tracks ------------------------------------------------------------
# Two rows, not three. The smoothed track is drawn ON TOP OF the DNase outline rather than
# beside it, because the claim being made is that the AREA is unchanged while the structure
# is gone. Side-by-side panels with a shared y axis make the smoothed track look like less
# signal, which is the opposite of what the control does; overlaying shows equal area
# directly instead of asserting it in a caption.
# BINNED AT 25 bp FOR DISPLAY. Raw 5-prime tracks put many reads on single bases, so on a
# linear axis a smoothed track's HEIGHT falls ~50x while its AREA is unchanged -- which makes
# the control look like a loss of signal, the exact opposite of what it does. Summing into
# 25 bp bins puts the two on the same footing: bin sums are directly comparable and equal
# area becomes equal height. The smoothing itself is applied at base resolution (above), as
# 0.36 does genome-wide; only the drawing is binned.
BIN = 25


def binned(v):
    return v.reshape(-1, BIN).sum(1)


x = (np.arange(-HALF, HALF, BIN) + BIN / 2)
atac_b, dnase_b, smooth_b = binned(atac), binned(dnase), binned(smooth)

ax = fig.add_subplot(gs[0, 0])
ax.fill_between(x, 0, atac_b, color="#2166AC", lw=0, alpha=0.85, step="mid")
ax.set_xlim(-HALF, HALF); ax.set_ylim(0, atac_b.max() * 1.05)
ax.set_yticks([]); ax.set_xticks([])
ax.spines[["top", "right", "left"]].set_visible(False)
ax.set_title("Real ATAC", fontsize=7, loc="left", pad=2)
ax.text(0.99, 0.9, f"{atac.sum():,.0f} reads in this window", transform=ax.transAxes,
        ha="right", va="top", fontsize=5.4, color="#555555")
add_panel_label(ax, "a")

ax = fig.add_subplot(gs[1, 0])
ax.fill_between(x, 0, dnase_b, color="#404040", lw=0, alpha=0.35, step="mid",
                label="real DNase")
ax.plot(x, smooth_b, color=SCHEMATIC_ACCENT, lw=1.6,
        label=f"smoothed {W} bp (the control)")
ax.set_xlim(-HALF, HALF); ax.set_ylim(0, max(dnase_b.max(), smooth_b.max()) * 1.05)
ax.set_yticks([])
ax.spines[["top", "right", "left"]].set_visible(False)
ax.set_xlabel("position relative to element centre (bp)", fontsize=6)
ax.set_title("Real DNase, and the control built from it", fontsize=7, loc="left", pad=2)
ax.legend(loc="upper right", frameon=False, fontsize=5.4, handlelength=1.3,
          labelspacing=0.25)
# Over a 1 kb window the two totals differ slightly because smoothing exchanges signal
# across the window edge. Conservation is a genome-wide property and 0.36 asserts it there;
# claiming exactness for one window would be wrong.
ax.text(0.01, 0.03, f"{dnase.sum():,.0f} vs {smooth.sum():,.0f} reads here\n"
                   f"(equal genome-wide; 0.36 asserts it)",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=5.4, color="#555555")

# --- b: the branch ------------------------------------------------------------
ax = fig.add_subplot(gs[:, 1])
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
add_panel_label(ax, "b")

ax.text(50, 97, "Why does real DNase beat real ATAC\nas the model's input?",
        ha="center", va="top", fontsize=7.5, color="#1a1a1a", weight="medium")
ax.text(50, 82, "Two properties differ at once, and they have\n"
                "opposite implications for the converter.",
        ha="center", va="top", fontsize=6, color="#555555")

ax.annotate("", xy=(50, 70), xytext=(50, 76),
            arrowprops=dict(arrowstyle="-|>", color=SCHEMATIC_EDGE, lw=1.1))
ax.text(50, 68, f"Feed the model DNase smoothed at {W} bp:\n"
                "magnitude kept, shape removed",
        ha="center", va="top", fontsize=6.5, color="#1a1a1a")

for x0, verdict, col, then in (
        (26, "advantage\nSURVIVES", "#8c2d04",
         "shape never mattered.\nThe converter's strength is\nworth nothing here; only\nbetter counts could help."),
        (74, "advantage\nCOLLAPSES", "#238b45",
         "shape is what DNase was\nwinning on, which is exactly\nwhat the converter is\ngood at.")):
    ax.annotate("", xy=(x0, 48), xytext=(50, 56),
                arrowprops=dict(arrowstyle="-|>", color=SCHEMATIC_EDGE, lw=1.0))
    ax.text(x0, 46, verdict, ha="center", va="top", fontsize=6.5, color=col,
            weight="medium")
    ax.text(x0, 33, then, ha="center", va="top", fontsize=5.8, color="#333333")

ax.text(50, 8, f"{W} bp is chosen, not arbitrary: it is where observed ATAC and\n"
               "DNase already agree about shape (r = 0.92 at 250 bp against 0.29\n"
               "at 1 bp), so smoothing there removes close to exactly the\n"
               "information ATAC does not already carry.",
        ha="center", va="top", fontsize=5.2, color="#666666")

out = save_fig(fig, f"{P}/figures/fig18_control_schematic")
print("wrote", out)
