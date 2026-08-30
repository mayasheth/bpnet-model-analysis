#!/usr/bin/env python
"""Compare the two paired-end 5'-end target definitions for TeloHAEC.

TeloHAEC H3K27ac is paired-end; K562 and GM12878 are single-end. A plain `genomecov -5`
on PE data counts TWO 5' ends per fragment, one at each end, which is not the quantity the
SE tracks measure. Two candidate targets:

  r1    read 1 only (-f 64): one 5' end per fragment, matching the SE convention.
  both  both mates: all reads, but signal deposited at both fragment ends.

THE COMPARISON HAS A BUILT-IN CONFOUND. `r1` has roughly half the reads by construction, so
its inter-replicate ceiling is depressed by Poisson noise regardless of whether it is the
better-defined quantity. Ranking the two on raw ceiling would favour `both` trivially. So:

  1. Meta-profile SHAPE is the primary evidence, normalised to unit mass so depth cancels.
     If counting both mates smears the bimodal H3K27ac shoulders, that is the same defect
     that disqualified the 250 bp fragment-extended track.
  2. The ceiling is reported BOTH raw and DEPTH-MATCHED: `both` counts are binomially
     thinned to `r1`'s total before the correlation, which removes the read-count advantage
     and leaves only the question of whether the quantity is more reproducible.

Ceiling correction is the project standard: sqrt(2r/(1+r)) -- Spearman-Brown for the merged
two-replicate target, then sqrt because a model predicts the expected signal, not a draw.
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import pearsonr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from nature_style import apply_rcparams
except Exception:
    apply_rcparams = None
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NP_COLS = ["chr", "start", "end", "name", "score", "strand",
           "signalValue", "pValue", "qValue", "summit"]
SEED = 0


def load_elements(path):
    df = pd.read_csv(path, sep="\t", header=None, names=NP_COLS)
    mid = (df["end"] - df["start"]) // 2
    df["center"] = df["start"] + df["summit"].where(df["summit"] >= 0, mid)
    return df


def stranded_values(paths, chrom, s, e):
    """Summed plus+minus coverage over [s,e), zeros outside data."""
    tot = None
    for p in paths:
        bw = BW[p]
        if chrom not in bw.chroms() or s < 0 or e > bw.chroms()[chrom]:
            return None
        v = np.nan_to_num(np.array(bw.values(chrom, s, e), dtype=np.float64))
        tot = v if tot is None else tot + v
    return tot


BW = {}


def open_bws(paths):
    for p in paths:
        if p not in BW:
            BW[p] = pyBigWig.open(p)


def meta_profile(paths, els, flank, binsize, n_sample, rng):
    open_bws(paths)
    idx = rng.choice(len(els), size=min(n_sample, len(els)), replace=False)
    nb = (2 * flank) // binsize
    acc = np.zeros(nb)
    n = 0
    for i in idx:
        r = els.iloc[i]
        s, e = int(r["center"]) - flank, int(r["center"]) + flank
        v = stranded_values(paths, r["chr"], s, e)
        if v is None:
            continue
        acc += v.reshape(nb, binsize).sum(axis=1)
        n += 1
    return acc / max(n, 1), n


def window_counts(paths, els, hw):
    open_bws(paths)
    out = np.full(len(els), np.nan)
    for i, r in enumerate(els.itertuples()):
        v = stranded_values(paths, r.chr, int(r.center) - hw, int(r.center) + hw)
        if v is not None:
            out[i] = v.sum()
    return out


def ceiling(c1, c2, label, rng, thin_to=None):
    """Corrected inter-replicate ceiling; optionally depth-match by binomial thinning."""
    ok = np.isfinite(c1) & np.isfinite(c2)
    a, b = c1[ok], c2[ok]
    note = ""
    if thin_to is not None:
        tot = a.sum() + b.sum()
        if tot > thin_to:
            p = thin_to / tot
            a = rng.binomial(np.round(a).astype(np.int64), p).astype(float)
            b = rng.binomial(np.round(b).astype(np.int64), p).astype(float)
            note = f"thinned p={p:.3f}"
    la, lb = np.log1p(a), np.log1p(b)
    rows = []
    for stratum, mask in [("all", np.ones(len(la), bool)),
                          ("top_quintile", (la + lb) >= np.quantile(la + lb, 0.8))]:
        r = pearsonr(la[mask], lb[mask])[0]
        sb = 2 * r / (1 + r)
        rows.append({"label": label, "stratum": stratum, "n": int(mask.sum()),
                     "raw_r": round(float(r), 4),
                     "spearman_brown": round(float(sb), 4),
                     "corrected_ceiling": round(float(np.sqrt(sb)), 4),
                     "total_counts": int(a.sum() + b.sum()), "note": note})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-dir", required=True)
    ap.add_argument("--elements", required=True)
    ap.add_argument("--condition", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--half-window", type=int, default=500)
    ap.add_argument("--flank", type=int, default=3000)
    ap.add_argument("--bin-size", type=int, default=50)
    ap.add_argument("--n-sample", type=int, default=20000)
    a = ap.parse_args()
    if apply_rcparams:
        apply_rcparams()
    rng = np.random.default_rng(SEED)
    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(a.figdir, exist_ok=True)

    els = load_elements(a.elements)
    print(f"{a.condition}: {len(els):,} elements")

    D, C = a.panel_dir, a.condition
    prof, ceil_rows = {}, []
    counts = {}
    for variant in ("r1", "both"):
        merged = [f"{D}/TeloHAEC_{C}_h3k27ac_{variant}_5p_{s}.bw" for s in ("plus", "minus")]
        reps = [[f"{D}/TeloHAEC_{C}_h3k27ac_rep{k}_{variant}_5p_{s}.bw"
                 for s in ("plus", "minus")] for k in (1, 2)]
        for p in merged + reps[0] + reps[1]:
            if not os.path.exists(p):
                raise SystemExit(f"missing {p}")
        prof[variant], nprof = meta_profile(merged, els, a.flank, a.bin_size,
                                            a.n_sample, np.random.default_rng(SEED))
        print(f"  {variant}: meta-profile over {nprof:,} elements", flush=True)
        counts[variant] = (window_counts(reps[0], els, a.half_window),
                           window_counts(reps[1], els, a.half_window))

    r1_total = float(np.nansum(counts["r1"][0]) + np.nansum(counts["r1"][1]))
    for variant in ("r1", "both"):
        c1, c2 = counts[variant]
        ceil_rows += ceiling(c1, c2, variant, rng)
        if variant == "both":
            ceil_rows += ceiling(c1, c2, "both_depthmatched", rng, thin_to=r1_total)

    cdf = pd.DataFrame(ceil_rows)
    p1 = f"{a.outdir}/telohaec_{C}_pe_variant_ceiling.tsv"
    cdf.to_csv(p1, sep="\t", index=False)
    print("\nWrote", p1)
    print(cdf.to_string(index=False))

    # profile shapes, each normalised to unit mass so depth cancels
    nb = len(prof["r1"])
    x = (np.arange(nb) - nb / 2) * a.bin_size + a.bin_size / 2
    pdf = pd.DataFrame({"distance": x,
                        "r1": prof["r1"] / prof["r1"].sum(),
                        "both": prof["both"] / prof["both"].sum()})
    p2 = f"{a.outdir}/telohaec_{C}_pe_variant_profile.tsv"
    pdf.to_csv(p2, sep="\t", index=False)
    print("Wrote", p2)

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.4))
    ax = axes[0]
    for v, col in (("r1", "#c0392b"), ("both", "#2c3e50")):
        ax.plot(x, pdf[v], color=col, lw=1.2, label=v)
    ax.set_xlabel("Distance from element center (bp)")
    ax.set_ylabel("Normalised mean 5' counts")
    ax.set_title(f"TeloHAEC {C}: PE 5' target shape", fontsize=9)
    ax.legend(frameon=False, fontsize=7)

    ax = axes[1]
    sub = cdf[cdf["stratum"] == "top_quintile"]
    ax.bar(range(len(sub)), sub["corrected_ceiling"],
           color=["#c0392b", "#2c3e50", "#7f8c8d"])
    ax.set_xticks(range(len(sub)))
    ax.set_xticklabels(sub["label"], fontsize=7, rotation=15, ha="right")
    ax.set_ylabel("Corrected ceiling (top quintile)")
    ax.set_title("Depth-matched is the fair comparison", fontsize=9)
    for ax in axes:
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.figdir}/telohaec_{C}_pe_variant.{ext}", dpi=300,
                    bbox_inches="tight")
    print(f"Wrote {a.figdir}/telohaec_{C}_pe_variant.{{pdf,png}}")


if __name__ == "__main__":
    main()
