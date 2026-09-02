#!/usr/bin/env python3
"""Inter-replicate ceiling on PROFILE SHAPE as a function of bin size.

WHY. `0.3.replicate_ceiling_by_window.py` bounds how well any model can predict COUNTS.
Nothing bounds how well any model can predict the profile SHAPE, and the profile head is
trained at 1 bp resolution with an MNLL term that carries real gradient
(`loss = profile_loss + count_loss_weight * count_loss`). If two replicates of the same
experiment cannot agree on the 1 bp profile, that term is fitting Poisson noise and is
spending trunk capacity for nothing. This measures the ceiling directly, with no training,
and the bin size at which it saturates is the resolution worth predicting -- rather than
picking 20-50 bp by intuition.

METHOD. For a sample of elements, take the +/-500 bp window (the training out_window) from
each replicate, bin at increasing widths, and correlate the two replicates' binned profiles
*within* each element across positions. Correlating within an element removes the count
signal, which `0.3` already covers, so what is left is shape alone. Reported as the mean
over elements, split by signal quintile because low-count elements are pure noise and would
drag any average toward zero.

Also reports the ceiling on the SUMMED counts of the same windows as an internal check
against `0.3`, and does the analysis per strand as well as strand-summed: 5' ends are
strand-asymmetric around a nucleosome, so the stranded task the model is actually given can
have a different ceiling from the unstranded shape.

Usage:
  0.25.profile_ceiling_by_binsize.py --rep1-plus ... --rep1-minus ... \
      --rep2-plus ... --rep2-minus ... --elements ... --label k562 [--n 20000]
"""
import argparse, json, os
import numpy as np
import pandas as pd
import pyBigWig

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
BINS = [1, 5, 10, 25, 50, 100, 250]
# 500 bp is excluded on purpose: a 1000 bp window gives only 2 bins, below the
# 3-bin floor in shape_corr, so a within-element correlation is not defined.
HW = 500

ap = argparse.ArgumentParser()
ap.add_argument("--rep1-plus", required=True); ap.add_argument("--rep1-minus", required=True)
ap.add_argument("--rep2-plus", required=True); ap.add_argument("--rep2-minus", required=True)
ap.add_argument("--elements", required=True); ap.add_argument("--label", required=True)
ap.add_argument("--n", type=int, default=20000)
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()

els = pd.read_csv(a.elements, sep="\t", header=None,
                  names=["chr", "start", "end", "name", "score", "strand",
                         "signal", "pval", "qval", "summit"])
rng = np.random.RandomState(a.seed)
idx = rng.choice(len(els), size=min(a.n, len(els)), replace=False)
els = els.iloc[np.sort(idx)].reset_index(drop=True)
print(f"{a.label}: {len(els):,} sampled elements", flush=True)

bws = {k: pyBigWig.open(v) for k, v in
       dict(r1p=a.rep1_plus, r1m=a.rep1_minus, r2p=a.rep2_plus, r2m=a.rep2_minus).items()}
sizes = bws["r1p"].chroms()

rows = []
for _, r in els.iterrows():
    c = r["chr"]
    if c not in sizes:
        continue
    ctr = int(r["start"]) + int(r["summit"])
    s, e = ctr - HW, ctr + HW
    if s < 0 or e > sizes[c]:
        continue
    v = {}
    ok = True
    for k, bw in bws.items():
        x = bw.values(c, s, e, numpy=True)
        if x is None or len(x) != 2 * HW:
            ok = False
            break
        v[k] = np.nan_to_num(x, nan=0.0)
    if ok:
        rows.append(v)
for bw in bws.values():
    bw.close()
print(f"extracted {len(rows):,} windows", flush=True)

r1p = np.stack([r["r1p"] for r in rows]); r1m = np.stack([r["r1m"] for r in rows])
r2p = np.stack([r["r2p"] for r in rows]); r2m = np.stack([r["r2m"] for r in rows])
del rows

tot = r1p.sum(1) + r1m.sum(1) + r2p.sum(1) + r2m.sum(1)
q = np.quantile(tot, [0.2, 0.4, 0.6, 0.8])
quint = np.digitize(tot, q)  # 0..4, 4 = top quintile


def shape_corr(x, y, b):
    """Mean within-element Pearson r between two binned profiles.

    Binning by reshape-and-sum; correlating along the position axis inside each element
    removes the element's total, so this is shape only. Elements with no variation in
    either replicate at this bin size carry no shape information and are excluded rather
    than counted as r = 0, which would confound "unreproducible" with "flat".
    """
    n, L = x.shape
    xb = x.reshape(n, L // b, b).sum(2)
    yb = y.reshape(n, L // b, b).sum(2)
    if xb.shape[1] < 3:
        return np.nan, 0
    xc = xb - xb.mean(1, keepdims=True)
    yc = yb - yb.mean(1, keepdims=True)
    num = (xc * yc).sum(1)
    den = np.sqrt((xc ** 2).sum(1) * (yc ** 2).sum(1))
    good = den > 0
    if good.sum() == 0:
        return np.nan, 0
    return float((num[good] / den[good]).mean()), int(good.sum())


def pooled_ceiling(r):
    """Ceiling on predicting the POOLED profile, from replicate-vs-replicate agreement.

    Raw replicate agreement understates what a model can reach. If each replicate is
    signal + independent noise, then r(rep1, rep2) = var(sig)/(var(sig)+var(noise)), while a
    smooth predictor is scored against the pooled track, whose noise variance is halved.
    The Spearman-Brown reliability of a 2-replicate mean is 2r/(1+r), and correlation with
    it goes as the square root -- the same correction 0.3 applies to the count ceiling. So a
    1 bp replicate agreement of 0.02 does NOT mean a model can only reach 0.02.
    """
    if not np.isfinite(r) or r <= 0:
        return np.nan
    return float(np.sqrt(2 * r / (1 + r)))


out = []
for b in BINS:
    for strat, mask in (("all", np.ones(len(tot), bool)), ("topq", quint == 4)):
        r_uns, n_uns = shape_corr((r1p + r1m)[mask], (r2p + r2m)[mask], b)
        r_pls, _ = shape_corr(r1p[mask], r2p[mask], b)
        r_min, _ = shape_corr(r1m[mask], r2m[mask], b)
        out.append({"bin_bp": b, "stratum": strat, "n_elements": n_uns,
                    "shape_r_unstranded": r_uns,
                    "shape_r_plus": r_pls, "shape_r_minus": r_min,
                    "ceiling_unstranded": pooled_ceiling(r_uns),
                    "ceiling_plus": pooled_ceiling(r_pls),
                    "ceiling_minus": pooled_ceiling(r_min),
                    "n_bins": (2 * HW) // b})
        print(f"  bin={b:>4} {strat:<5} n={n_uns:>6,} rep-vs-rep r={r_uns:.4f}  "
              f"-> pooled ceiling {pooled_ceiling(r_uns):.4f}  "
              f"(plus {pooled_ceiling(r_pls):.4f} minus {pooled_ceiling(r_min):.4f})",
              flush=True)

df = pd.DataFrame(out)
p1 = f"{P}/results/profile_ceiling_binsize_{a.label}.tsv"
df.round(4).to_csv(p1, sep="\t", index=False)
print("\nWrote", p1)

# Internal check against 0.3: the count ceiling on this same +/-1000 bp window.
c1 = np.log1p(r1p.sum(1) + r1m.sum(1))
c2 = np.log1p(r2p.sum(1) + r2m.sum(1))
for strat, mask in (("all", np.ones(len(tot), bool)), ("topq", quint == 4)):
    r = np.corrcoef(c1[mask], c2[mask])[0, 1]
    print(f"count ceiling ({strat}, 1000 bp window): r={r:.4f}  "
          f"Spearman-Brown corrected sqrt(2r/(1+r))={np.sqrt(2*r/(1+r)):.4f}")
