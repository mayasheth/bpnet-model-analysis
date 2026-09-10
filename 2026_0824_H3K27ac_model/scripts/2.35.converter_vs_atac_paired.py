#!/usr/bin/env python
"""The converter's headline number: how much more DNase-like is the painted track than ATAC?

Everything else measured so far answers a neighbouring question. 2.31 gives fraction of the
DNase inter-replicate ceiling, which says how close the converter gets to DNase. 2.15 gives
the sequence margin over an ATAC-only model, which says whether sequence is contributing.
Neither is the deployment decision. The downstream model currently reads the observed ATAC
track; swapping in a converted track is only worth doing if the converted track is more
DNase-like than the one being replaced. That difference, paired within fold, is what this
computes.

PAIRED WITHIN FOLD, because fold sd on this project runs 0.041-0.046, several times most of
the effects being measured, and an unpaired comparison of two 5-point clouds resolves
nothing. Same reasoning as 4.17 and 2.20.

INPUTS. `converter_profile_<label>_<cell>.tsv` from 2.31 and
`atac_vs_dnase_profile_baseline.tsv` from 0.35. Both compute shape with the identical
`shape_corr` on the identical held-out windows per fold, which is the only reason they can
be differenced at all.

Usage: 2.35.converter_vs_atac_paired.py --label converter5f [--label ataconly5f]
"""
import argparse
import numpy as np
import pandas as pd
from scipy.stats import ttest_rel, t as tdist

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
TC = tdist.ppf(0.975, df=4)

ap = argparse.ArgumentParser()
ap.add_argument("--label", action="append", required=True)
ap.add_argument("--cells", nargs="+", default=["k562", "gm12878"])
ap.add_argument("--bins", type=int, nargs="+", default=None)
a = ap.parse_args()

base = pd.read_csv(f"{P}/results/atac_vs_dnase_profile_baseline.tsv", sep="\t")
if "fold" not in base.columns:
    raise SystemExit("baseline table has no fold column; rerun 0.35 (per-fold version)")

rows = []
for label in a.label:
    for cell in a.cells:
        f = f"{P}/results/converter_profile_{label}_{cell}.tsv"
        try:
            m = pd.read_csv(f, sep="\t")
        except FileNotFoundError:
            print(f"  (skipping {label} on {cell}: {f} not found)")
            continue
        b = base[base["cell"] == cell]
        j = m.merge(b[["fold", "bin_bp", "stratum", "shape_r_observed_atac_vs_dnase"]],
                    on=["fold", "bin_bp", "stratum"], how="inner")
        if len(j) != len(m):
            raise SystemExit(f"{label}/{cell}: {len(m)} model rows matched {len(j)} baseline "
                             f"rows; the two tables do not cover the same folds and bins")
        for (b_, st), g in j.groupby(["bin_bp", "stratum"]):
            g = g.sort_values("fold")
            d = (g["shape_r_model"] - g["shape_r_observed_atac_vs_dnase"]).to_numpy(float)
            k = len(d)
            mu = d.mean()
            half = TC * d.std(ddof=1) / np.sqrt(k) if k > 1 else np.nan
            pv = (ttest_rel(g["shape_r_model"], g["shape_r_observed_atac_vs_dnase"]).pvalue
                  if k > 1 else np.nan)
            rows.append({"label": label, "cell": cell, "bin_bp": b_, "stratum": st,
                         "n_folds": k,
                         "atac": g["shape_r_observed_atac_vs_dnase"].mean(),
                         "model": g["shape_r_model"].mean(),
                         "ceiling": g["ceiling"].mean(),
                         "delta": mu, "ci_lo": mu - half, "ci_hi": mu + half, "p": pv})

df = pd.DataFrame(rows)
out = f"{P}/results/converter_vs_atac_paired.tsv"
df.round(4).to_csv(out, sep="\t", index=False)
print("Wrote", out)

show = df if a.bins is None else df[df["bin_bp"].isin(a.bins)]
print("\nShape correlation with observed DNase. delta = converter - observed ATAC, paired "
      "within fold.")
print(f"\n{'label':>13}{'cell':>9}{'bin':>5}{'strat':>6}{'ATAC':>8}{'model':>8}"
      f"{'ceil':>8}{'delta':>9}{'95% CI':>20}{'p':>9}")
for _, r in show.sort_values(["label", "cell", "stratum", "bin_bp"]).iterrows():
    ci = (f"[{r['ci_lo']:+.3f}, {r['ci_hi']:+.3f}]" if np.isfinite(r["ci_lo"]) else "-")
    pv = f"{r['p']:.4f}" if np.isfinite(r["p"]) else "-"
    print(f"{r['label']:>13}{r['cell']:>9}{int(r['bin_bp']):>5}{r['stratum']:>6}"
          f"{r['atac']:>8.3f}{r['model']:>8.3f}{r['ceiling']:>8.3f}{r['delta']:>+9.3f}"
          f"{ci:>20}{pv:>9}")
print("\nThe 1 bp top-quintile row is the deployment-relevant one: it is the resolution the")
print("downstream accessibility branch reads, on the elements that carry signal.")
