#!/usr/bin/env python3
"""How much does a predicted activity track know about observed p300 BEYOND ATAC?

WHY. rE2G gains +0.12 AUPRC from observed p300 and nothing from predicted p300 (F-029), and
the oracle decomposition puts most of ABC's loss in how the tested elements are ranked
(F-028). One explanation covers both: the predictions mostly re-encode the ATAC input, so
their ranking of elements carries little that ATAC does not already supply. This measures
that directly, as the partial rank correlation between predicted and observed p300 with
ATAC partialled out:

    rho(pred, obs | ATAC) = corr(resid(rank(pred) ~ rank(ATAC)), resid(rank(obs) ~ rank(ATAC)))

computed within ABC's promoter and non-promoter classes, since ABC ranks within class.

WHICH REGIONS. All 153,545 ABC candidate regions -- the headline, and the only version fit
for comparing models, because it never touches the benchmark. The CRISPR-tested regions are
reported separately as a DIAGNOSTIC; selecting models on them would leak the benchmark into
model choice.

Values are ABC's own per-region numbers from each arm's EnhancerList.txt (H3K27ac.RPM is the
activity slot, ATAC.RPM the accessibility), so this sees exactly what ABC ranked.

Usage: 4.35.partial_corr_given_atac.py [--out-tsv T]
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ABC = "/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction/results"
CR = ("/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v3/CRISPR_comparison/"
      "workflow/results/2026_0922_p300_multicell/expt_pred_merged_annot.txt.gz")
OBS = f"{ABC}/2026_0905_p300_activity/p300obs_k562"
ARMS = [
    ("observed H3K27ac (reference)", f"{ABC}/2026_0721_h3k27ac_counting_comparison/"
                                     "K562_ATAC_H3K27ac_element"),
    ("predicted p300, K562 in-cell", f"{ABC}/2026_0905_p300_activity/p300pred_k562_multimodal"),
    ("predicted p300, multi-cell", f"{ABC}/2026_0922_p300_multicell/p300pred_multicell_holdK562"),
    ("predicted p300, GM12878 target norm",
     f"{ABC}/2026_0922_p300_multicell/p300pred_gm12878_corrected_tgtnorm"),
    ("predicted p300, GM12878 own stats",
     f"{ABC}/2026_0919_p300_corrected/p300pred_gm12878_corrected"),
]
KEY = ["chr", "start", "end"]

ap = argparse.ArgumentParser()
ap.add_argument("--out-tsv", default=None)
a = ap.parse_args()


def load(arm_dir):
    return pd.read_csv(f"{arm_dir}/Neighborhoods/EnhancerList.txt", sep="\t",
                       usecols=KEY + ["class", "H3K27ac.RPM", "ATAC.RPM"])


def partial_rho(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    Z = np.column_stack([np.ones_like(rz), rz])
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


obs = load(OBS).rename(columns={"H3K27ac.RPM": "obs", "ATAC.RPM": "atac"})

# tested regions: candidate regions overlapping any CRISPR-tested element (diagnostic only)
cr = pd.read_csv(CR, sep="\t", usecols=["chrom", "chromStart", "chromEnd"]).drop_duplicates()
obs["tested"] = False
for c, g in obs.groupby("chr"):
    e = cr[cr.chrom == c]
    if e.empty:
        continue
    o = np.argsort(g["start"].values)
    s, en, idx = g["start"].values[o], g["end"].values[o], g.index.values[o]
    for es, ee in zip(e.chromStart.values, e.chromEnd.values):
        lo, hi = np.searchsorted(en, es, side="right"), np.searchsorted(s, ee, side="left")
        if hi > lo:
            obs.loc[idx[lo:hi], "tested"] = True

rows = []
for name, d in ARMS:
    p = load(d).rename(columns={"H3K27ac.RPM": "pred", "ATAC.RPM": "atac_arm"})
    m = obs.merge(p, on=KEY, suffixes=("", "_p"), validate="one_to_one")
    assert len(m) == len(obs), f"{name}: region sets differ ({len(m)} vs {len(obs)})"
    assert np.allclose(m["atac"], m["atac_arm"]), f"{name}: ATAC input differs from obs arm"
    m["prom"] = m["class"].isin(["tss", "promoter"])
    for subset, sel in (("all candidates", np.ones(len(m), bool)),
                        ("tested (diagnostic)", m["tested"].values)):
        for cls, cm in (("nonpromoter", ~m["prom"].values), ("promoter", m["prom"].values)):
            k = sel & cm
            x, y, z = m.loc[k, "pred"].values, m.loc[k, "obs"].values, m.loc[k, "atac"].values
            rows.append(dict(arm=name, regions=subset, cls=cls, n=int(k.sum()),
                             rho_pred_obs=spearmanr(x, y)[0],
                             rho_atac_obs=spearmanr(z, y)[0],
                             rho_pred_atac=spearmanr(x, z)[0],
                             partial_rho_given_atac=partial_rho(x, y, z)))

r = pd.DataFrame(rows)
pd.set_option("display.width", 200)
for subset in ("all candidates", "tested (diagnostic)"):
    print(f"\n=== {subset} ===")
    print(r[r.regions == subset].drop(columns="regions").to_string(
        index=False, float_format=lambda v: f"{v:.3f}"))
if a.out_tsv:
    r.to_csv(a.out_tsv, sep="\t", index=False)
    print("\nwrote", a.out_tsv)
