#!/usr/bin/env python3
"""Why doesn't better H3K27ac help ABC? Three mechanistic checks.

Observed H3K27ac beats ATAC alone by +0.062 AUPRC, and predicted H3K27ac recovers none of
that on precision-at-min-sensitivity. Either the prediction is bad specifically where the
benchmark looks, or it is accurate but loses the DISCRIMINATION that ABC's rank-based qnorm
needs, or the benchmark cannot resolve the difference. These are separable.

  A. ACCURACY WHERE IT MATTERS. Predicted vs observed H3K27ac per ABC region, genome-wide,
     then restricted to regions the CRISPR benchmark actually tests, then to regions of
     CRISPR-positive pairs. A top-quintile r of 0.69 over all elements says nothing about
     the few thousand elements the benchmark scores.

  B. DYNAMIC-RANGE COMPRESSION. ABC's activity depends on qnorm RANKS, so a regression model
     that shrinks extremes can be well correlated yet unable to order the top elements
     against each other. Compares spread within the tested regions, before and after qnorm.

  C. DISCORDANCE ON POSITIVES. Among the 466 regulated pairs, which does the observed arm
     score above its threshold while the predicted arm does not, and what do those elements
     look like? This localises the loss to specific pairs rather than a summary statistic.
"""
import gzip
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

D = "/oak/stanford/groups/engreitz/Users/sheth"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction/results"
OBS = f"{ABC}/2026_0721_h3k27ac_counting_comparison/K562_ATAC_H3K27ac_element/Neighborhoods/EnhancerList.txt"
PRED = f"{ABC}/2026_0903_predicted_activity/pred_k562_multimodal/Neighborhoods/EnhancerList.txt"
MERGED = (f"{D}/CRISPR_comparison_v3/CRISPR_comparison/workflow/results/"
          f"2026_0904_predicted_activity/expt_pred_merged_annot.txt.gz")
COLS = ["chr", "start", "end", "H3K27ac.RPM", "ATAC.RPM", "normalized_h3k27ac",
        "normalized_atac", "activity_base", "class"]


def load(path):
    """Index on `chr:start-end`, which is what the merged table's pred_elements holds.

    The EnhancerList `name` column is `class|chr:start-end` (e.g. `promoter|chr1:35831-36331`),
    so joining on it directly matches nothing and silently yields empty strata.
    """
    df = pd.read_csv(path, sep="\t", usecols=COLS)
    df.index = (df["chr"].astype(str) + ":" + df["start"].astype(str)
                + "-" + df["end"].astype(str))
    return df.drop(columns=["chr", "start", "end"])


obs = load(OBS)
prd = load(PRED)
common = obs.index.intersection(prd.index)
assert len(common) == len(obs) == len(prd), (len(common), len(obs), len(prd))
print(f"{len(common):,} ABC regions joined on name\n")

m = pd.read_csv(MERGED, sep="\t", compression="gzip",
                usecols=["pred_elements", "pred_uid", "pred_value", "Regulated",
                         "measuredGeneSymbol", "EffectSize"])
one = m[m["pred_uid"] == "K562_ATAC_only.ABC.Score"]
tested = set(one["pred_elements"].dropna())
positive = set(one.loc[one["Regulated"] == True, "pred_elements"].dropna())
print(f"{len(tested):,} regions carry a CRISPR-tested pair; "
      f"{len(positive):,} carry a regulated pair\n")

STRATA = [("all regions", common),
          ("CRISPR-tested regions", common.intersection(tested)),
          ("regions of regulated pairs", common.intersection(positive))]
for lab, idx in STRATA[1:]:
    assert len(idx) > 100, (
        f"{lab}: only {len(idx)} regions matched. The join key is wrong -- EnhancerList must "
        f"be indexed on chr:start-end to match the merged table's pred_elements.")

print("A. PREDICTED vs OBSERVED H3K27ac")
print(f"{'stratum':<30}{'n':>8}{'Spearman':>10}{'Pearson(log)':>14}")
for lab, idx in STRATA:
    o = obs.loc[idx, "H3K27ac.RPM"].to_numpy(float)
    p = prd.loc[idx, "H3K27ac.RPM"].to_numpy(float)
    keep = np.isfinite(o) & np.isfinite(p)
    if keep.sum() < 3:
        print(f"{lab:<30}{keep.sum():>8,}{'--':>10}{'--':>14}   (too few to correlate)")
        continue
    rs = spearmanr(p[keep], o[keep]).statistic
    rp = pearsonr(np.log1p(p[keep]), np.log1p(o[keep]))[0]
    print(f"{lab:<30}{keep.sum():>8,}{rs:>10.4f}{rp:>14.4f}")

print("\nB. DYNAMIC RANGE (activity_base, what ABC actually multiplies by contact)")
print(f"{'stratum':<30}{'source':<10}{'p50':>10}{'p90':>10}{'p99':>10}{'p99/p50':>9}{'sd(log)':>9}")
for lab, idx in STRATA:
    for src, df in (("observed", obs), ("predicted", prd)):
        v = df.loc[idx, "activity_base"].to_numpy(float)
        v = v[np.isfinite(v) & (v > 0)]
        if len(v) < 10:
            continue
        q50, q90, q99 = np.percentile(v, [50, 90, 99])
        print(f"{lab if src=='observed' else '':<30}{src:<10}{q50:>10.3f}{q90:>10.3f}"
              f"{q99:>10.3f}{q99/q50:>9.2f}{np.std(np.log1p(v)):>9.3f}")

print("\n   Within the tested regions, how well does each ordering separate the top decile?")
for lab, idx in STRATA[1:]:
    for src, df in (("observed", obs), ("predicted", prd)):
        v = df.loc[idx, "normalized_h3k27ac"].to_numpy(float)
        v = v[np.isfinite(v)]
        top = v >= np.percentile(v, 90)
        print(f"   {lab:<28}{src:<10}top-decile mean/median ratio "
              f"{v[top].mean()/max(np.median(v),1e-9):>7.2f}")

print("\nC. DISCORDANCE ON REGULATED PAIRS")
piv = m.pivot_table(index=["pred_elements", "measuredGeneSymbol"], columns="pred_uid",
                    values="pred_value", aggfunc="first")
reg = m[(m["pred_uid"] == "K562_ATAC_only.ABC.Score")].set_index(
    ["pred_elements", "measuredGeneSymbol"])["Regulated"]
piv = piv.join(reg.rename("Regulated"), how="inner")
OBS_C, PRD_C, FLOOR_C = ("K562_ATAC_H3K27ac_element.ABC.Score",
                         "pred_k562_multimodal.ABC.Score", "K562_ATAC_only.ABC.Score")
THR = {OBS_C: 0.011611, PRD_C: 0.011154, FLOOR_C: 0.008852}   # thresholdMinSens per arm
p = piv[piv["Regulated"] == True].dropna(subset=[OBS_C, PRD_C, FLOOR_C])
print(f"   {len(p):,} regulated pairs with a score from all three arms")
for a, b in ((OBS_C, PRD_C), (OBS_C, FLOOR_C), (PRD_C, FLOOR_C)):
    ha, hb = p[a] >= THR[a], p[b] >= THR[b]
    print(f"   {a.split('.')[0]:<28} caught {ha.sum():>4}   "
          f"{b.split('.')[0]:<24} caught {hb.sum():>4}   "
          f"only-first {int((ha & ~hb).sum()):>3}   only-second {int((~ha & hb).sum()):>3}")

lost = p[(p[OBS_C] >= THR[OBS_C]) & (p[PRD_C] < THR[PRD_C])]
print(f"\n   {len(lost)} pairs the observed arm catches and the predicted arm misses.")
if len(lost):
    e = lost.index.get_level_values(0)
    e = [x for x in e if x in obs.index]
    if e:
        print(f"   at those elements: observed H3K27ac.RPM median "
              f"{obs.loc[e,'H3K27ac.RPM'].median():.2f}, predicted "
              f"{prd.loc[e,'H3K27ac.RPM'].median():.2f}")
        print(f"   genome-wide medians: observed {obs['H3K27ac.RPM'].median():.2f}, "
              f"predicted {prd['H3K27ac.RPM'].median():.2f}")
