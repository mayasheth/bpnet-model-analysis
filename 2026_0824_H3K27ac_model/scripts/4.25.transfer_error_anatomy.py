#!/usr/bin/env python3
"""What a transferred ATAC-input model gets wrong on the CRISPR benchmark, element by element.

THE COMPARISON IS STRUCTURALLY CLEAN, WHICH IS WHY IT IS WORTH DOING HERE. The in-cell and
transferred arms are the same ABC run with the same ATAC tagAligns, the same candidate regions,
the same contact model and the same qnorm. The ONLY input that differs is the predicted p300
bigwig in the activity slot. So every difference in ABC.Score traces to the predicted track, and
no counting-path or region-set confound can enter.

WHAT IT REPORTS, in the order the questions arise:
  1. Rank anatomy. AUPRC is a ranking statistic, so the loss has to appear as positives falling
     or negatives rising. Reported as percentile shifts, not as a single threshold.
  2. Where in the accessibility range the damage sits. ATAC is identical between the arms, so
     stratifying by ATAC decile asks whether transfer fails uniformly or in a specific stratum.
  3. Dynamic range on the regulated-pair elements. F-004 named compressed spread as the
     mechanism behind the H3K27ac negative result; this checks whether transfer reintroduces it.
  4. Counts accuracy against OBSERVED p300, restricted to the CRISPR elements. The project's
     headline Pearson is genome-wide over held-out folds, which is not the population the
     benchmark is decided on. This computes the same quantity where it matters.

Usage:
  4.25.transfer_error_anatomy.py --merged <expt_pred_merged_annot.txt.gz> \\
      --incell-arm p300pred_k562_multimodal.ABC.Score \\
      --transfer-arm p300pred_gm12878_multimodal.ABC.Score \\
      --floor-arm K562_ATAC_only.ABC.Score \\
      --putative incell=<dir>/EnhancerPredictionsAllPutative.tsv.gz \\
      --putative transfer=<...> --putative observed=<...>
"""
import argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--merged", required=True)
ap.add_argument("--incell-arm", required=True)
ap.add_argument("--transfer-arm", required=True)
ap.add_argument("--floor-arm", required=True)
ap.add_argument("--putative", action="append", default=[],
                metavar="LABEL=PATH", help="repeatable; labels incell/transfer/observed")
ap.add_argument("--top-frac", type=float, default=0.10,
                help="the 'called' fraction used for the demotion/promotion counts")
ap.add_argument("--out-prefix", default=None)
a = ap.parse_args()

PUT = dict(kv.split("=", 1) for kv in a.putative)

# ---------------------------------------------------------------- labels and per-arm scores
cols = ["chrom", "chromStart", "chromEnd", "name", "measuredGeneSymbol", "Regulated",
        "EffectSize", "startTSS", "pred_elements", "pred_uid", "pred_value"]
d = pd.read_csv(a.merged, sep="\t", low_memory=False, usecols=cols)
d["pair"] = d["name"].astype(str) + "|" + d["measuredGeneSymbol"].astype(str)

arms = [a.incell_arm, a.transfer_arm, a.floor_arm]
wide = d[d["pred_uid"].isin(arms)].pivot_table(index="pair", columns="pred_uid",
                                               values="pred_value", aggfunc="max")
missing = [x for x in arms if x not in wide.columns]
if missing:
    raise SystemExit(f"arm(s) absent from {a.merged}: {missing}")

meta = (d.drop_duplicates("pair")
         .set_index("pair")[["chrom", "chromStart", "chromEnd", "measuredGeneSymbol",
                             "Regulated", "EffectSize", "startTSS", "pred_elements"]])
t = wide.join(meta, how="inner").dropna(subset=arms)
t["y"] = t["Regulated"].astype(str).str.upper().eq("TRUE").astype(int)
t["dist"] = (t["chromStart"] + t["chromEnd"]) / 2 - t["startTSS"]
t["absdist"] = t["dist"].abs()

IN, TX, FL = a.incell_arm, a.transfer_arm, a.floor_arm
print(f"{len(t):,} element-gene pairs, {int(t['y'].sum()):,} regulated "
      f"({100 * t['y'].mean():.2f}%)\n")


def auprc(y, s):
    o = np.argsort(-np.asarray(s, float), kind="mergesort")
    y = np.asarray(y)[o]
    tp = np.cumsum(y)
    prec = tp / np.arange(1, len(y) + 1)
    return float((prec * y).sum() / max(y.sum(), 1))


print("=" * 78)
print("1. RANK ANATOMY: the loss has to be positives falling or negatives rising")
print("=" * 78)
for nm, c in (("in-cell", IN), ("transferred", TX), ("ATAC-only floor", FL)):
    print(f"  AUPRC {nm:<18} {auprc(t['y'], t[c]):.4f}")

for c in arms:
    t[f"pct_{c}"] = t[c].rank(pct=True)

pos, neg = t[t["y"] == 1], t[t["y"] == 0]
print(f"\n  mean percentile of the {len(pos):,} POSITIVES")
for nm, c in (("in-cell", IN), ("transferred", TX), ("floor", FL)):
    print(f"    {nm:<14} {pos[f'pct_{c}'].mean():.4f}")
print(f"  mean percentile of the {len(neg):,} NEGATIVES")
for nm, c in (("in-cell", IN), ("transferred", TX), ("floor", FL)):
    print(f"    {nm:<14} {neg[f'pct_{c}'].mean():.4f}")

k = int(round(a.top_frac * len(t)))
thr_in = t[IN].nlargest(k).iloc[-1]
thr_tx = t[TX].nlargest(k).iloc[-1]
called_in, called_tx = t[IN] >= thr_in, t[TX] >= thr_tx
print(f"\n  at the top {100 * a.top_frac:.0f}% of pairs ({k:,} called by each arm)")
print(f"    positives called by in-cell only  {int(((t['y'] == 1) & called_in & ~called_tx).sum()):>5}")
print(f"    positives called by transfer only {int(((t['y'] == 1) & ~called_in & called_tx).sum()):>5}")
print(f"    positives called by both          {int(((t['y'] == 1) & called_in & called_tx).sum()):>5}")
print(f"    negatives called by in-cell only  {int(((t['y'] == 0) & called_in & ~called_tx).sum()):>5}")
print(f"    negatives called by transfer only {int(((t['y'] == 0) & ~called_in & called_tx).sum()):>5}")

t["dpct"] = t[f"pct_{TX}"] - t[f"pct_{IN}"]
print(f"\n  percentile shift on transfer (transferred - in-cell)")
print(f"    positives  mean {pos.index.map(t['dpct']).to_series().mean():+.4f}   "
      f"median {np.median(t.loc[pos.index, 'dpct']):+.4f}")
print(f"    negatives  mean {t.loc[neg.index, 'dpct'].mean():+.4f}   "
      f"median {np.median(t.loc[neg.index, 'dpct']):+.4f}")
print("  A negative shift on positives means transfer demotes real enhancers. A positive")
print("  shift on negatives means it promotes non-functional ones. Both cost AUPRC.")

# ------------------------------------------------------------------- activity components
if PUT:
    print("\n" + "=" * 78)
    print("2. TRACING IT TO THE PREDICTED TRACK (ATAC is identical between the arms)")
    print("=" * 78)
    t["ekey"] = t["pred_elements"].astype(str) + "|" + t["measuredGeneSymbol"].astype(str)
    want = set(t["ekey"])
    keep = ["chr", "start", "end", "TargetGene", "normalized_atac_enh",
            "normalized_h3k27ac_enh", "activity_base", "distance", "ABC.Score"]
    comp = {}
    for label, path in PUT.items():
        got = []
        for ch in pd.read_csv(path, sep="\t", usecols=keep, chunksize=500_000,
                              low_memory=False):
            ch["ekey"] = (ch["chr"].astype(str) + ":" + ch["start"].astype(str) + "-"
                          + ch["end"].astype(str) + "|" + ch["TargetGene"].astype(str))
            got.append(ch[ch["ekey"].isin(want)])
        c = pd.concat(got).drop_duplicates("ekey").set_index("ekey")
        comp[label] = c
        print(f"  {label:<10} matched {len(c):,} of {len(want):,} element-gene keys")

    e = t.set_index("ekey")
    for label, c in comp.items():
        e[f"act_{label}"] = c["normalized_h3k27ac_enh"]
        e[f"base_{label}"] = c["activity_base"]
    e["atac"] = comp[list(comp)[0]]["normalized_atac_enh"]
    e = e.dropna(subset=[f"act_{l}" for l in comp])
    print(f"  {len(e):,} pairs with components from every arm\n")

    ep, en = e[e["y"] == 1], e[e["y"] == 0]
    if {"incell", "transfer"} <= set(comp):
        print("  predicted p300 in the activity slot, median over pairs")
        print(f"    {'':<12}{'positives':>12}{'negatives':>12}{'pos/neg ratio':>15}")
        for label in ("incell", "transfer", "observed"):
            if label not in comp:
                continue
            mp, mn = ep[f"act_{label}"].median(), en[f"act_{label}"].median()
            print(f"    {label:<12}{mp:>12.4f}{mn:>12.4f}{mp / max(mn, 1e-9):>15.3f}")
        print("  The pos/neg ratio is the separation the activity term provides. A transferred")
        print("  arm can lose AUPRC by shrinking it from either side.")

        print(f"\n  dynamic range of the predicted track on these elements (p99/p50)")
        for label in ("incell", "transfer", "observed"):
            if label not in comp:
                continue
            v = e[f"act_{label}"]
            print(f"    {label:<12}{np.percentile(v, 99) / max(np.percentile(v, 50), 1e-9):>8.2f}x"
                  f"   on regulated-pair elements only "
                  f"{np.percentile(ep[f'act_{label}'], 99) / max(np.percentile(ep[f'act_{label}'], 50), 1e-9):>6.2f}x")

        print("\n  ATAC decile (identical in both arms): where does transfer lose rank?")
        # Built with an explicit loop rather than groupby.apply: this pandas is old enough
        # that apply's signature differs, and the loop is clearer about the per-decile rows.
        e["atac_dec"] = pd.qcut(e["atac"].rank(method="first"), 10, labels=False) + 1
        rows = []
        for dec, x in e.groupby("atac_dec"):
            xp, xn = x[x["y"] == 1], x[x["y"] == 0]
            rows.append({
                "atac_dec": int(dec), "n": len(x), "n_pos": int(x["y"].sum()),
                "d_pct_pos": xp["dpct"].mean() if len(xp) else np.nan,
                "d_pct_neg": xn["dpct"].mean(),
                "med_act_incell": x["act_incell"].median(),
                "med_act_transfer": x["act_transfer"].median(),
                "med_act_obs": x["act_observed"].median() if "act_observed" in x else np.nan,
                "tx_over_in": x["act_transfer"].median() / max(x["act_incell"].median(), 1e-9),
                "tx_over_obs": (x["act_transfer"].median() / max(x["act_observed"].median(), 1e-9)
                                if "act_observed" in x else np.nan),
                "in_over_obs": (x["act_incell"].median() / max(x["act_observed"].median(), 1e-9)
                                if "act_observed" in x else np.nan),
            })
        print(pd.DataFrame(rows).to_string(index=False,
                                           float_format=lambda v: f"{v:9.4f}"))

    if "observed" in comp:
        print("\n" + "=" * 78)
        print("3. COUNTS ACCURACY ON THE CRISPR ELEMENTS, WHICH IS NOT THE GENOME-WIDE NUMBER")
        print("=" * 78)
        obs = e["act_observed"]
        print(f"  {'':<12}{'Pearson':>10}{'Spearman':>10}{'Pearson(log)':>14}"
              f"{'  on regulated-pair elements only':>34}")
        for label in ("incell", "transfer"):
            v = e[f"act_{label}"]
            pr = np.corrcoef(obs, v)[0, 1]
            # Rank then Pearson. pandas' spearman imports scipy, which is broken against
            # the numpy installed here (numpy.random.bit_generator missing).
            sr = np.corrcoef(pd.Series(obs).rank(), pd.Series(v).rank())[0, 1]
            lr = np.corrcoef(np.log1p(obs), np.log1p(v))[0, 1]
            po = ep["act_observed"]
            pv = ep[f"act_{label}"]
            prp = np.corrcoef(po, pv)[0, 1]
            print(f"  {label:<12}{pr:>10.4f}{sr:>10.4f}{lr:>14.4f}"
                  f"{prp:>34.4f}")
        print("  If transfer's Pearson here is far below its genome-wide value, the benchmark")
        print("  is being decided on exactly the elements the model transfers worst on.")

if a.out_prefix:
    t.to_csv(f"{a.out_prefix}_pairs.tsv", sep="\t")
    print(f"\nwrote {a.out_prefix}_pairs.tsv")
