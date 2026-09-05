#!/usr/bin/env python3
"""Characterise BOTH tails of the H3K27ac prediction error, not just the false negatives.

The ABC diagnostic found the model ranks CRISPR-functional elements too low while catching
the same number of positives as the floor, so its deficit is in suppressing negatives. Two
error directions produce that, and they have opposite biology:

  UNDER-prediction  strong observed H3K27ac called near-background. The 12 threshold-level
                    false negatives sit here (observed RPM 22.8, predicted 0.7).
  OVER-prediction   weak observed H3K27ac called high. Accessible-but-unacetylated elements
                    -- CTCF sites, insulators, structural loops -- would land here, because
                    a model leaning on accessibility has no way to know they are not
                    acetylated. This is the direction that would spoil negative suppression.

Widened from 12 elements to both full tails so there is power to compare them, and reported
against a matched-accessibility background: elements differ in ATAC, and without matching any
contrast would just rediscover that.

Outputs per stratum: n, observed and predicted H3K27ac, ATAC, GC content, CpG ratio, element
class, and the ATAC/H3K27ac ratio that defines the accessible-but-inactive phenotype.
"""
import argparse
import os
import numpy as np
import pandas as pd
import pyfaidx

_ap = argparse.ArgumentParser()
_ap.add_argument("--emit-beds", default=None,
                 help="Directory to write one sorted BED per stratum, for overlap testing "
                      "against ENCODE peak calls (see 4.11).")
_args = _ap.parse_args()

D = "/oak/stanford/groups/engreitz/Users/sheth"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction/results"
OBS = f"{ABC}/2026_0721_h3k27ac_counting_comparison/K562_ATAC_H3K27ac_element/Neighborhoods/EnhancerList.txt"
PRED = f"{ABC}/2026_0903_predicted_activity/pred_k562_multimodal/Neighborhoods/EnhancerList.txt"
GEN = f"{D}/hg38_resources/hg38.fa"
OUT = f"{D}/EP300_BPNet/2026_0824_H3K27ac_model/results/prediction_error_strata.tsv"
COLS = ["chr", "start", "end", "H3K27ac.RPM", "ATAC.RPM", "class"]


def load(path, tag):
    df = pd.read_csv(path, sep="\t", usecols=COLS)
    df.index = df["chr"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
    return df.rename(columns={"H3K27ac.RPM": f"k27_{tag}"})


obs, prd = load(OBS, "obs"), load(PRED, "pred")
d = obs.join(prd[["k27_pred"]], how="inner")
assert len(d) == len(obs) == len(prd), (len(d), len(obs), len(prd))
print(f"{len(d):,} regions")

# Work in log space; the residual is what "the model got this wrong" means here.
d["lo"] = np.log1p(d["k27_obs"])
d["lp"] = np.log1p(d["k27_pred"])
d["latac"] = np.log1p(d["ATAC.RPM"])

# Residualise the error on accessibility, so a stratum is not just "more open".
A = np.vstack([np.ones(len(d)), d["latac"].to_numpy()]).T
for src in ("lo", "lp"):
    beta, *_ = np.linalg.lstsq(A, d[src].to_numpy(), rcond=None)
    d[src + "_r"] = d[src].to_numpy() - A @ beta
d["err"] = d["lp_r"] - d["lo_r"]          # >0 over-predicted, <0 under-predicted

q = d["err"].quantile([0.01, 0.05, 0.95, 0.99])
STRATA = [
    ("under-predicted 1%", d["err"] <= q[0.01]),
    ("under-predicted 5%", d["err"] <= q[0.05]),
    ("typical (middle 50%)", (d["err"] > d["err"].quantile(0.25)) & (d["err"] < d["err"].quantile(0.75))),
    ("over-predicted 5%", d["err"] >= q[0.95]),
    ("over-predicted 1%", d["err"] >= q[0.99]),
]

gen = pyfaidx.Fasta(GEN)


def seqstats(sub, n_max=4000):
    """GC and CpG observed/expected over the element, on a capped random sample."""
    idx = sub.index if len(sub) <= n_max else sub.sample(n_max, random_state=0).index
    gc, cpg = [], []
    for k in idx:
        c, se = k.split(":")
        s, e2 = se.split("-")
        try:
            seq = str(gen[c][int(s):int(e2)]).upper()
        except (KeyError, ValueError):
            continue
        if len(seq) < 50:
            continue
        g, cc = seq.count("G"), seq.count("C")
        gc.append((g + cc) / len(seq))
        exp = (g * cc) / len(seq)
        cpg.append(seq.count("CG") / exp if exp > 0 else np.nan)
    return (np.nanmean(gc) if gc else np.nan), (np.nanmean(cpg) if cpg else np.nan)


rows = []
print(f"\n{'stratum':<24}{'n':>7}{'k27_obs':>9}{'k27_pred':>10}{'ATAC':>8}"
      f"{'ATAC/k27':>10}{'GC':>7}{'CpG o/e':>9}{'%promoter':>10}")
for lab, mask in STRATA:
    sub = d[mask]
    gc, cpg = seqstats(sub)
    atac_k27 = (sub["ATAC.RPM"] / sub["k27_obs"].replace(0, np.nan)).median()
    r = dict(stratum=lab, n=len(sub),
             k27_obs=sub["k27_obs"].median(), k27_pred=sub["k27_pred"].median(),
             atac=sub["ATAC.RPM"].median(), atac_over_k27=atac_k27,
             gc=gc, cpg_oe=cpg,
             pct_promoter=100.0 * (sub["class"] == "promoter").mean(),
             pct_intergenic=100.0 * (sub["class"] == "intergenic").mean())
    rows.append(r)
    print(f"{lab:<24}{r['n']:>7,}{r['k27_obs']:>9.2f}{r['k27_pred']:>10.2f}"
          f"{r['atac']:>8.2f}{r['atac_over_k27']:>10.2f}{r['gc']:>7.3f}"
          f"{r['cpg_oe']:>9.3f}{r['pct_promoter']:>10.1f}")

if _args.emit_beds:
    os.makedirs(_args.emit_beds, exist_ok=True)
    for lab, mask in STRATA:
        sub = d[mask]
        bed = pd.DataFrame({"chr": sub["chr"], "start": sub["start"], "end": sub["end"]})
        bed = bed.sort_values(["chr", "start"])
        safe = lab.replace(" ", "_").replace("%", "pct").replace("(", "").replace(")", "")
        bed.to_csv(f"{_args.emit_beds}/stratum_{safe}.bed", sep="\t",
                   header=False, index=False)

pd.DataFrame(rows).round(4).to_csv(OUT, sep="\t", index=False)
print(f"\nwrote {OUT}")
print("\nATAC/k27 is the accessible-but-unacetylated signature: high in a stratum means those")
print("elements are open without being acetylated, which is what a CTCF-site explanation")
print("predicts for the over-predicted tail.")
