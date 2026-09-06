#!/usr/bin/env python3
"""Is CTCF actually the cause of the over-prediction, or a passenger?

4.13 established that the over-predicted stratum has 2.6x the CTCF RPKM of a typical
element, with IgG BELOW background so it is not an antibody-accessibility artifact. That is
a STRATUM-LEVEL MEAN. It does not establish that the elements being over-predicted are the
CTCF-bound ones: a 2.6x mean can come from a minority, and 4.11's peak overlap says only
36.9% of the over-predicted 5% carries a CTCF peak at all, so ~63% of the stratum has no
CTCF to blame.

Four element-level tests, each of which the mean enrichment cannot answer:

  1. DISTRIBUTION, not mean. Quartiles of CTCF signal within each stratum, plus the fraction
     of each stratum above a background-defined CTCF threshold. If the enrichment is carried
     by a minority, the mean is a misleading summary of the stratum.
  2. DOES ERROR TRACK CTCF AT MATCHED ACCESSIBILITY. CTCF sites are accessible and the model
     leans on accessibility, so any raw error-vs-CTCF association is confounded. Correlate
     within ATAC deciles.
  3. IS CTCF DOING ANYTHING BEYOND GC / CpG / PROMOTER STATUS. CTCF sites are GC- and
     CpG-rich and the over-predicted tail is promoter-enriched, so these are collinear.
     Partial correlation of error with CTCF given ATAC, GC, CpG and promoter class.
  4. THE DECISIVE ONE -- ATTRIBUTABLE FRACTION. Remove the CTCF-high elements from the
     over-predicted tail and ask how much over-prediction is left. If the remainder is still
     grossly over-predicted, CTCF explains a minority of the phenotype and the CTCF story is
     a passenger, however real the enrichment.

Reports IgG through every test as the background control, since a CTCF column that behaves
like IgG is reporting antibody accessibility rather than binding.
"""
import os
import sys
import numpy as np
import pandas as pd
import pyfaidx

D = "/oak/stanford/groups/engreitz/Users/sheth"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction/results"
OBS = f"{ABC}/2026_0721_h3k27ac_counting_comparison/K562_ATAC_H3K27ac_element/Neighborhoods/EnhancerList.txt"
PRED = f"{ABC}/2026_0903_predicted_activity/pred_k562_multimodal/Neighborhoods/EnhancerList.txt"
SEQARM = f"{ABC}/2026_0903_predicted_activity/pred_k562_sequence/Neighborhoods/EnhancerList.txt"
GEN = f"{D}/hg38_resources/hg38.fa"
COUNTS = sys.argv[1]          # bedtools multicov output
MARKS = sys.argv[2].split(",")  # column names for the multicov count columns, in order
TOTALS = {k: float(v) for k, v in
          (p.split("=") for p in sys.argv[3].split(","))}  # mark=mapped_reads


def load(path, tag):
    df = pd.read_csv(path, sep="\t", usecols=["chr", "start", "end", "H3K27ac.RPM",
                                              "ATAC.RPM", "class"])
    df.index = df["chr"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
    return df.rename(columns={"H3K27ac.RPM": tag})


d = load(OBS, "obs")
d = d.join(load(PRED, "pred")[["pred"]], how="inner")
d = d.join(load(SEQARM, "seq")[["seq"]], how="inner")

# --- per-element mark signal, as RPKM so widths and libraries are comparable -------------
cnt = pd.read_csv(COUNTS, sep="\t", header=None)
ncol = cnt.shape[1]
assert ncol == 3 + len(MARKS), f"multicov has {ncol} cols, expected {3 + len(MARKS)}"
cnt.columns = ["chr", "start", "end"] + MARKS
cnt.index = cnt["chr"] + ":" + cnt["start"].astype(str) + "-" + cnt["end"].astype(str)
kb = (cnt["end"] - cnt["start"]) / 1000.0
for m in MARKS:
    cnt[m + "_rpkm"] = cnt[m] / (TOTALS[m] / 1e6) / kb
d = d.join(cnt[[m + "_rpkm" for m in MARKS]], how="inner")
print(f"{len(d):,} regions with predictions and per-element mark counts\n")

# --- GC and CpG per element ---------------------------------------------------------------
gen = pyfaidx.Fasta(GEN)
gc, cpg = np.full(len(d), np.nan), np.full(len(d), np.nan)
for i, (c, s, e) in enumerate(zip(d["chr"], d["start"], d["end"])):
    try:
        seq = str(gen[c][int(s):int(e)]).upper()
    except (KeyError, ValueError):
        continue
    if len(seq) < 50:
        continue
    g, cc = seq.count("G"), seq.count("C")
    gc[i] = (g + cc) / len(seq)
    cpg[i] = (seq.count("CG") * len(seq) / (g * cc)) if g and cc else np.nan
d["gc"], d["cpg_oe"] = gc, cpg
d["promoter"] = (d["class"] == "promoter").astype(float)

# --- prediction error, defined exactly as in 4.10 ----------------------------------------
for c in ["obs", "pred", "seq"]:
    d["l_" + c] = np.log1p(d[c])
d["latac"] = np.log1p(d["ATAC.RPM"])
A = np.vstack([np.ones(len(d)), d["latac"].to_numpy()]).T
for c in ["obs", "pred", "seq"]:
    beta, *_ = np.linalg.lstsq(A, d["l_" + c].to_numpy(), rcond=None)
    d["r_" + c] = d["l_" + c].to_numpy() - A @ beta
d["err"] = d["r_pred"] - d["r_obs"]

q = d["err"].quantile([0.01, 0.05, 0.25, 0.75, 0.95, 0.99])
STRATA = [("under-predicted 1%", d["err"] <= q[0.01]),
          ("under-predicted 5%", d["err"] <= q[0.05]),
          ("typical (middle 50%)", (d["err"] > q[0.25]) & (d["err"] < q[0.75])),
          ("over-predicted 5%", d["err"] >= q[0.95]),
          ("over-predicted 1%", d["err"] >= q[0.99])]

RP = {m: m + "_rpkm" for m in MARKS}
OUT_EL = f"{D}/EP300_BPNet/2026_0824_H3K27ac_model/results/error_strata_ctcf_elements.tsv"
CT, IG = RP["CTCF"], RP["IgG"]

# =========================================================================================
# Per-element table for the figure, so the plot and these tests share one computation.
_lab = pd.Series("other", index=d.index)
for _l, _m in STRATA:
    if "middle" in _l or "1%" in _l:
        _lab[_m] = _l
_out = d[["chr", "start", "end", "obs", "pred", "seq", "ATAC.RPM", "gc", "cpg_oe",
          "promoter", "err", RP["CTCF"], RP["IgG"]]].copy()
_out["stratum"] = _lab
_out["over5"] = STRATA[3][1].astype(int)
_out.round(5).to_csv(OUT_EL, sep="\t", index=False)
print(f"wrote {OUT_EL}\n")

print("TEST 1  distribution of CTCF signal within each stratum, not just the mean")
print("        A mean enrichment carried by a minority is a misleading stratum summary.")
typ = d[STRATA[2][1]]
# background threshold: the 95th percentile of IgG-scaled CTCF in the typical stratum
thr = np.nanpercentile(typ[CT], 90)
print(f"        CTCF-high threshold = 90th pct of the typical stratum = {thr:.2f} RPKM\n")
print(f"{'stratum':<24}{'n':>7}{'CTCF p25':>10}{'p50':>8}{'p75':>8}{'p90':>8}"
      f"{'mean':>8}{'%>thr':>8}{'IgG p50':>9}")
for lab, mask in STRATA:
    s = d[mask]
    print(f"{lab:<24}{len(s):>7,}"
          f"{s[CT].quantile(.25):>10.2f}{s[CT].median():>8.2f}{s[CT].quantile(.75):>8.2f}"
          f"{s[CT].quantile(.90):>8.2f}{s[CT].mean():>8.2f}"
          f"{100.0 * (s[CT] > thr).mean():>8.1f}{s[IG].median():>9.2f}")

# =========================================================================================
print("\nTEST 2  does the error track CTCF at MATCHED accessibility")
print("        Raw association is confounded: CTCF sites are open and the model reads ATAC.")
d["atac_dec"] = pd.qcut(d["latac"], 10, labels=False, duplicates="drop")
rows = []
for dec, s in d.groupby("atac_dec"):
    if len(s) < 200:
        continue
    rows.append((int(dec), len(s),
                 s["err"].corr(np.log1p(s[CT]), method="spearman"),
                 s["err"].corr(np.log1p(s[IG]), method="spearman"),
                 s["err"].corr(s["gc"], method="spearman")))
print(f"{'ATAC decile':<13}{'n':>8}{'rho(err,CTCF)':>15}{'rho(err,IgG)':>14}{'rho(err,GC)':>13}")
for dec, n, rc, ri, rg in rows:
    print(f"{dec:<13}{n:>8,}{rc:>15.3f}{ri:>14.3f}{rg:>13.3f}")
allc = d["err"].corr(np.log1p(d[CT]), method="spearman")
print(f"\n        pooled rho(err, CTCF) = {allc:.3f}   "
      f"median within-decile = {np.median([r[2] for r in rows]):.3f}")

# =========================================================================================
print("\nTEST 3  partial correlation of error with CTCF given ATAC, GC, CpG and promoter")
print("        These are collinear with CTCF, so an unadjusted association proves little.")
cov = ["latac", "gc", "cpg_oe", "promoter"]
sub = d[["err", CT, IG] + cov].dropna()
X = np.column_stack([np.ones(len(sub))] + [sub[c].to_numpy() for c in cov])


def partial(y_name):
    y = sub["err"].to_numpy()
    x = np.log1p(sub[y_name].to_numpy())
    by, *_ = np.linalg.lstsq(X, y, rcond=None)
    bx, *_ = np.linalg.lstsq(X, x, rcond=None)
    ry, rx = y - X @ by, x - X @ bx
    return np.corrcoef(ry, rx)[0, 1]


print(f"        n = {len(sub):,}")
print(f"        partial rho(err, CTCF | ATAC, GC, CpG, promoter) = {partial(CT):+.3f}")
print(f"        partial rho(err, IgG  | same covariates)         = {partial(IG):+.3f}   <- control")
print(f"        raw     rho(err, CTCF)                           = {allc:+.3f}")

# =========================================================================================
print("\nTEST 4  attributable fraction: how much over-prediction survives removing CTCF-high")
print("        This is the decisive test. If the CTCF-low remainder is still grossly")
print("        over-predicted, CTCF explains a minority of the phenotype.")
over = d[STRATA[3][1]]
hi, lo = over[over[CT] > thr], over[over[CT] <= thr]
med = {c: d[c].median() for c in ["obs", "pred", "seq"]}
print(f"\n{'subset':<34}{'n':>8}{'% of tail':>11}{'obs fold':>10}{'pred fold':>11}"
      f"{'median err':>12}{'CTCF':>8}")
for lab, s in (("over-predicted 5%, all", over),
               ("  CTCF-high (> threshold)", hi),
               ("  CTCF-low", lo),
               ("typical (middle 50%)", typ)):
    print(f"{lab:<34}{len(s):>8,}{100.0 * len(s) / len(over):>11.1f}"
          f"{s['obs'].median() / med['obs']:>10.2f}{s['pred'].median() / med['pred']:>11.2f}"
          f"{s['err'].median():>12.3f}{s[CT].median():>8.2f}")

frac_hi = len(hi) / len(over)
err_hi, err_lo = hi["err"].median(), lo["err"].median()
print("\nVERDICT")
print(f"   CTCF-high elements are {100 * frac_hi:.0f}% of the over-predicted tail.")
if err_lo >= 0.7 * err_hi:
    print(f"   The CTCF-LOW remainder is over-predicted almost as badly ({err_lo:.3f} against")
    print(f"   {err_hi:.3f} for CTCF-high), so CTCF is largely a PASSENGER: the enrichment is")
    print("   real but it is not what causes the over-prediction. The phenotype is")
    print("   'accessible but unacetylated', of which CTCF sites are one instance.")
elif err_lo <= 0.4 * err_hi:
    print(f"   The CTCF-low remainder is much less over-predicted ({err_lo:.3f} against")
    print(f"   {err_hi:.3f}), so CTCF binding is doing real explanatory work.")
else:
    print(f"   Intermediate: CTCF-high {err_hi:.3f} against CTCF-low {err_lo:.3f}. CTCF")
    print("   contributes but does not account for the tail on its own.")
