#!/usr/bin/env python3
"""Gate for feature B: do per-element TOBIAS footprints (2.2, validated in 2.3 / F-048) track
observed K562 p300 beyond the ATAC scE2G already uses?

REGIONS. scE2G's candidate regions for the Xu run (the 155,134 TOBIAS ran on), minus promoters
(GENCODE TSS +/- 250 bp, as 4.60). Signals over the region midpoint +/- 250 bp: Xu ATAC (TOBIAS's
uncorrected insertion track, i.e. the multiome ATAC scE2G sees), ENCODE K562 ATAC and p300
(both strands) from config/p300_panel.json, as in 4.43 / 4.60.

FEATURES per element, per motif, from BINDetect's site table: n_sites, max motif match score
(TFBS_score), and the footprint score `K562_Xu_score`: fp_max and fp_sum over the element's sites.
Motifs: the erythroid set GATA1, GATA1::TAL1, NFE2, MAF::NFE2, TAL1::TCF3 and their sum
(`erythroid`); CTCF and HNF4A (not expressed in K562) as controls. Plus TF-agnostic `fp_any`, the
max of ScoreBigwig's footprint track over the window.

CONTROLS Z, as 4.60: Xu ATAC percentile and its square, GC, CpG o/e, log10 distance to TSS. All
signals percentile-ranked, so every number is a partial Spearman. 95% intervals: 200 bootstraps
over elements of the residual correlation.

TEST 1, all non-promoter elements (no site -> lowest rank): feature vs p300 beyond Z.
TEST 2, elements with >= 1 site of the motif: fp_max vs p300 beyond Z + n_sites + max motif score,
i.e. does the footprint add anything over motif presence and strength. Feature B passes if
TEST 2's interval excludes 0 for the erythroid set. Reference scale: observed H3K27ac 0.376 and
H3K4me1 0.542 partial rho with p300 in K562 (F-045, different region set).

Run in the tobias pixi env. Usage: 2.4.footprint_p300_gate.py [--out-prefix results/2026_1003_footprint_p300_gate]
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd
import pyBigWig
import pysam
from scipy.stats import rankdata

D = "/oak/stanford/groups/engreitz/Users/sheth"
TB = os.path.join(os.environ["SCRATCH"], "tobias_k562_xu")
PANEL = json.load(open(f"{D}/EP300_BPNet/2026_0824_H3K27ac_model/config/p300_panel.json"))["K562"]
TSS = f"{D}/scE2G_temp/scE2G/resources/genome_annotations/CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6"
GENOME = f"{D}/hg38_resources/hg38.fa"
MOTIFS = {"GATA1": "MA0035.5", "GATA1::TAL1": "MA0140.3", "NFE2": "MA0841.2",
          "MAF::NFE2": "MA0501.2", "TAL1::TCF3": "MA0091.2", "CTCF": "MA0139.2", "HNF4A": "MA0114.5"}
ERY = ["GATA1", "GATA1::TAL1", "NFE2", "MAF::NFE2", "TAL1::TCF3"]
HALF = 250
KEY = ["chr", "start", "end"]

ap = argparse.ArgumentParser()
ap.add_argument("--out-prefix", default="results/2026_1003_footprint_p300_gate")
ap.add_argument("--n-boot", type=int, default=200)
ap.add_argument("--z-encode-atac", action="store_true",
                help="also control for ENCODE K562 ATAC (percentile + square): soaks up Xu ATAC noise")
ap.add_argument("--score", choices=["tobias", "print"], default="tobias",
                help="print: PRINT (3.1) scores at the same BINDetect sites (fp_any stays TOBIAS)")
a = ap.parse_args()
if a.score == "print":
    a.out_prefix += ".print"
from footprint_scores import site_scores  # noqa: E402
assert os.path.exists(f"{TB}/tobias_complete.txt"), "TOBIAS 2.2 not complete"
if a.z_encode_atac:
    a.out_prefix += ".z_encode_atac"
print(f"Z: Xu ATAC{' + ENCODE ATAC' if a.z_encode_atac else ''} + GC + CpG o/e + log dist TSS", flush=True)


def overlaps(reg, iv):
    hit = np.zeros(len(reg), bool)
    for c, g in reg.groupby("chr"):
        e = iv[iv["chr"] == c].sort_values("start")
        if e.empty:
            continue
        es, ee = e["start"].values, np.maximum.accumulate(e["end"].values)
        k = np.searchsorted(es, g["end"].values, side="left") - 1
        ok = k >= 0
        hit[g.index.values[ok]] = ee[k[ok]] > g["start"].values[ok]
    return hit


def window(path, how="sum", absval=False):
    bw = pyBigWig.open(path)
    ch = bw.chroms()
    out = np.zeros(len(U))
    for i, (c, m) in enumerate(zip(U["chr"].values, U["mid"].values)):
        if c not in ch:
            continue
        v = np.nan_to_num(bw.values(c, max(0, m - HALF), min(m + HALF, ch[c]), numpy=True))
        v = np.abs(v) if absval else v
        out[i] = v.sum() if how == "sum" else v.max()
    bw.close()
    return out


# ---------------------------------------------------------------- regions and signals
U = pd.read_csv(f"{TB}/peaks.bed", sep="\t", header=None, names=KEY)
U["mid"] = (U["start"] + U["end"]) // 2
tss = pd.read_csv(TSS, sep="\t", header=None, usecols=[0, 1, 2], names=KEY)
U["prom"] = overlaps(U, tss)
tss["m"] = (tss["start"] + tss["end"]) // 2
dist = np.full(len(U), 1e7)
for c, g in U.groupby("chr"):
    t = np.sort(tss.loc[tss["chr"] == c, "m"].values)
    if len(t) > 1:
        k = np.clip(np.searchsorted(t, g["mid"].values), 1, len(t) - 1)
        dist[g.index.values] = np.minimum(np.abs(g["mid"].values - t[k - 1]),
                                          np.abs(g["mid"].values - t[k]))
U["log_dist_tss"] = np.log10(dist + 1)
fa = pysam.FastaFile(GENOME)
gc, cpg = np.empty(len(U)), np.empty(len(U))
for i, (ch, mid) in enumerate(zip(U["chr"].values, U["mid"].values)):
    s = fa.fetch(ch, max(0, mid - HALF), mid + HALF).upper()
    n = max(sum(s.count(b) for b in "ACGT"), 1)
    nc, ng = s.count("C"), s.count("G")
    gc[i], cpg[i] = (nc + ng) / n, s.count("CG") * n / max(nc * ng, 1)
U["gc"], U["cpg_oe"] = gc, cpg
U["atac_xu"] = window(f"{TB}/atacorrect/k562_xu_uncorrected.bw")
U["atac_encode"] = window(PANEL["accessibility_bw"])
U["p300"] = window(PANEL["signal_plus_bw"], absval=True) + window(PANEL["signal_minus_bw"], absval=True)
U["fp_any"] = window(f"{TB}/k562_xu_footprints.bw", how="max")
print(f"{len(U):,} regions, {int((~U['prom']).sum()):,} non-promoter", flush=True)

# ---------------------------------------------------------------- per-element motif features
ix = pd.Series(np.arange(len(U)), index=pd.MultiIndex.from_frame(U[KEY]))
assert ix.index.is_unique, "duplicate regions in peaks.bed"
for name, mid in MOTIFS.items():
    (ov,) = glob.glob(f"{TB}/bindetect/*_{mid}/*_{mid}_overview.txt")
    s = site_scores(pd.read_csv(ov, sep="\t"), a.score, name)
    s["r"] = ix.reindex(pd.MultiIndex.from_frame(
        s[["peak_chr", "peak_start", "peak_end"]].set_axis(KEY, axis=1))).values
    assert s["r"].notna().all(), f"{name}: sites outside the region set"
    agg = s.groupby("r").agg(n=("K562_Xu_score", "size"), motif_max=("TFBS_score", "max"),
                             fp_max=("K562_Xu_score", "max"), fp_sum=("K562_Xu_score", "sum"))
    for col in agg:
        v = np.zeros(len(U))
        v[agg.index.values.astype(int)] = agg[col].values
        U[f"{name}_{col}"] = v
U["erythroid_n"] = U[[f"{m}_n" for m in ERY]].sum(1)
U["erythroid_motif_max"] = U[[f"{m}_motif_max" for m in ERY]].max(1)
U["erythroid_fp_sum"] = U[[f"{m}_fp_sum" for m in ERY]].sum(1)
U["erythroid_fp_max"] = U[[f"{m}_fp_max" for m in ERY]].max(1)

# ---------------------------------------------------------------- measures
pct = lambda v: (rankdata(v) - 0.5) / len(v)
rng = np.random.default_rng(0)


def resid(y, Z):
    Q, _ = np.linalg.qr(Z)
    return y - Q @ (Q.T @ y)


def partial(x, y, Z):
    rx, ry = resid(x, Z), resid(y, Z)
    est = float(np.corrcoef(rx, ry)[0, 1])
    bs = []
    for _ in range(a.n_boot):
        i = rng.integers(0, len(rx), len(rx))
        bs.append(np.corrcoef(rx[i], ry[i])[0, 1])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return est, lo, hi


def Zmat(g, extra=()):
    at = pct(g["atac_xu"].values)
    cols = [np.ones(len(g)), at, at ** 2, g["gc"].values, g["cpg_oe"].values, g["log_dist_tss"].values]
    if a.z_encode_atac:
        ae = pct(g["atac_encode"].values)
        cols += [ae, ae ** 2]
    return np.column_stack(cols + [pct(g[c].values) for c in extra])


G = U[~U["prom"]].reset_index(drop=True)
Z = Zmat(G)
p3 = pct(G["p300"].values)
rows = []
for f in ["atac_encode", "fp_any"] + [f"{m}_{k}" for m in list(MOTIFS) + ["erythroid"]
                                      for k in ("n", "fp_max", "fp_sum")]:
    est, lo, hi = partial(pct(G[f].values), p3, Z)
    rows.append(dict(feature=f, n=len(G), frac_nonzero=float((G[f] > 0).mean()),
                     partial_rho=est, ci_lo=lo, ci_hi=hi))
T1 = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print("\n=== test 1: partial Spearman with p300 beyond Xu ATAC + composition, all non-promoter elements ===")
print(T1.round(3).to_string(index=False))

rows = []
for m in list(MOTIFS) + ["erythroid"]:
    g = G[G[f"{m}_n"] > 0].reset_index(drop=True)
    Zg = Zmat(g, extra=(f"{m}_n", f"{m}_motif_max"))
    est, lo, hi = partial(pct(g[f"{m}_fp_max"].values), pct(g["p300"].values), Zg)
    rows.append(dict(motif=m, n_elements=len(g), partial_rho_fp_max=est, ci_lo=lo, ci_hi=hi))
T2 = pd.DataFrame(rows)
print("\n=== test 2: fp_max beyond Xu ATAC + composition + motif count + motif strength, elements with the motif ===")
print(T2.round(3).to_string(index=False))

ery = T2.set_index("motif").loc["erythroid"]
print(f"\nGATE (erythroid fp_max beyond motif, interval excludes 0): "
      f"{'PASS' if ery['ci_lo'] > 0 else 'FAIL'}  {ery['partial_rho_fp_max']:.3f} "
      f"[{ery['ci_lo']:.3f}, {ery['ci_hi']:.3f}]")

os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
T1.to_csv(f"{a.out_prefix}.all_elements.tsv", sep="\t", index=False)
T2.to_csv(f"{a.out_prefix}.beyond_motif.tsv", sep="\t", index=False)
U.to_csv(f"{a.out_prefix}.element_features.tsv.gz", sep="\t", index=False)
print(f"\nwrote {a.out_prefix}.all_elements.tsv, .beyond_motif.tsv, .element_features.tsv.gz")
