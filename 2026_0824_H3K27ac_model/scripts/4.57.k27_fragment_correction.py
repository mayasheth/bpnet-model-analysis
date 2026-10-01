#!/usr/bin/env python3
"""Correct the transferred H3K27ac prediction with within-cell-ranked ATAC fragment features,
fitted in GM12878 and applied in K562, and paint the results for ABC.

WHY. F-039: element-level fragment features (4.51), ranked within each cell, predict H3K27ac
beyond ATAC (partial rho 0.48-0.52) and their weights carry K562 <-> GM12878 losing 0.01. F-005
fed fragment-size information as raw per-base input channels and the in-cell gain was null
transferred. This asks whether the ranked, element-level form gives a transferable gain in
the benchmark.

DESIGN. Deployment direction: GM12878-trained H3K27ac model applied to K562
(`predk27ac_gm12878_multimodal.bw`, ABC arm pred_gm12878_multimodal). The correction is fitted
in GM12878 on the MIRROR prediction, the K562-trained model applied to GM12878
(`predk27ac_k562_on_gm12878_multimodal.bw`), so in both cells the prediction entering the
correction is a transferred one and no in-cell H3K27ac leaks into K562. Target: observed
H3K27ac 5' counts over the candidate, ranked within the cell. Every predictor is a within-cell
percentile rank, so no library scale crosses cells.

    tx_repaint   transferred prediction alone, repainted here (control for the painting path)
    tx_base      OLS: prediction, ATAC and composition, NO fragment features (control: what
                 recalibrating the prediction in GM12878 alone buys)
    tx_frag      tx_base + the nine fragment features (the test)
    tx_frag_gbm  same predictors, gradient-boosted trees
    atac_frag    OLS: ATAC, composition and fragment features, no sequence model
    tx_breadth   tx_base + ATAC within +/- 1 kb of the centre (pct, pct^2), no fragment
                 lengths. 4.56 showed the H3K27ac-high group has accessibility spread over
                 +/- 1 kb rather than one sharp central peak, so this asks whether WHERE the
                 accessibility lies explains the gain without fragment LENGTH
    tx_breadth_frag  tx_base + breadth + fragment features

Predictors: prediction pct, pct^2, ATAC pct, pct^2, prediction x ATAC; GC and CpG o/e (+/- 250
bp), log10 TSS distance, promoter flag; fragment features from 4.51 ($SCRATCH/nucsig, union
region overlapping the candidate; median-imputed and flagged when absent). Each arm is painted
over the K562 ABC candidates (July K562_ATAC_only region set) as (pct score x 1000 + 1) /
width, so ABC's per-region sum is a positive monotone of the score; qnorm is rank-based.

Usage: 4.57.k27_fragment_correction.py [--out-prefix results/2026_0930_k27frag]
"""
import argparse
import os

import numpy as np
import pandas as pd
import pyBigWig
import pyfaidx
from scipy.stats import rankdata, spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor

D = "/oak/stanford/groups/engreitz/Users/sheth"
P = f"{D}/EP300_BPNet"
A = f"{P}/2026_0824_H3K27ac_model"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction"
PRED = f"{A}/data/abc_predicted"
CACHE = f"{A}/data/tf_activity_recovery"
GENOME = f"{D}/hg38_resources/hg38.fa"
TSS = (f"{D}/scE2G_temp/scE2G/resources/genome_annotations/"
       "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6")
SIZES = f"{ABC}/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv"
PK = "Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed"
NUC = f"{os.environ.get('SCRATCH', '/tmp')}/nucsig"
CELL = {
    "K562": dict(cands=f"{ABC}/results/2026_0721_h3k27ac_counting_comparison/K562_ATAC_only/{PK}",
                 pred=f"{PRED}/predk27ac_gm12878_multimodal.bw",
                 incell=f"{PRED}/predk27ac_k562_multimodal.bw",
                 atac=f"{P}/2026_0529_multimodal_p300_model/data/atac.bw",
                 k27=f"{A}/data/h3k27ac_5p_{{}}.bw"),
    "GM12878": dict(cands=f"{ABC}/results/2026_0924_p300_panel_obs/p300obs_GM12878/{PK}",
                    pred=f"{PRED}/predk27ac_k562_on_gm12878_multimodal.bw",
                    incell=None,
                    atac=f"{P}/2026_0606_GM12878_transferability/data/atac.bw",
                    k27=f"{A}/data/gm12878_h3k27ac_5p_{{}}.bw"),
}
KEY = ["chr", "start", "end"]
FEATS = ["nfr", "mono_flank", "mono_centre", "di_flank", "nfr_width", "positioning",
         "flank_ratio", "centre_occ", "short_frac"]

ap = argparse.ArgumentParser()
ap.add_argument("--out-prefix", default="results/2026_0930_k27frag")
a = ap.parse_args()
pct = lambda v: (rankdata(v) - 0.5) / len(v)


def bw_sum(path, g):
    bw = pyBigWig.open(path)
    L = bw.chroms()
    out = np.zeros(len(g))
    for i, (c, s, e) in enumerate(zip(g["chr"].values, g["start"].values, g["end"].values)):
        if c in L:
            v = bw.values(c, max(0, int(s)), min(int(e), L[c]), numpy=True)
            out[i] = np.nansum(np.abs(v))
    bw.close()
    return out


def map_to(reg, iv):
    """index of the interval in iv (sorted, non-overlapping per chr) overlapping each reg row"""
    out = np.full(len(reg), -1)
    for c, g in reg.groupby("chr"):
        e = iv[iv["chr"] == c].sort_values("start")
        if e.empty:
            continue
        es, ee, ei = e["start"].values, e["end"].values, e.index.values
        k = np.searchsorted(es, g["end"].values, side="left") - 1
        ok = (k >= 0)
        ok[ok] &= ee[k[ok]] > g["start"].values[ok]
        out[g.index.values[ok]] = ei[k[ok]]
    return out


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


U = pd.read_csv(f"{CACHE}/union_signals.tsv.gz", sep="\t", usecols=KEY)
tss = pd.read_csv(TSS, sep="\t", header=None, usecols=[0, 1, 2], names=KEY)
tss["m"] = (tss["start"] + tss["end"]) // 2
fa = pyfaidx.Fasta(GENOME, sequence_always_upper=True)

T = {}
for c, cf in CELL.items():
    g = pd.read_csv(cf["cands"], sep="\t", header=None, usecols=[0, 1, 2], names=KEY)
    g = g.sort_values(KEY).reset_index(drop=True)
    g["mid"] = (g["start"] + g["end"]) // 2
    g["prom"] = overlaps(g, tss).astype(float)
    dist = np.full(len(g), 1e7)
    for ch, gg in g.groupby("chr"):
        t = np.sort(tss.loc[tss["chr"] == ch, "m"].values)
        if len(t) > 1:
            k = np.clip(np.searchsorted(t, gg["mid"].values), 1, len(t) - 1)
            dist[gg.index.values] = np.minimum(np.abs(gg["mid"].values - t[k - 1]),
                                               np.abs(gg["mid"].values - t[k]))
    g["ldist"] = np.log10(dist + 1)
    gc, cpg = np.empty(len(g)), np.empty(len(g))
    for i, (ch, mid) in enumerate(zip(g["chr"].values, g["mid"].values)):
        s = str(fa[ch][max(0, mid - 250):mid + 250]) if ch in fa else ""
        n = max(sum(s.count(b) for b in "ACGT"), 1)
        nc, ng = s.count("C"), s.count("G")
        gc[i], cpg[i] = (nc + ng) / n, s.count("CG") * n / max(nc * ng, 1)
    g["gc"], g["cpg"] = gc, cpg
    g["pred"] = pct(bw_sum(cf["pred"], g))
    g["atac"] = pct(bw_sum(cf["atac"], g))
    g["atac1kb"] = pct(bw_sum(cf["atac"], pd.DataFrame({"chr": g["chr"], "start": g["mid"] - 1000,
                                                        "end": g["mid"] + 1000})))
    g["obs"] = pct(bw_sum(cf["k27"].format("plus"), g) + bw_sum(cf["k27"].format("minus"), g))
    if cf["incell"]:
        g["incell"] = pct(bw_sum(cf["incell"], g))
    f = pd.read_csv(f"{NUC}/{c}.tsv.gz", sep="\t").set_index("u")
    ui = map_to(g, U)
    has = np.isin(ui, f.index.values)
    g["has_frag"] = has.astype(float)
    for k in FEATS:
        v = np.full(len(g), np.nan)
        v[has] = f.loc[ui[has], k].values
        ok = np.isfinite(v)
        v[~ok] = np.nanmedian(v)
        g[k] = pct(v)
    T[c] = g
    print(f"{c}: {len(g):,} candidates, {has.mean():.1%} with fragment features; "
          f"Spearman with observed H3K27ac: transferred prediction {spearmanr(g['pred'], g['obs'])[0]:.3f}, "
          f"ATAC {spearmanr(g['atac'], g['obs'])[0]:.3f}", flush=True)


def design(g, use_pred, use_frag, use_breadth=False):
    cols = [np.ones(len(g)), g["atac"], g["atac"] ** 2, g["gc"], g["cpg"], g["ldist"], g["prom"]]
    if use_pred:
        cols += [g["pred"], g["pred"] ** 2, g["pred"] * g["atac"]]
    if use_breadth:
        cols += [g["atac1kb"], g["atac1kb"] ** 2, g["atac1kb"] * g["atac"]]
    if use_frag:
        cols += [g[k] for k in FEATS] + [g["has_frag"]]
    return np.column_stack(cols)


tr, te = T["GM12878"], T["K562"]
SC = {"tx_repaint": te["pred"].values}
for name, up, uf, ub in (("tx_base", True, False, False), ("tx_frag", True, True, False),
                         ("atac_frag", False, True, False), ("tx_breadth", True, False, True),
                         ("tx_breadth_frag", True, True, True)):
    w = np.linalg.lstsq(design(tr, up, uf, ub), tr["obs"].values, rcond=None)[0]
    SC[name] = design(te, up, uf, ub) @ w
Xtr = design(tr, True, True)[:, 1:]
gbm = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_leaf_nodes=31,
                                    l2_regularization=1.0, random_state=0).fit(Xtr, tr["obs"].values)
SC["tx_frag_gbm"] = gbm.predict(design(te, True, True)[:, 1:])

enh = te["prom"].values == 0
print("\n=== K562 Spearman with observed H3K27ac (weights fitted in GM12878) ===")
rows = []
refs = {"ATAC": te["atac"].values, "ATAC +/- 1 kb": te["atac1kb"].values, "in-cell K562 model (reference)": te["incell"].values}
for name, v in {**refs, **SC}.items():
    r_all, r_enh = spearmanr(v, te["obs"])[0], spearmanr(v[enh], te["obs"].values[enh])[0]
    rows.append(dict(score=name, spearman_all=r_all, spearman_nonpromoter=r_enh))
    print(f"  {name:<32} all {r_all:.3f}   non-promoter {r_enh:.3f}")
os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
pd.DataFrame(rows).to_csv(f"{a.out_prefix}.element_spearman.tsv", sep="\t", index=False)

sizes = dict(l.split()[:2] for l in open(SIZES))
order = [c for c in sizes if c in set(te["chr"])]
te_o = te.assign(_o=te["chr"].map({c: i for i, c in enumerate(order)})).sort_values(["_o", "start"])
assert (te_o.groupby("chr").apply(lambda x: (x["start"].values[1:] >= x["end"].values[:-1]).all())).all(), \
    "candidate regions overlap"
for name, v in SC.items():
    s = pd.Series(pct(v) * 1000 + 1, index=te.index).loc[te_o.index].values
    out = f"{PRED}/predk27frag_{name}.bw"
    bw = pyBigWig.open(out, "w")
    bw.addHeader([(c, int(sizes[c])) for c in order])
    width = (te_o["end"] - te_o["start"]).values
    bw.addEntries(te_o["chr"].tolist(), te_o["start"].tolist(), ends=te_o["end"].tolist(),
                  values=(s / width).astype(float).tolist())
    bw.close()
    print(f"wrote {out}")
te.assign(**SC).to_csv(f"{a.out_prefix}.k562_scores.tsv.gz", sep="\t", index=False)
