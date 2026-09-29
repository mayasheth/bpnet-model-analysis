#!/usr/bin/env python3
"""Which elements do the p300 models mis-rank, and what do those elements have in common?

WHY. Every intervention on the transferred p300 model has moved the CRISPR benchmark by at most
~0.01, while observed p300 sits ~0.09 higher (F-027, F-028). Most of that gap is how the tested
elements are ranked (F-028), and the transferred predictions are mostly ATAC re-encoded
(partial rho 0.28 against 0.44 in-cell, F-029). So some property of the elements determines
p300 in K562 and is invisible to the model. This screens for it.

LEADING HYPOTHESIS. The multi-cell panel (GM12878, HepG2, A549, MCF-7) has no erythroid cell
type, so the model never saw GATA1/TAL1 motifs recruit p300. At GATA-driven K562 enhancers it
can only fall back on ATAC. More generally: any TF active in K562 and weak in the panel.

RESIDUAL. For each ABC candidate region, within ABC's class (promoter / non-promoter), in
percentile points: r = rank(observed p300) - rank(predicted p300). r > 0 is UNDER-prediction.
Computed for the in-cell model, the multi-cell (K562 held out) model and the ranking-loss w=10
model. The transfer-specific residual is r_multicell - r_incell: where the transferred model
ranks an element lower than the in-cell model does.

FEATURES. JASPAR 2024 CORE vertebrate motifs (FIMO, p < 1e-4, on region centre +/- 250 bp),
ENCODE K562 peak overlaps (CTCF, EP300, H3K4me1, H3K27me3, H3K27ac), GC, CpG o/e, distance to
the nearest TSS.

TEST. For each feature, the OLS coefficient of r on the feature after partialling out ATAC rank
(linear + quadratic), GC, CpG o/e and log distance to TSS (Frisch-Waugh-Lovell, so one pass per
feature), in percentile points: "elements with this motif are under-predicted by X points at
matched accessibility and composition". Reported with t. With 132k regions almost everything is
nominally significant, so rank by effect size and read t only as a floor.

ALL CANDIDATE REGIONS are the screening population; the CRISPR-tested regions are summarised
separately and are too few (3,406) to screen 900 features on.

Usage: 4.41.p300_error_anatomy.py [--out-prefix results/2026_0925_p300_error_anatomy]
       [--motif-cache data/p300_error_anatomy/fimo_counts.npz]
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata

D = "/oak/stanford/groups/engreitz/Users/sheth"
P = f"{D}/EP300_BPNet"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction/results"
K562 = f"{D}/Data/ENCODE/K562"
GENOME = f"{D}/hg38_resources/hg38.fa"
JASPAR = f"{P}/reference/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"
TSS = (f"{D}/scE2G_temp/scE2G/resources/genome_annotations/"
       "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6")
CR = (f"{D}/CRISPR_comparison_v3/CRISPR_comparison/workflow/results/"
      "2026_0924_p300_rank/expt_pred_merged_annot.txt.gz")
OBS = f"{ABC}/2026_0905_p300_activity/p300obs_k562"
MODELS = {
    "incell": f"{ABC}/2026_0905_p300_activity/p300pred_k562_multimodal",
    "multicell": f"{ABC}/2026_0922_p300_multicell/p300pred_multicell_holdK562",
    "rank10": f"{ABC}/2026_0924_p300_rank/p300pred_multicell_holdK562_rank10",
}
PEAKS = {"CTCF": "ENCFF519CXF", "EP300": "ENCFF702XPO", "H3K4me1": "ENCFF135ZLM",
         "H3K27me3": "ENCFF323WOT", "H3K27ac": "ENCFF544LXB"}
KEY = ["chr", "start", "end"]
HALF = 250

ap = argparse.ArgumentParser()
ap.add_argument("--out-prefix", default="results/2026_0925_p300_error_anatomy")
ap.add_argument("--motif-cache", default="data/p300_error_anatomy/fimo_counts.npz")
ap.add_argument("--threshold", type=float, default=1e-4)
ap.add_argument("--extra", action="append", default=[], metavar="NAME=ABC_DIR",
                help="another prediction's ABC run to screen alongside (repeatable), e.g. a "
                     "Test 1 model; adds r_NAME and its coefficients to the tables")
ap.add_argument("--test-n", type=int, default=0,
                help="run on a random subset of N regions, write nothing (smoke test)")
a = ap.parse_args()
for e_ in a.extra:
    k_, v_ = e_.split("=", 1)
    MODELS[k_] = v_
EXTRA = [e_.split("=", 1)[0] for e_ in a.extra]


def load(d, col="H3K27ac.RPM"):
    return pd.read_csv(f"{d}/Neighborhoods/EnhancerList.txt", sep="\t",
                       usecols=KEY + ["class", col, "ATAC.RPM"])


def pct_within(df, col):
    out = np.empty(len(df))
    for sel in (df["prom"].values, ~df["prom"].values):
        out[sel] = 100 * (rankdata(df.loc[sel, col].values) - 0.5) / sel.sum()
    return out


def overlaps(reg, iv):
    """Boolean per region: overlaps any interval in iv (chr, start, end)."""
    hit = np.zeros(len(reg), bool)
    for c, g in reg.groupby("chr"):
        e = iv[iv["chr"] == c].sort_values("start")
        if e.empty:
            continue
        es, ee = e["start"].values, np.maximum.accumulate(e["end"].values)
        # any interval with start < region end and end > region start
        k = np.searchsorted(es, g["end"].values, side="left") - 1
        ok = k >= 0
        hit[g.index.values[ok]] = ee[k[ok]] > g["start"].values[ok]
    return hit


# ---------------------------------------------------------------- regions and residuals
reg = load(OBS).rename(columns={"H3K27ac.RPM": "obs_raw", "ATAC.RPM": "atac_raw"})
reg["prom"] = reg["class"].isin(["tss", "promoter"])
for name, d in MODELS.items():
    p = load(d).rename(columns={"H3K27ac.RPM": f"{name}_raw", "ATAC.RPM": "_a"})
    m = reg[KEY].merge(p[KEY + [f"{name}_raw", "_a"]], on=KEY, how="left", validate="one_to_one")
    assert m[f"{name}_raw"].notna().all() and np.allclose(m["_a"], reg["atac_raw"])
    reg[f"{name}_raw"] = m[f"{name}_raw"].values
reg["obs"] = pct_within(reg, "obs_raw")
reg["atac"] = pct_within(reg, "atac_raw")
for name in MODELS:
    reg[f"pred_{name}"] = pct_within(reg, f"{name}_raw")
    reg[f"r_{name}"] = reg["obs"] - reg[f"pred_{name}"]
reg["r_transfer"] = reg["r_multicell"] - reg["r_incell"]
reg["mid"] = (reg["start"] + reg["end"]) // 2
if a.test_n:
    reg = reg.sample(a.test_n, random_state=0).sort_values(KEY).reset_index(drop=True)
    a.motif_cache = "/nonexistent/skip.npz"
print(f"{len(reg):,} candidate regions ({int(reg.prom.sum()):,} promoter-class)", flush=True)

# ---------------------------------------------------------------- CRISPR-tested flags
cr = pd.read_csv(CR, sep="\t", usecols=["chrom", "chromStart", "chromEnd", "Regulated"])
cr["Regulated"] = cr["Regulated"].astype(str).str.upper().eq("TRUE")
iv = cr.rename(columns={"chrom": "chr", "chromStart": "start", "chromEnd": "end"})
reg["tested"] = overlaps(reg, iv.drop_duplicates(KEY))
reg["regulated"] = overlaps(reg, iv[iv["Regulated"]].drop_duplicates(KEY))

# ---------------------------------------------------------------- peaks, TSS distance
feats = {}
for mark, acc in PEAKS.items():
    pk = pd.read_csv(f"{K562}/{acc}.bed.gz", sep="\t", header=None, usecols=[0, 1, 2],
                     names=KEY)
    feats[f"peak:{mark}"] = overlaps(reg, pk).astype(float)
tss = pd.read_csv(TSS, sep="\t", header=None, usecols=[0, 1, 2], names=KEY)
tss["m"] = (tss["start"] + tss["end"]) // 2
dist = np.full(len(reg), 1e7)
for c, g in reg.groupby("chr"):
    t = np.sort(tss.loc[tss["chr"] == c, "m"].values)
    if not len(t):
        continue
    k = np.clip(np.searchsorted(t, g["mid"].values), 1, len(t) - 1)
    dist[g.index.values] = np.minimum(np.abs(g["mid"].values - t[k - 1]),
                                      np.abs(g["mid"].values - t[k]))
reg["log_dist_tss"] = np.log10(dist + 1)

# ---------------------------------------------------------------- sequences, GC, CpG, motifs
import pyfaidx  # noqa: E402

fa = pyfaidx.Fasta(GENOME, sequence_always_upper=True)
LUT = np.full(256, -1, np.int8)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i
idx = np.empty((len(reg), 2 * HALF), np.int8)
for i, (c, mid) in enumerate(zip(reg["chr"].values, reg["mid"].values)):
    s = str(fa[c][max(0, mid - HALF):mid + HALF]).encode()
    s = s.ljust(2 * HALF, b"N")
    idx[i] = LUT[np.frombuffer(s, np.uint8)]
valid = idx >= 0
gc = ((idx == 1) | (idx == 2)).sum(1) / np.maximum(valid.sum(1), 1)
cg = ((idx[:, :-1] == 1) & (idx[:, 1:] == 2)).sum(1)
nc, ng = (idx == 1).sum(1), (idx == 2).sum(1)
reg["gc"] = gc
reg["cpg_oe"] = cg * valid.sum(1) / np.maximum(nc * ng, 1)
print("sequences loaded", flush=True)

if os.path.exists(a.motif_cache):
    z = np.load(a.motif_cache, allow_pickle=True)
    counts, mnames = z["counts"], list(z["names"])
    assert counts.shape[0] == len(reg)
    print(f"motif counts from cache: {counts.shape}", flush=True)
else:
    from memelite import fimo
    from memelite.io import read_meme
    motifs = read_meme(JASPAR)
    mnames = list(motifs)
    X = np.zeros((len(reg), 4, 2 * HALF), np.float32)
    ii, jj = np.nonzero(valid)
    X[ii, idx[ii, jj], jj] = 1.0
    counts = np.zeros((len(reg), len(mnames)), np.uint8)
    B = 20000
    for s in range(0, len(reg), B):
        # return_counts=True sums over ALL sequences (one number per motif), so count hits
        # per sequence from the hit tables instead; sequence_name is the row within the batch.
        hits = fimo(motifs, X[s:s + B], threshold=a.threshold)
        blk = np.zeros((min(B, len(reg) - s), len(mnames)), np.int32)
        for j, h in enumerate(hits):
            if len(h):
                np.add.at(blk[:, j], h["sequence_name"].astype(int).values, 1)
        counts[s:s + B] = np.clip(blk, 0, 255)
        print(f"  fimo {min(s + B, len(reg)):,}/{len(reg):,}", flush=True)
    del X
    if not a.test_n:
        os.makedirs(os.path.dirname(a.motif_cache), exist_ok=True)
        np.savez_compressed(a.motif_cache, counts=counts, names=np.array(mnames, dtype=object))
for j, n in enumerate(mnames):
    feats[f"motif:{n}"] = (counts[:, j] > 0).astype(float)
print(f"{len(feats):,} features", flush=True)

# ---------------------------------------------------------------- per-class screen
RES = ["r_incell", "r_multicell", "r_rank10", "r_transfer"] + [f"r_{k}" for k in EXTRA]
rows = []
for cls, sel in (("nonpromoter", ~reg["prom"].values), ("promoter", reg["prom"].values)):
    g = reg[sel]
    Z = np.column_stack([np.ones(len(g)), g["atac"], g["atac"] ** 2, g["gc"], g["cpg_oe"],
                         g["log_dist_tss"]])
    Q, _ = np.linalg.qr(Z)
    R = {k: g[k].values - Q @ (Q.T @ g[k].values) for k in RES}
    df_ = len(g) - Z.shape[1] - 1
    for fname, f in feats.items():
        fv = f[sel]
        n1 = int(fv.sum())
        if n1 < (5 if a.test_n else 50) or n1 > len(g) - (5 if a.test_n else 50):
            continue
        ft = fv - Q @ (Q.T @ fv)
        ss = ft @ ft
        row = dict(cls=cls, feature=fname, n_with=n1, frac=n1 / len(g))
        for k in RES:
            b = (ft @ R[k]) / ss
            e = R[k] - b * ft
            se = np.sqrt((e @ e) / df_ / ss)
            row[f"coef_{k}"], row[f"t_{k}"] = b, b / se
        # raw mean residuals, no covariates, for reading alongside
        row["mean_r_multicell_with"] = g["r_multicell"].values[fv > 0].mean()
        row["mean_r_multicell_without"] = g["r_multicell"].values[fv == 0].mean()
        rows.append(row)
    print(f"screened {cls}", flush=True)
S = pd.DataFrame(rows)
if a.test_n:
    a.out_prefix = os.path.join(os.environ.get("SCRATCH", "/tmp"), "p300_error_anatomy_test")
os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
S.to_csv(f"{a.out_prefix}.features.tsv", sep="\t", index=False)

pd.set_option("display.width", 220)
fmt = lambda v: f"{v:+.2f}" if isinstance(v, float) else v
show = ["feature", "n_with", "coef_r_incell", "coef_r_multicell", "coef_r_rank10",
        "coef_r_transfer", "t_r_transfer"] + [f"coef_r_{k}" for k in EXTRA]
for cls in ("nonpromoter", "promoter"):
    s = S[S.cls == cls]
    print(f"\n=== {cls}: peaks and composition (percentile points, + = under-predicted) ===")
    print(s[s.feature.str.startswith("peak:")][show].to_string(index=False, float_format=fmt))
    print(f"\n=== {cls}: top 25 motifs the TRANSFERRED model under-predicts relative to in-cell ===")
    print(s[s.feature.str.startswith("motif:")].sort_values("coef_r_transfer", ascending=False)
          .head(25)[show].to_string(index=False, float_format=fmt))
    print(f"\n=== {cls}: top 15 motifs the transferred model under-predicts, absolute ===")
    print(s[s.feature.str.startswith("motif:")].sort_values("coef_r_multicell", ascending=False)
          .head(15)[show].to_string(index=False, float_format=fmt))
    print(f"\n=== {cls}: top 15 motifs the transferred model OVER-predicts ===")
    print(s[s.feature.str.startswith("motif:")].sort_values("coef_r_multicell")
          .head(15)[show].to_string(index=False, float_format=fmt))

# ---------------------------------------------------------------- tested elements
t = reg[reg["tested"]]
print(f"\n=== CRISPR-tested regions: {len(t):,} ({int(t.regulated.sum())} overlap a regulated "
      "pair) ===")
for grp, s in (("regulated", t[t.regulated]), ("not regulated", t[~t.regulated])):
    print(f"  {grp:<14} n={len(s):>5}  mean r: in-cell {s.r_incell.mean():+.1f}  "
          f"multi-cell {s.r_multicell.mean():+.1f}  rank10 {s.r_rank10.mean():+.1f}  "
          f"transfer-specific {s.r_transfer.mean():+.1f}")
reg[KEY + ["class", "tested", "regulated", "obs", "atac", "gc", "cpg_oe", "log_dist_tss"]
    + [f"pred_{m}" for m in MODELS] + RES].to_csv(f"{a.out_prefix}.regions.tsv.gz", sep="\t",
                                                  index=False)
print(f"\nwrote {a.out_prefix}.features.tsv and .regions.tsv.gz")
