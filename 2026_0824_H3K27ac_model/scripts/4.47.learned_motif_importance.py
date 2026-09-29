#!/usr/bin/env python3
"""Learn a motif's p300 importance in a cell from that cell's ATAC and sequence alone.

WHY. F-033: ATAC recovers the TF activities that matter for K562 (GATA, CEBP) but no hand-made
estimator is good enough. Raw within-cell accessibility association ranks CTCF first; the
specificity estimators remove CTCF but are relative to the reference panel, so a new cell type
similar to a training cell would have its own important TFs cancelled. This learns the mapping
instead, from features that describe how each motif behaves in the cell's OWN ATAC, with no
reference to the other cell types, and no motif identity (so it cannot memorise "CTCF is bad").

EXAMPLES. One per (cell type, motif): ~780 motifs x 5 panel cell types.
TARGET. 4.43's truth: the motif's coefficient for p300 beyond ATAC in that cell (standardised
within cell).
FEATURES, per (cell, motif), from the cell's ATAC and the sequence only:
    assoc        within-cell ATAC association (4.43's E_within)
    tail         coefficient of motif presence on "in the cell's top 10% most accessible"
    prom_frac    fraction of the motif's regions in this cell that are promoters
    prevalence   fraction of the cell's regions carrying the motif
    ctcf_cooc    fraction of the motif's regions that also carry a CTCF motif
    fp_core      ATAC at the motif core relative to its 50-100 bp flanks (footprint depth)
    fp_flank     ATAC at the flanks relative to 150-200 bp background (flank height)
    motif_gc     GC fraction of the motif's consensus (sequence only)
Cell-varying features are standardised within cell.

EVALUATION. Leave one cell type out: fit gradient boosting on four, predict the fifth. Scored per
cell (Spearman across motifs, top-20 overlap) against two baselines: raw assoc, and assoc with the
CTCF-family motifs removed by hand. For K562 also: the K562-minus-panel contrast of the prediction
against the true contrast, and against 4.41's transfer-specific residuals.

Usage: 4.47.learned_motif_importance.py [--out-prefix results/2026_0928_learned_motif_importance]
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import rankdata, spearmanr

D = "/oak/stanford/groups/engreitz/Users/sheth"
P = f"{D}/EP300_BPNet"
A = f"{P}/2026_0824_H3K27ac_model"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction/results"
CANDS = {"K562": f"{ABC}/2026_0905_p300_activity/p300obs_k562",
         **{c: f"{ABC}/2026_0924_p300_panel_obs/p300obs_{c.replace('-', '_')}"
            for c in ("GM12878", "HepG2", "A549", "MCF-7")}}
CELLS = list(CANDS)
PANEL = json.load(open(f"{A}/config/p300_panel.json"))
GENOME = f"{D}/hg38_resources/hg38.fa"
JASPAR = f"{P}/reference/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"
TSS = (f"{D}/scE2G_temp/scE2G/resources/genome_annotations/"
       "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6")
REC = f"{A}/results/2026_0926_tf_activity_recovery.motifs.tsv"
ANAT = f"{A}/results/2026_0925_p300_error_anatomy.features.tsv"
CACHE = f"{A}/data/tf_activity_recovery"
HITS = os.path.join(os.environ.get("SCRATCH", "/tmp"), "motif_hits")
KEY = ["chr", "start", "end"]
HALF = 250

ap = argparse.ArgumentParser()
ap.add_argument("--out-prefix", default="results/2026_0928_learned_motif_importance")
ap.add_argument("--n-fp", type=int, default=2000, help="hits per motif per cell for footprints")
a = ap.parse_args()
rng = np.random.default_rng(0)


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


def fwl_coefs(y, covs, F):
    Q, _ = np.linalg.qr(np.column_stack([np.ones(len(y))] + covs))
    yr = y - Q @ (Q.T @ y)
    Fr = F - Q @ (Q.T @ F)
    return (Fr * yr[:, None]).sum(0) / np.maximum((Fr * Fr).sum(0), 1e-12)


pct = lambda v: 100 * (rankdata(v) - 0.5) / len(v)

# ---------------------------------------------------------------- 4.43's union set and caches
S = pd.read_csv(f"{CACHE}/union_signals.tsv.gz", sep="\t")
U = S[KEY].copy()
U["mid"] = (U["start"] + U["end"]) // 2
z = np.load(f"{CACHE}/union_fimo_presence.npz", allow_pickle=True)
M_all, names_all = z["presence"], list(z["names"])
assert M_all.shape[0] == len(U), "presence cache does not match the union signal cache"
R = pd.read_csv(REC, sep="\t")
mnames = list(R["motif"])
M = M_all[:, [names_all.index(n) for n in mnames]]
print(f"{len(U):,} union regions, {len(mnames)} motifs", flush=True)
for c, d in CANDS.items():
    U[f"in_{c}"] = overlaps(U, pd.read_csv(f"{d}/Neighborhoods/EnhancerList.txt", sep="\t",
                                           usecols=KEY))
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

import pyfaidx  # noqa: E402

fa = pyfaidx.Fasta(GENOME, sequence_always_upper=True)
gc = np.empty(len(U))
cpg = np.empty(len(U))
for i, (c, mid) in enumerate(zip(U["chr"].values, U["mid"].values)):
    s = str(fa[c][max(0, mid - HALF):mid + HALF])
    n = max(sum(s.count(b) for b in "ACGT"), 1)
    nc, ng = s.count("C"), s.count("G")
    gc[i] = (nc + ng) / n
    cpg[i] = s.count("CG") * n / max(nc * ng, 1)
U["gc"], U["cpg_oe"] = gc, cpg
print("union rebuilt", flush=True)

ctcf_cols = [i for i, n in enumerate(mnames) if n.split(" ", 1)[1].upper().startswith("CTCF")]
has_ctcf = M[:, ctcf_cols].any(1)

# motif consensus GC, from the PWMs
from memelite.io import read_meme  # noqa: E402

pwm = read_meme(JASPAR)
motif_gc = np.array([float(pwm[n][1:3].sum() / pwm[n].sum()) for n in mnames])

# ---------------------------------------------------------------- per-cell features
feat_rows = []
for c in CELLS:
    sel = U[f"in_{c}"].values & ~U["prom"].values
    selall = U[f"in_{c}"].values
    g = U[sel]
    atac = pct(S.loc[sel, f"atac_{c}"].values)
    comp = [g["gc"].values, g["cpg_oe"].values, g["log_dist_tss"].values]
    F = M[sel].astype(float)
    assoc = fwl_coefs(atac, comp, F)
    tail = fwl_coefs((atac >= 90).astype(float), comp, F)
    Mc = M[selall]
    n_m = np.maximum(Mc.sum(0), 1)
    prom_frac = (Mc & U["prom"].values[selall, None]).sum(0) / n_m
    prevalence = Mc.mean(0)
    ctcf_cooc = (Mc & has_ctcf[selall, None]).sum(0) / n_m
    feat_rows.append(pd.DataFrame({"cell": c, "motif": mnames, "assoc": assoc, "tail": tail,
                                   "prom_frac": prom_frac, "prevalence": prevalence,
                                   "ctcf_cooc": ctcf_cooc, "motif_gc": motif_gc}))
    print(f"features: {c}", flush=True)
Fe = pd.concat(feat_rows, ignore_index=True)

# ---------------------------------------------------------------- footprints from base-level ATAC
# Genome-wide hits (4.44), restricted to each cell's candidate regions, up to n_fp per motif.
motif_index = {n: i for i, n in enumerate(names_all)}
want = np.array([motif_index[n] for n in mnames])
fp_core = {c: np.full(len(mnames), np.nan) for c in CELLS}
fp_flank = {c: np.full(len(mnames), np.nan) for c in CELLS}
files = sorted(glob.glob(f"{HITS}/chr*.npz"))
assert files, f"no motif hits in {HITS}"
genome_len = sum(int(np.load(f, allow_pickle=True)["chrom_len"]) for f in files)
samp = {c: [] for c in CELLS}                                  # (motif, chrom, centre) blocks
for f in files:
    h = np.load(f, allow_pickle=True)
    chrom, L = os.path.basename(f)[:-4], int(h["chrom_len"])
    lut = np.full(len(names_all), -1)
    lut[want] = np.arange(len(want))
    j = lut[h["motif"].astype(np.int64)]
    pos = (h["start"].astype(np.int64) + h["end"].astype(np.int64)) // 2
    ok = j >= 0
    j, pos = j[ok], pos[ok]
    quota = max(1, int(np.ceil(a.n_fp * L / genome_len * 1.5)))
    for c in CELLS:
        reg = U[(U["chr"] == chrom) & U[f"in_{c}"]].sort_values("start")
        if reg.empty:
            continue
        rs, re_ = reg["start"].values, reg["end"].values
        k = np.searchsorted(rs, pos, side="right") - 1
        inreg = (k >= 0) & (pos < re_[np.clip(k, 0, None)])
        jj, pp = j[inreg], pos[inreg]
        o = rng.permutation(len(jj))                          # random order, then
        jj, pp = jj[o], pp[o]
        o = np.argsort(jj, kind="stable")                     # group by motif, keep order
        jj, pp = jj[o], pp[o]
        starts = np.searchsorted(jj, jj, side="left")
        keep = (np.arange(len(jj)) - starts) < quota          # first `quota` per motif
        samp[c].append((jj[keep], np.full(int(keep.sum()), chrom, dtype=object), pp[keep]))
    print(f"footprint sampling: {chrom}", flush=True)
acc = {}
for c in CELLS:
    jj = np.concatenate([x[0] for x in samp[c]])
    cc = np.concatenate([x[1] for x in samp[c]])
    pp = np.concatenate([x[2] for x in samp[c]])
    acc[c] = [list(zip(cc[jj == m_], pp[jj == m_])) for m_ in range(len(mnames))]
for c in CELLS:
    bw = pyBigWig.open(PANEL[c]["accessibility_bw"])
    ch = bw.chroms()
    for m_ in range(len(mnames)):
        pts = acc[c][m_]
        if len(pts) > a.n_fp:
            pts = [pts[i] for i in rng.choice(len(pts), a.n_fp, replace=False)]
        prof = []
        for chrom, p_ in pts:
            if chrom not in ch or p_ - 200 < 0 or p_ + 200 > ch[chrom]:
                continue
            prof.append(np.nan_to_num(bw.values(chrom, p_ - 200, p_ + 200, numpy=True)))
        if len(prof) < 50:
            continue
        pr = np.mean(prof, 0)                                     # 400 bp, centre at 200
        core = pr[190:210].mean()
        flank = np.r_[pr[100:150], pr[250:300]].mean()
        bg = np.r_[pr[0:50], pr[350:400]].mean()
        fp_core[c][m_] = np.log((core + 1e-3) / (flank + 1e-3))
        fp_flank[c][m_] = np.log((flank + 1e-3) / (bg + 1e-3))
    bw.close()
    print(f"footprints: {c}", flush=True)
Fe["fp_core"] = np.concatenate([fp_core[c] for c in CELLS])
Fe["fp_flank"] = np.concatenate([fp_flank[c] for c in CELLS])

# ---------------------------------------------------------------- target, standardisation
Fe["truth"] = np.concatenate([R[f"truth_{c}"].values for c in CELLS])
CELLFEATS = ["assoc", "tail", "prom_frac", "prevalence", "ctcf_cooc", "fp_core", "fp_flank"]
FEATS = CELLFEATS + ["motif_gc"]
for col in CELLFEATS + ["truth"]:
    Fe[col] = Fe.groupby("cell")[col].transform(lambda v: (v - v.mean()) / v.std())

# ---------------------------------------------------------------- leave one cell out
from sklearn.ensemble import HistGradientBoostingRegressor  # noqa: E402

Fe["pred"] = np.nan
for c in CELLS:
    tr, te = Fe["cell"] != c, Fe["cell"] == c
    gb = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_leaf_nodes=15,
                                       min_samples_leaf=20, random_state=0)
    gb.fit(Fe.loc[tr, FEATS].values, Fe.loc[tr, "truth"].values)
    Fe.loc[te, "pred"] = gb.predict(Fe.loc[te, FEATS].values)

is_ctcf = np.isin(np.arange(len(mnames)), ctcf_cols)
top = lambda v, k=20: set(np.argsort(-v)[:k])
rows = []
for c in CELLS:
    s = Fe[Fe["cell"] == c]
    tv = s["truth"].values
    no_ctcf = s["assoc"].values.copy()
    no_ctcf[is_ctcf] = no_ctcf.min() - 1
    for nm, est in (("raw ATAC association", s["assoc"].values),
                    ("ATAC association, CTCF removed by hand", no_ctcf),
                    ("learned (leave this cell out)", s["pred"].values)):
        rows.append(dict(cell=c, estimate=nm, spearman=spearmanr(est, tv)[0],
                         top20_overlap=len(top(est) & top(tv)),
                         ctcf_best_rank=int(rankdata(-est)[is_ctcf].min())))
Rt = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print("\n=== estimate vs truth, per held-out cell ===")
print(Rt.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

panel = [c for c in CELLS if c != "K562"]
get = lambda col, c: Fe.loc[Fe["cell"] == c, col].values
dT = get("truth", "K562") - np.mean([get("truth", c) for c in panel], 0)
an = pd.read_csv(ANAT, sep="\t")
an = an[(an.cls == "nonpromoter") & an.feature.str.startswith("motif:")].copy()
rt = an.assign(motif=an.feature.str[6:]).set_index("motif")["coef_r_transfer"].reindex(mnames)
ok = rt.notna().values
print("\n=== K562 minus panel ===")
for col, nm in (("assoc", "raw ATAC association"), ("pred", "learned")):
    dE = get(col, "K562") - np.mean([get(col, c) for c in panel], 0)
    print(f"  {nm:<22} vs true contrast {spearmanr(dE, dT)[0]:+.3f} (top-20 shared "
          f"{len(top(dE) & top(dT))});  vs 4.41 transfer residual "
          f"{spearmanr(dE[ok], rt.values[ok])[0]:+.3f}")
gi = [i for i, n in enumerate(mnames) if "GATA" in n.upper() or "TRPS1" in n.upper()]
print("\nGATA family, mean standardised score:")
for col in ("truth", "assoc", "pred"):
    print(f"  {col:<6} " + ", ".join(f"{c} {get(col, c)[gi].mean():+.2f}" for c in CELLS))

gb = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_leaf_nodes=15,
                                   min_samples_leaf=20, random_state=0)
gb.fit(Fe[FEATS].values, Fe["truth"].values)
from sklearn.inspection import permutation_importance  # noqa: E402

pi = permutation_importance(gb, Fe[FEATS].values, Fe["truth"].values, n_repeats=10,
                            random_state=0)
print("\nfeature importance (permutation, all cells):")
for f_, v in sorted(zip(FEATS, pi.importances_mean), key=lambda x: -x[1]):
    print(f"  {f_:<11} {v:.3f}")

os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
wide = Fe.pivot(index="motif", columns="cell", values="pred").add_prefix("learned_").reset_index()
wide.to_csv(f"{a.out_prefix}.motifs.tsv", sep="\t", index=False)
Fe.to_csv(f"{a.out_prefix}.features.tsv.gz", sep="\t", index=False)
Rt.to_csv(f"{a.out_prefix}.scores.tsv", sep="\t", index=False)
print(f"\nwrote {a.out_prefix}.motifs.tsv (per-cell leave-one-out importances), .features, .scores")
