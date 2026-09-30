#!/usr/bin/env python3
"""Does the nucleosome arrangement around an element carry p300 signal beyond its accessibility,
and does it pick out the p300-important motifs better than accessibility does?

WHY. p300 acetylates the nucleosomes flanking an enhancer's nucleosome-free region. Active
enhancers are expected to have a wider nucleosome-free region and well-positioned flanking
nucleosomes, and paired-end ATAC fragment lengths show both. F-005 found fragment-size input
channels gain in-cell and are null transferred (library-specific), but here each cell's
signature is only used within that cell, standardised there.

STAGE features (one cell per call): for each ABC candidate region of the cell (4.43/4.47's union
set, restricted to that cell's candidates), paired-end fragments with MAPQ >= 30 whose midpoint
lies within +/- 1 kb of the region centre, from the cell's ATAC BAMs:
    n_total      all fragments, midpoint within +/- 450 bp
    nfr          fragments < 120 bp, midpoint within +/- 100 bp (the open centre)
    mono_flank   150-250 bp fragments (one nucleosome), midpoint 100-450 bp either side
    mono_centre  150-250 bp fragments, midpoint within +/- 100 bp (a nucleosome over the centre)
    di_flank     300-450 bp fragments, midpoint 150-600 bp either side
    nfr_width    distance between the strongest mono-nucleosome dyad left and right of centre
                 (10 bp bins, 3-bin smoothing; >= 5 mono fragments each side, else NaN)
    positioning  fraction of flank mono midpoints within +/- 40 bp of those two dyads
Derived: flank_ratio = log((mono_flank+1)/(nfr+1)), centre_occ = log((mono_centre+1)/(nfr+1)),
short_frac = nfr / n_total.

STAGE analyze (all five cells):
 1. Per cell, non-promoter candidates: partial Spearman of each feature with observed p300,
    controlling ATAC (rank, rank^2), GC, CpG o/e, log distance to TSS. A joint score (OLS of
    p300 residual on all features, weights fitted on the OTHER four cells, so deployable) gives
    the total nucleosome signal beyond ATAC; K562 reference points: the in-cell model's partial
    rho 0.44 and the transferred multi-cell model's 0.28-0.29 (F-029, F-031).
 2. Motif association with that joint nucleosome score (beyond ATAC and composition) as a TF-list
    estimator, alone and added to raw ATAC association, scored against 4.43's truth leave one
    cell out, as 4.48/4.50 did. These are proxies; a better score does not guarantee a benchmark
    gain (F-030, F-036).

Usage: 4.51.nucleosome_signature.py features --cell K562
       4.51.nucleosome_signature.py analyze [--out-prefix results/2026_0930_nucleosome_signature]
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

D = "/oak/stanford/groups/engreitz/Users/sheth"
P = f"{D}/EP300_BPNet"
A = f"{P}/2026_0824_H3K27ac_model"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction/results"
CANDS = {"K562": f"{ABC}/2026_0905_p300_activity/p300obs_k562",
         **{c: f"{ABC}/2026_0924_p300_panel_obs/p300obs_{c.replace('-', '_')}"
            for c in ("GM12878", "HepG2", "A549", "MCF-7")}}
CELLS = list(CANDS)
SCR = os.environ.get("SCRATCH", "/tmp")
BAMS = {"K562": sorted(glob.glob(f"{SCR}/atac_pe/ENCFF*.pe.bam")),
        "GM12878": [f"{D}/Data/ENCODE/GM12878/ATAC/{x}.bam"
                    for x in ("ENCFF440GRZ", "ENCFF962FMH", "ENCFF981FXV")],
        **{c: sorted(glob.glob(f"{D}/Data/ENCODE/{c}/ATAC/ENCFF*.bam"))
           for c in ("HepG2", "A549", "MCF-7")}}
GENOME = f"{D}/hg38_resources/hg38.fa"
TSS = (f"{D}/scE2G_temp/scE2G/resources/genome_annotations/"
       "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6")
REC = f"{A}/results/2026_0926_tf_activity_recovery.motifs.tsv"
CACHE = f"{A}/data/tf_activity_recovery"
OUT = f"{SCR}/nucsig"
KEY = ["chr", "start", "end"]
FEATS = ["n_total", "nfr", "mono_flank", "mono_centre", "di_flank", "nfr_width", "positioning",
         "flank_ratio", "centre_occ", "short_frac"]

ap = argparse.ArgumentParser()
ap.add_argument("stage", choices=("features", "analyze"))
ap.add_argument("--cell")
ap.add_argument("--out-prefix", default="results/2026_0930_nucleosome_signature")
a = ap.parse_args()


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


S = pd.read_csv(f"{CACHE}/union_signals.tsv.gz", sep="\t")
U = S[KEY].copy()
U["mid"] = (U["start"] + U["end"]) // 2

# ======================================================================== features
if a.stage == "features":
    import pysam

    c = a.cell
    os.makedirs(OUT, exist_ok=True)
    idx_dir = f"{SCR}/bam_index"
    os.makedirs(idx_dir, exist_ok=True)
    inc = overlaps(U, pd.read_csv(f"{CANDS[c]}/Neighborhoods/EnhancerList.txt", sep="\t",
                                  usecols=KEY))
    R = U[inc].reset_index()
    bams = []
    for b in BAMS[c]:
        bai = b + ".bai" if os.path.exists(b + ".bai") else f"{idx_dir}/{os.path.basename(b)}.bai"
        if not os.path.exists(bai):
            pysam.index(b, bai)
            print(f"indexed {b} -> {bai}", flush=True)
        bams.append(pysam.AlignmentFile(b, index_filename=bai))
    print(f"{c}: {len(R):,} candidate regions, {len(bams)} BAMs", flush=True)
    out = np.full((len(R), len(FEATS)), np.nan)
    bins = np.arange(-450, 451, 10)
    for i, (ch, mid) in enumerate(zip(R["chr"].values, R["mid"].values)):
        rel, fl = [], []
        for bam in bams:
            for r in bam.fetch(ch, max(0, mid - 1000), mid + 1000):
                if (r.is_secondary or r.is_supplementary or not r.is_proper_pair
                        or r.mapping_quality < 30 or r.template_length <= 0):
                    continue
                t = r.template_length
                rel.append(r.reference_start + t // 2 - mid)
                fl.append(t)
        rel, fl = np.asarray(rel), np.asarray(fl)
        ar = np.abs(rel)
        mono = (fl >= 150) & (fl <= 250)
        n_total = int((ar <= 450).sum())
        nfr = int(((fl < 120) & (ar <= 100)).sum())
        mfl = mono & (ar > 100) & (ar <= 450)
        mono_flank = int(mfl.sum())
        mono_centre = int((mono & (ar <= 100)).sum())
        di_flank = int(((fl >= 300) & (fl <= 450) & (ar > 150) & (ar <= 600)).sum())
        width = posn = np.nan
        L, Rt = rel[mfl & (rel < 0)], rel[mfl & (rel > 0)]
        if len(L) >= 5 and len(Rt) >= 5:
            h, _ = np.histogram(rel[mfl], bins)
            hs = np.convolve(h, np.ones(3), "same")
            ctr = (bins[:-1] + bins[1:]) / 2
            lft, rgt = ctr < 0, ctr > 0
            dl = ctr[lft][np.argmax(hs[lft])]
            dr = ctr[rgt][np.argmax(hs[rgt])]
            width = dr - dl
            posn = ((np.abs(L - dl) <= 40).sum() + (np.abs(Rt - dr) <= 40).sum()) / mono_flank
        out[i] = [n_total, nfr, mono_flank, mono_centre, di_flank, width, posn,
                  np.log((mono_flank + 1) / (nfr + 1)), np.log((mono_centre + 1) / (nfr + 1)),
                  nfr / n_total if n_total else np.nan]
        if i % 20000 == 0:
            print(f"  {i:,}/{len(R):,}", flush=True)
    df = pd.DataFrame(out, columns=FEATS)
    df.insert(0, "u", R["index"].values)
    df.to_csv(f"{OUT}/{c.replace('-', '_')}.tsv.gz", sep="\t", index=False)
    print(f"wrote {OUT}/{c.replace('-', '_')}.tsv.gz")
    sys.exit(0)

# ======================================================================== analyze
tss = pd.read_csv(TSS, sep="\t", header=None, usecols=[0, 1, 2], names=KEY)
U["prom"] = overlaps(U, tss)
tss["m"] = (tss["start"] + tss["end"]) // 2
dist = np.full(len(U), 1e7)
for ch, g in U.groupby("chr"):
    t = np.sort(tss.loc[tss["chr"] == ch, "m"].values)
    if len(t) > 1:
        k = np.clip(np.searchsorted(t, g["mid"].values), 1, len(t) - 1)
        dist[g.index.values] = np.minimum(np.abs(g["mid"].values - t[k - 1]),
                                          np.abs(g["mid"].values - t[k]))
U["ldist"] = np.log10(dist + 1)
import pyfaidx  # noqa: E402

fa = pyfaidx.Fasta(GENOME, sequence_always_upper=True)
gc, cpg = np.empty(len(U)), np.empty(len(U))
for i, (ch, mid) in enumerate(zip(U["chr"].values, U["mid"].values)):
    s = str(fa[ch][max(0, mid - 250):mid + 250])
    n = max(sum(s.count(b) for b in "ACGT"), 1)
    nc, ng = s.count("C"), s.count("G")
    gc[i], cpg[i] = (nc + ng) / n, s.count("CG") * n / max(nc * ng, 1)
U["gc"], U["cpg"] = gc, cpg
z = np.load(f"{CACHE}/union_fimo_presence.npz", allow_pickle=True)
M_all, names_all = z["presence"], list(z["names"])
Rm = pd.read_csv(REC, sep="\t")
mnames = list(Rm["motif"])
M = M_all[:, [names_all.index(n) for n in mnames]]
ctcf = np.array([n.split(" ", 1)[1].upper().startswith("CTCF") for n in mnames])
pct = lambda v: (rankdata(v) - 0.5) / len(v)


def resid(y, Z):
    Q, _ = np.linalg.qr(Z)
    return y - Q @ (Q.T @ y)


DATA = {}
print("\n=== 1. partial Spearman with observed p300, beyond ATAC and composition "
      "(non-promoter candidates) ===")
rows = []
for c in CELLS:
    f = pd.read_csv(f"{OUT}/{c.replace('-', '_')}.tsv.gz", sep="\t")
    g = U.loc[f["u"].values].reset_index(drop=True)
    keep = ~g["prom"].values & (f["n_total"].values >= 10)
    f, g = f[keep].reset_index(drop=True), g[keep].reset_index(drop=True)
    atac = pct(S.loc[f["u"].values, f"atac_{c}"].values)
    y = pct(S.loc[f["u"].values, f"p300_{c}"].values)
    Z = np.column_stack([np.ones(len(g)), atac, atac ** 2, g["gc"], g["cpg"], g["ldist"]])
    yr = resid(y, Z)
    X = {}
    for k in FEATS[1:]:
        v = f[k].values.astype(float)
        v = np.where(np.isfinite(v), v, np.nanmedian(v))
        X[k] = resid(pct(v), Z)
        rows.append(dict(cell=c, feature=k, partial_rho=np.corrcoef(X[k], yr)[0, 1]))
    DATA[c] = dict(u=f["u"].values, Z=Z, yr=yr, X=np.column_stack([X[k] for k in FEATS[1:]]),
                   n=len(f))
    print(f"  {c}: {len(f):,} regions", flush=True)
PR = pd.DataFrame(rows).pivot(index="feature", columns="cell", values="partial_rho")[CELLS]
print(PR.round(3).to_string())

print("\njoint nucleosome score (weights from the other four cells), partial rho with p300:")
SCORE = {}
for c in CELLS:
    oth = [d for d in CELLS if d != c]
    Xo = np.vstack([DATA[d]["X"] for d in oth])
    yo = np.concatenate([DATA[d]["yr"] for d in oth])
    w = np.linalg.lstsq(np.column_stack([np.ones(len(yo)), Xo]), yo, rcond=None)[0]
    sc = DATA[c]["X"] @ w[1:]
    SCORE[c] = sc
    rin = np.corrcoef(DATA[c]["X"] @ np.linalg.lstsq(DATA[c]["X"], DATA[c]["yr"], rcond=None)[0],
                      DATA[c]["yr"])[0, 1]
    print(f"  {c:<8} leave-one-cell-out {np.corrcoef(sc, DATA[c]['yr'])[0, 1]:+.3f}   "
          f"(fitted in-cell, optimistic: {rin:+.3f})")
print("  K562 references: in-cell sequence+ATAC model 0.44; transferred multi-cell 0.28-0.29")

print("\n=== 2. motifs associated with the nucleosome score, as a TF-list estimator ===")
EST = {k: {} for k in ("raw ATAC association", "nucleosome score association",
                       "ATAC + nucleosome (z sum)")}
T = {}
zs = lambda v: (v - np.nanmean(v)) / np.nanstd(v)
for c in CELLS:
    d = DATA[c]
    F = M[d["u"]].astype(np.float32)
    Fr = F - (lambda Q: Q @ (Q.T @ F))(np.linalg.qr(d["Z"])[0])
    Zc = d["Z"][:, [0, 3, 4, 5]]                  # composition only, for the ATAC association
    Fc = F - (lambda Q: Q @ (Q.T @ F))(np.linalg.qr(Zc)[0])
    atac = d["Z"][:, 1]
    ar = resid(atac, Zc)
    b_atac = (Fc * ar[:, None]).sum(0) / np.maximum((Fc * Fc).sum(0), 1e-12)
    sr = resid(SCORE[c], d["Z"])
    b_nuc = (Fr * sr[:, None]).sum(0) / np.maximum((Fr * Fr).sum(0), 1e-12)
    few = F.sum(0) < 200
    b_atac[few], b_nuc[few] = np.nan, np.nan
    EST["raw ATAC association"][c] = zs(b_atac)
    EST["nucleosome score association"][c] = zs(b_nuc)
    EST["ATAC + nucleosome (z sum)"][c] = zs(zs(b_atac) + zs(b_nuc))
    T[c] = zs(Rm[f"truth_{c}"].values)
top = lambda v, k=20: set(np.argsort(-np.nan_to_num(v, nan=-9))[:k])
rows = []
for name, E in EST.items():
    for c in CELLS:
        prior = np.mean([T[d] for d in CELLS if d != c], 0)
        e, t = E[c], T[c]
        ok = np.isfinite(e) & np.isfinite(t)
        rows.append(dict(estimate=name, cell=c, spearman=spearmanr(e[ok], t[ok])[0],
                         dep_spearman=spearmanr((e - prior)[ok], (t - prior)[ok])[0],
                         top20=len(top(e) & top(t)),
                         ctcf_best_rank=int(rankdata(-np.nan_to_num(e, nan=-9))[ctcf].min())))
Rt = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print(Rt.groupby("estimate", sort=False)[["spearman", "dep_spearman", "top20"]].mean()
      .round(3).to_string())
print(Rt.pivot(index="estimate", columns="cell", values="spearman").reindex(EST)[CELLS]
      .round(3).to_string())
gi = [i for i, n in enumerate(mnames) if "GATA" in n.upper() or "TRPS1" in n.upper()]
for name, E in EST.items():
    e = E["K562"]
    print(f"  K562 {name:<30} GATA {np.nanmean(e[gi]):+.2f}  top10: "
          + ", ".join(mnames[i].split(" ", 1)[1] for i in np.argsort(-np.nan_to_num(e, nan=-9))[:10]))
print(f"  K562 truth{'':<26}GATA {T['K562'][gi].mean():+.2f}")

os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
PR.to_csv(f"{a.out_prefix}.partial_rho.tsv", sep="\t")
Rt.to_csv(f"{a.out_prefix}.tf_list_scores.tsv", sep="\t", index=False)
wide = pd.DataFrame({"motif": mnames})
for name, E in EST.items():
    for c in CELLS:
        wide[f"{name.split()[0]}_{c}"] = E[c]
wide.to_csv(f"{a.out_prefix}.motifs.tsv", sep="\t", index=False)
print(f"\nwrote {a.out_prefix}.partial_rho.tsv, .tf_list_scores.tsv, .motifs.tsv")
