#!/usr/bin/env python3
"""Are nucleosome-spanning ATAC fragments slightly longer next to H3K27ac-high (or p300-high)
peaks than next to H3K27ac-low peaks?

HYPOTHESIS (Maya, 2026-09-30). Acetylated nucleosomes flanking an active element protect, or are
spaced by, a little more DNA, so the one- and two-nucleosome fragment peaks in the flanks sit
~10-15 bp longer than in the flanks of inactive accessible sites. fig35 (4.56) could not see
this: its window included the peak itself and it compared counts, not the exact peak position.

PEAKS. K562 ATAC MACS2 peaks from the July ABC run (one row per distinct peak, strongest summit).
Confident peaks: -log10 q >= 2. Kept: non-promoter (summit > 1 kb from a GENCODE v43 TSS),
ISOLATED (no other confident peak within --isolation bp of either peak bound).

CLASSIFICATION (percentiles among the kept peaks):
    ATAC      5' ATAC bigwig over summit +/- 250 (strong = top third, weak = bottom third)
    H3K27ac   5' ChIP bigwig over summit +/- 1 kb (low = bottom 40%, high = top 20%)
    p300      5' ChIP bigwig over summit +/- 250 (same cut-offs)
    CTCF      motif: FIMO MA0139.2 / MA1929.2 / MA1930.2 hit, -log10 p >= 5, within 100 bp of
              the summit; ChIP: ENCFF216XRV reads within 250 bp (CTCF-bound = motif + top 25%
              ChIP; CTCF-free = no motif and ChIP below median)
Groups, for MARK = H3K27ac and again for MARK = p300:
    CTCF, MARK-low            CTCF-bound, MARK low
    weak ATAC, MARK-low       CTCF-free
    strong ATAC, MARK-low     CTCF-free
    strong ATAC, MARK-high    CTCF-free   (ATAC-matched comparator for the line above)
    any ATAC, MARK-high       CTCF-free

FRAGMENTS. Paired-end, proper pair, MAPQ >= 30, not duplicate, from the three K562 ATAC BAMs.
A fragment counts only if it lies ENTIRELY inside a flank: [peak start - W, peak start - 50) or
[peak end + 50, peak end + W), W = --window. The peak itself (MACS2 bounds + 50 bp) is excluded.

DISTANCE. Each flank fragment is also binned by the distance of its inner end from the peak
bound (50-150, 150-300, 300-500, 500-750, 750-1000, 1000-1500, 1500-2000 bp), because the
acetylated nucleosomes are the ones next to the open site and a wide window dilutes them.

STATISTICS. Pooled 1-bp length histogram per group. Mono-nucleosome peak = mode of the
5-bp-Gaussian-smoothed histogram in 130-260 bp; di-nucleosome peak = mode in 280-460 bp;
repeat length = di - mono. Mean length of 140-260 bp fragments. 95% intervals from 300
bootstrap resamples of peaks. Flank GC reported per group (fragment length depends on GC).

Usage: 4.59.flank_fragment_lengths.py [--window 1000] [--isolation 2000] [--n 4000]
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyBigWig  # noqa: E402
import pyfaidx  # noqa: E402
import pysam  # noqa: E402
from scipy.ndimage import gaussian_filter1d  # noqa: E402

D = "/oak/stanford/groups/engreitz/Users/sheth"
P = f"{D}/EP300_BPNet"
A = f"{P}/2026_0824_H3K27ac_model"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction"
SCR = os.environ.get("SCRATCH", "/tmp")
PEAKS = f"{ABC}/results/2026_0721_h3k27ac_counting_comparison/K562_ATAC_only/Peaks/macs2_peaks.narrowPeak"
BAMS = sorted(glob.glob(f"{SCR}/atac_pe/ENCFF*.pe.bam"))
CTCF_BAM = f"{D}/Data/ENCODE/K562/ENCFF216XRV.filtered.sorted.bam"
ATAC_BW = [f"{P}/2026_0529_multimodal_p300_model/data/atac.bw"]
K27_BW = [f"{A}/data/h3k27ac_5p_plus.bw", f"{A}/data/h3k27ac_5p_minus.bw"]
P300_BW = [f"{P}/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig",
           f"{P}/2025_0703_retrain_p300_model/data/ENCSR000EGE_minus.bigWig"]
GENOME = f"{D}/hg38_resources/hg38.fa"
TSS = (f"{D}/scE2G_temp/scE2G/resources/genome_annotations/"
       "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6")
CTCF_MOTIFS = ("MA0139.2 CTCF", "MA1929.2 CTCF", "MA1930.2 CTCF")
L = np.arange(1, 1001)
DBINS = np.array([50, 150, 300, 500, 750, 1000, 1500, 2000])
MW = (L >= 140) & (L <= 260)

ap = argparse.ArgumentParser()
ap.add_argument("--window", type=int, default=1000)
ap.add_argument("--isolation", type=int, default=2000)
ap.add_argument("--n", type=int, default=4000, help="max peaks per group")
ap.add_argument("--out-prefix", default=f"{A}/results/2026_0930_flank_fragment_lengths")
ap.add_argument("--outdir", default=f"{A}/figures")
a = ap.parse_args()
rng = np.random.default_rng(0)

# ---------------------------------------------------------------- peaks
pk = pd.read_csv(PEAKS, sep="\t", header=None,
                 names=["chr", "start", "end", "name", "score", "strand", "fe", "p", "q", "summit"])
pk = pk[pk["chr"].str.match(r"^chr[0-9XY]+$")]
pk = pk.sort_values("p", ascending=False).drop_duplicates(["chr", "start", "end"])
conf = pk[pk["q"] >= 2].sort_values(["chr", "start"]).reset_index(drop=True)
conf["smt"] = conf["start"] + conf["summit"]
gap_l = conf["start"] - conf.groupby("chr")["end"].shift(1)
gap_r = conf.groupby("chr")["start"].shift(-1) - conf["end"]
iso = (gap_l.fillna(1e9) >= a.isolation) & (gap_r.fillna(1e9) >= a.isolation)
tss = pd.read_csv(TSS, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "s", "e"])
tss["m"] = (tss["s"] + tss["e"]) // 2
far = np.ones(len(conf), bool)
for ch, g in conf.groupby("chr"):
    t = np.sort(tss.loc[tss["chr"] == ch, "m"].values)
    if len(t):
        k = np.clip(np.searchsorted(t, g["smt"].values), 1, len(t) - 1)
        d = np.minimum(np.abs(g["smt"].values - t[k - 1]), np.abs(g["smt"].values - t[k]))
        far[g.index.values] = d > 1000
K = conf[iso.values & far].reset_index(drop=True)
print(f"{len(pk):,} distinct MACS2 peaks; {len(conf):,} with q < 0.01; {int(iso.sum()):,} isolated "
      f"(>= {a.isolation} bp from another); {len(K):,} isolated and non-promoter", flush=True)


def bw_sum(paths, half):
    out = np.zeros(len(K))
    for p in paths:
        bw = pyBigWig.open(p)
        Ls = bw.chroms()
        for i, (c, m) in enumerate(zip(K["chr"].values, K["smt"].values)):
            if c in Ls:
                v = bw.values(c, max(0, m - half), min(m + half, Ls[c]), numpy=True)
                out[i] += np.nansum(np.abs(v))
        bw.close()
    return out


pct = lambda v: pd.Series(v).rank(pct=True).values
K["atac"] = pct(bw_sum(ATAC_BW, 250))
K["k27"] = pct(bw_sum(K27_BW, 1000))
K["p300"] = pct(bw_sum(P300_BW, 250))
cb = pysam.AlignmentFile(CTCF_BAM)
K["ctcf_chip"] = pct(np.array([cb.count(c, max(0, m - 250), m + 250)
                               for c, m in zip(K["chr"].values, K["smt"].values)]))
mot = np.zeros(len(K), bool)
for ch, g in K.groupby("chr"):
    f = f"{SCR}/motif_hits/{ch}.npz"
    if not os.path.exists(f):
        continue
    z = np.load(f, allow_pickle=True)
    ids = [i for i, n in enumerate(z["names"]) if n in CTCF_MOTIFS]
    sel = np.isin(z["motif"], ids) & (z["neglog10p"].astype(float) >= 5)
    ctr = np.sort(((z["start"][sel] + z["end"][sel]) // 2).astype(np.int64))
    if len(ctr) == 0:
        continue
    s = g["smt"].values
    k = np.clip(np.searchsorted(ctr, s), 1, len(ctr) - 1)
    d = np.minimum(np.abs(s - ctr[k - 1]), np.abs(s - ctr[k]))
    mot[g.index.values] = d <= 100
K["ctcf_motif"] = mot
bound = K["ctcf_motif"] & (K["ctcf_chip"] >= 0.75)
free = ~K["ctcf_motif"] & (K["ctcf_chip"] < 0.5)
print(f"CTCF motif at {mot.mean():.1%} of kept peaks; CTCF-bound {bound.mean():.1%}, "
      f"CTCF-free {free.mean():.1%}", flush=True)

GROUPS = {}
for mark in ("k27", "p300"):
    lo, hi = K[mark] <= 0.4, K[mark] >= 0.8
    strong, weak = K["atac"] >= 2 / 3, K["atac"] <= 1 / 3
    GROUPS[mark] = {"CTCF-bound, mark low": bound & lo,
                    "weak ATAC, mark low": free & weak & lo,
                    "strong ATAC, mark low": free & strong & lo,
                    "strong ATAC, mark high": free & strong & hi,
                    "any ATAC, mark high": free & hi}

# ---------------------------------------------------------------- fragments
bams = [pysam.AlignmentFile(b) for b in BAMS]
fa = pyfaidx.Fasta(GENOME, sequence_always_upper=True)
W, PAD = a.window, 50


def peak_hist(i):
    """per-BAM length histograms: fragments entirely in the flanks, and every fragment with
    its midpoint inside [peak start - W, peak end + W] (peak included, the standard view)"""
    c, s, e = K.at[i, "chr"], int(K.at[i, "start"]), int(K.at[i, "end"])
    hf = np.zeros((len(bams), len(L)))
    ha = np.zeros((len(bams), len(L)))
    hd = np.zeros((len(DBINS) - 1, len(L)))
    for b, bam in enumerate(bams):
        for r in bam.fetch(c, max(0, s - W), e + W):
            if (r.is_secondary or r.is_supplementary or r.is_duplicate or not r.is_proper_pair
                    or r.mapping_quality < 30 or r.template_length <= 0):
                continue
            t = r.template_length
            if t > L[-1]:
                continue
            f0, f1 = r.reference_start, r.reference_start + t
            if s - W <= f0 + t // 2 < e + W:
                ha[b, t - 1] += 1
            if (f0 >= s - W and f1 <= s - PAD) or (f0 >= e + PAD and f1 <= e + W):
                hf[b, t - 1] += 1
                d = s - f1 if f1 <= s - PAD else f0 - e
                k = np.searchsorted(DBINS, d, side="right") - 1
                if 0 <= k < len(DBINS) - 1:
                    hd[k, t - 1] += 1
    seq = str(fa[c][max(0, s - W):s - PAD]) + str(fa[c][e + PAD:e + W])
    gc = (seq.count("G") + seq.count("C")) / max(sum(seq.count(b) for b in "ACGT"), 1)
    # per distance bin: pooled histogram, and (n, sum of lengths) of 140-260 bp fragments
    return hf, ha, gc, hd, np.stack([hd[:, MW].sum(1), (hd[:, MW] * L[MW]).sum(1)], 1)


def modes(h):
    sm = gaussian_filter1d(h.astype(float), 5)
    m1 = L[(L >= 130) & (L <= 260)][np.argmax(sm[(L >= 130) & (L <= 260)])]
    m2 = L[(L >= 280) & (L <= 460)][np.argmax(sm[(L >= 280) & (L <= 460)])]
    w = (L >= 140) & (L <= 260)
    return m1, m2, m2 - m1, (h[w] * L[w]).sum() / max(h[w].sum(), 1)


cache, rows, H, HB, DR, HDG = {}, [], {}, {}, [], {}
for mark, gs in GROUPS.items():
    for gname, sel in gs.items():
        idx = np.flatnonzero(sel.values)
        idx = rng.choice(idx, min(a.n, len(idx)), replace=False)
        per = []
        for i in idx:
            if i not in cache:
                cache[i] = peak_hist(i)
            per.append(cache[i])
        hs = np.array([p[0].sum(0) for p in per])
        gcs = np.array([p[2] for p in per])
        HB[(mark, gname)] = (np.sum([p[0] for p in per], 0), np.sum([p[1] for p in per], 0))
        HDG[(mark, gname)] = np.sum([p[3] for p in per], 0)
        ns = np.array([p[4] for p in per])            # peaks x bins x (n, sum)
        for kb in range(len(DBINS) - 1):
            if DBINS[kb + 1] > W:
                break
            n_, s_ = ns[:, kb, 0], ns[:, kb, 1]
            bt = []
            for _ in range(300):
                j_ = rng.integers(0, len(n_), len(n_))
                bt.append(s_[j_].sum() / max(n_[j_].sum(), 1))
            DR.append(dict(mark=mark, group=gname, d0=DBINS[kb], d1=DBINS[kb + 1],
                           n_mono=int(n_.sum()), mono_mean=s_.sum() / max(n_.sum(), 1),
                           lo=np.percentile(bt, 2.5), hi=np.percentile(bt, 97.5)))
        tot = hs.sum(0)
        m1, m2, nrl, mean_mono = modes(tot)
        bs = np.array([modes(hs[rng.integers(0, len(hs), len(hs))].sum(0)) for _ in range(300)])
        lo_, hi_ = np.percentile(bs, 2.5, 0), np.percentile(bs, 97.5, 0)
        H[(mark, gname)] = tot
        rows.append(dict(mark=mark, group=gname, n_peaks=len(idx), n_fragments=int(tot.sum()),
                         mean_atac_pct=K["atac"].values[idx].mean(),
                         mean_mark_pct=K[mark].values[idx].mean(), flank_gc=gcs.mean(),
                         mono_mode=m1, mono_lo=lo_[0], mono_hi=hi_[0],
                         di_mode=m2, di_lo=lo_[1], di_hi=hi_[1],
                         repeat=nrl, repeat_lo=lo_[2], repeat_hi=hi_[2],
                         mono_mean=mean_mono, mono_mean_lo=lo_[3], mono_mean_hi=hi_[3]))
        r = rows[-1]
        print(f"{mark:<5} {gname:<24} peaks {r['n_peaks']:>5,} frags {r['n_fragments']:>9,} "
              f"ATAC {r['mean_atac_pct']:.2f} GC {r['flank_gc']:.3f} | mono {m1} [{lo_[0]:.0f},{hi_[0]:.0f}] "
              f"di {m2} [{lo_[1]:.0f},{hi_[1]:.0f}] repeat {nrl} | mean 140-260 {mean_mono:.1f} "
              f"[{lo_[3]:.1f},{hi_[3]:.1f}]", flush=True)
R = pd.DataFrame(rows)
R.to_csv(f"{a.out_prefix}.summary.tsv", sep="\t", index=False)
pd.DataFrame({f"{m}|{g}": h for (m, g), h in H.items()}, index=L).to_csv(
    f"{a.out_prefix}.hist.tsv.gz", sep="\t")
REPS = [os.path.basename(b).split(".")[0] for b in BAMS]
pd.DataFrame({f"{m}|{g}|{kind}|{REPS[b]}": HB[(m, g)][k][b] for (m, g) in HB
              for k, kind in enumerate(("flank", "all")) for b in range(len(REPS))},
             index=L).to_csv(f"{a.out_prefix}.hist_by_library.tsv.gz", sep="\t")

# ---------------------------------------------------------------- figure
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 160, "savefig.dpi": 160, "savefig.bbox": "tight",
})
COLS = {"CTCF-bound, mark low": "#1D9E75", "weak ATAC, mark low": "#B4B2A9",
        "strong ATAC, mark low": "#5F5E5A", "strong ATAC, mark high": "#D85A30",
        "any ATAC, mark high": "#F0997B"}
LAB = {"k27": "H3K27ac", "p300": "p300"}
fig, ax = plt.subplots(2, 4, figsize=(13, 6), gridspec_kw={"width_ratios": [1.3, 1, 1, 1.1]})
for i, mark in enumerate(GROUPS):
    for gname in GROUPS[mark]:
        h = H[(mark, gname)]
        lab = gname.replace("mark", LAB[mark])
        dens = gaussian_filter1d(h / h.sum(), 2)
        ax[i, 0].plot(L, dens, color=COLS[gname], lw=1.2, label=lab)
        for j, (x0, x1) in ((1, (120, 280)), (2, (270, 470))):
            w = (L >= x0) & (L <= x1)
            ax[i, j].plot(L[w], gaussian_filter1d(h[w] / h[w].sum(), 2), color=COLS[gname], lw=1.2)
    ax[i, 0].legend(frameon=False, fontsize=7)
    ax[i, 0].set_ylabel(f"groups by {LAB[mark]}\nfraction of flank fragments")
    r = R[R["mark"] == mark].reset_index(drop=True)
    y = np.arange(len(r))[::-1]
    for k_, (col, lo_, hi_, mk) in enumerate((("mono_mode", "mono_lo", "mono_hi", "o"),
                                              ("mono_mean", "mono_mean_lo", "mono_mean_hi", "s"))):
        ax[i, 3].errorbar(r[col], y + (0.15 if k_ == 0 else -0.15),
                          xerr=[r[col] - r[lo_], r[hi_] - r[col]], fmt=mk, ms=4, lw=1,
                          color="0.2" if k_ == 0 else "0.6",
                          label="mode" if k_ == 0 else "mean of 140-260 bp")
    ax[i, 3].set_yticks(y)
    ax[i, 3].set_yticklabels([g.replace("mark", LAB[mark]) for g in r["group"]], fontsize=7)
    ax[i, 3].legend(frameon=False, fontsize=7, loc="lower right")
ax[0, 0].set_title("a  flank fragment lengths")
ax[0, 1].set_title("b  one nucleosome (normalised)")
ax[0, 2].set_title("c  two nucleosomes (normalised)")
ax[0, 3].set_title("d  one-nucleosome length, 95% CI")
for j, lab in enumerate(("fragment length (bp)",) * 3 + ("bp",)):
    ax[1, j].set_xlabel(lab)
fig.suptitle(f"K562 isolated ATAC peaks: fragments lying entirely in the flanks "
             f"(peak bounds + 50 bp to {W} bp), peak itself excluded", fontsize=10)
fig.tight_layout()
name = f"fig36_flank_fragment_lengths_w{W}"
for ext in ("png", "pdf"):
    fig.savefig(f"{a.outdir}/{name}.{ext}")
print(f"wrote {a.outdir}/{name}.png and {a.out_prefix}.summary.tsv")


# ---------------------------------------------------------------- per-library, report 1 style
RCOL = ["#0096a0", "#792374", "#ee6a00"]
SUB, MONO = (1, 99), (180, 247)
edges = np.arange(0, 805, 5)
for k, kind, desc in ((0, "flank", f"fragments lying entirely in the flanks (peak + 50 bp to {W} bp)"),
                      (1, "all", f"all fragments with midpoint within the peak +/- {W} bp")):
    fig, ax = plt.subplots(2, 5, figsize=(16, 5.6), sharex=True, sharey=True)
    for i, mark in enumerate(GROUPS):
        for j, gname in enumerate(GROUPS[mark]):
            h3 = HB[(mark, gname)][k]
            for b in range(len(REPS)):
                c = np.add.reduceat(h3[b][:800], np.arange(0, 800, 5))
                ax[i, j].plot(edges[:-1] + 2.5, c / c.sum() / 5, color=RCOL[b], lw=1.1,
                              label=REPS[b])
            ax[i, j].axvspan(*SUB, color="#0096a0", alpha=0.08)
            ax[i, j].axvspan(*MONO, color="#792374", alpha=0.08)
            n = R[(R["mark"] == mark) & (R["group"] == gname)].iloc[0]
            ax[i, j].set_title(f"{gname.replace('mark', LAB[mark])}\n{int(n['n_peaks']):,} peaks",
                               fontsize=8)
        ax[i, 0].set_ylabel(f"groups by {LAB[mark]}\ndensity")
    ax[0, 4].legend(frameon=False, fontsize=7)
    for j in range(5):
        ax[1, j].set_xlabel("fragment length (bp)")
    fig.suptitle(f"K562 isolated ATAC peaks, three ENCSR868FGK libraries: {desc}; "
                 "shaded sub (1-99 bp) and mono (180-247 bp) channels as in report 1", fontsize=10)
    fig.tight_layout()
    name = f"fig37_fragment_lengths_by_library_{kind}_w{W}"
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.outdir}/{name}.{ext}")
    plt.close(fig)
    print(f"wrote {a.outdir}/{name}.png")


# ---------------------------------------------------------------- length vs distance from the peak
DR = pd.DataFrame(DR)
DR.to_csv(f"{a.out_prefix}.by_distance.tsv", sep="\t", index=False)
print(DR.pivot_table(index=["mark", "group"], columns="d0", values="mono_mean").round(1).to_string())
fig, ax = plt.subplots(2, 2, figsize=(10, 6), gridspec_kw={"width_ratios": [1.3, 1]})
for i, mark in enumerate(GROUPS):
    for gname in GROUPS[mark]:
        r = DR[(DR["mark"] == mark) & (DR["group"] == gname)]
        x = (r["d0"] + r["d1"]) / 2
        ax[i, 0].plot(x, r["mono_mean"], "o-", ms=3, color=COLS[gname],
                      label=gname.replace("mark", LAB[mark]))
        ax[i, 0].fill_between(x, r["lo"], r["hi"], color=COLS[gname], alpha=0.2, lw=0)
        h = HDG[(mark, gname)][0]
        w = (L >= 120) & (L <= 280)
        ax[i, 1].plot(L[w], gaussian_filter1d(h[w] / h[w].sum(), 2), color=COLS[gname], lw=1.2)
    ax[i, 0].set_ylabel(f"groups by {LAB[mark]}\nmean length of 140-260 bp fragments (bp)")
    ax[i, 0].legend(frameon=False, fontsize=7)
    ax[i, 1].set_ylabel("fraction (normalised in window)")
ax[0, 0].set_title("a  one-nucleosome fragment length vs distance from the peak, 95% CI")
ax[0, 1].set_title(f"b  nearest flank only ({DBINS[0]}-{DBINS[1]} bp from the peak)")
ax[1, 0].set_xlabel("distance of the fragment's inner end from the peak bound (bp)")
ax[1, 1].set_xlabel("fragment length (bp)")
fig.tight_layout()
name = f"fig38_flank_length_by_distance_w{W}"
for ext in ("png", "pdf"):
    fig.savefig(f"{a.outdir}/{name}.{ext}")
print(f"wrote {a.outdir}/{name}.png")
