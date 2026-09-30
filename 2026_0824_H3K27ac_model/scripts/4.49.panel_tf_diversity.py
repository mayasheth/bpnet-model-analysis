#!/usr/bin/env python3
"""Which cell types to add to the p300 training panel: pick for coverage of active TF families.

WHY. The transferred model learns TF effects from the DNA of its training cell types and cannot
tell "this motif matters" from "this TF is active here" (F-035, F-037). A TF family active in
none of the training cell types is one the model has never seen switch p300 on, which is what
happened with GATA in K562 (F-032). So the panel should contain, for as many TF families as
possible, at least one cell type where the family is active, and others where it is not.

ACTIVITY WITHOUT p300. For each candidate (every ENCODE biosample with released human EP300
ChIP-seq and ATAC or DNase), the estimator that kept the cell-specific signal best in F-034:
within-cell association of accessibility with motif presence at non-promoter peaks, partialling
out GC, CpG o/e and log distance to TSS (Frisch-Waugh-Lovell), standardised within the cell.
Accessibility = within-cell percentile of the ENCODE peak signalValue; regions = summit +/- 250
bp, overlapping regions resolved to the stronger. Motif presence from the genome-wide FIMO hits
(4.44). CTCF-family motifs excluded.

FAMILIES. JASPAR has near-duplicate matrices (eight GATA), which would count one TF eight times.
Motifs are clustered by co-occurrence of their hits across all candidates' regions (average
linkage on 1 - Pearson r of presence, cut at r = 0.5); a family's score is its members' mean.

DEPARTURE. Family score minus its mean across all candidates: what is specific to the cell.

SPLIT-HALF RELIABILITY (added after the first run, 2026-09-29, which ranked noisy zinc-finger
and GC-rich families above K562's GATA/TAL1). Everything is computed three times: all
chromosomes, odd chromosomes only, even only. A family counts as active in a cell only as far as
BOTH halves agree: reliable = max(0, min(departure_odd, departure_even)). Noise does not
replicate across independent halves of the genome; a TF programme does. Candidates with fewer
than --min-regions regions are dropped (a few tissues had 6-22k IDR peaks).

SELECTION. Facility location, greedy (submodular, so greedy is within 1 - 1/e of the best):
coverage(S) = sum over families of max over c in S of max(0, departure[f, c]). Run (a) from the
current training panel, adding cell types one at a time, and (b) from scratch. The held-out
cell type is never a candidate; its own top families are used to check what each panel covers.

Usage: 4.49.panel_tf_diversity.py --manifest config/panel_candidates.tsv --peak-dir DIR
       [--holdout K562] [--base GM12878,HepG2,A549,MCF-7] [--k 6]
       [--out-prefix results/2026_0929_panel_tf_diversity]
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata, spearmanr

D = "/oak/stanford/groups/engreitz/Users/sheth"
A = f"{D}/EP300_BPNet/2026_0824_H3K27ac_model"
GENOME = f"{D}/hg38_resources/hg38.fa"
TSS = (f"{D}/scE2G_temp/scE2G/resources/genome_annotations/"
       "CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6")
REC = f"{A}/results/2026_0926_tf_activity_recovery.motifs.tsv"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
HALF = 250

ap = argparse.ArgumentParser()
ap.add_argument("--manifest", required=True)
ap.add_argument("--peak-dir", required=True)
ap.add_argument("--hits-dir", default=os.path.join(os.environ.get("SCRATCH", "/tmp"), "motif_hits"))
ap.add_argument("--holdout", default="K562")
ap.add_argument("--base", default="GM12878,HepG2,A549,MCF-7")
ap.add_argument("--k", type=int, default=6)
ap.add_argument("--max-regions", type=int, default=300000)
ap.add_argument("--min-regions", type=int, default=20000)
ap.add_argument("--min-hits", type=int, default=100)
ap.add_argument("--family-r", type=float, default=0.5)
ap.add_argument("--out-prefix", default="results/2026_0929_panel_tf_diversity")
a = ap.parse_args()

man = pd.read_csv(a.manifest, sep="\t")
CELLS = list(man["cell"])
base = [c for c in a.base.split(",") if c]
assert a.holdout in CELLS and all(c in CELLS for c in base)


# ---------------------------------------------------------------- regions per cell
def load_regions(path):
    p = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2, 6, 9],
                    names=["chr", "start", "end", "signal", "summit"])
    p = p[p["chr"].isin(CHROMS)]
    ctr = np.where(p["summit"] >= 0, p["start"] + p["summit"], (p["start"] + p["end"]) // 2)
    r = pd.DataFrame({"chr": p["chr"].values, "c": ctr.astype(np.int64),
                      "signal": p["signal"].values.astype(float)})
    r = r.sort_values("signal", ascending=False).drop_duplicates(["chr", "c"])
    keep = []
    for ch, g in r.groupby("chr", sort=False):          # resolve overlaps to the stronger
        g = g.sort_values("c")
        cs, sg = g["c"].values, g["signal"].values
        ok = np.ones(len(g), bool)
        last = -1
        for i in range(len(g)):
            if last >= 0 and cs[i] - cs[last] < 2 * HALF:
                if sg[i] > sg[last]:
                    ok[last] = False
                    last = i
                else:
                    ok[i] = False
            else:
                last = i
        keep.append(g[ok])
    r = pd.concat(keep).sort_values("signal", ascending=False).head(a.max_regions)
    return r.sort_values(["chr", "c"]).reset_index(drop=True)


REG = {}
for _, row in man.iterrows():
    f = glob.glob(f"{a.peak_dir}/{row['accession']}*.bed.gz")
    assert f, f"missing peaks for {row['cell']}: {row['accession']}"
    REG[row["cell"]] = load_regions(f[0])
    print(f"{row['cell']:<30} {row['assay']:<10} {len(REG[row['cell']]):>7,} regions", flush=True)
small = [c for c in CELLS if len(REG[c]) < a.min_regions]
assert a.holdout not in small and not set(base) & set(small), f"too few regions: {small}"
if small:
    print(f"dropped (fewer than {a.min_regions:,} regions): {', '.join(small)}")
CELLS = [c for c in CELLS if c not in small]
REG = {c: REG[c] for c in CELLS}

# ---------------------------------------------------------------- TSS, GC, CpG
tss = pd.read_csv(TSS, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "s", "e"])
tss["m"] = (tss["s"] + tss["e"]) // 2
import pyfaidx  # noqa: E402

fa = pyfaidx.Fasta(GENOME, sequence_always_upper=True)
names0 = None
PRES = {c: None for c in CELLS}
COV = {c: {k: np.empty(len(REG[c])) for k in ("gc", "cpg", "dist", "prom")} for c in CELLS}
IDX = {c: {ch: np.nonzero(REG[c]["chr"].values == ch)[0] for ch in CHROMS} for c in CELLS}
for ch in CHROMS:
    seq = np.frombuffer(str(fa[ch][:]).encode(), np.uint8)
    isgc = np.isin(seq, np.frombuffer(b"GC", np.uint8))
    isc, isg = seq == ord("C"), seq == ord("G")
    isacgt = np.isin(seq, np.frombuffer(b"ACGT", np.uint8))
    cg = np.zeros(len(seq), bool)
    cg[:-1] = isc[:-1] & isg[1:]
    cum = {k: np.concatenate([[0], np.cumsum(v, dtype=np.int32)]) for k, v in
           (("gc", isgc), ("c", isc), ("g", isg), ("n", isacgt), ("cg", cg))}
    t = np.sort(tss.loc[tss["chr"] == ch, "m"].values)
    z = np.load(f"{a.hits_dir}/{ch}.npz", allow_pickle=True)
    if names0 is None:
        names0 = list(z["names"])
        for c in CELLS:
            PRES[c] = np.zeros((len(REG[c]), len(names0)), np.uint8)
    hm = ((z["start"].astype(np.int64) + z["end"].astype(np.int64)) // 2)
    hmot = z["motif"].astype(np.int64)
    for c in CELLS:
        ii = IDX[c][ch]
        if not len(ii):
            continue
        ctr = REG[c]["c"].values[ii]
        s, e = np.clip(ctr - HALF, 0, len(seq)), np.clip(ctr + HALF, 0, len(seq))
        n = np.maximum(cum["n"][e] - cum["n"][s], 1)
        nc, ng = cum["c"][e] - cum["c"][s], cum["g"][e] - cum["g"][s]
        COV[c]["gc"][ii] = (cum["gc"][e] - cum["gc"][s]) / n
        COV[c]["cpg"][ii] = (cum["cg"][e] - cum["cg"][s]) * n / np.maximum(nc * ng, 1)
        if len(t) > 1:
            k = np.clip(np.searchsorted(t, ctr), 1, len(t) - 1)
            dd = np.minimum(np.abs(ctr - t[k - 1]), np.abs(ctr - t[k]))
        else:
            dd = np.full(len(ctr), 1e7)
        COV[c]["dist"][ii] = np.log10(dd + 1)
        COV[c]["prom"][ii] = dd <= 500
        # regions are non-overlapping and sorted, so each hit midpoint falls in at most one
        k = np.searchsorted(s, hm, side="right") - 1
        ok = (k >= 0) & (hm < e[np.clip(k, 0, None)])
        PRES[c][ii[k[ok]], hmot[ok]] = 1
    print(f"  {ch}: motif presence and composition done", flush=True)

# ---------------------------------------------------------------- association per cell
ctcf = np.array([n.split(" ", 1)[1].upper().startswith("CTCF") for n in names0])
ODD = {f"chr{i}" for i in range(1, 23, 2)}


def assoc(c, rows):
    sel = rows & ~COV[c]["prom"].astype(bool)
    y = rankdata(REG[c]["signal"].values[sel]) / sel.sum()
    Q, _ = np.linalg.qr(np.column_stack([np.ones(sel.sum()), COV[c]["gc"][sel],
                                         COV[c]["cpg"][sel], COV[c]["dist"][sel]]))
    F = PRES[c][sel].astype(np.float32)
    nh = F.sum(0)
    yr = y - Q @ (Q.T @ y)
    Fr = F - Q @ (Q.T @ F)
    b = (Fr * yr[:, None]).sum(0) / np.maximum((Fr * Fr).sum(0), 1e-12)
    b[(nh < a.min_hits) | ctcf] = np.nan
    return (b - np.nanmean(b)) / np.nanstd(b), int(sel.sum())


Z, Z1, Z2 = (pd.DataFrame(index=names0, columns=CELLS, dtype=float) for _ in range(3))
for c in CELLS:
    odd = REG[c]["chr"].isin(ODD).values
    Z[c], n = assoc(c, np.ones(len(odd), bool))
    Z1[c], _ = assoc(c, odd)
    Z2[c], _ = assoc(c, ~odd)
    print(f"association: {c} ({n:,} non-promoter regions, "
          f"{int(Z[c].notna().sum())} motifs scored)", flush=True)

# sanity: against 4.43's within-cell association on ABC candidate regions, where it exists
if os.path.exists(REC):
    R = pd.read_csv(REC, sep="\t").set_index("motif")
    for c in CELLS:
        col = f"within_{c}"
        if col in R:
            j = Z[c].dropna().index.intersection(R.index)
            print(f"  check {c}: Spearman with 4.43 within-cell association "
                  f"{spearmanr(Z.loc[j, c], R.loc[j, col])[0]:.3f}")

# ---------------------------------------------------------------- motif families
rng = np.random.default_rng(0)
samp = np.vstack([PRES[c][rng.choice(len(PRES[c]), min(20000, len(PRES[c])), replace=False)]
                  for c in CELLS]).astype(np.float32)
use = np.nonzero(~ctcf & (samp.sum(0) >= 50))[0]
r = np.corrcoef(samp[:, use].T)
r = np.nan_to_num(r)
dist = np.clip(1 - r, 0, 2)
np.fill_diagonal(dist, 0)
lab = fcluster(linkage(squareform(dist, checks=False), "average"), 1 - a.family_r, "distance")
fam = pd.Series(lab, index=np.array(names0)[use])
short = lambda n: n.split(" ", 1)[1]
famname = {}
for f_, members in fam.groupby(fam):
    ms = [short(m) for m in members.index]
    famname[f_] = ms[0] + (f" (+{len(ms) - 1})" if len(ms) > 1 else "")


def family_departure(Zm):
    fz = Zm.loc[fam.index].groupby(fam.values).mean()
    fz.index = [famname[i] for i in fz.index]
    fz = fz.fillna(0)
    return fz.sub(fz.mean(axis=1), axis=0)


DEP, DEP1, DEP2 = family_departure(Z), family_departure(Z1), family_departure(Z2)
FZ = DEP
REL = np.minimum(DEP1, DEP2).clip(lower=0)
print("split-half agreement of family departures (Spearman, odd vs even chromosomes):")
print("  " + ", ".join(f"{c} {spearmanr(DEP1[c], DEP2[c])[0]:.2f}" for c in CELLS))
print(f"\n{len(use)} motifs in {len(FZ)} families (co-occurrence r >= {a.family_r})")
gata = [i for i in FZ.index if "GATA" in i.upper()]
print("GATA families, departure:", "; ".join(
    f"{g}: " + ", ".join(f"{c} {DEP.loc[g, c]:+.2f}" for c in CELLS if abs(DEP.loc[g, c]) > 1)
    for g in gata))

print("\n=== top 5 RELIABLE cell-specific families per candidate (min of the two halves) ===")
for c in CELLS:
    s = REL[c].sort_values(ascending=False).head(5)
    print(f"  {c:<30} " + ", ".join(f"{i} {v:+.1f}" for i, v in s.items()))
rk = REL[a.holdout].rank(ascending=False)
print(f"\nVALIDATION: {a.holdout} rank of GATA families by reliable departure: "
      + ", ".join(f"{g} {int(rk[g])}" for g in gata) + f" (of {len(rk)})")

# ---------------------------------------------------------------- selection
P = REL
pool = [c for c in CELLS if c != a.holdout]


def cover(S):
    return float(P[S].max(axis=1).sum()) if S else 0.0


def greedy(start, k):
    S, out = list(start), []
    for _ in range(k):
        cand = [c for c in pool if c not in S]
        if not cand:
            break
        g = {c: cover(S + [c]) - cover(S) for c in cand}
        best = max(g, key=g.get)
        gain_f = (P[S + [best]].max(axis=1) - (P[S].max(axis=1) if S else 0)).sort_values(
            ascending=False)
        out.append((best, g[best], ", ".join(f"{i} +{v:.1f}" for i, v in gain_f.head(4).items())))
        S.append(best)
    return S, out


tops = REL[a.holdout].sort_values(ascending=False)
tops = tops[tops > 0.5].head(10)


def holdout_cover(S):
    """Fraction of the held-out cell's top departures reached by some panel cell (>= half)."""
    got = [(i, P.loc[i, S].max() if S else 0.0) for i in tops.index]
    return np.mean([g >= 0.5 * tops[i] for i, g in got]), got


rows = []
for label, start, k in (("extend current panel", base, a.k), ("from scratch", [], a.k + len(base))):
    S, out = greedy(start, k)
    print(f"\n=== greedy, {label} (start: {', '.join(start) or 'none'}; "
          f"coverage {cover(list(start)):.1f}) ===")
    cur = list(start)
    for i, (c, g, why) in enumerate(out, 1):
        cur.append(c)
        hc, _ = holdout_cover(cur)
        print(f"  {i}. {c:<30} +{g:5.1f} -> {cover(cur):6.1f}; {a.holdout} top families "
              f"covered {100 * hc:3.0f}%  | adds {why}")
        rows.append(dict(run=label, step=i, cell=c, gain=g, coverage=cover(cur),
                         holdout_covered=hc, adds=why))
hc, got = holdout_cover(base)
print(f"\n{a.holdout} top families (departure) and best coverage by the CURRENT panel:")
for i, g in got:
    print(f"  {i:<30} {tops[i]:+.2f}   best panel cell {g:+.2f}")

os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
Z.to_csv(f"{a.out_prefix}.motif_z.tsv", sep="\t", index_label="motif")
DEP.to_csv(f"{a.out_prefix}.family_departure.tsv", sep="\t", index_label="family")
REL.to_csv(f"{a.out_prefix}.family_reliable.tsv", sep="\t", index_label="family")
pd.DataFrame(rows).to_csv(f"{a.out_prefix}.selection.tsv", sep="\t", index=False)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

order = [a.holdout] + base + [r_["cell"] for r_ in rows if r_["run"] == "extend current panel"]
order += [c for c in CELLS if c not in order]
topf = P.max(axis=1).sort_values(ascending=False).head(40).index
M = DEP.loc[topf, order]
M = M.iloc[np.lexsort((-M.max(axis=1).values, M.values.argmax(axis=1)))]
fig, ax = plt.subplots(figsize=(0.42 * len(order) + 3, 0.26 * len(M) + 2))
v = np.nanpercentile(np.abs(M.values), 98)
im = ax.imshow(M.values, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
ax.set_xticks(range(len(order)))
ax.set_xticklabels(order, rotation=60, ha="right", fontsize=8)
ax.set_yticks(range(len(M)))
ax.set_yticklabels(M.index, fontsize=7)
nb = 1 + len(base)
ax.axvline(0.5, color="k", lw=0.8)
ax.axvline(nb - 0.5, color="k", lw=0.8)
ax.set_title(f"TF-family activity specific to each cell (ATAC/DNase motif association minus "
             f"mean)\n{a.holdout} | current panel | greedy additions in order | rest", fontsize=9)
fig.colorbar(im, ax=ax, fraction=0.03, label="departure (sd units)")
fig.tight_layout()
fig.savefig(f"{a.out_prefix}.heatmap.png", dpi=150)
print(f"\nwrote {a.out_prefix}.motif_z.tsv, .family_departure.tsv, .selection.tsv, .heatmap.png")
