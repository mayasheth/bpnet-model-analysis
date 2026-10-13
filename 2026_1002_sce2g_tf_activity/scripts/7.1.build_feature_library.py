#!/usr/bin/env python3
"""7.1 Deployable element feature library over the Xu K562 scE2G candidates, for the p300 screen
(7.2, plan in plans/7_feature_screen.md). Sequence and Xu multiome only; no ENCODE K562 assay.

UCSC hg38 conservation (phyloP100way), RepeatMasker and CpG islands from $SCRATCH/ucsc_hg38.
Elements, target (p300), controls (atac_xu, log_dist_tss, prom), GC, CpG o/e and footprint
reference columns come from the 2.4 tables. Windows are the element midpoint +/- 250 bp unless
named otherwise.

WRITES ($SCRATCH/feature_screen/, too large for git)
    library.elements.tsv.gz   one row per element: keys, target, controls, non-motif features
    library.motifs.npz        n (int16) and score_max (float32), elements x 879 JASPAR motifs, + names;
                              sites with their midpoint within the element midpoint +/- 250 bp

Run in the tobias pixi env. Usage: 7.1.build_feature_library.py [--chroms chr21,chr22] [--out-dir DIR]
"""
import argparse
import glob
import os
import re
import shutil
import subprocess

import numpy as np
import pandas as pd
import pyBigWig
import pysam

D = "/oak/stanford/groups/engreitz/Users/sheth"
A2 = f"{D}/EP300_BPNet/2026_1002_sce2g_tf_activity"
SC = os.environ["SCRATCH"]
TB = f"{SC}/tobias_k562_xu"
GENOME = f"{D}/hg38_resources/hg38.fa"
UMAP = f"{D}/hg38_resources/k36.Umap.MultiTrackMappability.bw"
UCSC = f"{SC}/ucsc_hg38"  # hg38.phyloP100way.bw, rmsk.txt.gz, cpgIslandExt.txt.gz (log/ucsc_download)
REPEAT_CLASSES = ["SINE", "LINE", "LTR", "DNA", "Simple_repeat", "Low_complexity", "Satellite"]
FRAG = "/oak/stanford/groups/engreitz/Projects/scE2G/data/K562_260122/K562_Xu/atac_fragments.tsv.gz"
MAT = f"{SC}/featA/atac_matrix"
HALF, CENTRE, NEIGH = 250, 50, 50_000
KEY = ["chr", "start", "end"]

ap = argparse.ArgumentParser()
ap.add_argument("--chroms", default="", help="comma list; default all (smoke tests use one)")
ap.add_argument("--out-dir", default=None, help="default $SCRATCH/feature_screen, or .../smoke_<chroms>")
a = ap.parse_args()
a.out_dir = a.out_dir or (f"{SC}/feature_screen/smoke_{a.chroms.replace(',', '_')}" if a.chroms
                          else f"{SC}/feature_screen")
os.makedirs(a.out_dir, exist_ok=True)

# ---------------------------------------------------------------- elements, target, controls
base = ["mid", "prom", "log_dist_tss", "gc", "cpg_oe", "atac_xu", "p300", "fp_any"]
E = pd.read_csv(f"{A2}/results/2026_1003_footprint_p300_gate.print.element_features.tsv.gz", sep="\t",
                usecols=KEY + base + ["erythroid_fp_max", "CTCF_fp_max"]).rename(
    columns={"erythroid_fp_max": "ref_ery_fp_print", "CTCF_fp_max": "ref_ctcf_fp_print", "fp_any": "ref_fp_any_tobias"})
T = pd.read_csv(f"{A2}/results/2026_1003_footprint_p300_gate.element_features.tsv.gz", sep="\t",
                usecols=KEY + ["erythroid_fp_max"]).rename(columns={"erythroid_fp_max": "ref_ery_fp_tobias"})
E = E.merge(T, on=KEY, how="left", validate="1:1")
assert not E.duplicated(KEY).any()
if a.chroms:
    E = E[E["chr"].isin(a.chroms.split(","))]
E = E.sort_values(KEY).reset_index(drop=True)
print(f"{len(E):,} elements on {E['chr'].nunique()} chromosomes", flush=True)
mid = E["mid"].values

# ---------------------------------------------------------------- composition
DINUC = [x + y for x in "ACGT" for y in "ACGT"]
fa = pysam.FastaFile(GENOME)
comp = np.zeros((len(E), len(DINUC) + 3), np.float32)
run = re.compile(r"A+|T+")
for i, (c, m) in enumerate(zip(E["chr"].values, mid)):
    s = fa.fetch(c, max(0, m - HALF), m + HALF).upper()
    pairs = [s[j:j + 2] for j in range(len(s) - 1)]
    cnt = pd.Series(pairs).value_counts()
    tot = max(sum(cnt.get(d, 0) for d in DINUC), 1)
    comp[i, :16] = [cnt.get(d, 0) / tot for d in DINUC]
    runs = [len(r) for r in run.findall(s)] or [0]
    comp[i, 16:] = cnt.get("CG", 0), max(runs), sum(r >= 5 for r in runs)
for k, name in enumerate([f"dinuc_{d}" for d in DINUC] + ["cpg_n", "at_run_max", "polyAT5_n"]):
    E[name] = comp[:, k]
print("composition done", flush=True)

# ---------------------------------------------------------------- element shape
E["width"] = E["end"] - E["start"]
bw = pyBigWig.open(UMAP)
E["mappability"] = [np.nan_to_num(bw.stats(c, int(max(0, m - HALF)), int(m + HALF), type="mean")[0] or 0.0)
                    for c, m in zip(E["chr"].values, mid)]
bw.close()

# ---------------------------------------------------------------- conservation, repeats, CpG islands (UCSC)
bw = pyBigWig.open(f"{UCSC}/hg38.phyloP100way.bw")
cons = np.zeros((len(E), 2))
for i, (c, m) in enumerate(zip(E["chr"].values, mid)):
    v = bw.values(c, int(max(0, m - HALF)), int(m + HALF), numpy=True)
    v = v[np.isfinite(v)]
    cons[i] = (v.mean(), (v > 2).mean()) if len(v) else (0.0, 0.0)
bw.close()
E["phylop_mean"], E["phylop_frac_gt2"] = cons[:, 0], cons[:, 1]


def covered_bp(iv, c, lo, hi):
    """bp of [lo, hi) covered by the merged intervals of iv on chromosome c."""
    g = iv[iv["chr"] == c].sort_values("s")
    if g.empty:
        return np.zeros(len(lo))
    s_, e_ = g["s"].values, np.maximum.accumulate(g["e"].values)
    keep = np.r_[True, s_[1:] > e_[:-1]]  # merge overlapping intervals
    grp = np.cumsum(keep) - 1
    ms = s_[keep]
    me = pd.Series(e_).groupby(grp).max().values
    cum = np.r_[0, np.cumsum(me - ms)]

    def upto(x):  # covered bp in (-inf, x)
        k = np.searchsorted(ms, x, side="right")
        part = np.clip(x - ms[np.maximum(k - 1, 0)], 0, me[np.maximum(k - 1, 0)] - ms[np.maximum(k - 1, 0)])
        return cum[np.maximum(k - 1, 0)] + np.where(k > 0, part, 0)
    return upto(hi) - upto(lo)


R = pd.read_csv(f"{UCSC}/rmsk.txt.gz", sep="\t", header=None, usecols=[5, 6, 7, 11], names=["chr", "s", "e", "cls"])
G = pd.read_csv(f"{UCSC}/cpgIslandExt.txt.gz", sep="\t", header=None, usecols=[1, 2, 3], names=["chr", "s", "e"])
w = 2 * HALF
for name, iv in [("repeat_any", R)] + [(f"repeat_{k.lower()}", R[R["cls"] == k]) for k in REPEAT_CLASSES] + [("cgi", G)]:
    E[f"{name}_frac"] = 0.0
    for c, g in E.groupby("chr"):
        E.loc[g.index, f"{name}_frac"] = covered_bp(iv, c, g["mid"].values - HALF, g["mid"].values + HALF) / w
E["cgi_log_dist"] = np.nan
for c, g in E.groupby("chr"):
    gm = np.sort(((G.loc[G["chr"] == c, "s"] + G.loc[G["chr"] == c, "e"]) // 2).values)
    if len(gm):
        k = np.clip(np.searchsorted(gm, g["mid"].values), 1, max(len(gm) - 1, 1))
        d = np.minimum(np.abs(g["mid"].values - gm[k - 1]), np.abs(g["mid"].values - gm[np.minimum(k, len(gm) - 1)]))
        E.loc[g.index, "cgi_log_dist"] = np.log10(d + 1)
assert (E.filter(like="_frac").max() <= 1 + 1e-9).all()
print("UCSC: " + ", ".join(f"{c} {E[c].mean():.3f}" for c in ["phylop_mean", "repeat_any_frac", "cgi_frac"]) + " (means)", flush=True)

# ---------------------------------------------------------------- neighbourhood
E["neigh_n50k"] = 0
E["neigh_atac50k"] = 0.0
for c, g in E.groupby("chr"):
    m, at = g["mid"].values, g["atac_xu"].values
    cs = np.concatenate([[0], np.cumsum(at)])
    lo, hi = np.searchsorted(m, m - NEIGH), np.searchsorted(m, m + NEIGH, side="right")
    E.loc[g.index, "neigh_n50k"] = hi - lo - 1
    E.loc[g.index, "neigh_atac50k"] = cs[hi] - cs[lo] - at

# ---------------------------------------------------------------- cells (Kendall matrix, peaks x cells)
nr, nc = (int(x) for x in open(f"{MAT}.dim.txt").read().split())
i_ = np.fromfile(f"{MAT}.i.int32", np.int32)
x_ = np.fromfile(f"{MAT}.x.float64", np.float64)
nnz = np.bincount(i_, minlength=nr)
tot = np.bincount(i_, weights=x_, minlength=nr)
rows = pd.Series(open(f"{MAT}.rows.txt").read().split()).str.split("-", expand=True)
K = pd.DataFrame({"chr": rows[0], "kmid": (rows[1].astype(int) + rows[2].astype(int)) // 2,
                  "frac": nnz / nc, "mean_pos": np.where(nnz > 0, tot / np.maximum(nnz, 1), 0)})
cells = np.full((len(E), 2), np.nan)
for c, g in E.groupby("chr"):
    k = K[K["chr"] == c].sort_values("kmid")
    km = k["kmid"].values
    lo, hi = np.searchsorted(km, g["start"].values), np.searchsorted(km, g["end"].values)
    for j, (l, h) in enumerate(zip(lo, hi)):
        if h > l:
            sel = k.iloc[l:h]
            b = sel["frac"].values.argmax()
            cells[g.index[j]] = sel["frac"].values[b], sel["mean_pos"].values[b]
E["cells_frac_open"], E["cells_mean_count"] = cells[:, 0], cells[:, 1]
print(f"cells: {np.isfinite(cells[:, 0]).mean():.4f} of elements matched a Kendall peak", flush=True)
cols_bc = set(open(f"{MAT}.cols.txt").read().split())

# ---------------------------------------------------------------- fragments (insertion-weighted)
TABIX = shutil.which("tabix")
fr = np.full((len(E), 4), np.nan)
for c, g in E.groupby("chr"):
    if TABIX:
        p = subprocess.Popen([TABIX, FRAG, c], stdout=subprocess.PIPE)
        F = pd.read_csv(p.stdout, sep="\t", header=None, usecols=[1, 2, 3], names=["s", "e", "bc"])
        assert p.wait() == 0, f"tabix failed on {c}"
    else:
        F = pd.DataFrame([ln.split("\t")[1:4] for ln in pysam.TabixFile(FRAG).fetch(c)], columns=["s", "e", "bc"])
        F[["s", "e"]] = F[["s", "e"]].astype(int)
    F = F[F["bc"].isin(cols_bc)]
    L = (F["e"] - F["s"]).values
    pos = np.concatenate([F["s"].values, F["e"].values - 1])
    L2 = np.concatenate([L, L])
    o = np.argsort(pos, kind="stable")
    pos, L2 = pos[o], L2[o]
    cs = {k: np.concatenate([[0], np.cumsum(v)]) for k, v in
          (("short", L2 < 100), ("mono", (L2 >= 180) & (L2 <= 247)), ("logL", np.log2(L2)))}
    m = g["mid"].values
    lo, hi = np.searchsorted(pos, m - HALF), np.searchsorted(pos, m + HALF)
    clo, chi = np.searchsorted(pos, m - CENTRE), np.searchsorted(pos, m + CENTRE)
    n = hi - lo
    ok = n > 0
    nn = np.maximum(n, 1)
    vals = np.column_stack([(cs["short"][hi] - cs["short"][lo]) / nn, (cs["mono"][hi] - cs["mono"][lo]) / nn,
                            (cs["logL"][hi] - cs["logL"][lo]) / nn, (chi - clo) / nn])
    vals[~ok] = np.nan
    fr[g.index.values] = vals
    print(f"  {c}: {len(F):,} fragments from called cells, {ok.mean():.4f} of elements with insertions", flush=True)
for k, name in enumerate(["frag_short_frac", "frag_mono_frac", "frag_mean_log2_len", "ins_centre_frac"]):
    E[name] = fr[:, k]

E.to_csv(f"{a.out_dir}/library.elements.tsv.gz", sep="\t", index=False, float_format="%.6g")
print(f"wrote {a.out_dir}/library.elements.tsv.gz ({E.shape[1]} columns)", flush=True)

# ---------------------------------------------------------------- motifs (BINDetect sites per peak)
ix = pd.Series(np.arange(len(E)), index=pd.MultiIndex.from_frame(E[KEY]))
files = sorted(glob.glob(f"{TB}/bindetect/*/*_overview.txt"))
assert len(files) == 879, f"{len(files)} BINDetect motifs, expected 879"
names, N, S = [], np.zeros((len(E), len(files)), np.int16), np.zeros((len(E), len(files)), np.float32)
for j, f in enumerate(files):
    s = pd.read_csv(f, sep="\t", usecols=["TFBS_name", "TFBS_start", "TFBS_end", "TFBS_score", "peak_chr",
                                         "peak_start", "peak_end"])
    names.append(s["TFBS_name"].iloc[0] if len(s) else os.path.basename(os.path.dirname(f)))
    r = ix.reindex(pd.MultiIndex.from_frame(s[["peak_chr", "peak_start", "peak_end"]].set_axis(KEY, axis=1))).values
    if not a.chroms:
        assert np.isfinite(r).all(), f"{f}: sites outside the element set"
    keep = np.isfinite(r)
    # same window as every other feature: site midpoint within the element midpoint +/- HALF
    smid = (s["TFBS_start"].values + s["TFBS_end"].values) // 2
    keep[keep] = np.abs(smid[keep] - mid[r[keep].astype(int)]) <= HALF
    agg = s[keep].assign(r=r[keep].astype(int)).groupby("r")["TFBS_score"].agg(["size", "max"])
    N[agg.index.values, j] = np.minimum(agg["size"].values, np.iinfo(np.int16).max)
    S[agg.index.values, j] = agg["max"].values
    if j % 100 == 0:
        print(f"  motifs {j}/{len(files)}", flush=True)
assert len(set(names)) == len(names), "duplicate motif names"
np.savez_compressed(f"{a.out_dir}/library.motifs.npz", n=N, score_max=S, names=np.array(names),
                    chr=E["chr"].values, start=E["start"].values, end=E["end"].values)
print(f"wrote {a.out_dir}/library.motifs.npz: {len(names)} motifs; median motifs per element "
      f"{int(np.median((N > 0).sum(1)))}", flush=True)
