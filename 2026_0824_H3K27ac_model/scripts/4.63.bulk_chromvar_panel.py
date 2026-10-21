#!/usr/bin/env python3
"""Bulk chromVAR (Schep et al. 2017) TF-activity vector per cell, from each cell's own ATAC.

WHY. The multi-cell p300 model needs a label-free conditioning vector (FiLM) saying which TFs
are active in the target cell. At inference the target is K562, so the vector must come from
ATAC alone: no p300, no CRISPR. One deviation z-score per JASPAR 2024 motif per cell.

INPUTS (ATAC only; check_inputs() asserts every track is an ATAC 5' insertion bigwig).
- 5' insertion bigwigs: the seven panel cells in data/p300_panel/*_atac_5p.bw and K562
  2026_0529_multimodal_p300_model/data/atac_5p.bw (the track the accs5p models used), so all
  eight cells use the same track type (integer Tn5 5' counts per base).
- Motif hits: $SCRATCH/motif_hits/{chrom}.npz from 4.44 (JASPAR 2024 CORE vertebrate
  non-redundant, FIMO p < 1e-4, both strands, N-free).
- hg38 FASTA for GC (same file 4.44 scans).

ELEMENT SET. Union of ENCODE ATAC-seq IDR thresholded peaks of the eight cells
(config/panel_candidates.tsv rows with assay ATAC-seq; files already in $SCRATCH/panel_candidates
from 4.49). Why this and not K562 + GM12878 DNase candidates: a blood-only element set would
miss liver/breast/lung-specific elements, which is where HNF4A, FOXA, ESR1, GRHL2 and TEAD
signal lives, and would bias the vector toward K562/GM12878. Each cell contributes its own
peaks equally, peaks are called from ATAC only, and the set is fixed once built (written to
--elements and reused). Build: 500 bp around each peak summit, pooled, overlapping windows
merged, then one 500 bp window centred on each merged cluster's midpoint. Clusters are
disjoint, so the recentred windows do not overlap. chr1-22 and chrX only.

METHOD.
1. X[e, c] = sum of 5' insertions in window e for cell c. Elements with zero total dropped.
2. E[e, c] = (row total e / grand total) * column total c.
3. Annotation a_m[e] = 1 if a hit of motif m has its midpoint inside window e.
   raw[m, c] = (a_m.X[:, c] - a_m.E[:, c]) / a_m.E[:, c].
4. Background: elements binned on GC x log10 mean count (quantile bins); for each element,
   --n-bg background elements sampled uniformly from its own bin (with replacement), as in
   chromVAR's getBackgroundPeaks. Iteration i replaces every annotated element by its i-th
   background element and recomputes the deviation. z = (raw - mean_bg) / sd_bg.
   The `deviation` column is chromVAR's bias-corrected deviation, raw - mean_bg.
   With eight bulk samples the expectation is the panel average, so z is relative to the panel.

OUTPUTS. results/2026_1020_bulk_chromvar_panel.tsv (motif, motif_name, cell, deviation, z,
n_elements) and ..._z_wide.tsv for seed --seeds[0]; extra seeds go to ..._seed{s}_z_wide.tsv;
per-cell Spearman of z between seeds in ..._seed_stability.tsv. Counts, GC and annotation are
cached in --cache so reruns skip the bigwig pass.
"""
import argparse
import os
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr

R = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
A = f"{R}/2026_0824_H3K27ac_model"
SCR = os.environ.get("SCRATCH", "/tmp")
GENOME = "/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
PANEL = ["A549", "GM12878", "HepG2", "MCF_7", "stomach", "transverse_colon", "ovary"]
BIGWIGS = {c: f"{A}/data/p300_panel/{c}_atac_5p.bw" for c in PANEL}
BIGWIGS["K562"] = f"{R}/2026_0529_multimodal_p300_model/data/atac_5p.bw"
CELLS = ["K562"] + PANEL
MANIFEST_CELL = {"K562": "K562", "GM12878": "GM12878", "HepG2": "HepG2", "A549": "A549",
                 "MCF-7": "MCF_7", "ovary": "ovary", "stomach": "stomach",
                 "transverse colon": "transverse_colon"}
HALF = 250


def check_inputs(bigwigs, manifest):
    """Fail if any signal track is not an ATAC 5' bigwig or any peak file is not ATAC."""
    for c, p in bigwigs.items():
        b = os.path.basename(p).lower()
        if not b.endswith("atac_5p.bw") or "p300" in b or "ep300" in b or "crispr" in b:
            raise SystemExit(f"non-ATAC signal track for {c}: {p}")
        if not os.path.exists(p):
            raise SystemExit(f"missing: {p}")
    if not (manifest["assay"] == "ATAC-seq").all():
        raise SystemExit("element manifest has non-ATAC rows")


def build_elements(manifest, peak_dir):
    rows = []
    for _, r in manifest.iterrows():
        f = f"{peak_dir}/{r['accession']}.bed.gz"
        d = pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 9])
        d = d[d[0].isin(CHROMS)]
        s = d[1].values + d[9].values
        rows.append(pd.DataFrame({"chrom": d[0].values, "start": s - HALF, "end": s + HALF}))
        print(f"  {r['cell']}: {len(d):,} peaks ({r['accession']})", flush=True)
    w = pd.concat(rows).sort_values(["chrom", "start"])
    out = []
    for ch, g in w.groupby("chrom", sort=False):
        st, en = g["start"].values, g["end"].values
        run_end = np.maximum.accumulate(en)
        new = np.r_[True, st[1:] >= run_end[:-1]]
        cs = np.minimum.reduceat(st, np.flatnonzero(new))
        ce = np.maximum.reduceat(en, np.flatnonzero(new))
        mid = (cs + ce) // 2
        out.append(pd.DataFrame({"chrom": ch, "start": mid - HALF, "end": mid + HALF}))
    el = pd.concat(out)
    el = el[el["start"] >= 0]
    el["chrom"] = pd.Categorical(el["chrom"], CHROMS)
    return el.sort_values(["chrom", "start"]).reset_index(drop=True)


def count_cell(args):
    import pyBigWig
    cell, path, el = args
    bw = pyBigWig.open(path)
    chroms = bw.chroms()
    out = np.zeros(len(el), np.float64)
    for ch, g in el.groupby("chrom", observed=True):
        if ch not in chroms:
            print(f"  {cell}: {ch} absent from bigwig", flush=True)
            continue
        v = bw.values(ch, 0, chroms[ch], numpy=True)
        v[np.isnan(v)] = 0
        cs = np.concatenate([[0.0], np.cumsum(v, dtype=np.float64)])
        del v
        e = np.minimum(g["end"].values, chroms[ch])
        out[g.index.values] = cs[e] - cs[g["start"].values]
    bw.close()
    print(f"  counted {cell}: total {out.sum():,.0f}", flush=True)
    return out


def gc_content(el):
    import pyfaidx
    fa = pyfaidx.Fasta(GENOME, sequence_always_upper=True)
    gc = np.zeros(len(el))
    for ch, g in el.groupby("chrom", observed=True):
        s = np.frombuffer(str(fa[ch][:]).encode(), np.uint8)
        isgc = (s == ord("G")) | (s == ord("C"))
        isacgt = isgc | (s == ord("A")) | (s == ord("T"))
        c1 = np.concatenate([[0], np.cumsum(isgc)])
        c2 = np.concatenate([[0], np.cumsum(isacgt)])
        st, en = g["start"].values, g["end"].values
        gc[g.index.values] = (c1[en] - c1[st]) / np.maximum(c2[en] - c2[st], 1)
    return gc


def annotate(el, hit_dir):
    rows, cols, names = [], [], None
    for ch, g in el.groupby("chrom", observed=True):
        z = np.load(f"{hit_dir}/{ch}.npz", allow_pickle=True)
        nm = list(z["names"])
        if names is None:
            names = nm
        elif nm != names:
            raise SystemExit(f"motif order differs in {ch}.npz")
        mid = (z["start"].astype(np.int64) + z["end"].astype(np.int64)) // 2
        st, en = g["start"].values, g["end"].values
        j = np.searchsorted(st, mid, side="right") - 1
        ok = (j >= 0) & (mid < en[np.maximum(j, 0)])
        rows.append(g.index.values[j[ok]])
        cols.append(z["motif"][ok].astype(np.int64))
        print(f"  {ch}: {ok.sum():,} of {len(mid):,} hits in elements", flush=True)
    r, c = np.concatenate(rows), np.concatenate(cols)
    M = sp.csr_matrix((np.ones(len(r), np.float32), (r, c)), shape=(len(el), len(names)))
    M.data[:] = 1.0
    return M, names


def background(gc, logmean, n_bg, n_bins, seed):
    rng = np.random.default_rng(seed)
    gb = pd.qcut(gc, n_bins, labels=False, duplicates="drop")
    ab = pd.qcut(logmean, n_bins, labels=False, duplicates="drop")
    b = gb * (ab.max() + 1) + ab
    bg = np.empty((len(gc), n_bg), np.int64)
    for k in np.unique(b):
        idx = np.flatnonzero(b == k)
        bg[idx] = rng.choice(idx, size=(len(idx), n_bg), replace=True)
    return bg


def deviations(X, E, Ann, bg):
    num = Ann.T @ (X - E)
    den = Ann.T @ E
    raw = num / den
    D = X - E
    bgd = np.empty((bg.shape[1],) + raw.shape)
    for i in range(bg.shape[1]):
        bgd[i] = (Ann.T @ D[bg[:, i]]) / (Ann.T @ E[bg[:, i]])
    mu, sd = bgd.mean(0), bgd.std(0, ddof=1)
    return raw - mu, (raw - mu) / sd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default=f"{A}/config/panel_candidates.tsv")
    ap.add_argument("--peak-dir", default=f"{SCR}/panel_candidates")
    ap.add_argument("--hit-dir", default=f"{SCR}/motif_hits")
    ap.add_argument("--cache", default=f"{SCR}/bulk_chromvar/counts_gc_annot.npz")
    ap.add_argument("--elements", default=f"{A}/data/bulk_chromvar/elements_atac_idr_union_500bp.bed.gz")
    ap.add_argument("--out-prefix", default=f"{A}/results/2026_1020_bulk_chromvar_panel")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1])
    ap.add_argument("--n-bg", type=int, default=50)
    ap.add_argument("--n-bins", type=int, default=20)
    ap.add_argument("--procs", type=int, default=8)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.cache), exist_ok=True)
    os.makedirs(os.path.dirname(a.elements), exist_ok=True)

    man = pd.read_csv(a.manifest, sep="\t")
    man = man[man["cell"].isin(MANIFEST_CELL)].copy()
    check_inputs(BIGWIGS, man)
    if sorted(MANIFEST_CELL[c] for c in man["cell"]) != sorted(CELLS):
        raise SystemExit(f"manifest cells {list(man['cell'])} do not match {CELLS}")
    t0 = time.time()

    if os.path.exists(a.cache):
        z = np.load(a.cache, allow_pickle=True)
        if list(z["cells"]) != CELLS:
            raise SystemExit("cache built for other cells")
        el = pd.read_csv(a.elements, sep="\t", header=None, names=["chrom", "start", "end"])
        X, gc, names = z["X"], z["gc"], list(z["names"])
        Ann = sp.csr_matrix((z["a_data"], z["a_ind"], z["a_ptr"]), shape=tuple(z["a_shape"]))
        print(f"loaded cache {a.cache}", flush=True)
    else:
        print("building element set", flush=True)
        el = build_elements(man, a.peak_dir)
        print(f"{len(el):,} elements", flush=True)
        with Pool(min(a.procs, len(CELLS))) as p:
            cols = p.map(count_cell, [(c, BIGWIGS[c], el) for c in CELLS])
        X = np.column_stack(cols)
        keep = X.sum(1) > 0
        print(f"dropping {(~keep).sum():,} zero-count elements", flush=True)
        el, X = el[keep].reset_index(drop=True), X[keep]
        print(f"counts done ({time.time() - t0:.0f} s); GC", flush=True)
        gc = gc_content(el)
        print(f"annotating ({time.time() - t0:.0f} s)", flush=True)
        Ann, names = annotate(el, a.hit_dir)
        el.to_csv(a.elements, sep="\t", header=False, index=False, compression="gzip")
        tmp = a.cache[:-4] + ".tmp.npz"
        np.savez_compressed(tmp, X=X, gc=gc, names=np.array(names, dtype=object),
                            cells=np.array(CELLS, dtype=object), a_data=Ann.data,
                            a_ind=Ann.indices, a_ptr=Ann.indptr, a_shape=np.array(Ann.shape))
        os.replace(tmp, a.cache)
        print(f"wrote {a.elements} and {a.cache}", flush=True)

    print(f"elements {X.shape[0]:,}; motifs {len(names)}; column totals in elements: "
          + ", ".join(f"{c} {t:,.0f}" for c, t in zip(CELLS, X.sum(0))), flush=True)
    np.seterr(invalid="ignore")
    E = X.sum(1, keepdims=True) / X.sum() * X.sum(0, keepdims=True)
    logmean = np.log10(X.mean(1) + 1e-3)
    n_el = np.asarray(Ann.sum(0)).ravel().astype(int)
    ids = [n.split(" ", 1)[0] for n in names]
    tf = [n.split(" ", 1)[1] if " " in n else n for n in names]

    zs = {}
    for s in a.seeds:
        bg = background(gc, logmean, a.n_bg, a.n_bins, s)
        dev, zz = deviations(X, E, Ann, bg)
        zs[s] = zz
        wide = pd.DataFrame(zz, columns=CELLS)
        wide.insert(0, "motif_name", tf)
        wide.insert(0, "motif", ids)
        if s == a.seeds[0]:
            long = pd.DataFrame({"motif": np.repeat(ids, len(CELLS)),
                                 "motif_name": np.repeat(tf, len(CELLS)),
                                 "cell": np.tile(CELLS, len(ids)),
                                 "deviation": dev.ravel(), "z": zz.ravel(),
                                 "n_elements": np.repeat(n_el, len(CELLS))})
            long.to_csv(f"{a.out_prefix}.tsv", sep="\t", index=False, float_format="%.5g")
            wide.to_csv(f"{a.out_prefix}_z_wide.tsv", sep="\t", index=False, float_format="%.5g")
            main_z = wide
        else:
            wide.to_csv(f"{a.out_prefix}_seed{s}_z_wide.tsv", sep="\t", index=False,
                        float_format="%.5g")
        print(f"seed {s} done ({time.time() - t0:.0f} s)", flush=True)

    if len(a.seeds) > 1:
        s0, s1 = a.seeds[:2]
        stab = pd.DataFrame({"cell": CELLS, "spearman_z": [
            spearmanr(zs[s0][:, k], zs[s1][:, k], nan_policy="omit")[0] for k in range(len(CELLS))]})
        stab.to_csv(f"{a.out_prefix}_seed_stability.tsv", sep="\t", index=False)
        print(f"\nseed stability (Spearman z, seed {s0} vs {s1}):\n{stab.to_string(index=False)}")

    print(f"\nmotifs with 0 annotated elements (z = NaN, ranked last): {(n_el == 0).sum()}: "
          + ", ".join(np.array(tf)[n_el == 0]))
    print(f"motifs with < 100 annotated elements: {(n_el < 100).sum()}; "
          f"median n_elements {int(np.median(n_el)):,}")
    for c in CELLS:
        top = main_z.sort_values(c, ascending=False).head(10)
        print(f"\ntop 10 {c}: " + ", ".join(f"{n} {v:.1f}" for n, v in zip(top["motif_name"], top[c])))

    controls = {"K562": ["GATA1", "GATA2", "GATA1::TAL1"],
                "HepG2": ["HNF4A", "HNF1A", "FOXA1", "FOXA2"],
                "MCF_7": ["ESR1", "GRHL2", "FOXA1"],
                "GM12878": ["SPI1", "IRF4", "IRF8", "RUNX1", "EBF1"],
                "A549": ["FOS::JUN", "FOSL1::JUN", "TEAD1", "TEAD4"]}
    families = {"K562": r"^GATA", "HepG2": r"^(HNF4|HNF1|FOXA)", "MCF_7": r"^(ESR|GRHL|FOXA)",
                "GM12878": r"^(SPI1|IRF|RUNX|EBF)", "A549": r"^(FOS|JUN|FOSL|TEAD)"}
    print("\ncontrol ranks (1 = highest z of", len(names), "motifs)")
    for c, tfs in controls.items():
        rk = main_z[c].rank(ascending=False, method="min", na_option="bottom").astype(int)
        up = main_z["motif_name"].str.upper()
        parts = []
        for t in tfs:
            hit = up == t.upper()
            parts.append(f"{t} {rk[hit].min() if hit.any() else 'absent'}")
        fam = up.str.contains(families[c], regex=True)
        parts.append(f"family best {rk[fam].min()} (n={fam.sum()})")
        print(f"  {c}: " + "; ".join(parts))
    print(f"\ndone ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
