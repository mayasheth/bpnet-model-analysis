#!/usr/bin/env python3
"""9.2 eQTL enrichment and linking recall for each arm, pooled over matched clusters (plan 9).

Definitions follow eQTLEnrichment (gather_aggregate_enrichment_recall.R), on 9.1's reductions:
  enrichment(t) = (variant locations in predicted elements / variant locations)
                  / (background SNPs in predicted elements / background SNPs)
  recall(t)     = (variant-eGene pairs whose variant is in an element linked to the eGene / pairs)
with "predicted" = scE2G quantile >= t, counts summed over the matched (biosample, eQTL cell type)
pairs of a group before the ratios, as the pipeline's all-matches curve.

READOUTS (pre-registered, decisions [2026-10-12]):
  R1  enrichment at recall 0.10 (log-enrichment interpolated linearly in recall)
  R2  mean enrichment over recall 0.02 to RMAX (50 evenly spaced points); RMAX = 0.20 unless the
      n6 curve of the primary group does not reach 0.20, then 0.9 x its maximum
Intervals: 1,000 bootstrap resamples of eGenes, the same draws for every arm (paired); a variant
location with several eGenes is assigned to its first eGene for the enrichment counts.

GROUPS: primary = every matched pair except K562 (K562 is the training cell); also K562 alone,
each sample (PBMC, BMMC22, GM12878) and each eQTL cell type.

Usage: 9.2.eqtl_transfer_gather.py [--out-dir $SCRATCH/eqtl_transfer] [--results results/2026_1012_eqtl_transfer]
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

SC = os.environ["SCRATCH"]
ap = argparse.ArgumentParser()
ap.add_argument("--out-dir", default=f"{SC}/eqtl_transfer")
ap.add_argument("--results", default="results/2026_1012_eqtl_transfer")
ap.add_argument("--n-boot", type=int, default=1000)
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--extra-arms", nargs="*", default=[], help="arms written by 9.3, merged in")
a = ap.parse_args()
od = a.out_dir

mt = pd.read_csv(f"{od}/matches.tsv", sep="\t")
done = sorted(set(mt["biosample"]))  # biosample names contain dots, so files are looked up by name
missing = [b for b in done if not os.path.exists(f"{od}/clusters/{b}.variants.tsv.gz")]
missing += [f"{b}.{x}" for x in a.extra_arms for b in done if not os.path.exists(f"{od}/clusters/{b}.{x}.variants.tsv.gz")]
assert not missing, f"clusters not scored: {missing}"
n_bg = int(open(f"{od}/bg_count.txt").read())
# rows differing only in pip (same variant, eGene, cell type from two studies) collapse here; the
# metrics count distinct variants and variant-eGene pairs, so this changes no number
VK = ["chr", "start", "end", "varID_hg38", "gene_hgnc", "tissue", "biosample"]
V = pd.concat([pd.read_csv(f"{od}/clusters/{b}.variants.tsv.gz", sep="\t") for b in done])
V = V.drop(columns="pip").groupby(VK, as_index=False).max()
BG = {b: pd.read_csv(f"{od}/clusters/{b}.bg.tsv", sep="\t") for b in done}
for x in a.extra_arms:
    X = pd.concat([pd.read_csv(f"{od}/clusters/{b}.{x}.variants.tsv.gz", sep="\t") for b in done])
    X = X.drop(columns="pip").groupby(VK, as_index=False).max()
    n0 = len(V)
    V = V.merge(X, on=VK, how="left", validate="1:1")
    assert len(V) == n0 and V[f"elem_{x}"].notna().all(), f"{x}: variant rows do not match 9.1's"
    for b in done:
        xb = pd.read_csv(f"{od}/clusters/{b}.{x}.bg.tsv", sep="\t")
        assert np.allclose(xb["threshold"], BG[b]["threshold"])
        BG[b][x] = xb[x].values
GRID = BG[next(iter(BG))]["threshold"].values
ARMS = [c.replace("elem_", "") for c in V.columns if c.startswith("elem_")]
V = V.merge(mt[["biosample", "tissue"]], on=["biosample", "tissue"], how="inner")  # matched pairs only
V["pair"] = V["biosample"] + "|" + V["tissue"]
V["loc"] = V["chr"] + ":" + V["start"].astype(str)

# rows for linking: distinct (pair, loc, eGene); rows for enrichment: distinct (pair, loc), first eGene
L = V.groupby(["pair", "biosample", "tissue", "loc", "gene_hgnc"], as_index=False)[[f"link_{k}" for k in ARMS]].max()
En = (V.sort_values("gene_hgnc").groupby(["pair", "biosample", "tissue", "loc"], as_index=False)
      .agg(gene_hgnc=("gene_hgnc", "first"), **{f"elem_{k}": (f"elem_{k}", "max") for k in ARMS}))
genes = np.array(sorted(set(L["gene_hgnc"]) | set(En["gene_hgnc"])))
gi = {g: i for i, g in enumerate(genes)}
L["g"] = L["gene_hgnc"].map(gi)
En["g"] = En["gene_hgnc"].map(gi)


def kidx(v):
    """Number of grid thresholds <= v (so the row counts at thresholds 0..k-1); -1 values count nowhere."""
    return np.searchsorted(GRID, v, side="right")


def counts(k, w):
    """Weighted count of rows with value >= GRID[t] for every t."""
    h = np.bincount(k, weights=w, minlength=len(GRID) + 1)
    return np.cumsum(h[::-1])[::-1][1:]


def prepare(pairs):
    """Cut a group's rows once: gene index and grid position per row and arm, and summed bg counts."""
    l_, e_ = L[L["pair"].isin(pairs)], En[En["pair"].isin(pairs)]
    return dict(n=len(pairs), gl=l_["g"].values, ge=e_["g"].values,
                kl={k: kidx(l_[f"link_{k}"].values) for k in ARMS}, ke={k: kidx(e_[f"elem_{k}"].values) for k in ARMS},
                x2={k: sum(BG[p.split("|")[0]][k].values for p in pairs) / (n_bg * len(pairs)) for k in ARMS})


def curves(pg, W):
    """Recall and enrichment per arm on the grid for gene weights W (len(genes))."""
    wl, we = W[pg["gl"]], W[pg["ge"]]
    out = {}
    for k in ARMS:
        rec = counts(pg["kl"][k], wl) / wl.sum()
        x1 = counts(pg["ke"][k], we) / we.sum()
        with np.errstate(divide="ignore", invalid="ignore"):
            out[k] = (rec, np.where(pg["x2"][k] > 0, x1 / pg["x2"][k], np.nan))
    return out


def readouts(rec, enr, rmax):
    o = np.argsort(rec)
    r, le = rec[o], np.log(enr[o])
    ok = np.isfinite(le) & (enr[o] > 0)
    r, le = r[ok], le[ok]
    r, idx = np.unique(r, return_index=True)
    le = le[idx]
    f = lambda x: float(np.exp(np.interp(x, r, le))) if r[0] <= x <= r[-1] else np.nan  # noqa: E731
    xs = np.linspace(0.02, rmax, 50)
    vals = [f(x) for x in xs]
    return f(0.10), (float(np.mean(vals)) if np.all(np.isfinite(vals)) else np.nan), float(r[-1])


pairs_all = sorted(set(L["pair"]))
groups = {"primary (all matches except K562)": [p for p in pairs_all if not p.startswith("K562|")],
          "K562 (training cell)": [p for p in pairs_all if p.startswith("K562|")],
          "PBMC": [p for p in pairs_all if p.startswith("PBMC_")],
          "BMMC22": [p for p in pairs_all if p.startswith("BMMC22_")],
          "GM12878": [p for p in pairs_all if p.startswith("GM12878|")]}
for t in sorted(set(L["tissue"])):
    groups[f"cell type: {t}"] = [p for p in groups["primary (all matches except K562)"] if p.endswith(f"|{t}")]
groups = {g: p for g, p in groups.items() if p}

ones = np.ones(len(genes))
PG = {g: prepare(p) for g, p in groups.items()}
c0 = curves(PG["primary (all matches except K562)"], ones)
n6max = float(np.nanmax(c0["n6"][0]))
RMAX = 0.20 if n6max >= 0.20 else 0.9 * n6max
print(f"primary n6 maximum recall {n6max:.3f} -> RMAX {RMAX:.3f}", flush=True)
pd.DataFrame({"threshold": GRID, **{f"{k}_{m}": c0[k][i] for k in ARMS for i, m in enumerate(["recall", "enrichment"])}}
             ).to_csv(f"{a.results}.primary_curves.tsv", sep="\t", index=False, float_format="%.5g")

rng = np.random.default_rng(a.seed)
Wb = rng.multinomial(len(genes), np.full(len(genes), 1 / len(genes)), size=a.n_boot).astype(float)
rows = []
for gname, pairs in groups.items():
    pt = curves(PG[gname], ones)
    est = {k: readouts(*pt[k], RMAX) for k in ARMS}
    boot = {k: [] for k in ARMS}
    for b in range(a.n_boot):
        cb = curves(PG[gname], Wb[b])
        for k in ARMS:
            boot[k].append(readouts(*cb[k], RMAX)[:2])
    nl = int(L[L["pair"].isin(pairs)].shape[0])
    for k in ARMS:
        bk = np.array(boot[k])
        row = dict(group=gname, n_pairs=len(pairs), n_variant_gene=nl, arm=k, enr_at_r10=est[k][0], enr_mean=est[k][1],
                   max_recall=est[k][2])
        for j, rname in enumerate(["enr_at_r10", "enr_mean"]):
            d = bk[:, j] - np.array(boot["n6"])[:, j]
            d = d[np.isfinite(d)]
            row[f"{rname}_vs_n6"] = est[k][j] - est["n6"][j]
            row[f"{rname}_vs_n6_lo"] = np.quantile(d, 0.025) if len(d) else np.nan
            row[f"{rname}_vs_n6_hi"] = np.quantile(d, 0.975) if len(d) else np.nan
            row[f"{rname}_n_boot_finite"] = len(d)
        rows.append(row)
    print(f"  {gname}: {len(pairs)} pairs, {nl:,} variant-gene pairs", flush=True)
R = pd.DataFrame(rows)
R.to_csv(f"{a.results}.readouts.tsv", sep="\t", index=False, float_format="%.5g")
pd.set_option("display.width", 250)
show = ["group", "arm", "n_pairs", "n_variant_gene", "max_recall", "enr_at_r10", "enr_at_r10_vs_n6", "enr_at_r10_vs_n6_lo",
        "enr_at_r10_vs_n6_hi", "enr_mean", "enr_mean_vs_n6", "enr_mean_vs_n6_lo", "enr_mean_vs_n6_hi"]
print(R[show].round(3).to_string(index=False))

# pre-registered verdict
p = R[(R["group"] == "primary (all matches except K562)")].set_index("arm")
for arm in ["n6_three", "n6_cgi"]:
    up = [p.loc[arm, f"{r}_vs_n6_lo"] > 0 for r in ["enr_at_r10", "enr_mean"]]
    down = [p.loc[arm, f"{r}_vs_n6_hi"] < 0 for r in ["enr_at_r10", "enr_mean"]]
    v = "PASS" if all(up) else ("WORSE" if any(down) else ("PARTIAL" if any(up) else "FAIL"))
    print(f"{arm} vs n6, primary: R1 {p.loc[arm, 'enr_at_r10_vs_n6']:+.3f} [{p.loc[arm, 'enr_at_r10_vs_n6_lo']:+.3f}, "
          f"{p.loc[arm, 'enr_at_r10_vs_n6_hi']:+.3f}], R2 {p.loc[arm, 'enr_mean_vs_n6']:+.3f} "
          f"[{p.loc[arm, 'enr_mean_vs_n6_lo']:+.3f}, {p.loc[arm, 'enr_mean_vs_n6_hi']:+.3f}] -> {v}")
json.dump(dict(rmax=RMAX, n6_max_recall_primary=n6max, n_genes=len(genes), n_boot=a.n_boot),
          open(f"{a.results}.meta.json", "w"), indent=1)
print(f"wrote {a.results}.readouts.tsv, .primary_curves.tsv, .meta.json")
