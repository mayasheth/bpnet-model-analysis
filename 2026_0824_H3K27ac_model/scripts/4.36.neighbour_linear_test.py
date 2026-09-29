#!/usr/bin/env python3
"""Do an element's neighbours say anything about its p300 that its own prediction does not?

WHY. A graph network over the sequence+accessibility predictions, with power-law edges within
1 Mb, is only worth building if neighbour information carries signal the per-element model
misses. This is the one-layer linear version of that graph: for each element, the
power-law-weighted mean and sum of its neighbours' predictions and ATAC, fitted with a linear
model. If this does not move the metric, a graph network over the same node values has little
to find; if it does, the graph has something to learn.

MODELS, fitted per ABC class (promoter / non-promoter), all in within-class percentile ranks
because ABC keeps only ranks:
    M0  the prediction as ABC saw it (no fit)
    M1  obs ~ pred + atac                          (refit without neighbours: the comparator)
    M2  obs ~ pred + atac + neighbour terms        (the test)
    M3a gradient boosting on pred + atac           (nonlinear comparator)
    M3  gradient boosting on the M2 features       (nonlinear, with neighbours)
M2 against M1, and M3 against M3a, are the contrasts, since each pair is fitted the same way. Neighbour terms, with
w_ij = ABC's power law (gamma 1.0242, scale 5.9595, clipped at 5 kb), j within 1 Mb, j != i:
    mean_w(pred_j), mean_w(atac_j), log sum_w(pred_j), log sum_w(atac_j)

CROSS-VALIDATION. Five chromosome groups; every element is predicted by a model fitted on the
other four, so the scored track is fully held out. Edges never cross a chromosome.

LABELS. The fit uses K562 observed p300. For the multi-cell (transfer) model this is NOT a
deployable track: it is an upper bound on what these neighbour terms can give. A deployable
version would fit the few coefficients in the training cell types.

METRICS, per class, over all candidate regions (never the tested elements alone):
    partial_rho   rank correlation with observed p300 given ATAC (as in 4.35)
    local_rho     correlation of (x_i - mean of x_j within 50 kb) between track and observed,
                  i.e. how well the track ranks elements WITHIN a locus. Smoothing raises
                  partial_rho while lowering this; a real gain should not lower it.

GRAPH. --graph same (the first version) builds one graph per ABC class, so an enhancer's
neighbours are only enhancers. --graph cross builds one graph over all regions and splits
every neighbour term by the NEIGHBOUR's class (promoter / non-promoter), so an enhancer sees
the promoters near it and vice versa; the models are still fitted and scored within class,
which is equivalent to a class indicator interacting with every feature. cross also adds
    M4  boosting pooled over both classes with an is_promoter indicator
which lets the 21k promoters borrow strength from the 132k non-promoters.

Usage: 4.36.neighbour_linear_test.py [--graph same|cross] [--out-tsv T] [--no-gbm]
       [--write-tracks DIR]

--write-tracks writes each arm's held-out boosting tracks (M3, M3a) per candidate region, for
4.37 to paint into bigwigs and run through ABC.
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ABC = "/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction/results"
OBS = f"{ABC}/2026_0905_p300_activity/p300obs_k562"
ARMS = [
    ("K562 in-cell", f"{ABC}/2026_0905_p300_activity/p300pred_k562_multimodal"),
    ("multi-cell, K562 held out", f"{ABC}/2026_0922_p300_multicell/p300pred_multicell_holdK562"),
    ("ATAC only (pred = ATAC)", None),
]
KEY = ["chr", "start", "end"]
GAMMA, SCALE, MIN_D, MAX_D, LOCAL_D = 1.024238616787792, 5.9594510043736655, 5000, 1_000_000, 50_000

ap = argparse.ArgumentParser()
ap.add_argument("--out-tsv", default=None)
ap.add_argument("--no-gbm", action="store_true")
ap.add_argument("--write-tracks", default=None, metavar="DIR")
ap.add_argument("--graph", choices=("same", "cross"), default="same")
a = ap.parse_args()


def load(d):
    return pd.read_csv(f"{d}/Neighborhoods/EnhancerList.txt", sep="\t",
                       usecols=KEY + ["class", "H3K27ac.RPM", "ATAC.RPM"])


def pct(v):
    return (rankdata(v) - 0.5) / len(v)


def resid(y, z):
    Z = np.column_stack([np.ones_like(z), z])
    return y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]


def partial_rho(x, y, z):
    return float(np.corrcoef(resid(rankdata(x), rankdata(z)), resid(rankdata(y), rankdata(z)))[0, 1])


def graphs(g):
    """Sparse neighbour matrices for elements sorted by chr, mid: W = ABC power law within
    1 Mb, L = unweighted mean within 50 kb (row-normalised). Self excluded. Edges stay within
    a chromosome."""
    from scipy import sparse
    ri, ci, wv, li, lj = [], [], [], [], []
    for _, c in g.groupby("chr", sort=False):
        idx, mid = c.index.values, c["mid"].values
        for k in range(len(mid)):
            lo, hi = np.searchsorted(mid, [mid[k] - MAX_D, mid[k] + MAX_D], side="left")[0], \
                np.searchsorted(mid, mid[k] + MAX_D, side="right")
            j = np.r_[lo:k, k + 1:hi]
            d = np.abs(mid[j] - mid[k]).astype(float)
            ri.append(np.full(len(j), idx[k])); ci.append(idx[j])
            wv.append(np.exp(SCALE - GAMMA * np.log(np.clip(d, MIN_D, None) + 1)))
            near = j[d <= LOCAL_D]
            li.append(np.full(len(near), idx[k])); lj.append(idx[near])
    n = len(g)
    W = sparse.csr_matrix((np.concatenate(wv), (np.concatenate(ri), np.concatenate(ci))), (n, n))
    L = sparse.csr_matrix((np.ones(sum(map(len, li))), (np.concatenate(li), np.concatenate(lj))),
                          (n, n))
    deg = np.asarray(L.sum(1)).ravel()
    L = sparse.diags(np.where(deg > 0, 1 / np.maximum(deg, 1), 0)) @ L
    return W, L, deg > 0


def neighbour_feats(W, v):
    s = W @ v
    return s / np.asarray(W.sum(1)).ravel().clip(1e-12), np.log1p(s)


def local_rho(track, obs, L, has):
    """Within-locus contrast: each value minus the mean of its neighbours within 50 kb."""
    t, o = track[has] - (L @ track)[has], obs[has] - (L @ obs)[has]
    return float(np.corrcoef(t, o)[0, 1])


obs = load(OBS).rename(columns={"H3K27ac.RPM": "obs_raw", "ATAC.RPM": "atac_raw"})
obs["mid"] = (obs["start"] + obs["end"]) // 2
obs["prom"] = obs["class"].isin(["tss", "promoter"])
CHR_ORDER = sorted(obs["chr"].unique())
grp = {c: i % 5 for i, c in enumerate(CHR_ORDER)}          # 5 chromosome groups
obs["fold"] = obs["chr"].map(grp)

G = {}
for cls, sel in (("nonpromoter", ~obs["prom"]), ("promoter", obs["prom"])):
    g0 = obs[sel].sort_values(["chr", "mid"]).reset_index(drop=True)
    G[cls] = (g0, *graphs(g0))
    print(f"graph {cls}: {len(g0):,} nodes, {G[cls][1].nnz:,} edges within 1 Mb, "
          f"median {int(np.median(np.diff(G[cls][1].indptr)))} neighbours", flush=True)

if a.graph == "cross":
    from scipy import sparse
    gall = obs.sort_values(["chr", "mid"]).reset_index(drop=True)
    Wall, _, _ = graphs(gall)
    isp = gall["prom"].values
    WP = Wall @ sparse.diags(isp.astype(float))          # edges to promoter neighbours
    WE = Wall @ sparse.diags((~isp).astype(float))       # edges to non-promoter neighbours
    print(f"cross graph: {len(gall):,} nodes, {Wall.nnz:,} edges; "
          f"{WP.nnz:,} to promoters, {WE.nnz:,} to non-promoters", flush=True)


def within_class_pct(df, col):
    out = np.empty(len(df))
    for sel in (df["prom"].values, ~df["prom"].values):
        out[sel] = pct(df.loc[sel, col].values)
    return out


rows, coefs, out_tracks = [], [], {}
for arm, d in ARMS:
    if d is None:
        m = obs.assign(pred_raw=obs["atac_raw"])
    else:
        p = load(d).rename(columns={"H3K27ac.RPM": "pred_raw", "ATAC.RPM": "atac_arm"})
        m = obs.merge(p[KEY + ["pred_raw", "atac_arm"]], on=KEY, validate="one_to_one")
        assert len(m) == len(obs) and np.allclose(m["atac_raw"], m["atac_arm"])
    if a.graph == "cross":
        ga = gall[KEY + ["prom"]].merge(m[KEY + ["pred_raw", "atac_raw", "obs_raw", "fold"]],
                                        on=KEY, how="left", validate="one_to_one")
        for c in ("pred", "atac", "obs"):
            ga[c] = within_class_pct(ga, f"{c}_raw")
        for c in ("pred", "atac"):
            for side, Wx in (("P", WP), ("E", WE)):
                ga[f"{c}_nmean{side}"], ga[f"{c}_nlogsum{side}"] = neighbour_feats(Wx, ga[c].values)
        NEIGH = [f"{c}_{k}{side}" for c in ("pred", "atac") for k in ("nmean", "nlogsum")
                 for side in ("P", "E")]
        pooled = {}
        if not a.no_gbm:
            from sklearn.ensemble import HistGradientBoostingRegressor
            Xa = ga[["pred", "atac"] + NEIGH].assign(is_prom=ga["prom"].astype(float))
            for name, feats in (("M4a pooled boosting + indicator, no neighbours",
                                 ["pred", "atac", "is_prom"]),
                                ("M4 pooled boosting + indicator, + neighbours",
                                 ["pred", "atac", "is_prom"] + NEIGH)):
                t = np.zeros(len(ga))
                for f in range(5):
                    tr, te = ga["fold"].values != f, ga["fold"].values == f
                    gb = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05,
                                                       max_leaf_nodes=31, random_state=0)
                    gb.fit(Xa.loc[tr, feats].values, ga.loc[tr, "obs"].values)
                    t[te] = gb.predict(Xa.loc[te, feats].values)
                pooled[name] = t
        ga_tr = ga[KEY + NEIGH].assign(**pooled)
    for cls in ("nonpromoter", "promoter"):
        g0, W, L, has = G[cls]
        g = g0[KEY].merge(m[KEY + ["pred_raw", "atac_raw", "obs_raw", "fold"]], on=KEY,
                          how="left", validate="one_to_one")
        assert g["pred_raw"].notna().all()
        g["pred"], g["atac"], g["obs"] = pct(g["pred_raw"]), pct(g["atac_raw"]), pct(g["obs_raw"])
        base = ["pred", "atac"]
        if a.graph == "cross":
            nf = g[KEY].merge(ga_tr, on=KEY, how="left", validate="one_to_one")
            neigh = NEIGH
        else:
            nf = pd.DataFrame(index=g.index)
            for c in ("pred", "atac"):
                nf[f"{c}_nmean"], nf[f"{c}_nlogsum"] = neighbour_feats(W, g[c].values)
            neigh = ["pred_nmean", "atac_nmean", "pred_nlogsum", "atac_nlogsum"]
        X = pd.concat([g[base], nf[neigh]], axis=1)
        tracks = {"M0 prediction as is": g["pred"].values}
        for name, feats in (("M1 refit, no neighbours", base), ("M2 + neighbour terms", base + neigh)):
            t = np.zeros(len(g))
            for f in range(5):
                tr, te = g["fold"].values != f, g["fold"].values == f
                A = np.column_stack([np.ones(tr.sum()), X.loc[tr, feats].values])
                b = np.linalg.lstsq(A, g.loc[tr, "obs"].values, rcond=None)[0]
                t[te] = np.column_stack([np.ones(te.sum()), X.loc[te, feats].values]) @ b
            A = np.column_stack([np.ones(len(g)), (X[feats] - X[feats].mean()) / X[feats].std()])
            b = np.linalg.lstsq(A, g["obs"].values, rcond=None)[0]
            coefs.append(dict(arm=arm, cls=cls, model=name,
                              **{f"std_coef_{k}": v for k, v in zip(feats, b[1:])}))
            tracks[name] = t
        if not a.no_gbm:
            from sklearn.ensemble import HistGradientBoostingRegressor
            # M3a is M3's comparator: boosting can reshape ATAC nonlinearly, which alone raises
            # a partial correlation that removes ATAC only linearly in ranks.
            for name, feats in (("M3a boosting, no neighbours", base),
                                ("M3 boosting, + neighbour terms", base + neigh)):
                t = np.zeros(len(g))
                for f in range(5):
                    tr, te = g["fold"].values != f, g["fold"].values == f
                    gb = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05,
                                                       max_leaf_nodes=31, random_state=0)
                    gb.fit(X.loc[tr, feats].values, g.loc[tr, "obs"].values)
                    t[te] = gb.predict(X.loc[te, feats].values)
                tracks[name] = t
        if a.graph == "cross" and not a.no_gbm:
            for name in pooled:
                tracks[name] = nf[name].values
        if a.write_tracks and not a.no_gbm:
            out_tracks.setdefault(arm, []).append(g[KEY].assign(
                cls=cls, boost_nbr=tracks["M3 boosting, + neighbour terms"],
                boost_nonbr=tracks["M3a boosting, no neighbours"]))
        for name, t in tracks.items():
            rows.append(dict(arm=arm, cls=cls, model=name, n=len(g),
                             rho_obs=float(np.corrcoef(rankdata(t), g["obs"])[0, 1]),
                             partial_rho=partial_rho(t, g["obs"].values, g["atac"].values),
                             local_rho=local_rho(rankdata(t) / len(t), g["obs"].values, L, has)))
        print(f"done: {arm} / {cls}", flush=True)

r = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print("\n" + r.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
print("\nstandardised coefficients (full-data fit, for reading only):")
print(pd.DataFrame(coefs).to_string(index=False, float_format=lambda v: f"{v:+.3f}"))
if a.write_tracks:
    import os
    os.makedirs(a.write_tracks, exist_ok=True)
    tag = {"K562 in-cell": "incell", "multi-cell, K562 held out": "multicell",
           "ATAC only (pred = ATAC)": "atac"}
    for arm, parts in out_tracks.items():
        f = f"{a.write_tracks}/nbr_tracks_{tag[arm]}.tsv.gz"
        pd.concat(parts).to_csv(f, sep="\t", index=False)
        print("wrote", f)
if a.out_tsv:
    r.to_csv(a.out_tsv, sep="\t", index=False)
    print("\nwrote", a.out_tsv)
