#!/usr/bin/env python3
"""9.0 Prepare the eQTL transfer test of the F-059 features (plan 9): models, features, variants.

WHY. F-059's three features were chosen and confirmed on the K562 CRISPR benchmark only. Fine-mapped
eQTLs in matched immune cell types (eQTL Catalogue v7) test them in cells the model never saw. This
script does everything that is the same for every cluster, once; 9.1 scores each cluster and 9.2
computes the metrics.

STEPS
1. Motif sites. scE2G clusters have no BINDetect sites, so GATA counts come from the genome-wide
   FIMO scan (4.44: JASPAR 2024, p < 1e-4, both strands) for training and application alike: sites
   of GATA6 MA1104.3, Gata3 MA0037.5 and TRPS1 MA1970.2, counted when the site midpoint lies within
   the element midpoint +/- 250 bp (7.1's window); GATA cluster = max over the three (8.2). IKZF2
   MA2326.1 is 6 bp wide, so FIMO at p < 1e-4 finds no site (first run, 2026-10-12); it is counted
   as exact consensus matches (AGGAAG, CTTCCT, overlapping) in the same window.
2. K562 training features. The three features are rebuilt on the 155,134 Xu elements with the
   functions 9.1 applies to clusters, encoded 1 + percentile over the cell's own elements (8.2), and
   mapped to the CRISPR pairs exactly as 8.2. Checks: CGI distance equals 7.1's library value, its
   pair-level percentile equals 8.2's candidates, and FIMO counts are reported against 7.1's
   BINDetect counts.
3. Models. scE2G's logistic regression (training_params.pkl, log(|x| + 0.01)) fitted on ALL K562
   pairs (the deployed form), per arm: n6; n6 + CGI distance; n6 + the three (FIMO counts).
   Replica LOCO AUPRC on K562 is reported for each arm and for the BINDetect version (descriptive).
4. Variants. eQTL Catalogue v7 fine tissues, PIP >= 0.5, distal noncoding (partition ABC /
   AllPeaks / Other / OtherIntron), eGene in the GENCODE v43 gene universe, matched cell types
   only; 1000G background SNPs filtered to distal noncoding the same way (as eQTLEnrichment).
5. Matching. Biosample -> eQTL cell types from the fine-grained config (config_predictions_hd_fine).

Usage: 9.0.eqtl_transfer_prep.py [--out-dir $SCRATCH/eqtl_transfer]
"""
import argparse
import importlib.util
import json
import os
import re

import numpy as np
import pandas as pd
import pysam
from scipy.stats import rankdata, spearmanr
from sklearn.linear_model import LogisticRegression

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("replica", os.path.join(HERE, "8.0.sce2g_replica.py"))
rep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rep)

SC = os.environ["SCRATCH"]
D = "/oak/stanford/groups/engreitz/Users/sheth"
EQ = f"{D}/eQTLEnrichment-integrated/eQTLEnrichment"
VARIANTS = f"{D}/Data/eQTL_Catalogue_v7/results/2024_1003_fine_tissues/eQTL_catalogue_v7.processed.PIP0.5.tsv.gz"
BG = f"{D}/hg38_resources/backgroundSNPs/all.bg.SNPs.hg38.baseline.v1.1.bed.sorted"
PARTITION = f"{EQ}/resources/genome_annotation/PartitionCombined.bed"
UNIVERSE = f"{EQ}/resources/genome_annotation/CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed6"
MATCH = f"{EQ}/config/mine/config_predictions_hd_fine.tsv"
CGI = f"{SC}/ucsc_hg38/cpgIslandExt.txt.gz"
HITS = f"{SC}/motif_hits"
RUN = f"{SC}/sce2g_runs/2026_1007_k562_xu_discovery/K562_Xu_pl"
GATA = ["GATA6_MA1104.3", "Gata3_MA0037.5", "TRPS1_MA1970.2"]
IKZF2 = "IKZF2_MA2326.1"
IKZF2_KMERS = ["AGGAAG", "CTTCCT"]
GENOME = f"{D}/hg38_resources/hg38.fa"
N6 = ["numTSSEnhGene", "normalizedATAC_prom", "numNearbyEnhancers", "ubiqExpressed", "numCandidateEnhGene", "ARC.E2G.Score"]
ARMS = {"n6": N6, "n6_cgi": N6 + ["cgi_p1"], "n6_three": N6 + ["cgi_p1", "gata_p1", "ikzf2_p1"],
        # activity-only positive control (F-064 caveat): ARC-E2G -> ABC, -0.050 AUPRC on K562 CRISPR (8.1)
        "n6_abc": [f for f in N6 if f != "ARC.E2G.Score"] + ["ABC.Score"]}
DISTAL = {"ABC", "AllPeaks", "Other", "OtherIntron"}
EK = ["chr", "start", "end"]
KEY = ["chr", "start", "end", "TargetGene", "name"]
HALF = 250

ap = argparse.ArgumentParser()
ap.add_argument("--out-dir", default=f"{SC}/eqtl_transfer")
a = ap.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
pct1 = lambda v: 1 + rankdata(np.nan_to_num(np.asarray(v, float), nan=-np.inf)) / len(v)  # noqa: E731


# ---------------------------------------------------------------- shared feature functions (9.1 imports these)
def cgi_log_dist(chrs, mids, G):
    """log10(1 + bp) from each midpoint to the nearest CGI midpoint on its chromosome (7.1)."""
    out = np.full(len(mids), np.nan)
    for c in np.unique(chrs):
        gm = np.sort(((G.loc[G["chr"] == c, "s"] + G.loc[G["chr"] == c, "e"]) // 2).values)
        i = np.where(chrs == c)[0]
        if len(gm):
            k = np.clip(np.searchsorted(gm, mids[i]), 1, max(len(gm) - 1, 1))
            d = np.minimum(np.abs(mids[i] - gm[k - 1]), np.abs(mids[i] - gm[np.minimum(k, len(gm) - 1)]))
            out[i] = np.log10(d + 1)
    return out


def site_counts(chrs, mids, S):
    """Sites of one motif (table chr, mid) with midpoint within element midpoint +/- HALF."""
    out = np.zeros(len(mids), np.int32)
    for c in np.unique(chrs):
        p = np.sort(S.loc[S["chr"] == c, "mid"].values)
        i = np.where(chrs == c)[0]
        out[i] = np.searchsorted(p, mids[i] + HALF, side="right") - np.searchsorted(p, mids[i] - HALF, side="left")
    return out


def kmer_counts(chrs, mids, kmers, fa):
    """Overlapping occurrences of any of `kmers` in the element midpoint +/- HALF."""
    pats = [re.compile(f"(?={k})") for k in kmers]
    out = np.zeros(len(mids), np.int32)
    for i, (c, m) in enumerate(zip(chrs, mids)):
        s = fa.fetch(c, max(0, int(m) - HALF), int(m) + HALF).upper()
        out[i] = sum(len(p.findall(s)) for p in pats)
    return out


def element_features(chrs, mids, G, H, fa):
    gata = np.max(np.column_stack([site_counts(chrs, mids, H[H["motif"] == m]) for m in GATA]), axis=1)
    # IKZF2 MA2326.1 is 6 bp wide, so no site can reach FIMO p < 1e-4 (a perfect match is ~2.4e-4);
    # count exact consensus matches on both strands instead (amendment, decisions [2026-10-12])
    ikz = kmer_counts(chrs, mids, IKZF2_KMERS, fa)
    cgi = cgi_log_dist(chrs, mids, G)
    return pd.DataFrame({"cgi_log_dist": cgi, "gata_n": gata, "ikzf2_n": ikz,
                         "cgi_p1": pct1(cgi), "gata_p1": pct1(gata), "ikzf2_p1": pct1(ikz)})


if __name__ == "__main__":
    # ------------------------------------------------------------ 1. motif sites
    hit_file = f"{a.out_dir}/motif_sites_gata_ikzf2.tsv.gz"
    if not os.path.exists(hit_file):
        rows = []
        for f in sorted(os.listdir(HITS)):
            z = np.load(f"{HITS}/{f}", allow_pickle=True)
            # 4.44 names motifs "MA1104.3 GATA6"; 7.1 / 8.2 use "GATA6_MA1104.3"
            names = ["_".join(str(x).split(" ", 1)[::-1]) for x in z["names"]]
            want = {names.index(m): m for m in GATA + [IKZF2]}
            assert len(want) == 4, f"{f}: motifs missing"
            sel = np.isin(z["motif"], list(want))
            mid = (z["start"][sel].astype(np.int64) + z["end"][sel]) // 2
            rows.append(pd.DataFrame({"chr": f.replace(".npz", ""), "motif": [want[i] for i in z["motif"][sel]], "mid": mid}))
            print(f"  {f}: {sel.sum():,} sites", flush=True)
        pd.concat(rows).to_csv(hit_file, sep="\t", index=False)
    H = pd.read_csv(hit_file, sep="\t")
    print(f"motif sites: {H['motif'].value_counts().to_dict()}", flush=True)
    G = pd.read_csv(CGI, sep="\t", header=None, usecols=[1, 2, 3], names=["chr", "s", "e"])

    # ------------------------------------------------------------ 2. K562 training features
    E = pd.read_csv(f"{SC}/feature_screen/library.elements.tsv.gz", sep="\t", usecols=EK + ["mid", "cgi_log_dist", "atac_xu"])
    F = element_features(E["chr"].values, E["mid"].values, G, H, pysam.FastaFile(GENOME))
    ok = np.isclose(F["cgi_log_dist"], E["cgi_log_dist"], equal_nan=True)
    assert ok.all(), f"{(~ok).sum()} Xu elements: CGI distance differs from 7.1"
    M = np.load(f"{SC}/feature_screen/library.motifs.npz", allow_pickle=True)
    mn = [str(x) for x in M["names"]]
    bind_gata = M["n"][:, [mn.index(m) for m in GATA]].max(1)
    bind_ikz = M["n"][:, mn.index(IKZF2)]
    print("FIMO vs BINDetect counts over Xu elements: GATA cluster Spearman "
          f"{spearmanr(F['gata_n'], bind_gata)[0]:.3f} (any site {np.mean(F['gata_n'] > 0):.3f} vs {np.mean(bind_gata > 0):.3f}); "
          f"IKZF2 {spearmanr(F['ikzf2_n'], bind_ikz)[0]:.3f} (any {np.mean(F['ikzf2_n'] > 0):.3f} vs {np.mean(bind_ikz > 0):.3f})", flush=True)
    for f_ in ["gata_n", "ikzf2_n"]:
        assert (F[f_] > 0).mean() > 0.05, f"{f_}: sites in {(F[f_] > 0).mean():.3f} of Xu elements, scan is broken"
    E = pd.concat([E, F.drop(columns="cgi_log_dist")], axis=1)

    df, miss, params = rep.load_run(RUN, "multiome_powerlaw_n6")
    el = df[EK].drop_duplicates().reset_index(drop=True)
    best = np.full(len(el), -1)
    for c, g in el.groupby("chr"):
        u_ = E.index[E["chr"] == c].values
        if len(u_) == 0:
            continue
        u_ = u_[np.argsort(E.loc[u_, "start"].values)]
        us, ue, at = E.loc[u_, "start"].values, E.loc[u_, "end"].values, E.loc[u_, "atac_xu"].values
        lo = np.searchsorted(ue, g["start"].values, side="right")
        hi = np.searchsorted(us, g["end"].values, side="left")
        for i, l, h in zip(g.index.values, lo, hi):
            if h > l:
                best[i] = u_[l + np.argmax(at[l:h])]
    P1 = ["cgi_p1", "gata_p1", "ikzf2_p1"]
    el[P1] = np.where(best[:, None] >= 0, E[P1].values[np.maximum(best, 0)], 1.0)
    df = df.merge(el, on=EK, how="left", validate="m:1")
    cand = pd.read_csv(f"{SC}/crispr_discovery/candidates.tsv.gz", sep="\t").drop(columns=["Regulated"]).drop_duplicates(KEY)
    BIND = ["cgi_log_dist", "motif:GATA6_MA1104.3+2", "motif:IKZF2_MA2326.1"]
    df = df.merge(cand[KEY + BIND], on=KEY, how="left", validate="1:1")
    bad = ~np.isclose(df["cgi_p1"], df["cgi_log_dist"])
    assert not bad.any(), f"{bad.sum()} pairs: CGI percentile differs from 8.2's candidates"
    df["ubiqExpressed"] = df["ubiqExpressed"].astype(float)

    # ------------------------------------------------------------ 3. models
    y = df["Regulated"].astype(bool).values
    models, rows = {}, []
    for arm, feats in {**ARMS, "n6_three_bindetect": N6 + BIND}.items():
        X = np.log(np.abs(df[feats].values.astype(float)) + rep.EPS)
        m = LogisticRegression(**params).fit(X, y.astype(int))
        if arm in ARMS:
            models[arm] = dict(features=feats, coef=m.coef_[0].tolist(), intercept=float(m.intercept_[0]), eps=rep.EPS)
        au = rep.auprc(y, rep.loco(df, feats, params), miss.sum(), len(miss))
        rows.append(dict(arm=arm, k562_loco_auprc=au, **{f"coef_{f}": c for f, c in zip(feats, m.coef_[0])}))
    R = pd.DataFrame(rows)
    R["vs_n6"] = R["k562_loco_auprc"] - R.loc[R["arm"] == "n6", "k562_loco_auprc"].iloc[0]
    print(R[["arm", "k562_loco_auprc", "vs_n6"]].round(4).to_string(index=False), flush=True)
    with open(f"{a.out_dir}/models.json", "w") as fh:
        json.dump(models, fh, indent=1)
    R.to_csv(f"{a.out_dir}/k562_arms.tsv", sep="\t", index=False, float_format="%.5g")

    # ------------------------------------------------------------ 5. matching
    mt = pd.read_csv(MATCH, sep="\t")
    mt = mt[["biosample", "scE2G_multiome", "GTExTissue"]].rename(columns={"scE2G_multiome": "pred_file"})
    mt = mt.assign(tissue=mt["GTExTissue"].str.split(",")).explode("tissue").drop(columns="GTExTissue")
    mt["tissue"] = mt["tissue"].str.strip()
    mt.to_csv(f"{a.out_dir}/matches.tsv", sep="\t", index=False)
    print(f"{mt['biosample'].nunique()} biosamples, {len(mt)} matched (biosample, eQTL cell type) pairs", flush=True)

    # ------------------------------------------------------------ 4. variants and background SNPs
    part = pd.read_csv(PARTITION, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chr", "s", "e", "cat"])
    part = part[part["cat"].isin(DISTAL)].sort_values(["chr", "s"])

    def in_partition(chrs, starts, ends):
        keep = np.zeros(len(starts), bool)
        for c in np.unique(chrs):
            p = part[part["chr"] == c]
            if p.empty:
                continue
            ps, pe = p["s"].values, np.maximum.accumulate(p["e"].values)
            i = np.where(chrs == c)[0]
            k = np.searchsorted(ps, ends[i], side="left") - 1
            keep[i[k >= 0]] = pe[k[k >= 0]] > starts[i[k >= 0]]
        return keep

    V = pd.read_csv(VARIANTS, sep="\t", usecols=["chr", "start", "end", "varID_hg38", "gene_hgnc", "tissue", "pip"])
    n0 = len(V)
    V = V[V["pip"] >= 0.5].drop_duplicates()
    V = V[in_partition(V["chr"].values, V["start"].values, V["end"].values)]
    uni = set(pd.read_csv(UNIVERSE, sep="\t", header=None, usecols=[3])[3])
    V = V[V["gene_hgnc"].isin(uni) & V["tissue"].isin(set(mt["tissue"]))]
    V.to_csv(f"{a.out_dir}/variants.tsv.gz", sep="\t", index=False)
    print(f"variants: {n0:,} -> {len(V):,} rows (PIP >= 0.5, distal noncoding, gene universe, matched cell types); "
          f"per cell type: {V.drop_duplicates(['chr', 'start', 'end', 'tissue'])['tissue'].value_counts().to_dict()}", flush=True)

    B = pd.read_csv(BG, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"])
    B = B[in_partition(B["chr"].values, B["start"].values, B["end"].values)]
    np.savez_compressed(f"{a.out_dir}/bg_snps.npz", **{c: g["start"].values.astype(np.int64)
                                                      for c, g in B.groupby("chr")})
    with open(f"{a.out_dir}/bg_count.txt", "w") as fh:
        fh.write(f"{len(B)}\n")
    print(f"background SNPs, distal noncoding: {len(B):,}", flush=True)
    print("done", flush=True)
