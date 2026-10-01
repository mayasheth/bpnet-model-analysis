#!/usr/bin/env python3
"""Paint K562 activity tracks that combine p300 and H3K27ac, observed and predicted.

WHY. p300 and H3K27ac both mark activity but carry different information beyond ATAC: p300 is
high where the central site is sharper and more open, H3K27ac where the element sits in a broad
accessible domain (F-039, fig35). If the observed combination beats either observed track in
ABC or rE2G, that is the ceiling a two-target model could aim for; if not, a joint model is
not worth training. The predicted combinations need no new training: they combine tracks that
already exist and are deployment-legal in K562.

    obs_k27          observed H3K27ac, 5' ChIP bigwigs summed over the candidate
    obs_p300         observed p300 (ENCSR000EGE), same
    obs_p300k27      sqrt(q(p300) x q(H3K27ac))
    pred_p300k27     sqrt(q(predicted p300, multi-cell MSE-only, K562 held out) x
                          q(corrected transferred H3K27ac, 4.57 tx_breadth_frag))
    pred_p300ch_k27  same with the assoc_noctcf TF-activity-channel p300 model (F-038)
    pred_p300_k27tx  predicted p300 (MSE-only) x uncorrected transferred H3K27ac (tx_repaint)

q() maps each track's within-cell rank onto ABC's own H3K27ac quantile reference
(EnhancersQNormRef.K562.txt, enh_class "any"), which is what ABC's run_qnorm does to anything in
its H3K27ac column, so the product weighs the two tracks the way ABC would. Painted over the
July K562 candidate set as value / width (ABC sums over the region; qnorm is rank-based).

Usage: 4.58.combine_p300_h3k27ac.py
"""
import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import rankdata

D = "/oak/stanford/groups/engreitz/Users/sheth"
P = f"{D}/EP300_BPNet"
A = f"{P}/2026_0824_H3K27ac_model"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction"
PRED = f"{A}/data/abc_predicted"
CANDS = f"{ABC}/results/2026_0721_h3k27ac_counting_comparison/K562_ATAC_only/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed"
SIZES = f"{ABC}/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv"
QREF = f"{ABC}/reference/EnhancersQNormRef.K562.txt"
TRACKS = {
    "k27": [f"{A}/data/h3k27ac_5p_plus.bw", f"{A}/data/h3k27ac_5p_minus.bw"],
    "p300": [f"{P}/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig",
             f"{P}/2025_0703_retrain_p300_model/data/ENCSR000EGE_minus.bigWig"],
    "pp300": [f"{PRED}/predp300mc_multicell_holdK562.bw"],
    "pp300ch": [f"{PRED}/predp300mc_multicell_holdK562_motifact_assoc_noctcf.bw"],
    "pk27": [f"{PRED}/predk27frag_tx_breadth_frag.bw"],
    "pk27tx": [f"{PRED}/predk27frag_tx_repaint.bw"],
}

g = pd.read_csv(CANDS, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"])
sizes = dict(l.split()[:2] for l in open(SIZES))
order = [c for c in sizes if c in set(g["chr"])]
g = g.assign(_o=g["chr"].map({c: i for i, c in enumerate(order)})).sort_values(["_o", "start"])
g = g.reset_index(drop=True)


def bw_sum(paths):
    out = np.zeros(len(g))
    for p in paths:
        bw = pyBigWig.open(p)
        L = bw.chroms()
        for i, (c, s, e) in enumerate(zip(g["chr"].values, g["start"].values, g["end"].values)):
            if c in L:
                v = bw.values(c, int(s), min(int(e), L[c]), numpy=True)
                out[i] += np.nansum(np.abs(v))
        bw.close()
    return out


ref = pd.read_csv(QREF, sep="\t")
ref = ref[ref["enh_class"] == "any"].sort_values("quantile")
eps = ref.loc[ref["H3K27ac.RPM"] > 0, "H3K27ac.RPM"].min() / 2


def q(v):
    pc = (rankdata(v) - 0.5) / len(v)
    return np.interp(pc, ref["quantile"], ref["H3K27ac.RPM"]) + eps


V = {k: bw_sum(p) for k, p in TRACKS.items()}
for k, v in V.items():
    print(f"{k:<8} zero regions {np.mean(v == 0):.1%}", flush=True)
Q = {k: q(v) for k, v in V.items()}
rho = lambda x, y: np.corrcoef(rankdata(x), rankdata(y))[0, 1]
print(f"Spearman observed p300 vs observed H3K27ac: {rho(V['p300'], V['k27']):.3f}; "
      f"predicted p300 vs corrected predicted H3K27ac: {rho(V['pp300'], V['pk27']):.3f}")
OUT = {
    "obs_k27": V["k27"] + eps,
    "obs_p300": V["p300"] + eps,
    "obs_p300k27": np.sqrt(Q["p300"] * Q["k27"]),
    "pred_p300k27": np.sqrt(Q["pp300"] * Q["pk27"]),
    "pred_p300ch_k27": np.sqrt(Q["pp300ch"] * Q["pk27"]),
    "pred_p300_k27tx": np.sqrt(Q["pp300"] * Q["pk27tx"]),
}
width = (g["end"] - g["start"]).values
for name, v in OUT.items():
    out = f"{PRED}/predk27frag_{name}.bw"
    bw = pyBigWig.open(out, "w")
    bw.addHeader([(c, int(sizes[c])) for c in order])
    bw.addEntries(g["chr"].tolist(), g["start"].tolist(), ends=g["end"].tolist(),
                  values=(v / width).astype(float).tolist())
    bw.close()
    print(f"wrote {out}; Spearman with observed H3K27ac {rho(v, V['k27']):.3f}, "
          f"with observed p300 {rho(v, V['p300']):.3f}", flush=True)
