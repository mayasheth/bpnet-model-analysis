#!/usr/bin/env python3
"""GATA oracle: give only GATA-motif elements another track's within-class RANKS.

WHAT THIS IS FOR. 4.41 found the transferred p300 model under-ranks GATA-motif elements by
~9 percentile points relative to the in-cell model (and over-ranks CEBP/ESR elements). This
asks how much of the benchmark gap between them that accounts for. Three tracks, all built on
the multi-cell (K562 held out) prediction:
    gata_incell   GATA-motif regions take the IN-CELL model's within-class ranks
    gata_obs      GATA-motif regions take OBSERVED p300's within-class ranks (upper bound)
    ctrl_incell   the same NUMBER of non-GATA regions, matched on class and ATAC decile, take
                  the in-cell ranks: the control for "replacing any N regions helps"
    ctrl_obs      the same matched non-GATA regions take OBSERVED ranks: gata_obs's control
                  (added 2026-09-26; --only ctrl_obs rebuilds just this track)
gata_incell - multicell is what knowing K562's GATA biology would buy; minus ctrl_incell is
the part specific to GATA.

Ranks transfer as in 4.31: a selected region gets the multi-cell track's value at the source
track's within-class quantile; every other region is unchanged.

GATA-motif = a FIMO hit (4.41 cache, p < 1e-4) for any of the GATA-family motifs that 4.41
ranked top: GATA1, GATA1::TAL1, GATA2, Gata3, GATA4, GATA5, GATA6, TRPS1.

Usage: 4.42.build_subset_oracle.py [--out-dir data/abc_predicted]
"""
import argparse

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import rankdata, spearmanr

D = "/oak/stanford/groups/engreitz/Users/sheth"
A = f"{D}/EP300_BPNet/2026_0824_H3K27ac_model"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction/results"
PRED_BW = f"{A}/data/abc_predicted/predp300mc_multicell_holdK562.bw"
SRC = {"multicell": f"{ABC}/2026_0922_p300_multicell/p300pred_multicell_holdK562",
       "incell": f"{ABC}/2026_0905_p300_activity/p300pred_k562_multimodal",
       "obs": f"{ABC}/2026_0905_p300_activity/p300obs_k562"}
CACHE = f"{A}/data/p300_error_anatomy/fimo_counts.npz"
GATA = ["MA0035.5 GATA1", "MA0140.3 GATA1::TAL1", "MA0036.4 GATA2", "MA0037.5 Gata3",
        "MA0482.3 GATA4", "MA0766.3 GATA5", "MA1104.3 GATA6", "MA1970.2 TRPS1"]
KEY = ["chr", "start", "end"]

ap = argparse.ArgumentParser()
ap.add_argument("--out-dir", default=f"{A}/data/abc_predicted")
ap.add_argument("--only", default=None, help="build only this track")
a = ap.parse_args()


def load(d, count=False):
    x = pd.read_csv(f"{d}/Neighborhoods/EnhancerList.txt", sep="\t")
    cols = KEY + ["class", "H3K27ac.RPM", "ATAC.RPM"]
    if count:
        rc = [c for c in x.columns if c.startswith("H3K27ac.") and c.endswith(".readCount")]
        assert len(rc) == 1, rc
        return x[cols + rc].rename(columns={rc[0]: "readCount"})
    return x[cols]


# Region order must match 4.41's, which is the observed arm's EnhancerList order.
reg = load(SRC["obs"]).rename(columns={"H3K27ac.RPM": "obs"})
mc = load(SRC["multicell"], count=True)
inc = load(SRC["incell"]).rename(columns={"H3K27ac.RPM": "incell"})
reg = reg.merge(mc[KEY + ["readCount"]], on=KEY, how="left", validate="one_to_one")
reg = reg.merge(inc[KEY + ["incell"]], on=KEY, how="left", validate="one_to_one")
assert reg[["readCount", "incell"]].notna().all().all()
reg["prom"] = reg["class"].isin(["tss", "promoter"])

z = np.load(CACHE, allow_pickle=True)
names = list(z["names"])
assert z["counts"].shape[0] == len(reg), "motif cache rows do not match the region set"
cols = [names.index(g) for g in GATA]
reg["gata"] = (z["counts"][:, cols] > 0).any(1)
print(f"{len(reg):,} regions; GATA-motif: {int(reg.gata.sum()):,} "
      f"({100 * reg.gata.mean():.1f}%; non-promoter {int((reg.gata & ~reg.prom).sum()):,}, "
      f"promoter {int((reg.gata & reg.prom).sum()):,})")

# ATAC-decile-matched control, per class
rng = np.random.default_rng(0)
reg["ctrl"] = False
for prom in (True, False):
    cls = reg.prom == prom
    dec = pd.qcut(rankdata(reg.loc[cls, "ATAC.RPM"]), 10, labels=False)
    dec = pd.Series(dec, index=reg.index[cls])
    for q in range(10):
        pool = dec.index[(dec == q) & ~reg.loc[dec.index, "gata"]]
        need = int(((dec == q) & reg.loc[dec.index, "gata"]).sum())
        reg.loc[rng.choice(pool, size=min(need, len(pool)), replace=False), "ctrl"] = True
print(f"matched control: {int(reg.ctrl.sum()):,} non-GATA regions")


def build(mask, src_col, out):
    target = reg["readCount"].astype(float).copy()
    for prom in (True, False):
        cls = reg.prom == prom
        q = reg.loc[cls, src_col].rank(method="first") / cls.sum()
        vals = np.sort(reg.loc[cls, "readCount"].values)
        t = cls & mask
        pos = np.clip(np.ceil(q[t] * len(vals)).astype(int) - 1, 0, len(vals) - 1)
        target[t] = vals[pos]
    new = {(r.chr, int(r.start), int(r.end)): v / (r.end - r.start)
           for r, v in zip(reg[mask].itertuples(), target[mask])}
    src, dst = pyBigWig.open(PRED_BW), pyBigWig.open(out, "w")
    ch = src.chroms()
    dst.addHeader([(c, ch[c]) for c in sorted(ch)])
    hit = 0
    for c in sorted(ch):
        iv = src.intervals(c) or []
        if not iv:
            continue
        v = []
        for s_, e_, val in iv:
            k = (c, s_, e_)
            v.append(float(new[k]) if k in new else float(val))
            hit += k in new
        dst.addEntries([c] * len(iv), [x[0] for x in iv], ends=[x[1] for x in iv], values=v)
    dst.close()
    src.close()
    assert hit == int(mask.sum()), f"repainted {hit} of {int(mask.sum())}"
    for prom in (True, False):
        s = mask & (reg.prom == prom)
        print(f"  {out.split('/')[-1]:<44} {'promoter' if prom else 'nonpromoter':<12} "
              f"rho(new, source) on selected {spearmanr(target[s], reg.loc[s, src_col])[0]:.4f}")


TRACKS = {"gata_incell": ("gata", "incell"), "gata_obs": ("gata", "obs"),
          "ctrl_incell": ("ctrl", "incell"), "ctrl_obs": ("ctrl", "obs")}
for tag, (mask, src) in TRACKS.items():
    if a.only and tag != a.only:
        continue
    build(reg[mask], src, f"{a.out_dir}/predp300gata_{tag}.bw")
