#!/usr/bin/env python3
"""8.9 Near-CGI encodings of CpG-island distance for scE2G, on the 8.0 replica.

WHY. F-059's CGI-distance feature is encoded as 1 + percentile of log distance to the nearest CpG
island over all Xu elements, which assumes "farther is better" over the whole range. F-060 found the
cell-general part is the near end only: the nearest decile has less p300 than ATAC predicts in all
five panel cells, while the far-end trend differs by cell (rises in K562, HepG2, A549; flat in
MCF-7; falls in GM12878). If an encoding that only marks the near end keeps most of the K562 CRISPR
gain, that is the version to carry into a test in another cell (Maya, 2026-10-09).

ENCODINGS, all label-free (thresholds are percentiles over all 155,134 Xu elements, from 8.2's
encoding; p = percentile, nearest = 0):
  pct      1 + p (the F-059 feature)
  near10   1 + [p > 0.10]   PRIMARY: flags the nearest decile, as F-060 found
  near20   1 + [p > 0.20]
  cap10    1 + min(p, 0.10) / 0.10   graded inside the nearest decile, flat beyond
  cap20    1 + min(p, 0.20) / 0.20
Models: n6, n6 + each encoding, n6 + near10 + GATA + IKZF2, and n6 + the Phase 4 three (pct).
LOCO + paired bootstrap (2,000) on the replica, CRISPR pairs that miss every candidate at score 0.
A pre-specified comparison of five encodings, not a search, so no nested selection.

READING, fixed before the run:
- near10 "keeps most of the gain" if its gain over n6 is >= 75% of pct's. Then use near10 for the
  external test, since it encodes only the part F-060 found in every cell.
- Otherwise the far end also carries K562 CRISPR signal; keep pct for the external test, and expect
  the far-end part to transfer less (F-060).
- The others (near20, cap10, cap20) are descriptive; picking the best of them after the fact
  would be selection on these pairs.

Usage: 8.9.near_cgi_encoding.py [--out results/2026_1009_near_cgi_encoding.tsv]
"""
import argparse
import importlib.util
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("replica", os.path.join(HERE, "8.0.sce2g_replica.py"))
rep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rep)

SC = os.environ["SCRATCH"]
KEY = ["chr", "start", "end", "TargetGene", "name"]
ap = argparse.ArgumentParser()
ap.add_argument("--run", default=f"{SC}/sce2g_runs/2026_1007_k562_xu_discovery/K562_Xu_pl")
ap.add_argument("--n-boot", type=int, default=2000)
ap.add_argument("--out", default="results/2026_1009_near_cgi_encoding.tsv")
a = ap.parse_args()

# raw distances, only to report what the thresholds mean in bp
E = pd.read_csv(f"{SC}/feature_screen/library.elements.tsv.gz", sep="\t", usecols=["cgi_log_dist"])
q = E["cgi_log_dist"].quantile([0.1, 0.2, 0.5])
print(f"{len(E):,} Xu elements; distance to nearest CGI midpoint at the 10th / 20th / 50th percentile: "
      + " / ".join(f"{10 ** v - 1:,.0f} bp" for v in q), flush=True)

df, miss, params = rep.load_run(a.run, "multiome_powerlaw_n6")
cand = pd.read_csv(f"{SC}/crispr_discovery/candidates.tsv.gz", sep="\t").drop(columns=["Regulated"]).drop_duplicates(KEY)
MOT = ["motif:GATA6_MA1104.3+2", "motif:IKZF2_MA2326.1"]
df = df.merge(cand[KEY + ["cgi_log_dist"] + MOT], on=KEY, how="left", validate="1:1")
assert df[["cgi_log_dist"] + MOT].notna().all().all()
assert df["cgi_log_dist"].between(1, 2).all(), "expected 8.2's 1 + percentile encoding"
p = df["cgi_log_dist"].values - 1
df["pct"] = 1 + p
df["near10"] = 1 + (p > 0.10)
df["near20"] = 1 + (p > 0.20)
df["cap10"] = 1 + np.minimum(p, 0.10) / 0.10
df["cap20"] = 1 + np.minimum(p, 0.20) / 0.20
y = df["Regulated"].astype(bool).values
print(f"{len(df):,} pairs, {y.sum()} regulated (+ {len(miss):,} missing, {miss.sum()} regulated); "
      f"nearest decile holds {(p <= 0.10).mean():.3f} of pairs, regulated rate {y[p <= 0.10].mean():.3f} "
      f"vs {y[p > 0.10].mean():.3f}", flush=True)

N6 = ["numTSSEnhGene", "normalizedATAC_prom", "numNearbyEnhancers", "ubiqExpressed", "numCandidateEnhGene", "ARC.E2G.Score"]
ENC = ["pct", "near10", "near20", "cap10", "cap20"]
MODELS = {"n6": N6, **{f"n6 + {e}": N6 + [e] for e in ENC},
          "n6 + near10 + GATA + IKZF2": N6 + ["near10"] + MOT,
          "n6 + pct + GATA + IKZF2 (Phase 4)": N6 + ["pct"] + MOT}
S = {m: rep.loco(df, f, params) for m, f in MODELS.items()}
ys = np.r_[y, miss]
Sx = {m: np.r_[s, np.zeros(len(miss))] for m, s in S.items()}
rng = np.random.default_rng(0)
nb = len(ys)
REFS = ["n6", "n6 + pct", "n6 + pct + GATA + IKZF2 (Phase 4)"]
boot = {(m, r): [] for m in MODELS for r in REFS}
for _ in range(a.n_boot):
    ii = rng.integers(0, nb, nb)
    au = {m: rep.auprc(ys[ii], Sx[m][ii]) for m in MODELS}
    for m in MODELS:
        for r in REFS:
            boot[(m, r)].append(au[m] - au[r])
A = {m: rep.auprc(y, S[m], miss.sum(), len(miss)) for m in MODELS}
rows = []
for m in MODELS:
    row = dict(model=m, auprc=A[m])
    for r, tag in zip(REFS, ["n6", "pct", "phase4"]):
        b = np.array(boot[(m, r)])
        row.update({f"vs_{tag}": A[m] - A[r], f"lo_{tag}": np.quantile(b, 0.025), f"hi_{tag}": np.quantile(b, 0.975)})
    rows.append(row)
R = pd.DataFrame(rows)
pd.set_option("display.width", 220)
print(R.round(4).to_string(index=False))

g_pct, g_near = A["n6 + pct"] - A["n6"], A["n6 + near10"] - A["n6"]
share = g_near / g_pct
print(f"\nnear10 keeps {share:.0%} of pct's gain ({g_near:+.4f} of {g_pct:+.4f}) -> "
      + ("near10 for the external test" if share >= 0.75 else "keep pct; far end carries K562 signal too"))
R.to_csv(a.out, sep="\t", index=False, float_format="%.5g")
print(f"wrote {a.out}")
