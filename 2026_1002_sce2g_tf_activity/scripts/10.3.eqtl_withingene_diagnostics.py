#!/usr/bin/env python3
"""10.3 Diagnostics for plan 10's failed calibration (F-066): unblinded arms only (n6, n6_abc,
v3_published) plus distance-only and random references, E1/E2 point estimates, and how differently
n6 and n6_abc rank one gene's candidates. Never reads the blinded arms (n6_three, n6_cgi).

Usage (from 2026_1002_sce2g_tf_activity): srun ... python scripts/10.3.eqtl_withingene_diagnostics.py
"""
import importlib.util
import os
import re

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

od = f"{os.environ['SCRATCH']}/eqtl_transfer"
src = open("scripts/10.2.eqtl_withingene_gather.py").read()
ns = {"np": np, "pd": pd}
for name in ["ap_steps", "weighted_ap", "unit_ap"]:
    exec(re.search(rf"^def {name}\(.*?(?=^\S)", src, re.S | re.M).group(0), ns)

mt = pd.read_csv(f"{od}/matches.tsv", sep="\t")
D = pd.concat([pd.read_csv(f"{od}/withingene/{b}.tsv.gz", sep="\t").assign(biosample=b) for b in sorted(set(mt["biosample"]))])
D = D.merge(mt[["biosample", "tissue"]].drop_duplicates(), on=["biosample", "tissue"])
D["unit"] = pd.factorize(D["biosample"] + D["tissue"] + D["TargetGene"])[0]
D["neg_distance"] = -D["distance"].astype(float)
D["random"] = np.random.default_rng(0).random(len(D))
# blinded arms are not read here
arms = ["n6", "n6_abc", "v3_published", "neg_distance", "random"]
for grp, g in [("primary", D[D["biosample"] != "K562"]), ("K562", D[D["biosample"] == "K562"])]:
    print(f"== {grp}: {len(g):,} candidates, positive rate {g['positive'].mean():.4f}")
    for k in arms:
        o, st = ns["ap_steps"](g[k].values.astype(float))
        e1 = ns["weighted_ap"](g["positive"].values[o].astype(float), np.ones(len(g)), st, st[-1] + 1)
        e2 = ns["unit_ap"](g, k).mean()
        print(f"  {k:13s} E1 {e1:.4f}  E2 {e2:.4f}")
    rs = g.groupby("unit").apply(lambda x: spearmanr(x["n6"], x["n6_abc"])[0] if len(x) > 2 else np.nan)
    rd = g.groupby("unit").apply(lambda x: spearmanr(x["n6"], -x["distance"])[0] if len(x) > 2 else np.nan)
    top_same = g.loc[g.groupby("unit")["n6"].idxmax()].set_index("unit")["eid"].eq(g.loc[g.groupby("unit")["n6_abc"].idxmax()].set_index("unit")["eid"]).mean()
    print(f"  within-gene Spearman n6 vs n6_abc median {rs.median():.3f}; n6 vs -distance median {rd.median():.3f}; same top element {top_same:.3f}")
    pos = g[g["positive"] == 1]
    print(f"  positives: median distance {pos['distance'].median():,.0f} bp vs all {g['distance'].median():,.0f}; "
          f"median n6 quantile {pos['n6'].median():.3f} vs {g['n6'].median():.3f}; n6_abc {pos['n6_abc'].median():.3f} vs {g['n6_abc'].median():.3f}")
