#!/usr/bin/env python3
"""8.5 Phase 3 readout (plans/8_crispr_feature_discovery.md): each procedure's real nested gain over the
base (8.3) against its gains when the candidate vectors are shuffled among pairs within TSS-distance bins
(8.4; labels and base untouched, so the null gain is the optimism of searching the library). A procedure
counts only if its real gain beats the 95th percentile of its null gains (empirical p = (1 + #null >= real)
/ (1 + n_null)). The null base AUPRC must equal the real one; checked below.

Usage: 8.5.compare_null.py [--dir $SCRATCH/crispr_discovery] [--out results/2026_1006_nested_null.tsv]
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default=f"{os.environ['SCRATCH']}/crispr_discovery")
ap.add_argument("--out", default="results/2026_1006_nested_null.tsv")
a = ap.parse_args()

real = pd.read_csv(f"{a.dir}/nested_real.auprc.tsv", sep="\t").set_index("model")
files = sorted(glob.glob(f"{a.dir}/nested_fperm*.auprc.tsv"))
null = pd.concat([pd.read_csv(f, sep="\t") for f in files])
nb = null.query('model == "base"')['auprc']
print(f"{len(files)} null runs; null base AUPRC {nb.min():.4f}-{nb.max():.4f} (real {real.loc['base', 'auprc']:.4f})")
assert np.allclose(nb, real.loc["base", "auprc"], atol=1e-6), "the null changed the base model"

rows = []
for m in sorted(set(null["model"]) - {"base"}):
    g = null.loc[null["model"] == m, "vs_base"].values
    r = real.loc[m, "vs_base"]
    rows.append(dict(model=m, real_gain=r, real_lo=real.loc[m, "lo"], real_hi=real.loc[m, "hi"], n_null=len(g),
                     null_mean=g.mean(), null_q95=np.quantile(g, 0.95), null_max=g.max(),
                     p_emp=(1 + (g >= r).sum()) / (1 + len(g)), real_minus_null_mean=r - g.mean()))
R = pd.DataFrame(rows).sort_values("p_emp")
R.to_csv(a.out, sep="\t", index=False, float_format="%.5g")
print(R.round(4).to_string(index=False))
print(f"wrote {a.out}")
