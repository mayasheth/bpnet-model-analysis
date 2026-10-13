#!/usr/bin/env python3
"""10.2 Within-gene eQTL readouts per arm, pooled over matched clusters, with the activity calibration (plan 10).

Pre-registered in decisions [2026-10-13] (plan 10). Reads 10.1's candidates.
  E1  pooled AUPRC over every candidate of a group (step-wise average precision, ties as one step)
  E2  mean over (matched pair, eGene) of the within-gene average precision
Intervals: 1,000 paired bootstrap resamples of eGenes (same draws for every arm), as 9.2.

CALIBRATION (read first): n6_abc - n6 interval below 0, and SCALE x |n6_abc - n6| against the test arm's
interval half-width h: >= 1.43 h calibrated (80% power), >= h weak, else not calibrated. SCALE is the
arm's K562 CRISPR gain over the ABC swap's loss (0.021 / 0.050 for n6_three, 0.023 / 0.050 for n6_cgi).
With --blind ARM ..., those arms' differences are written only as half-widths (no estimate, no bounds)
to <results>.calibration.tsv; without it, the full table goes to <results>.readouts.tsv.

Usage: 10.2.eqtl_withingene_gather.py --blind n6_three n6_cgi [--results results/2026_1013_eqtl_withingene]
"""
import argparse
import os

import numpy as np
import pandas as pd

SC = os.environ["SCRATCH"]
SCALE = {"n6_three": 0.021 / 0.050, "n6_cgi": 0.023 / 0.050}
ap = argparse.ArgumentParser()
ap.add_argument("--out-dir", default=f"{SC}/eqtl_transfer")
ap.add_argument("--results", default="results/2026_1013_eqtl_withingene")
ap.add_argument("--n-boot", type=int, default=1000)
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--blind", nargs="*", default=[])
a = ap.parse_args()
od = a.out_dir

mt = pd.read_csv(f"{od}/matches.tsv", sep="\t")
done = sorted(set(mt["biosample"]))
missing = [b for b in done if not os.path.exists(f"{od}/withingene/{b}.tsv.gz")]
assert not missing, f"clusters not scored: {missing}"
D = pd.concat([pd.read_csv(f"{od}/withingene/{b}.tsv.gz", sep="\t").assign(biosample=b) for b in done], ignore_index=True)
D = D.merge(mt[["biosample", "tissue"]].drop_duplicates(), on=["biosample", "tissue"], how="inner")
ARMS = [c for c in D.columns if c not in {"tissue", "eid", "TargetGene", "distance", "positive", "biosample"}]
assert "n6" in ARMS and "n6_abc" in ARMS, ARMS
D["pair"] = D["biosample"] + "|" + D["tissue"]
D["unit"] = pd.factorize(D["pair"] + "|" + D["TargetGene"])[0]
genes = np.array(sorted(set(D["TargetGene"])))
D["g"] = pd.Categorical(D["TargetGene"], categories=genes).codes
print(f"{len(D):,} candidates, {D['positive'].sum():,} positive, {D['unit'].nunique():,} (pair, eGene) units, {len(genes):,} eGenes", flush=True)


def ap_steps(score):
    """Order rows by score (descending) and give each row its tie step, for weighted AP."""
    o = np.argsort(-score, kind="stable")
    s = score[o]
    step = np.r_[0, np.cumsum(s[1:] != s[:-1])]
    return o, step


def weighted_ap(pos_w, all_w, step, nstep):
    tp = np.bincount(step, weights=pos_w, minlength=nstep)
    n = np.bincount(step, weights=all_w, minlength=nstep)
    TP, N = np.cumsum(tp), np.cumsum(n)
    tot = TP[-1]
    if tot <= 0:
        return np.nan
    with np.errstate(divide="ignore", invalid="ignore"):
        prec = np.where(N > 0, TP / N, 0)
    return float(np.sum(tp * prec) / tot)


def unit_ap(g, k):
    """Within-unit AP for arm k, every unit (one matched pair x eGene) at once."""
    x = g.sort_values(["unit", k], ascending=[True, False])
    u, s, p = x["unit"].values, x[k].values, x["positive"].values.astype(float)
    new_unit = np.r_[True, u[1:] != u[:-1]]
    new_step = new_unit | np.r_[True, s[1:] != s[:-1]]
    step = np.cumsum(new_step) - 1
    tp = np.bincount(step, weights=p)
    n = np.bincount(step)
    su = u[new_step]
    # cumulative within unit
    first = np.r_[True, su[1:] != su[:-1]]
    grp = np.cumsum(first) - 1
    TP = np.cumsum(tp) - np.repeat((np.cumsum(tp) - tp)[first], np.bincount(grp))
    N = np.cumsum(n) - np.repeat((np.cumsum(n) - n)[first], np.bincount(grp))
    num = np.bincount(grp, weights=tp * TP / N)
    den = np.bincount(grp, weights=tp)
    return pd.Series(num / den, index=su[first])


def group_estimates(g, W_list):
    """E1 and E2 per arm for each gene-weight vector in W_list."""
    gi = g["g"].values
    pos = g["positive"].values.astype(float)
    ug = g.drop_duplicates("unit").set_index("unit")["g"]
    out = {k: [] for k in ARMS}
    for k in ARMS:
        o, step = ap_steps(g[k].values.astype(float))
        nstep = step[-1] + 1
        uap = unit_ap(g, k)
        ugi = ug.loc[uap.index].values
        for W in W_list:
            w = W[gi][o]
            e1 = weighted_ap(pos[o] * w, w, step, nstep)
            uw = W[ugi]
            e2 = float(np.sum(uap.values * uw) / uw.sum())
            out[k].append((e1, e2))
    return {k: np.array(v) for k, v in out.items()}


pairs_all = sorted(set(D["pair"]))
groups = {"primary (all matches except K562)": [p for p in pairs_all if not p.startswith("K562|")],
          "K562 (training cell)": [p for p in pairs_all if p.startswith("K562|")],
          "PBMC": [p for p in pairs_all if p.startswith("PBMC_")],
          "BMMC22": [p for p in pairs_all if p.startswith("BMMC22_")],
          "GM12878": [p for p in pairs_all if p.startswith("GM12878|")]}
groups = {g: p for g, p in groups.items() if p}

rng = np.random.default_rng(a.seed)
Wb = rng.multinomial(len(genes), np.full(len(genes), 1 / len(genes)), size=a.n_boot).astype(float)
ones = np.ones(len(genes))
rows = []
for gname, pairs in groups.items():
    g = D[D["pair"].isin(pairs)].reset_index(drop=True)
    est = group_estimates(g, [ones])
    boot = group_estimates(g, list(Wb))
    for k in ARMS:
        row = dict(group=gname, arm=k, n_pairs=len(pairs), n_candidates=len(g), n_positive=int(g["positive"].sum()),
                   n_units=g["unit"].nunique(), E1=est[k][0, 0], E2=est[k][0, 1])
        for j, r in enumerate(["E1", "E2"]):
            d = boot[k][:, j] - boot["n6"][:, j]
            d = d[np.isfinite(d)]
            lo, hi = (np.quantile(d, 0.025), np.quantile(d, 0.975)) if len(d) else (np.nan, np.nan)
            row.update({f"{r}_vs_n6": est[k][0, j] - est["n6"][0, j], f"{r}_vs_n6_lo": lo, f"{r}_vs_n6_hi": hi,
                        f"{r}_vs_n6_halfwidth": (hi - lo) / 2})
        rows.append(row)
    print(f"  {gname}: {len(pairs)} pairs, {len(g):,} candidates, {int(g['positive'].sum()):,} positive", flush=True)
R = pd.DataFrame(rows)

# calibration on the primary group
p = R[R["group"] == "primary (all matches except K562)"].set_index("arm")
cal = []
for r in ["E1", "E2"]:
    d_abc, hi_abc = p.loc["n6_abc", f"{r}_vs_n6"], p.loc["n6_abc", f"{r}_vs_n6_hi"]
    for arm, s in SCALE.items():
        if arm not in p.index:
            continue
        h = p.loc[arm, f"{r}_vs_n6_halfwidth"]
        pred = s * abs(d_abc)
        verdict = ("calibrated" if pred >= 1.43 * h else "weak" if pred >= h else "not calibrated") if hi_abc < 0 else "not calibrated (ABC interval not below 0)"
        cal.append(dict(readout=r, test_arm=arm, n6=p.loc["n6", r], abc_vs_n6=d_abc, abc_lo=p.loc["n6_abc", f"{r}_vs_n6_lo"],
                        abc_hi=hi_abc, predicted_gain=pred, test_halfwidth=h, ratio=pred / h, verdict=verdict))
C = pd.DataFrame(cal)
pd.set_option("display.width", 250)
print(C.round(4).to_string(index=False))

if a.blind:
    hide = [c for c in R.columns if c.endswith(("_vs_n6", "_vs_n6_lo", "_vs_n6_hi")) or c in ("E1", "E2")]
    R.loc[R["arm"].isin(a.blind), hide] = np.nan
    R.to_csv(f"{a.results}.calibration_arms.tsv", sep="\t", index=False, float_format="%.5g")
    C.to_csv(f"{a.results}.calibration.tsv", sep="\t", index=False, float_format="%.5g")
    print(f"blinded {a.blind}; wrote {a.results}.calibration.tsv, .calibration_arms.tsv")
else:
    R.to_csv(f"{a.results}.readouts.tsv", sep="\t", index=False, float_format="%.5g")
    C.to_csv(f"{a.results}.calibration.tsv", sep="\t", index=False, float_format="%.5g")
    show = ["group", "arm", "n_units", "E1", "E1_vs_n6", "E1_vs_n6_lo", "E1_vs_n6_hi", "E2", "E2_vs_n6", "E2_vs_n6_lo", "E2_vs_n6_hi"]
    print(R[show].round(4).to_string(index=False))
    for arm in SCALE:
        ok = C[(C["test_arm"] == arm) & (C["verdict"] == "calibrated")]["readout"].tolist()
        if not ok:
            print(f"{arm}: no calibrated readout, no verdict")
            continue
        up = [p.loc[arm, f"{r}_vs_n6_lo"] > 0 for r in ok]
        down = [p.loc[arm, f"{r}_vs_n6_hi"] < 0 for r in ok]
        v = "PASS" if all(up) else ("WORSE" if any(down) else "FAIL")
        print(f"{arm} vs n6, primary, on {ok}: " + "; ".join(
            f"{r} {p.loc[arm, f'{r}_vs_n6']:+.4f} [{p.loc[arm, f'{r}_vs_n6_lo']:+.4f}, {p.loc[arm, f'{r}_vs_n6_hi']:+.4f}]" for r in ok) + f" -> {v}")
    print(f"wrote {a.results}.readouts.tsv, .calibration.tsv")
