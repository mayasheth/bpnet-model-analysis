#!/usr/bin/env python
"""Score the ATAC->DNase converter's PROFILE against the inter-replicate ceiling.

THE METRIC THIS PROJECT USUALLY REPORTS WOULD BE THE WRONG ONE HERE. 2.15's
`profile_pearson` comes from bpnetlite's `calculate_performance_measures` with
`kernel_sigma=7, kernel_width=81`, i.e. both profiles are Gaussian-SMOOTHED before
correlating. The DNase ceiling that unblocked this whole line of work (0.848 top quintile,
K562) is a RAW within-element correlation at 1 bp from 0.25. Reading a smoothed model score
against a raw ceiling would say the converter is closer to the ceiling than it is, and it
is precisely the class of apples-to-oranges comparison report 2 catalogues. So this script
scores the model with 0.25's `shape_corr`, unchanged, and computes the ceiling on the SAME
held-out windows in the same run -- no element-set or estimator difference is left to argue
about.

WHY THE BIN-SIZE CURVE AND NOT ONE NUMBER. The converter exists to feed the downstream
model, whose accessibility branch reads a base-resolution track over 2,114 bp. The question
is not "is the profile good" but "at what resolution is it trustworthy", because that is
what decides whether the painted track carries the structure the branch uses or is a smooth
bump that a flat fill would have matched. A converter that only works at 50 bp is a
different, weaker claim than one that works at 1 bp, and the curve says which we have.

SHAPE ONLY, COUNTS SEPARATELY. `shape_corr` centres each element before correlating, so the
element's total is divided out and any positive rescaling of the prediction is invisible.
That is deliberate: the counts head is scored by 2.15 and is not what is at issue.

CEILING CAVEAT, STATED BECAUSE IT RUNS THE WRONG WAY. The ceiling here comes from the two
ENCSR000EOT replicates (142M reads between them), Spearman-Brown corrected to the 2-replicate
pool. The converter's target is 0.34's 3-BAM pool at 301M reads, which is less noisy and so
has a HIGHER true ceiling than the number printed here. The gap to the ceiling reported by
this script is therefore an UNDERSTATEMENT of the gap to the reachable ceiling -- it flatters
the model, so do not read "close to the ceiling" as "nearly saturated".

Usage:
  2.31.converter_profile_eval.py --model-dir DIR --fold N [--fold N ...] \
      --label NAME [--cell k562|gm12878] [--mode multimodal|atac|sequence]
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import pyBigWig
import torch

R = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/scripts"
D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
P = f"{D}/2026_0824_H3K27ac_model"
sys.path.insert(0, R)
from train_multimodal_bpnet import extract_windows, load_peaks, normalize_accessibility

GEN = f"/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa"
FOLDS = f"{D}/reference/hg38_five_folds.json"
HW = 500
OUT_W = 2 * HW
BINS = [1, 5, 10, 25, 50, 100, 250]
dev = "cuda" if torch.cuda.is_available() else "cpu"

# Per cell type: the ATAC input the converter reads, the pooled stranded DNase it predicts,
# the two per-replicate stranded tracks the ceiling comes from (0.33), and the element set.
CELLS = {
    "k562": dict(
        acc=f"{D}/2026_0529_multimodal_p300_model/data/atac_5p.bw",
        tgt_plus=f"{P}/data/k562_dnase_5p_plus.bw",
        tgt_minus=f"{P}/data/k562_dnase_5p_minus.bw",
        rep=[(f"{P}/data/k562_dnase_rep1_5p_plus.bw", f"{P}/data/k562_dnase_rep1_5p_minus.bw"),
             (f"{P}/data/k562_dnase_rep2_5p_plus.bw", f"{P}/data/k562_dnase_rep2_5p_minus.bw")],
        elements=f"{D}/reference/K562_ATAC_candidate_elements.narrowPeak"),
    "gm12878": dict(
        acc=f"{D}/2026_0606_GM12878_transferability/data/atac_5p.bw",
        tgt_plus=f"{P}/data/gm12878_dnase_5p_plus.bw",
        tgt_minus=f"{P}/data/gm12878_dnase_5p_minus.bw",
        rep=[(f"{P}/data/gm12878_dnase_rep1_5p_plus.bw", f"{P}/data/gm12878_dnase_rep1_5p_minus.bw"),
             (f"{P}/data/gm12878_dnase_rep2_5p_plus.bw", f"{P}/data/gm12878_dnase_rep2_5p_minus.bw")],
        elements=f"{D}/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak"),
}

ap = argparse.ArgumentParser()
ap.add_argument("--model-dir", required=True)
ap.add_argument("--fold", type=int, action="append", required=True)
ap.add_argument("--label", required=True)
ap.add_argument("--cell", default="k562", choices=sorted(CELLS))
ap.add_argument("--mode", default="multimodal", choices=["multimodal", "atac", "sequence"])
ap.add_argument("--out", default=None)
a = ap.parse_args()
C = CELLS[a.cell]


def shape_corr(x, y, b):
    """Mean within-element Pearson r between two binned profiles.

    Copied verbatim from 0.25 rather than imported, because 0.25 is a script with
    argparse at import time. If one changes, the other must -- that is the point of the
    comparison, so any edit here without the matching edit there invalidates the ceiling
    column beside it.
    """
    n, L = x.shape
    xb = x.reshape(n, L // b, b).sum(2)
    yb = y.reshape(n, L // b, b).sum(2)
    if xb.shape[1] < 3:
        return np.nan, 0
    xc = xb - xb.mean(1, keepdims=True)
    yc = yb - yb.mean(1, keepdims=True)
    num = (xc * yc).sum(1)
    den = np.sqrt((xc ** 2).sum(1) * (yc ** 2).sum(1))
    good = den > 0
    if good.sum() == 0:
        return np.nan, 0
    return float((num[good] / den[good]).mean()), int(good.sum())


def pooled_ceiling(r):
    """Spearman-Brown: replicate-vs-replicate agreement -> agreement with the 2-rep pool."""
    if not np.isfinite(r) or r <= 0:
        return np.nan
    return float(np.sqrt(2 * r / (1 + r)))


def rep_windows(els, paths):
    """Per-replicate stranded windows on exactly the element rows the model was scored on."""
    bws = [pyBigWig.open(p) for p in paths]
    sizes = bws[0].chroms()
    out = [[] for _ in paths]
    keep = []
    for i, r in els.iterrows():
        c = r["chr"]
        ctr = int(r["start"]) + int(r["summit"])
        s, e = ctr - HW, ctr + HW
        if c not in sizes or s < 0 or e > sizes[c]:
            continue
        vals, ok = [], True
        for bw in bws:
            x = bw.values(c, s, e, numpy=True)
            if x is None or len(x) != OUT_W:
                ok = False
                break
            vals.append(np.nan_to_num(x, nan=0.0))
        if ok:
            keep.append(i)
            for j, v in enumerate(vals):
                out[j].append(v)
    for bw in bws:
        bw.close()
    return [np.stack(o) if o else np.zeros((0, OUT_W)) for o in out], keep


rows = []
folds_json = json.load(open(FOLDS))
for fold in a.fold:
    mp = f"{a.model_dir}/fold{fold}/multimodal_bpnet.torch"
    marker = f"{a.model_dir}/fold{fold}/training_complete.json"
    if not os.path.exists(marker):
        raise SystemExit(f"error: {marker} missing; refusing to score an unfinished fold.")
    m = torch.load(mp, map_location="cpu", weights_only=False)
    if not hasattr(m, "mode"):
        m.mode = a.mode
    in_w = OUT_W + 2 * m.trimming

    els = load_peaks(C["elements"], folds_json[str(fold)]["val"])
    seqs, sigs, accs, valid = extract_windows(
        els, GEN if a.mode != "atac" else None, C["tgt_plus"], C["tgt_minus"],
        C["acc"] if a.mode != "sequence" else None, in_w, OUT_W, 0, is_peak=True)
    els_kept = els[valid].reset_index(drop=True)

    x = accs
    if a.mode in ("multimodal", "atac"):
        st = json.load(open(f"{a.model_dir}/fold{fold}/acc_normalization_stats.json"))
        x = normalize_accessibility(accs, mean=st["acc_mean"], std=st["acc_std"])[0]
    X = (np.concatenate([seqs, x], axis=1) if a.mode == "multimodal"
         else seqs if a.mode == "sequence" else x).astype(np.float32)

    m = m.to(dev).eval()
    profs = []
    with torch.no_grad():
        for i in range(0, len(X), 64):
            xb = torch.from_numpy(X[i:i + 64]).to(dev)
            pr, _ = m(xb)
            sh = pr.shape
            # Softmax over the flattened (strand, position) axes, matching how the profile
            # head is trained. Shape only, so the per-element scale is irrelevant here.
            p = torch.softmax(pr.reshape(sh[0], -1), dim=-1).reshape(sh)
            profs.append(p.cpu().numpy())
    m.to("cpu")
    pred = np.concatenate(profs)
    del X, profs, seqs, accs, x

    # Observed pooled track, and the model prediction, both collapsed over strand for the
    # unstranded comparison and kept separate for the per-strand one.
    obs_p, obs_m = sigs[:, 0, :], sigs[:, 1, :]
    prd_p, prd_m = pred[:, 0, :], pred[:, 1, :]
    tot = obs_p.sum(1) + obs_m.sum(1)
    q = np.quantile(tot, [0.2, 0.4, 0.6, 0.8])
    quint = np.digitize(tot, q)

    # Ceiling on the SAME windows. rep_windows re-filters (a replicate track can be short
    # at a contig edge where the pooled one is not), so its rows are a subset and the
    # stratum mask has to be subset the same way before the two are compared.
    reps, keep = rep_windows(els_kept, [p for pair in C["rep"] for p in pair])
    r1p, r1m, r2p, r2m = reps
    kmask = np.zeros(len(els_kept), bool)
    kmask[keep] = True

    for b in BINS:
        for strat, mask in (("all", np.ones(len(tot), bool)), ("topq", quint == 4)):
            r_model, n_model = shape_corr((obs_p + obs_m)[mask], (prd_p + prd_m)[mask], b)
            mk = mask[kmask]
            if mk.sum() > 0:
                r_rep, n_rep = shape_corr((r1p + r1m)[mk], (r2p + r2m)[mk], b)
            else:
                r_rep, n_rep = np.nan, 0
            ceil = pooled_ceiling(r_rep)
            frac = r_model / ceil if np.isfinite(ceil) and ceil > 0 else np.nan
            rows.append({
                "label": a.label, "cell": a.cell, "fold": fold, "bin_bp": b,
                "stratum": strat, "n_elements": n_model,
                "shape_r_model": r_model,
                "shape_r_replicate": r_rep, "n_ceiling_elements": n_rep,
                "ceiling": ceil, "frac_of_ceiling": frac})
    print(f"fold{fold}: {len(tot):,} windows, ceiling on {kmask.sum():,}", flush=True)
    del sigs, pred

df = pd.DataFrame(rows)
out = a.out or f"{P}/results/converter_profile_{a.label}_{a.cell}.tsv"
df.round(4).to_csv(out, sep="\t", index=False)
print("\nWrote", out)

g = (df.groupby(["bin_bp", "stratum"])[["shape_r_model", "ceiling", "frac_of_ceiling"]]
       .mean().reset_index())
print(f"\n{a.label} on {a.cell}, mean over {len(a.fold)} fold(s)")
print(f"{'bin':>6}{'stratum':>9}{'model':>10}{'ceiling':>10}{'% of ceiling':>14}")
for _, r in g.iterrows():
    pc = 100 * r["frac_of_ceiling"] if np.isfinite(r["frac_of_ceiling"]) else np.nan
    print(f"{int(r['bin_bp']):>6}{r['stratum']:>9}{r['shape_r_model']:>10.4f}"
          f"{r['ceiling']:>10.4f}{pc:>13.1f}%")
print("\nThe 1 bp top-quintile row is the one that decides whether the painted track carries")
print("base-resolution structure. Ceiling here is from 2 replicates (142M reads) while the")
print("target is the 3-BAM pool (301M), so the true ceiling is higher and % of ceiling is")
print("optimistic -- see the docstring.")
