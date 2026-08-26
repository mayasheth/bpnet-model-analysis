#!/usr/bin/env python
"""
Evaluate what a model adds BEYOND accessibility.

Overall correlation is the wrong headline for this project. The deliverable is a model
that takes ATAC + sequence and predicts activity, and ATAC alone already reaches 0.543 on
active elements — so a joint model scoring 0.668 looks strong while saying little about
what the sequence contributed. The quantity of interest is the DEPARTURE from what
accessibility alone would predict.

Baseline is the ATAC-only MODEL's held-out prediction, not the raw ATAC track, so the
residual is what accessibility genuinely cannot explain rather than what a linear read of
the track cannot explain.

Metrics, per model:
  residual_pearson   r(observed - atac_pred, model_pred - atac_pred).
                     Does the model's departure from the ATAC expectation track the true
                     departure? This is the headline number.
  incremental_r2     R^2(model) - R^2(atac_only): variance explained beyond ATAC.
  overall_pearson    r(observed, model_pred), for continuity with earlier reporting.
  Stratified by |true residual| quintile — elements whose H3K27ac most departs from their
  ATAC expectation are the biologically interesting ones (accessible-but-unacetylated, and
  acetylated-beyond-accessibility), and are where a useful model must actually work.

Also writes the extreme-residual elements to TSV as a candidate list.

Usage:
  python 2.4.evaluate_residual.py --models-dir models --atac-config atac_hw500_clw1000 \
      --configs multimodal_hw500_clw1000 sequence_hw500_clw1000 \
      --elements <narrowPeak> --genome <hg38.fa> --signal-bw <target.bw> \
      --accessibility-bw <atac.bw> --fold-json <folds.json> \
      --outdir results --figdir figures
"""

import argparse
import gc
import json
import os
import re
import sys

import numpy as np
import pandas as pd
import torch
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_SCRIPTS = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/scripts"
sys.path.insert(0, REPO_SCRIPTS)
from train_multimodal_bpnet import (  # noqa: E402
    extract_windows, load_peaks, normalize_accessibility)

DIR_RE = re.compile(r"^(?P<mode>sequence|multimodal|atac)_hw(?P<hw>\d+)_clw(?P<clw>\d+)$")
TRIMMING = 557
COLOR = {"sequence": "#0096a0", "multimodal": "#792374", "atac": "#e96a00"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--models-dir", default="models")
    p.add_argument("--atac-config", required=True,
                   help="Config dir of the ATAC-only model that defines the baseline")
    p.add_argument("--configs", nargs="+", required=True,
                   help="Config dirs to evaluate against that baseline")
    p.add_argument("--elements", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--signal-bw", required=True)
    p.add_argument("--accessibility-bw", required=True)
    p.add_argument("--fold-json", required=True)
    p.add_argument("--outdir", default="results")
    p.add_argument("--figdir", default="figures")
    p.add_argument("--batch-size", type=int, default=256)
    p.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    return p.parse_args()


def predict_logcounts(model, X, batch_size, device):
    model = model.to(device).eval()
    out = []
    with torch.no_grad():
        for i in range(0, len(X), batch_size):
            _, lc = model(torch.from_numpy(X[i:i + batch_size]).to(device))
            out.append(lc.squeeze(-1).cpu().numpy())
    model.to("cpu")
    return np.concatenate(out)


def run_fold(cfg_dir, mode, hw, fold, args, want_observed):
    """Predicted logcounts for one fold (and observed, if asked)."""
    mp = os.path.join(cfg_dir, f"fold{fold}", "multimodal_bpnet.torch")
    if not os.path.exists(mp):
        return None
    out_window = 2 * hw
    in_window = out_window + 2 * TRIMMING
    with open(args.fold_json) as f:
        val_chroms = json.load(f)[str(fold)]["val"]
    els = load_peaks(args.elements, val_chroms)
    if len(els) == 0:
        return None

    acc_bw = args.accessibility_bw if mode in ("multimodal", "atac") else None
    genome = args.genome if mode != "atac" else None
    seqs, sigs, accs, _ = extract_windows(
        els, genome, args.signal_bw, None, acc_bw, in_window, out_window, 0,
        is_peak=True)
    observed = np.log1p(sigs.sum(axis=(1, 2))) if want_observed else None

    if acc_bw is not None:
        sp = os.path.join(cfg_dir, f"fold{fold}", "acc_normalization_stats.json")
        with open(sp) as f:
            st = json.load(f)
        accs = normalize_accessibility(accs, mean=st["acc_mean"], std=st["acc_std"])[0]

    X = (np.concatenate([seqs, accs], axis=1) if mode == "multimodal"
         else seqs if mode == "sequence" else accs)
    model = torch.load(mp, map_location="cpu", weights_only=False)
    if not hasattr(model, "mode"):
        model.mode = mode
    pred = predict_logcounts(model, X.astype(np.float32), args.batch_size, args.device)
    coords = (els["chr"].to_numpy(), els["center"].to_numpy()
              if "center" in els else (els["start"] + els["summit"]).to_numpy())
    del seqs, sigs, accs, X
    gc.collect()
    return pred, observed, coords


def r2(obs, pred):
    """Variance of obs explained by pred, allowing an affine rescaling."""
    return float(pearsonr(obs, pred)[0] ** 2)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.figdir, exist_ok=True)

    m = DIR_RE.match(args.atac_config)
    if not m or m["mode"] != "atac":
        raise SystemExit(f"--atac-config must be an atac_* config dir, got "
                         f"{args.atac_config}")
    hw = int(m["hw"])

    # baseline: ATAC-only predictions and the observed target, pooled over folds
    base_pred, observed, chrom_l, ctr_l, folds = [], [], [], [], []
    for fold in range(5):
        r = run_fold(os.path.join(args.models_dir, args.atac_config), "atac", hw,
                     fold, args, want_observed=True)
        if r is None:
            continue
        base_pred.append(r[0]); observed.append(r[1])
        chrom_l.append(r[2][0]); ctr_l.append(r[2][1])
        folds.append(fold)
        print(f"  baseline {args.atac_config} fold{fold}: n={len(r[0])}")
    if not base_pred:
        raise SystemExit("no ATAC-only folds found")
    atac_pred = np.concatenate(base_pred)
    obs = np.concatenate(observed)
    chroms = np.concatenate(chrom_l)
    centers = np.concatenate(ctr_l)
    true_resid = obs - atac_pred
    print(f"baseline pooled: {len(obs):,} elements over folds {folds}")

    # |true residual| quintiles: where accessibility is most wrong
    aq = pd.qcut(np.abs(true_resid).argsort().argsort(), 5, labels=False)

    rows = []
    preds = {"atac": atac_pred}
    for cfg in args.configs:
        mm = DIR_RE.match(cfg)
        if not mm:
            print(f"skip unrecognized config {cfg}")
            continue
        if int(mm["hw"]) != hw:
            raise SystemExit(f"{cfg} has hw={mm['hw']} but baseline hw={hw}; "
                             "windows must match for the residual to be meaningful")
        ps = []
        for fold in folds:
            r = run_fold(os.path.join(args.models_dir, cfg), mm["mode"], hw, fold,
                         args, want_observed=False)
            if r is None:
                raise SystemExit(f"{cfg} missing fold{fold}, which the baseline has")
            ps.append(r[0])
            print(f"  {cfg} fold{fold}: n={len(r[0])}")
        pred = np.concatenate(ps)
        preds[mm["mode"]] = pred
        model_resid = pred - atac_pred

        base_r2 = r2(obs, atac_pred)
        for name, mask in ([("all", np.ones(len(obs), bool))]
                           + [(f"abs_resid_q{i+1}", aq == i) for i in range(5)]):
            if mask.sum() < 10:
                continue
            o, tr, mr, pr = obs[mask], true_resid[mask], model_resid[mask], pred[mask]
            row = {"config": cfg, "mode": mm["mode"], "stratum": name,
                   "n": int(mask.sum()),
                   "overall_pearson": round(float(pearsonr(o, pr)[0]), 4)}
            if np.std(tr) > 1e-9 and np.std(mr) > 1e-9:
                row["residual_pearson"] = round(float(pearsonr(tr, mr)[0]), 4)
                row["residual_spearman"] = round(float(spearmanr(tr, mr)[0]), 4)
            else:
                row["residual_pearson"] = np.nan
                row["residual_spearman"] = np.nan
            row["incremental_r2"] = round(r2(o, pr) - r2(o, atac_pred[mask]), 4)
            rows.append(row)
        print(f"{cfg}: residual_pearson(all)="
              f"{[r for r in rows if r['config']==cfg and r['stratum']=='all'][0]['residual_pearson']}")

    res = pd.DataFrame(rows)
    p1 = os.path.join(args.outdir, "residual_evaluation.tsv")
    res.to_csv(p1, sep="\t", index=False)
    print(f"\nWrote {p1}")
    print(res.to_string(index=False))

    # extreme-residual elements as a candidate list
    order = np.argsort(true_resid)
    k = 2000
    ex = pd.DataFrame({
        "chr": np.concatenate([chroms[order[:k]], chroms[order[-k:]]]),
        "center": np.concatenate([centers[order[:k]], centers[order[-k:]]]),
        "observed_log": np.concatenate([obs[order[:k]], obs[order[-k:]]]),
        "atac_pred_log": np.concatenate([atac_pred[order[:k]], atac_pred[order[-k:]]]),
        "true_residual": np.concatenate([true_resid[order[:k]], true_resid[order[-k:]]]),
        "class": ["accessible_not_acetylated"] * k + ["acetylated_beyond_accessibility"] * k,
    })
    for mode, pr in preds.items():
        if mode == "atac":
            continue
        ex[f"{mode}_residual"] = np.concatenate([(pr - atac_pred)[order[:k]],
                                                (pr - atac_pred)[order[-k:]]])
    p2 = os.path.join(args.outdir, "extreme_residual_elements.tsv")
    ex.to_csv(p2, sep="\t", index=False)
    print(f"Wrote {p2}")

    # --- Figure -----------------------------------------------------------
    sub = res[res["stratum"] == "all"]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
    ax = axes[0]
    x = np.arange(len(sub))
    ax.bar(x - 0.2, sub["overall_pearson"], width=0.38, color="#bbbbbb",
           label="Overall r (observed vs pred)")
    ax.bar(x + 0.2, sub["residual_pearson"], width=0.38,
           color=[COLOR[m] for m in sub["mode"]], label="Residual r (beyond ATAC)")
    ax.set_xticks(x)
    ax.set_xticklabels(sub["mode"], fontsize=8)
    ax.set_ylabel("Pearson r")
    ax.set_title("Overall correlation overstates what sequence adds", fontsize=9)
    ax.legend(frameon=False, fontsize=7)

    ax = axes[1]
    for mode in sub["mode"]:
        s = res[(res["mode"] == mode) & res["stratum"].str.startswith("abs_resid_q")]
        ax.plot(range(1, len(s) + 1), s["residual_pearson"], marker="o", ms=4,
                color=COLOR[mode], label=mode)
    ax.set_xlabel("|true residual| quintile (5 = ATAC most wrong)")
    ax.set_ylabel("Residual Pearson r")
    ax.set_title("Where accessibility fails, does sequence help?", fontsize=9)
    ax.legend(frameon=False, fontsize=7)

    for ax in axes:
        ax.grid(False)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
        for s_ in ("left", "bottom"):
            ax.spines[s_].set_color("black")
        ax.tick_params(colors="black")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        p = os.path.join(args.figdir, f"residual_evaluation.{ext}")
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(f"Wrote {p}")


if __name__ == "__main__":
    main()
