#!/usr/bin/env python
"""
Stratified evaluation of H3K27ac models on held-out elements.

bpnetlite's reported count Pearson is computed over ALL validation elements, and most
of the 150k DNase candidate elements carry almost no H3K27ac. So that number is
dominated by the easy dead-vs-active contrast, which accessibility alone resolves — it
overstates skill on the elements that matter. This script recomputes correlations
stratified by observed signal, so "how well do we predict H3K27ac among elements that
actually have H3K27ac" gets a separate answer.

Predictions are assembled across folds into a genome-wide cross-validated set: each
fold contributes its own held-out chromosomes, which are disjoint, so concatenating is
not leakage.

Extraction and accessibility normalization are imported from train_multimodal_bpnet so
evaluation and training cannot drift apart.

Ceilings come from 0.3.replicate_ceiling_by_window.py with two corrections applied:
Spearman-Brown for the 2-replicate merge, then sqrt because a model predicts the
expected signal rather than a noisy draw of it.

Usage:
  python 2.2.evaluate_stratified.py --models-dir models \
      --elements  <candidate elements narrowPeak> \
      --genome    <hg38.fa> \
      --signal-bw <H3K27ac coverage bw> \
      --accessibility-bw <atac.bw> \
      --fold-json <hg38_five_folds.json> \
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
TRIMMING = 557  # n_layers=8
MODE_LABEL = {"sequence": "Sequence only", "multimodal": "Sequence + ATAC",
              "atac": "ATAC only"}
MODE_COLOR = {"sequence": "#0096a0", "multimodal": "#792374", "atac": "#e96a00"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--models-dir", default="models")
    p.add_argument("--elements", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--signal-bw", required=True, help="H3K27ac coverage BigWig (target)")
    p.add_argument("--accessibility-bw", required=True)
    p.add_argument("--fold-json", required=True)
    p.add_argument("--outdir", default="results")
    p.add_argument("--figdir", default="figures")
    p.add_argument("--ceiling", default="results/replicate_ceiling_by_window.tsv")
    p.add_argument("--batch-size", type=int, default=256)
    p.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    p.add_argument("--configs", nargs="*", default=None,
                   help="Limit to these config dir names (default: all found)")
    return p.parse_args()


def predict_logcounts(model, X, batch_size, device):
    model = model.to(device).eval()
    out = []
    with torch.no_grad():
        for i in range(0, len(X), batch_size):
            xb = torch.from_numpy(X[i:i + batch_size]).to(device)
            _, logcounts = model(xb)
            out.append(logcounts.squeeze(-1).cpu().numpy())
    model.to("cpu")
    return np.concatenate(out)


def eval_one_fold(cfg_dir, mode, hw, fold, args):
    """Return (observed_log, predicted_log) for one fold's held-out chromosomes."""
    model_path = os.path.join(cfg_dir, f"fold{fold}", "multimodal_bpnet.torch")
    if not os.path.exists(model_path):
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
        els, genome, args.signal_bw, None, acc_bw,
        in_window, out_window, 0, is_peak=True)

    # observed target: what the counts head was trained against
    observed = np.log1p(sigs[:, 0, :].sum(axis=1))

    if acc_bw is not None:
        stats_path = os.path.join(cfg_dir, f"fold{fold}", "acc_normalization_stats.json")
        if not os.path.exists(stats_path):
            raise SystemExit(f"error: missing {stats_path}; needed to match training "
                             "accessibility normalization")
        with open(stats_path) as f:
            st = json.load(f)
        accs = normalize_accessibility(accs, mean=st["acc_mean"], std=st["acc_std"])[0]

    if mode == "multimodal":
        X = np.concatenate([seqs, accs], axis=1)
    elif mode == "sequence":
        X = seqs
    else:
        X = accs

    model = torch.load(model_path, map_location="cpu", weights_only=False)
    if not hasattr(model, "mode"):
        model.mode = mode   # checkpoints pickled before `mode` was an attribute
    pred = predict_logcounts(model, X.astype(np.float32), args.batch_size, args.device)

    del seqs, sigs, accs, X
    gc.collect()
    return observed, pred


def strat_rows(obs, pred, label_prefix):
    """Correlations overall and by observed-signal quintile."""
    rows = []
    q = pd.qcut(obs.argsort().argsort(), 5, labels=False)  # rank-based, tie-safe
    subsets = [("all", np.ones(len(obs), bool))]
    subsets += [(f"quintile{i+1}", q == i) for i in range(5)]
    subsets += [("top_quintile", q == 4), ("top_two_quintiles", q >= 3)]
    for name, mask in subsets:
        if mask.sum() < 10:
            continue
        o, p = obs[mask], pred[mask]
        if np.std(o) < 1e-9 or np.std(p) < 1e-9:
            continue
        rows.append({**label_prefix, "stratum": name, "n": int(mask.sum()),
                     "pearson": round(float(pearsonr(o, p)[0]), 4),
                     "spearman": round(float(spearmanr(o, p)[0]), 4)})
    return rows


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.figdir, exist_ok=True)

    configs = sorted(d for d in os.listdir(args.models_dir) if DIR_RE.match(d))
    if args.configs:
        configs = [c for c in configs if c in args.configs]
    if not configs:
        raise SystemExit("no matching config directories found")

    all_rows = []
    for cfg in configs:
        m = DIR_RE.match(cfg)
        mode, hw = m["mode"], int(m["hw"])
        cfg_dir = os.path.join(args.models_dir, cfg)
        obs_all, pred_all, folds_used = [], [], []
        for fold in range(5):
            r = eval_one_fold(cfg_dir, mode, hw, fold, args)
            if r is None:
                continue
            obs_all.append(r[0])
            pred_all.append(r[1])
            folds_used.append(fold)
            print(f"  {cfg} fold{fold}: n={len(r[0])}")
        if not obs_all:
            print(f"{cfg}: no folds available yet, skipping")
            continue
        obs = np.concatenate(obs_all)
        pred = np.concatenate(pred_all)
        print(f"{cfg}: {len(folds_used)} folds, {len(obs):,} elements pooled")
        all_rows += strat_rows(obs, pred, {
            "config": cfg, "mode": mode, "half_window": hw,
            "n_folds": len(folds_used)})

    res = pd.DataFrame(all_rows)
    p1 = os.path.join(args.outdir, "stratified_evaluation.tsv")
    res.to_csv(p1, sep="\t", index=False)
    print(f"\nWrote {p1}")

    # ceilings, corrected the same way as 2.1
    ceil_all, ceil_top = {}, {}
    if os.path.exists(args.ceiling):
        c = pd.read_csv(args.ceiling, sep="\t")
        for _, r in c.iterrows():
            for col, dest in (("pearson_all", ceil_all),
                              ("pearson_top_quintile", ceil_top)):
                rel = 2 * r[col] / (1 + r[col])
                dest[int(r["half_window"])] = rel ** 0.5

    show = res[res["stratum"].isin(["all", "top_quintile"])].copy()
    show["ceiling"] = [
        (ceil_all if s == "all" else ceil_top).get(hw, np.nan)
        for s, hw in zip(show["stratum"], show["half_window"])]
    show["frac_of_ceiling"] = (show["pearson"] / show["ceiling"]).round(3)
    print("\n" + show.to_string(index=False))

    # --- Figure: all vs top quintile, per config ---------------------------
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8), sharey=True)
    for ax, stratum in zip(axes, ("all", "top_quintile")):
        sub = show[show["stratum"] == stratum].sort_values(["half_window", "mode"])
        x = np.arange(len(sub))
        ax.bar(x, sub["pearson"], color=[MODE_COLOR[m] for m in sub["mode"]], width=0.62)
        for hw in sorted(sub["half_window"].unique()):
            cd = (ceil_all if stratum == "all" else ceil_top).get(hw)
            if cd:
                ax.axhline(cd, color="black", lw=0.8, ls="--")
                ax.text(0, cd, f" ceiling ±{hw}", fontsize=6, va="bottom")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{MODE_LABEL[m]}\n±{h}"
                            for m, h in zip(sub["mode"], sub["half_window"])],
                           fontsize=6.5)
        ax.set_title("All elements" if stratum == "all"
                     else "Top signal quintile", fontsize=9)
        ax.grid(False)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color("black")
        ax.tick_params(colors="black")
    axes[0].set_ylabel("Pearson r (log counts)")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        p = os.path.join(args.figdir, f"h3k27ac_stratified_eval.{ext}")
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(f"Wrote {p}")


if __name__ == "__main__":
    main()
