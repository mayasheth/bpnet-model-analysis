#!/usr/bin/env python
"""
Aggregate H3K27ac training runs into one comparison table.

Walks models/{mode}_hw{W}_clw{C}/fold{N}/multimodal_bpnet.log and summarises validation
count correlation per config, against the inter-replicate ceiling for that window.

Caveat on the numbers: bpnetlite computes count Pearson over ALL validation elements,
and most of the 150k DNase candidate elements carry little H3K27ac. So a large part of
the correlation is the easy dead-vs-active contrast, which accessibility alone resolves.
Treat these as relative comparisons between configs, not as absolute skill; a
stratified (top-quintile) evaluation is needed for that.

Reports two numbers per run, because they answer different questions:
  checkpoint_*  metric at the last epoch bpnetlite actually saved (best validation
                loss). This is the model you would use, so it is the honest number.
  best_*        maximum over all epochs. Optimistic — it selects on the metric being
                reported — and shown only to reveal how much was left on the table.

Safe to run mid-grid; incomplete configs are reported with the fold count found.

Usage:
  python 2.1.aggregate_results.py --models-dir models --outdir results --figdir figures \
      [--ceiling results/replicate_ceiling_by_window.tsv]
"""

import argparse
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR_RE = re.compile(r"^(?P<mode>sequence|multimodal|atac)_hw(?P<hw>\d+)_clw(?P<clw>\d+)$")
MODE_ORDER = ["sequence", "multimodal", "atac"]
MODE_LABEL = {"sequence": "Sequence only", "multimodal": "Sequence + ATAC",
              "atac": "ATAC only"}
MODE_COLOR = {"sequence": "#0096a0", "multimodal": "#792374", "atac": "#e96a00"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--models-dir", default="models")
    p.add_argument("--outdir", default="results")
    p.add_argument("--figdir", default="figures")
    p.add_argument("--ceiling", default="results/replicate_ceiling_by_window.tsv",
                   help="Output of 0.3.replicate_ceiling_by_window.py; skipped if absent")
    return p.parse_args()


def read_run(path):
    """Parse one epoch log; return checkpoint and best metrics."""
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    cp_col = "Validation Count Pearson"
    pp_col = "Validation Profile Pearson"
    saved = df[df["Saved?"].astype(str).str.strip().str.lower() == "true"]
    if saved.empty:
        return None
    last = saved.iloc[-1]
    return {
        "n_epochs": len(df),
        "checkpoint_epoch": int(last["Epoch"]),
        "checkpoint_count_pearson": float(last[cp_col]),
        "checkpoint_profile_pearson": float(last[pp_col]),
        "best_count_pearson": float(df[cp_col].max()),
        "final_val_mnll": float(last["Validation MNLL"]),
    }


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.figdir, exist_ok=True)

    rows = []
    if not os.path.isdir(args.models_dir):
        raise SystemExit(f"error: no such directory: {args.models_dir}")
    for cfg in sorted(os.listdir(args.models_dir)):
        m = DIR_RE.match(cfg)
        if not m:
            continue
        for fold_dir in sorted(os.listdir(os.path.join(args.models_dir, cfg))):
            log = os.path.join(args.models_dir, cfg, fold_dir, "multimodal_bpnet.log")
            if not os.path.exists(log):
                continue
            rec = read_run(log)
            if rec is None:
                continue
            rec.update(mode=m["mode"], half_window=int(m["hw"]),
                       count_loss_weight=int(m["clw"]),
                       fold=int(fold_dir.replace("fold", "")))
            rows.append(rec)

    if not rows:
        raise SystemExit("No completed runs found yet.")

    per_fold = pd.DataFrame(rows).sort_values(
        ["half_window", "mode", "fold"]).reset_index(drop=True)
    cols = ["mode", "half_window", "count_loss_weight", "fold", "n_epochs",
            "checkpoint_epoch", "checkpoint_count_pearson",
            "checkpoint_profile_pearson", "best_count_pearson", "final_val_mnll"]
    per_fold = per_fold[cols]
    p1 = os.path.join(args.outdir, "training_results_per_fold.tsv")
    per_fold.to_csv(p1, sep="\t", index=False, float_format="%.4f")
    print(f"Wrote {p1}  ({len(per_fold)} runs)")

    summ = (per_fold.groupby(["mode", "half_window", "count_loss_weight"])
            .agg(n_folds=("fold", "count"),
                 count_pearson_mean=("checkpoint_count_pearson", "mean"),
                 count_pearson_sd=("checkpoint_count_pearson", "std"),
                 best_count_pearson_mean=("best_count_pearson", "mean"),
                 profile_pearson_mean=("checkpoint_profile_pearson", "mean"))
            .reset_index())

    ceiling = None
    if os.path.exists(args.ceiling):
        c = pd.read_csv(args.ceiling, sep="\t")
        # 0.3 reports r(rep1, rep2) — the reliability of a SINGLE replicate. The
        # models predict the MERGED track, which pools both replicates and is less
        # noisy, so single-replicate r is not the relevant bound. Spearman-Brown
        # converts it to the reliability of the 2-replicate average:
        #     r_merged = 2r / (1 + r)
        # Models can legitimately exceed the single-replicate number; exceeding this
        # corrected one is a signal to go looking for a problem.
        raw = dict(zip(c["half_window"], c["pearson_all"]))
        # Two corrections to r(rep1, rep2) before it is a bound on model performance:
        #   1. Spearman-Brown: the target is the 2-replicate MERGE, which is less noisy
        #      than one replicate.        rel = 2r / (1 + r)
        #   2. A model predicts the expected signal, not a noisy draw of it. The maximum
        #      correlation between a perfect predictor of the true value and a noisy
        #      observation of it is sqrt(reliability), not reliability.
        # Without (2) the "ceiling" is exceeded by any decent model, which is what
        # happened here before this was fixed.
        rel = {hw: 2 * r / (1 + r) for hw, r in raw.items()}
        ceiling = {hw: v ** 0.5 for hw, v in rel.items()}
        summ["r_rep1_rep2"] = summ["half_window"].map(raw)
        summ["reliability_merged"] = summ["half_window"].map(rel)
        summ["ceiling"] = summ["half_window"].map(ceiling)
        summ["frac_of_ceiling"] = summ["count_pearson_mean"] / summ["ceiling"]
    else:
        print(f"note: {args.ceiling} not found; skipping ceiling columns")

    summ["mode"] = pd.Categorical(summ["mode"], MODE_ORDER, ordered=True)
    summ = summ.sort_values(["half_window", "mode"]).reset_index(drop=True)
    p2 = os.path.join(args.outdir, "training_results_summary.tsv")
    summ.to_csv(p2, sep="\t", index=False, float_format="%.4f")
    print(f"Wrote {p2}\n")
    print(summ.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    incomplete = summ[summ["n_folds"] < 5]
    if len(incomplete):
        print(f"\nNOTE: {len(incomplete)} config(s) have fewer than 5 folds — "
              "grid still running, numbers will move.")

    # --- Figure -----------------------------------------------------------
    windows = sorted(summ["half_window"].unique())
    fig, axes = plt.subplots(1, len(windows), figsize=(4.2 * len(windows), 3.8),
                             sharey=True, squeeze=False)
    for ax, hw in zip(axes[0], windows):
        sub = summ[summ["half_window"] == hw]
        x = np.arange(len(sub))
        ax.bar(x, sub["count_pearson_mean"],
               yerr=sub["count_pearson_sd"].fillna(0), capsize=3,
               color=[MODE_COLOR[m] for m in sub["mode"]], width=0.62)
        if ceiling and hw in ceiling:
            ax.axhline(ceiling[hw], color="black", lw=0.9, ls="--")
            ax.text(len(sub) - 0.5, ceiling[hw], "  ceiling", va="bottom",
                    ha="right", fontsize=7)
        ax.set_xticks(x)
        ax.set_xticklabels([MODE_LABEL[m] for m in sub["mode"]], fontsize=7,
                           rotation=20, ha="right")
        ax.set_title(f"Window ±{hw} bp", fontsize=9)
        ax.grid(False)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color("black")
        ax.tick_params(colors="black")
    axes[0][0].set_ylabel("Validation count Pearson r")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        p = os.path.join(args.figdir, f"h3k27ac_model_comparison.{ext}")
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(f"Wrote {p}")


if __name__ == "__main__":
    main()
