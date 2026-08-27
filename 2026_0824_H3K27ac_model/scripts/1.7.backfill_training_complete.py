#!/usr/bin/env python
"""
Backfill training_complete.json for folds that finished before the trainer wrote it.

A run on the preemptible `owners` partition can be evicted mid-fit, leaving a best-so-far
checkpoint that is indistinguishable from a finished one. The trainer now writes
training_complete.json as its last action and the evaluators refuse folds without it — but
folds trained earlier have no marker even though they finished.

Completion is inferred from the epoch log, which is decisive without needing SLURM
records (those are reset when a job is requeued):
  finished  <=> the run reached --max-epochs, OR at least `early_stopping` epochs elapsed
               after the last improvement (which is what stopped it)
  otherwise the log simply stops mid-fit, which is what preemption looks like.

Usage:
  python 1.7.backfill_training_complete.py --model-dir models/residual5pFIXED_hw500_clw10 \
      [--early-stopping 10] [--max-epochs 100] [--dry-run]
"""

import argparse
import json
import os

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--model-dir", required=True)
    p.add_argument("--early-stopping", type=int, default=10)
    p.add_argument("--max-epochs", type=int, default=100)
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def looks_finished(log_path, early_stopping, max_epochs):
    df = pd.read_csv(log_path, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    n_epochs = len(df)
    saved = df["Saved?"].astype(str).str.strip().str.lower() == "true"
    if not saved.any():
        return False, n_epochs, None
    last_improve = int(df.loc[saved, "Epoch"].iloc[-1])
    last_epoch = int(df["Epoch"].iloc[-1])
    since = last_epoch - last_improve
    if n_epochs >= max_epochs:
        return True, n_epochs, since
    return since >= early_stopping, n_epochs, since


def main():
    args = parse_args()
    n_done = n_skip = n_incomplete = 0
    for fold in sorted(os.listdir(args.model_dir)):
        d = os.path.join(args.model_dir, fold)
        if not (fold.startswith("fold") and os.path.isdir(d)):
            continue
        log = os.path.join(d, "multimodal_bpnet.log")
        marker = os.path.join(d, "training_complete.json")
        if os.path.exists(marker):
            print(f"  {fold}: marker already present")
            n_skip += 1
            continue
        if not os.path.exists(log):
            print(f"  {fold}: no log, cannot judge")
            n_incomplete += 1
            continue
        ok, n_epochs, since = looks_finished(log, args.early_stopping, args.max_epochs)
        if ok:
            print(f"  {fold}: FINISHED ({n_epochs} epochs, {since} since last "
                  f"improvement >= {args.early_stopping})")
            if not args.dry_run:
                with open(marker, "w") as f:
                    json.dump({"_backfilled": True, "n_epochs": n_epochs,
                               "epochs_since_last_improvement": since,
                               "_note": "inferred from the epoch log by "
                                        "1.7.backfill_training_complete.py"}, f, indent=2)
            n_done += 1
        else:
            print(f"  {fold}: INCOMPLETE ({n_epochs} epochs, only {since} since last "
                  f"improvement) — preempted or still running, no marker written")
            n_incomplete += 1
    print(f"{'Would mark' if args.dry_run else 'Marked'} {n_done} finished, "
          f"{n_skip} already had markers, {n_incomplete} judged incomplete")


if __name__ == "__main__":
    main()
