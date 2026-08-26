#!/usr/bin/env python
"""
Backfill training_target.json for models trained before the trainer recorded it.

The offset-model guard in train_multimodal_bpnet.py compares the offset model's target
against the current run's. Models trained before 2026-08-26 have no record, so the guard
refuses them. This writes the record for model directories whose training configuration
is known unambiguously from the submit script that produced them.

Only use this where the configuration is genuinely known. If unsure, retrain.

Usage:
  python 1.5.backfill_training_target.py --model-dir models/atac5p_hw500_clw10 \
      --target fiveprime --mode atac --out-window 1000 [--dry-run]
"""

import argparse
import json
import os

PROJ = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
REPO = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"

TARGETS = {
    "fiveprime": {
        "signal_plus_bw": f"{PROJ}/data/h3k27ac_5p_plus.bw",
        "signal_minus_bw": f"{PROJ}/data/h3k27ac_5p_minus.bw",
    },
    # Retained only so pre-existing models can be labelled honestly. Do NOT train new
    # models on this target — it smears each read over 250 bp, which breaks MNLL's
    # read-count assumption and bleeds signal between neighbouring elements.
    "fragment250": {
        "signal_plus_bw": "/oak/stanford/groups/engreitz/Users/sheth/Data/share/IGV/ENCSR000AKP_coverage.bw",
        "signal_minus_bw": None,
    },
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--model-dir", required=True)
    p.add_argument("--target", required=True, choices=sorted(TARGETS))
    p.add_argument("--mode", required=True, choices=["sequence", "multimodal", "atac"])
    p.add_argument("--out-window", type=int, default=1000)
    p.add_argument("--n-layers", type=int, default=8)
    p.add_argument("--count-loss-weight", type=float, required=True)
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    tgt = TARGETS[args.target]
    rec = {
        "signal_plus_bw": tgt["signal_plus_bw"],
        "signal_minus_bw": tgt["signal_minus_bw"],
        "accessibility_bw": (None if args.mode == "sequence" else
                             f"{REPO}/2026_0529_multimodal_p300_model/data/atac.bw"),
        "peaks": f"{REPO}/reference/K562_DNase_candidate_elements.narrowPeak",
        "mode": args.mode,
        "in_window": args.out_window + 2 * (47 + sum(2 ** i for i in range(1, args.n_layers + 1))),
        "out_window": args.out_window,
        "n_layers": args.n_layers,
        "count_loss_weight": args.count_loss_weight,
        "count_offset_model": None,
        "_backfilled": True,
        "_backfill_note": "written by 1.5.backfill_training_target.py, not by the trainer",
    }
    n = 0
    for fold in sorted(os.listdir(args.model_dir)):
        d = os.path.join(args.model_dir, fold)
        if not (fold.startswith("fold") and os.path.isdir(d)):
            continue
        if not os.path.exists(os.path.join(d, "multimodal_bpnet.torch")):
            print(f"  skip {fold}: no checkpoint")
            continue
        out = os.path.join(d, "training_target.json")
        if os.path.exists(out):
            print(f"  skip {fold}: record already exists")
            continue
        print(f"  {'would write' if args.dry_run else 'wrote'} {out}")
        if not args.dry_run:
            with open(out, "w") as f:
                json.dump(rec, f, indent=2)
        n += 1
    print(f"{'Would backfill' if args.dry_run else 'Backfilled'} {n} fold(s) "
          f"in {args.model_dir} as target={args.target}")


if __name__ == "__main__":
    main()
