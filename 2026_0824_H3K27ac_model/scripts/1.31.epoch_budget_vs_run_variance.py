#!/usr/bin/env python3
"""Isolate the epoch-budget effect from run-to-run variance using only the training logs.

WHY THIS EXISTS. F-019 claims 100 epochs beats early stopping by +0.0086 AUPRC downstream.
But 1.29 is a SEPARATE training run from the 1.11 baseline, not a continuation, so that
number conflates the epoch budget with run variance at one run per condition. The handover
proposed scoring the same run's epoch-35 and epoch-100 checkpoints, but only two checkpoints
per fold were saved (best-validation and final), so that exact test is unavailable.

WHAT IS AVAILABLE, AND IT IS BETTER THAN NOTHING. Both runs logged validation metrics every
epoch. Comparing the two runs AT THE SAME EPOCH measures run variance directly, because the
epoch budget cannot affect a metric recorded before either run ended. That gives a yardstick
to hold the best-checkpoint difference against.
"""
import csv, statistics

P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model/models"
RUNS = {"baseline": "multimodal5p_accs5p_hw500_clw10",
        "ep100":    "multimodal5p_accs5p_ep100_hw500_clw10"}


def load(run, fold):
    with open(f"{P}/{RUNS[run]}/fold{fold}/multimodal_bpnet.log") as fh:
        r = csv.reader(fh, delimiter="\t")
        next(r)
        out = {}
        for row in r:
            if len(row) < 11:
                continue
            out[int(row[0])] = {"mnll": float(row[6]), "prof_r": float(row[7]),
                                "cnt_r": float(row[8]), "saved": row[10] == "True"}
        return out


def best(d):
    """The LAST saved epoch, which is what multimodal_bpnet.torch actually holds.

    Not the minimum-MNLL saved epoch. bpnetlite saves on the total validation loss, which
    is MNLL plus the weighted count MSE, so a later save can have a higher MNLL: baseline
    fold 1 saves at both 44 and 45 with MNLL 188.3003 then 188.3092. Picking the MNLL
    minimum would report weights that were overwritten and never scored.
    """
    ep = max(k for k, v in d.items() if v["saved"])
    return ep, d[ep]["mnll"]


print("RUN VARIANCE AT MATCHED EPOCHS")
print("Same config, same epoch, different run. All of this difference is run variance.\n")
hdr = ("fold", "epoch", "base MNLL", "ep100 MNLL", "delta", "base cnt_r", "ep100 cnt_r", "delta")
print(f"{hdr[0]:<6}{hdr[1]:>6}{hdr[2]:>12}{hdr[3]:>12}{hdr[4]:>10}{hdr[5]:>12}{hdr[6]:>13}{hdr[7]:>10}")
vm, vc = [], []
for fold in range(5):
    b, e = load("baseline", fold), load("ep100", fold)
    for ep in (20, 25, 30):
        if ep in b and ep in e:
            dm = e[ep]["mnll"] - b[ep]["mnll"]
            dc = e[ep]["cnt_r"] - b[ep]["cnt_r"]
            vm.append(dm)
            vc.append(dc)
            print(f"{fold:<6}{ep:>6}{b[ep]['mnll']:>12.4f}{e[ep]['mnll']:>12.4f}{dm:>+10.4f}"
                  f"{b[ep]['cnt_r']:>12.5f}{e[ep]['cnt_r']:>13.5f}{dc:>+10.5f}")
print(f"\n  MNLL          mean {statistics.mean(vm):+.4f}  sd {statistics.stdev(vm):.4f}  "
      f"max|d| {max(abs(x) for x in vm):.4f}")
print(f"  count Pearson mean {statistics.mean(vc):+.5f}  sd {statistics.stdev(vc):.5f}  "
      f"max|d| {max(abs(x) for x in vc):.5f}")

print("\n\nBEST-CHECKPOINT DIFFERENCE, WHICH IS WHAT F-019 ACTUALLY SCORED")
hdr = ("fold", "base ep", "ep100 ep", "base MNLL", "ep100 MNLL", "gain", "base cnt_r", "ep100 cnt_r", "gain")
print(f"{hdr[0]:<6}{hdr[1]:>9}{hdr[2]:>10}{hdr[3]:>12}{hdr[4]:>12}{hdr[5]:>10}{hdr[6]:>12}{hdr[7]:>13}{hdr[8]:>10}")
gm, gc, later = [], [], 0
for fold in range(5):
    b, e = load("baseline", fold), load("ep100", fold)
    be, bm = best(b)
    ee, em = best(e)
    dm, dc = em - bm, e[ee]["cnt_r"] - b[be]["cnt_r"]
    gm.append(dm)
    gc.append(dc)
    later += ee > be
    print(f"{fold:<6}{be:>9}{ee:>10}{bm:>12.4f}{em:>12.4f}{dm:>+10.4f}"
          f"{b[be]['cnt_r']:>12.5f}{e[ee]['cnt_r']:>13.5f}{dc:>+10.5f}")
print(f"\n  MNLL          mean {statistics.mean(gm):+.4f}  (negative = the longer run is better)")
print(f"  count Pearson mean {statistics.mean(gc):+.5f}")
print(f"  the longer run's minimum landed later in {later} of 5 folds")
print(f"\n  Run-variance sd at matched epochs was MNLL {statistics.stdev(vm):.4f} / "
      f"cnt_r {statistics.stdev(vc):.5f}.")
print("  If the best-checkpoint gain sits inside that, the budget is not resolved by this run pair.")

print("\n\nFINAL EPOCH vs BEST VALIDATION, WITHIN THE ep100 RUN")
print("Same run, so zero run variance. This asks whether training past the validation")
print("minimum is itself useful, which is the part a second seed cannot answer.\n")
hdr = ("fold", "best ep", "best MNLL", "ep99 MNLL", "delta", "best cnt_r", "ep99 cnt_r", "delta")
print(f"{hdr[0]:<6}{hdr[1]:>9}{hdr[2]:>12}{hdr[3]:>12}{hdr[4]:>10}{hdr[5]:>12}{hdr[6]:>12}{hdr[7]:>10}")
for fold in range(5):
    e = load("ep100", fold)
    ee, em = best(e)
    last = max(e)
    print(f"{fold:<6}{ee:>9}{em:>12.4f}{e[last]['mnll']:>12.4f}{e[last]['mnll'] - em:>+10.4f}"
          f"{e[ee]['cnt_r']:>12.5f}{e[last]['cnt_r']:>12.5f}"
          f"{e[last]['cnt_r'] - e[ee]['cnt_r']:>+10.5f}")
