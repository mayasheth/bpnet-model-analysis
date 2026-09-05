#!/usr/bin/env python3
"""Can the sequence branch see the CpG-island/CTCF signature, or is accessibility drowning it?

The multimodal model over-predicts H3K27ac at highly accessible, GC/CpG-rich, promoter-
enriched elements. CpG islands and CTCF motifs are strong, local, learnable sequence signals
and the first conv layer is 21 bp -- wide enough for a 19 bp CTCF motif -- so capacity is not
the obvious explanation. Two hypotheses:

  A. SEQUENCE CANNOT SEE IT. The sequence-only model over-predicts these elements too.
     Then the fix is representational: give the model explicit GC / CpG-density / motif
     channels, or more capacity.
  B. ACCESSIBILITY DROWNS IT OUT. The sequence-only model gets them RIGHT while the
     multimodal model does not. Then sequence already carries the signal, the accessibility
     channel overrides it, and the fix is about the loss or the architecture's balance --
     adding sequence features would not help because the information is already there and
     being ignored.

These make opposite predictions, so one comparison separates them. Also reports what the
ATAC-only model does, which bounds how much of the error is attributable to accessibility.
"""
import numpy as np
import pandas as pd

D = "/oak/stanford/groups/engreitz/Users/sheth"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction/results"
OBS = f"{ABC}/2026_0721_h3k27ac_counting_comparison/K562_ATAC_H3K27ac_element/Neighborhoods/EnhancerList.txt"
ARMS = {m: f"{ABC}/2026_0903_predicted_activity/pred_k562_{m}/Neighborhoods/EnhancerList.txt"
        for m in ("sequence", "atac", "multimodal")}
COLS = ["chr", "start", "end", "H3K27ac.RPM", "ATAC.RPM", "class"]


def load(path, tag):
    df = pd.read_csv(path, sep="\t", usecols=COLS)
    df.index = df["chr"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
    return df.rename(columns={"H3K27ac.RPM": tag})


d = load(OBS, "obs")
for m, p in ARMS.items():
    d = d.join(load(p, m)[[m]], how="inner")
print(f"{len(d):,} regions with all four tracks\n")

for c in ["obs", "sequence", "atac", "multimodal"]:
    d["l_" + c] = np.log1p(d[c])
d["latac"] = np.log1p(d["ATAC.RPM"])

# Residualise each prediction and the truth on accessibility, so the strata are not simply
# "more open" -- the same control 4.10 uses.
A = np.vstack([np.ones(len(d)), d["latac"].to_numpy()]).T
for c in ["obs", "sequence", "atac", "multimodal"]:
    beta, *_ = np.linalg.lstsq(A, d["l_" + c].to_numpy(), rcond=None)
    d["r_" + c] = d["l_" + c].to_numpy() - A @ beta

# Strata defined by the MULTIMODAL model's error, exactly as in 4.10
d["err_mm"] = d["r_multimodal"] - d["r_obs"]
q = d["err_mm"].quantile([0.01, 0.05, 0.25, 0.75, 0.95, 0.99])
STRATA = [("under-predicted 1%", d["err_mm"] <= q[0.01]),
          ("under-predicted 5%", d["err_mm"] <= q[0.05]),
          ("typical (middle 50%)", (d["err_mm"] > q[0.25]) & (d["err_mm"] < q[0.75])),
          ("over-predicted 5%", d["err_mm"] >= q[0.95]),
          ("over-predicted 1%", d["err_mm"] >= q[0.99])]

print("Median accessibility-residualised error by model (0 = right, >0 = over-predicts)")
print(f"{'stratum':<24}{'n':>7}{'multimodal':>12}{'sequence':>11}{'ATAC only':>11}")
for lab, mask in STRATA:
    s = d[mask]
    print(f"{lab:<24}{len(s):>7,}"
          f"{(s['r_multimodal']-s['r_obs']).median():>12.3f}"
          f"{(s['r_sequence']-s['r_obs']).median():>11.3f}"
          f"{(s['r_atac']-s['r_obs']).median():>11.3f}")

print("\nRaw medians in the tail where the multimodal model over-predicts most")
s = d[STRATA[-1][1]]
print(f"   observed H3K27ac.RPM {s['obs'].median():.2f}   ATAC.RPM {s['ATAC.RPM'].median():.2f}")
for m in ("sequence", "atac", "multimodal"):
    print(f"   predicted by {m:<11} {s[m].median():.2f}")
print(f"\n   genome-wide observed median {d['obs'].median():.2f}; "
      f"sequence {d['sequence'].median():.2f}, atac {d['atac'].median():.2f}, "
      f"multimodal {d['multimodal'].median():.2f}")

print("\nVERDICT")
mm = (s["r_multimodal"] - s["r_obs"]).median()
sq = (s["r_sequence"] - s["r_obs"]).median()
if sq < 0.5 * mm:
    print("   B: sequence-only is substantially closer to the truth here, so the signal IS in")
    print("   sequence and the accessibility channel is overriding it. Adding sequence-derived")
    print("   features would not help; the balance or the loss is what needs changing.")
else:
    print("   A: sequence-only over-predicts these elements too, so the sequence branch is not")
    print("   capturing the CpG/CTCF signature. Explicit GC / CpG-density / motif channels or")
    print("   more sequence capacity are the indicated fix.")
