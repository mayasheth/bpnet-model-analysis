#!/usr/bin/env python3
"""Emit the cell-type panel JSON that --cell-types-json consumes.

Generated rather than hand-written because the five cell types' tracks live in four
different project directories with four different naming conventions, and a typo in a path
would train silently against the wrong cell type's signal. Every path is checked to exist
before the file is written.

K562 and GM12878 point at their EXISTING tracks, not at rebuilds. That is deliberate: those
two are what every prior p300 result in this project was measured on, and swapping them for
rebuilds would make the new panel incomparable to F-021 and F-024. 0.45 rebuilds GM12878
separately, as a construction control for the three new cell types, and that rebuild is not
used for training.
"""
import argparse
import json
import os
import sys

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
DATA = "/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE"
PANEL_DIR = f"{D}/2026_0824_H3K27ac_model/data/p300_panel"
NEG = f"{D}/reference/genomewide_gc_stride_1000_flank_size_1057.gc.bed"

PANEL = {
    "K562": {
        "peaks": f"{D}/reference/ENCSR000EGE_peaks_inliers.narrowPeak",
        "signal_plus_bw": f"{D}/2025_0703_retrain_p300_model/data/ENCSR000EGE_plus.bigWig",
        "signal_minus_bw": f"{D}/2025_0703_retrain_p300_model/data/ENCSR000EGE_minus.bigWig",
        "accessibility_bw": f"{D}/2026_0529_multimodal_p300_model/data/atac.bw",
        "negatives": NEG,
    },
    "GM12878": {
        "peaks": f"{DATA}/GM12878/EP300/ENCFF926AKK.bed.gz",
        "signal_plus_bw": f"{D}/2026_0606_GM12878_transferability/data/ENCFF960OFK_plus.bw",
        "signal_minus_bw": f"{D}/2026_0606_GM12878_transferability/data/ENCFF941MGK_minus.bw",
        "accessibility_bw": f"{D}/2026_0606_GM12878_transferability/data/atac.bw",
        "negatives": NEG,
    },
}
for cell in ("A549", "HepG2", "MCF-7"):
    safe = cell.replace("-", "_")
    PANEL[cell] = {
        "peaks": f"{PANEL_DIR}/{safe}_ep300_peaks.narrowPeak",
        "signal_plus_bw": f"{PANEL_DIR}/{safe}_ep300_5p_plus.bw",
        "signal_minus_bw": f"{PANEL_DIR}/{safe}_ep300_5p_minus.bw",
        "accessibility_bw": f"{PANEL_DIR}/{safe}_atac.bw",
        "negatives": NEG,
    }

ap = argparse.ArgumentParser()
ap.add_argument("--out", default=f"{D}/2026_0824_H3K27ac_model/config/p300_panel.json")
ap.add_argument("--allow-missing", action="store_true",
                help="write the file even if tracks are absent, for a dry run")
a = ap.parse_args()

missing = [(c, k, v) for c, spec in PANEL.items() for k, v in spec.items()
           if not os.path.exists(v)]
for c, k, v in missing:
    print(f"MISSING  {c:<9} {k:<18} {v}", file=sys.stderr)
if missing and not a.allow_missing:
    raise SystemExit(f"\n{len(missing)} path(s) absent; build them first or pass "
                     f"--allow-missing")

os.makedirs(os.path.dirname(a.out), exist_ok=True)
with open(a.out, "w") as f:
    json.dump(PANEL, f, indent=2)
print(f"wrote {a.out} with {len(PANEL)} cell types: {', '.join(PANEL)}")
for c, spec in PANEL.items():
    print(f"  {c:<9} peaks={os.path.basename(spec['peaks'])}  "
          f"acc={os.path.basename(spec['accessibility_bw'])}")
