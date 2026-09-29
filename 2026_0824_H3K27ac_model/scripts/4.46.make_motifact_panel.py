#!/usr/bin/env python3
"""Panel JSON for Test 1: config/p300_panel.json plus each cell type's motif-activity track.

Usage: 4.46.make_motifact_panel.py --weighting oracle [--track-dir $SCRATCH/motif_activity]
"""
import argparse
import json
import os

A = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
ap = argparse.ArgumentParser()
ap.add_argument("--weighting", required=True)
ap.add_argument("--track-dir", default=os.path.join(os.environ.get("SCRATCH", "/tmp"),
                                                    "motif_activity"))
a = ap.parse_args()
panel = json.load(open(f"{A}/config/p300_panel.json"))
for cell, spec in panel.items():
    bw = f"{a.track_dir}/motifact_{a.weighting}_{cell.replace('-', '_')}.bw"
    assert os.path.exists(bw), f"missing {bw}"
    spec["motif_activity_bw"] = bw
out = f"{A}/config/p300_panel_motifact_{a.weighting}.json"
json.dump(panel, open(out, "w"), indent=2)
print("wrote", out)
