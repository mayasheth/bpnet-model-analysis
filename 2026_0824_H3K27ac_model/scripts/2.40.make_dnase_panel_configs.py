#!/usr/bin/env python3
"""Emit the three-way DNase-input transfer configs, one per TARGET cell type.

WHY A GENERATOR AND NOT THREE HAND-WRITTEN CONFIGS. Three targets x (one local arm + two
transferred arms + one floor) is twelve entries whose paths differ along three axes at once:
which cell type's model, which cell type's accessibility track, which cell type's H3K27ac.
Hand-writing that is how a transferred arm ends up scored against the SOURCE cell type's
H3K27ac, which would inflate it and look like a transfer success.

THE RULE, STATED ONCE AND APPLIED MECHANICALLY: an entry's model_dir comes from the SOURCE
cell type; its accessibility_bw, signal_*_bw and element set all come from the TARGET. A
transferred model is handed the target cell type's own version of the input it was trained
on, which is the deployment scenario, and is scored against the target's own observed mark.
This is the same convention 2.37 used for the 2x2.

THE FLOOR IS EACH TARGET'S OWN DNase-ONLY MODEL, and it exists for all three only because
K562's and GM12878's were trained on 2026-09-14. Every earlier accessibility-only model in
the project used ATAC, and THP-1 has no ATAC, so an ATAC floor would have meant a different
comparator in each column of a table whose subject is the DNase input. baseline is what
residual_pearson and incremental_r2 are defined against, so mixing assays there would make
those two columns incomparable across targets.

WHY THREE CELL TYPES AT ALL. F-009 (p300) and F-014 (DNase inputs) both report the same
asymmetry, K562-trained models keeping their advantage in GM12878 while GM12878-trained
models lose it in K562, and with two cell types "K562 is a good TRAINING cell type" and
"K562->GM12878 is a good PAIR" predict identical numbers. A third target separates them.

Usage: 2.40.make_dnase_panel_configs.py [--out-dir config] [--check]
  --check verifies that every fold of every arm left a completion marker before writing,
  which is the thing that fails when a fold was preempted. It cannot distinguish "still
  running" from "preempted", since both leave a best-so-far checkpoint and no marker, so it
  reports the observable state and leaves the cause to squeue. Run it ON SHERLOCK; /oak is
  not mounted locally and a local existence check reports everything missing.
"""
import argparse
import json
import os

R = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
P = f"{R}/2026_0824_H3K27ac_model"

# Per cell type: the DNase accessibility track, the H3K27ac target, the element set, the
# local multimodal DNase-input model, and the DNase-only floor.
CELLS = {
    "k562": {
        "acc": f"{P}/data/k562_dnase_5p.bw",
        "sig_plus": f"{P}/data/h3k27ac_5p_plus.bw",
        "sig_minus": f"{P}/data/h3k27ac_5p_minus.bw",
        "elements": f"{R}/reference/K562_DNase_candidate_elements.narrowPeak",
        "multimodal": f"{P}/models/multimodal5p_dnase_hw500_clw10",
        "floor": f"{P}/models/atac5p_dnase_hw500_clw10",
    },
    "gm12878": {
        "acc": f"{P}/data/gm12878_dnase_5p.bw",
        "sig_plus": f"{P}/data/gm12878_h3k27ac_5p_plus.bw",
        "sig_minus": f"{P}/data/gm12878_h3k27ac_5p_minus.bw",
        "elements": (f"{R}/2026_0606_GM12878_transferability/reference/"
                     f"GM12878_candidate_elements.narrowPeak"),
        "multimodal": f"{P}/models/gm12878_multimodal5p_dnase_hw500_clw10",
        "floor": f"{P}/models/gm12878_atac5p_dnase_hw500_clw10",
    },
    "thp1": {
        "acc": f"{P}/data/thp1_dnase_5p.bw",
        "sig_plus": f"{P}/data/thp1_h3k27ac_5p_plus.bw",
        "sig_minus": f"{P}/data/thp1_h3k27ac_5p_minus.bw",
        "elements": f"{R}/reference/THP1_DNase_candidate_elements.narrowPeak",
        "multimodal": f"{P}/models/thp1_multimodal5p_dnase_hw500_clw10",
        "floor": f"{P}/models/thp1_atac5p_dnase_hw500_clw10",
    },
}

ap = argparse.ArgumentParser()
ap.add_argument("--out-dir", default=f"{P}/config")
ap.add_argument("--check", action="store_true",
                help="require every fold0 checkpoint to exist; run on Sherlock")
a = ap.parse_args()


def entry(label, model_dir, target):
    t = CELLS[target]
    return {"label": label, "mode": "multimodal", "half_window": 500,
            "model_dir": model_dir,
            "accessibility_bw": t["acc"],
            "signal_plus_bw": t["sig_plus"], "signal_minus_bw": t["sig_minus"]}


missing, written = [], []
for target in CELLS:
    t = CELLS[target]
    sources = [c for c in CELLS if c != target]
    floor = {"label": "dnase_only_LOCAL", "mode": "atac", "half_window": 500,
             "model_dir": t["floor"],
             "accessibility_bw": t["acc"],
             "signal_plus_bw": t["sig_plus"], "signal_minus_bw": t["sig_minus"]}
    compare = [entry(f"{target}_LOCAL", t["multimodal"], target)]
    compare += [entry(f"{s}_to_{target}", CELLS[s]["multimodal"], target)
                for s in sorted(sources)]
    spec = {
        "_comment": (
            f"DNase-input H3K27ac panel, scored IN {target.upper()}. Baseline is "
            f"{target}'s own DNase-only model, so residual_pearson and incremental_r2 "
            f"measure what sequence adds beyond the accessibility that is locally "
            f"available. {target}_LOCAL is the model trained here; the other two are "
            f"transferred and are each handed {target}'s own DNase track and scored "
            f"against {target}'s own H3K27ac, which is the deployment scenario. Elements "
            f"are {target}'s own DNase-derived set. THP-1 has ONE DNase replicate, so it "
            f"has no DNase shape ceiling and cannot be a converter target; it is a "
            f"DNase-input panel member only."),
        "baseline": floor,
        "compare": compare,
        "_elements": t["elements"],
    }
    if a.check:
        # ALL FIVE FOLDS, not just fold0. 2.15 loads every fold, and a fold-3 preemption
        # breaks scoring exactly as badly as a fold-0 one while leaving fold0 looking fine.
        for e in [floor] + compare:
            for fold in range(5):
                ck = f'{e["model_dir"]}/fold{fold}/multimodal_bpnet.torch'
                mk = f'{e["model_dir"]}/fold{fold}/training_complete.json'
                if not os.path.exists(ck):
                    missing.append((e["model_dir"], fold, "no checkpoint"))
                elif not os.path.exists(mk):
                    # Deliberately does NOT claim preemption: a job still RUNNING also
                    # has a best-so-far checkpoint and no marker, and the two are
                    # indistinguishable from the filesystem. Check squeue for which.
                    missing.append((e["model_dir"], fold,
                                    "checkpoint present but no completion marker: still "
                                    "running, or preempted"))
        for k in ("acc", "sig_plus", "sig_minus", "elements"):
            if not os.path.exists(t[k]):
                missing.append((t[k], -1, f"missing {k} for {target}"))
    out = os.path.join(a.out_dir, f"dnase_panel_on_{target}_configs.json")
    written.append((out, spec))

if missing:
    # One unfinished model shows up in all three targets, as LOCAL in its own and as a
    # transferred arm in the other two, so report distinct problems rather than mentions.
    # Keyed on the MODEL PATH, not the label: `dnase_only_LOCAL` names a different model in
    # every target, so deduplicating by label silently merges three separate floors.
    uniq = sorted(set(missing))
    print(f"NOT WRITING. {len(uniq)} distinct problem(s):")
    for path, fold, why in uniq:
        where = os.path.basename(path) if fold < 0 else f"{os.path.basename(path)} fold{fold}"
        print(f"  {where}: {why}")
    raise SystemExit(1)

for out, spec in written:
    with open(out, "w") as f:
        json.dump(spec, f, indent=1)
    print(f"wrote {out}")
    print(f"  baseline {spec['baseline']['label']}  compare "
          f"{[c['label'] for c in spec['compare']]}")
    print(f"  elements {spec['_elements']}")
