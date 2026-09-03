#!/usr/bin/env python3
"""Build the ABC biosample table and results skeleton for the predicted-activity arms.

DESIGN. Only the activity term changes across arms, so every arm must score an IDENTICAL
region set. ABC derives candidate regions per biosample by running MACS2 on the accessibility
file, which would re-derive them nine times and risk drift. Instead each arm's `Peaks/`
directory is pre-populated from the completed July run, so Snakemake sees those outputs as
present and skips region calling.

ORDER MATTERS: Snakemake compares mtimes, so Peaks must be copied AFTER the predicted
bigwigs exist. Otherwise the newer bigwig makes Peaks look stale and MACS2 runs anyway --
on a bigwig, for the activity-only arms, which would fail.

THREE ARM FAMILIES.
  geomean arms   H3K27ac column = predicted bigwig, ATAC column = the real tagAligns.
                 ABC computes activity_base = sqrt(norm_h3k27ac * norm_atac).
  k27only arms   ATAC column = predicted bigwig, no H3K27ac column. compute_activity then
                 returns activity_base = normalized_atac, which IS the prediction. This is
                 how "predicted H3K27ac alone as activity" is expressed without patching ABC.
  baselines      already complete in the July run and reused as-is:
                 K562_ATAC_only and K562_ATAC_H3K27ac_element.

WHY THE SCALE OF THE PREDICTION DOES NOT MATTER: run_qnorm is rank-based, mapping each
region's within-sample quantile onto the K562 reference. Keep use_qnorm True.
"""
import argparse, os, shutil

D = "/oak/stanford/groups/engreitz/Users/sheth"
ABC = f"{D}/ABC_working/ABC-Enhancer-Gene-Prediction"
JULY = f"{ABC}/results/2026_0721_h3k27ac_counting_comparison"
PRED = f"{D}/EP300_BPNet/2026_0824_H3K27ac_model/data/abc_predicted"
DATA = f"{D}/Data/ENCODE/K562"

ATAC = ",".join(f"{DATA}/{a}.tn5.sorted.tagAlign.gz"
                for a in ("ENCFF077FBI", "ENCFF128WZG", "ENCFF534DCE"))

MODELS = [("k562_atac", "K562 ATAC-only model"),
          ("k562_sequence", "K562 sequence-only model"),
          ("k562_multimodal", "K562 multimodal model"),
          ("gm12878_atac", "GM12878 ATAC-only model, applied to K562"),
          ("gm12878_sequence", "GM12878 sequence-only model, applied to K562"),
          ("gm12878_multimodal", "GM12878 multimodal model, applied to K562")]
# "predicted K27ac alone as activity" only for the K562-trained models
K27ONLY = ["k562_atac", "k562_sequence", "k562_multimodal"]

COLS = ["biosample", "DHS", "ATAC", "H3K27ac", "default_accessibility_feature",
        "HiC_file", "HiC_type", "HiC_resolution", "alt_TSS", "alt_genes"]

ap = argparse.ArgumentParser()
ap.add_argument("--results-dir", default="results/2026_0903_predicted_activity")
ap.add_argument("--peaks-from", default=f"{JULY}/K562_ATAC_only/Peaks")
ap.add_argument("--copy-peaks", action="store_true",
                help="Also populate each arm's Peaks/ from the July run. Run this only "
                     "after every predicted bigwig exists, so Peaks is the newer file.")
a = ap.parse_args()

rows, missing = [], []


def add(name, atac, h3k27ac, feature):
    r = {c: "" for c in COLS}
    r.update(biosample=name, ATAC=atac, H3K27ac=h3k27ac,
             default_accessibility_feature=feature)
    rows.append(r)


for tag, _desc in MODELS:
    bw = f"{PRED}/predk27ac_{tag}.bw"
    if not os.path.exists(bw):
        missing.append(bw)
    add(f"pred_{tag}", ATAC, bw, "ATAC")
for tag in K27ONLY:
    bw = f"{PRED}/predk27ac_{tag}.bw"
    add(f"k27only_{tag}", bw, "", "ATAC")

out = f"{ABC}/config/mine/config_biosamples_predicted_activity.tsv"
with open(out, "w") as f:
    f.write("\t".join(COLS) + "\n")
    for r in rows:
        f.write("\t".join(r[c] for c in COLS) + "\n")
print(f"wrote {out} with {len(rows)} arms")
for r in rows:
    kind = "geomean" if r["H3K27ac"] else "k27ac-as-activity"
    print(f"  {r['biosample']:<28} {kind}")

cfg = f"{ABC}/config/mine/config_predicted_activity.yaml"
base = open(f"{ABC}/config/config.yaml").read()
base = base.replace('biosamplesTable: "config/config_biosamples_chr22.tsv"',
                    'biosamplesTable: "config/mine/config_biosamples_predicted_activity.tsv"')
base = base.replace('results_dir: "results/"', f'results_dir: "{a.results_dir}"')
with open(cfg, "w") as f:
    f.write(base)
print("wrote", cfg)

if missing:
    print(f"\nWARNING: {len(missing)} predicted bigwig(s) missing; "
          f"do not run ABC until they exist:")
    for m in missing:
        print("  ", m)

if a.copy_peaks:
    if missing:
        raise SystemExit("\nrefusing to copy Peaks while bigwigs are missing: Peaks must be "
                         "newer than every input or Snakemake will re-run region calling")
    src = a.peaks_from
    for r in rows:
        dst = f"{ABC}/{a.results_dir}/{r['biosample']}/Peaks"
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if os.path.exists(dst):
            print(f"  exists, skipping {dst}")
            continue
        shutil.copytree(src, dst)
        print(f"  populated {dst}")
    print("Peaks populated; MACS2 should be skipped for every arm")
else:
    print("\nPeaks NOT copied (pass --copy-peaks once the bigwigs exist)")
