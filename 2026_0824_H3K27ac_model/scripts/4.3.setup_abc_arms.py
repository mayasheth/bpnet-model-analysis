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
import argparse, os, shutil, time

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

def newest_input_mtime():
    """Latest mtime among every file named in the biosample table."""
    newest = 0.0
    for r in rows:
        for col in ("ATAC", "H3K27ac", "DHS"):
            for path in str(r[col]).split(","):
                if path and os.path.exists(path):
                    newest = max(newest, os.path.getmtime(path))
    return newest


# Peaks files in rule-chain order: each stage is stamped later than the one feeding it.
# Anything unmatched lands in the final bucket, which is safe -- a too-new file is never
# the thing that triggers a re-run of the rule that produced it.
PEAKS_ORDER = [
    lambda fs: [f for f in fs if f in ("macs2_peaks.xls", "macs2_summits.bed",
                                       "macs2_peaks.narrowPeak")],
    lambda fs: [f for f in fs if f == "macs2_peaks.narrowPeak.sorted"],
    lambda fs: [f for f in fs if f.endswith("Counts.bed")],
    lambda fs: [f for f in fs if f.endswith("candidateRegions.bed")],
    lambda fs: [f for f in fs if f in ("candidateRegions.qc.txt", "params.txt")],
    lambda fs: [f for f in fs
                if not any((f in ("macs2_peaks.xls", "macs2_summits.bed",
                                  "macs2_peaks.narrowPeak",
                                  "macs2_peaks.narrowPeak.sorted",
                                  "candidateRegions.qc.txt", "params.txt"),
                            f.endswith("Counts.bed"),
                            f.endswith("candidateRegions.bed")))],
]

ap = argparse.ArgumentParser()
ap.add_argument("--results-dir", default="results/2026_0903_predicted_activity")
ap.add_argument("--peaks-from", default=f"{JULY}/K562_ATAC_only/Peaks")
ap.add_argument("--force-peaks", action="store_true",
                help="Delete and re-copy each arm's Peaks/ even if present. REQUIRED after "
                     "regenerating the predicted bigwigs: the existing Peaks would then be "
                     "older than the new inputs, and Snakemake would re-run region calling, "
                     "giving each arm its own region set.")
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

    # `sort_narrowpeaks` takes a SECOND input beyond the narrowPeak: the per-run
    # results/<run>/tmp/<chrom sizes>.bed written by generate_chrom_sizes_bed_file. In a
    # fresh results directory that file is absent, so Snakemake generates it and then
    # re-runs everything downstream with reason "Input files updated by another job" --
    # which is what aborted the first two attempts, and is not an mtime problem at all.
    # Copy it too, stamped earlier than the Peaks chain.
    tmp_src = os.path.join(os.path.dirname(os.path.dirname(a.peaks_from)), "tmp")
    tmp_dst = f"{ABC}/{a.results_dir}/tmp"
    if a.force_peaks and os.path.isdir(tmp_dst):
        shutil.rmtree(tmp_dst)
    if os.path.isdir(tmp_src) and not os.path.exists(tmp_dst):
        shutil.copytree(tmp_src, tmp_dst)
        t0 = newest_input_mtime() + 30
        for root, _d, files in os.walk(tmp_dst):
            for fn in files:
                os.utime(os.path.join(root, fn), (t0, t0))
        print(f"  populated {tmp_dst} (per-run chrom-sizes bed)")
    elif not os.path.isdir(tmp_src):
        print(f"  WARNING: no tmp/ at {tmp_src}; generate_chrom_sizes_bed_file will run and "
              f"cascade into region calling")

    src = a.peaks_from
    for r in rows:
        dst = f"{ABC}/{a.results_dir}/{r['biosample']}/Peaks"
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if os.path.exists(dst):
            if not a.force_peaks:
                print(f"  exists, skipping {dst}")
                continue
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
        # copytree preserves mtimes, so the copies carry the ORIGINAL run's dates. Snakemake
        # compares mtimes, so July-dated Peaks against a September-dated predicted bigwig
        # reads as stale and re-runs region calling -- on a bigwig, for the activity-only
        # arms.
        #
        # Stamping everything with one timestamp is NOT enough, and the first attempt failed
        # exactly there: Snakemake requires an output to be STRICTLY NEWER than its input, so
        # equal mtimes still scheduled sort_narrowpeaks and make_candidate_regions. The files
        # must be staggered along the rule chain
        #   macs2_peaks.narrowPeak -> .sorted -> (+ Counts.bed) -> .candidateRegions.bed
        # with each stage later than the one feeding it, and the whole set later than every
        # external input.
        base = newest_input_mtime() + 60   # tmp/ is stamped at +30, so the chain follows it
        for step, names in enumerate(PEAKS_ORDER):
            t = base + 10 * step
            for fn in names(os.listdir(dst)):
                fp = os.path.join(dst, fn)
                if os.path.isfile(fp):
                    os.utime(fp, (t, t))
        os.utime(dst, (base + 10 * len(PEAKS_ORDER), ) * 2)
        print(f"  populated and staggered {dst}")
    # Verify the invariant rather than trusting it: every Peaks file must post-date every
    # input file named in the biosample table.
    newest_input = 0.0
    for r in rows:
        for col in ("ATAC", "H3K27ac", "DHS"):
            for path in str(r[col]).split(","):
                if path and os.path.exists(path):
                    newest_input = max(newest_input, os.path.getmtime(path))
    bad = []
    for r in rows:
        d = f"{ABC}/{a.results_dir}/{r['biosample']}/Peaks"
        mt = {fn: os.path.getmtime(os.path.join(d, fn)) for fn in os.listdir(d)}
        for fn, t in mt.items():
            if t < newest_input:
                bad.append(f"{r['biosample']}/{fn}: older than newest input")
        # each link of the rule chain must be strictly newer than the one before it
        chain = ["macs2_peaks.narrowPeak", "macs2_peaks.narrowPeak.sorted",
                 "macs2_peaks.narrowPeak.sorted.candidateRegions.bed"]
        for a_, b_ in zip(chain, chain[1:]):
            if a_ in mt and b_ in mt and not mt[b_] > mt[a_]:
                bad.append(f"{r['biosample']}: {b_} is not strictly newer than {a_}")
    tmp_dst = f"{ABC}/{a.results_dir}/tmp"
    if not os.path.isdir(tmp_dst) or not os.listdir(tmp_dst):
        bad.append("results/<run>/tmp is missing; generate_chrom_sizes_bed_file will run "
                   "and cascade into sort_narrowpeaks and make_candidate_regions")
    else:
        newest_tmp = max(os.path.getmtime(os.path.join(tmp_dst, f))
                         for f in os.listdir(tmp_dst)
                         if os.path.isfile(os.path.join(tmp_dst, f)))
        for r in rows:
            f = (f"{ABC}/{a.results_dir}/{r['biosample']}/Peaks/"
                 f"macs2_peaks.narrowPeak.sorted")
            if os.path.exists(f) and not os.path.getmtime(f) > newest_tmp:
                bad.append(f"{r['biosample']}: .sorted is not newer than results/<run>/tmp")
    if bad:
        raise SystemExit(f"ERROR: {len(bad)} problem(s); Snakemake would re-run region "
                         f"calling. First few: {bad[:3]}")
    print(f"Peaks populated, all files post-date every input (newest {newest_input:.0f}), "
          f"and the rule chain is strictly increasing; region calling will be skipped")
else:
    print("\nPeaks NOT copied (pass --copy-peaks once the bigwigs exist)")
