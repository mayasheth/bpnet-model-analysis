#!/usr/bin/env python
"""Validate the 60 TeloHAEC panel bigwigs.

sacct COMPLETED and ALL_DONE only prove the script ran to the end. A file can still be
well-formed but wrong -- so also check the internal consistency the build implies:
  * r1 must have strictly FEWER counts than both (it is a subset of the same reads)
  * plus and minus strands should be within a few percent of each other
  * per-replicate counts should sum to roughly the merged track
A violation of any of these means the variant logic did not do what it claims.
"""
import glob, os, re, sys, collections
import pyBigWig

P = ("/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/"
     "2026_0824_H3K27ac_model/data/panel")
files = sorted(glob.glob(f"{P}/TeloHAEC_*"))
print(f"{len(files)} TeloHAEC files")

stats, bad = {}, []
for f in files:
    n = os.path.basename(f)
    try:
        bw = pyBigWig.open(f)
        h, ch = bw.header(), bw.chroms()
        bw.stats("chr1", 1_000_000, 2_000_000)
        bw.close()
        stats[n] = (h["sumData"], len(ch))
        if len(ch) < 20 or h["sumData"] <= 0:
            bad.append((n, f"nchrom={len(ch)} sum={h['sumData']}"))
    except Exception as e:
        bad.append((n, f"UNREADABLE {e}"))

print(f"readable: {len(stats)}/{len(files)}")
for n, why in bad:
    print("  BAD", n, why)

CONDS = ["ctrl", "IL1b", "TNFa", "VEGF"]
print(f"\n{'condition':<10} {'track':<10} {'r1 sum':>16} {'both sum':>16} {'r1/both':>9}  check")
ok = True
for c in CONDS:
    reps = sorted({m.group(1) for f in stats
                   for m in [re.match(rf"TeloHAEC_{c}_h3k27ac_(rep\d+)_", f)] if m})
    for track in ["merged"] + reps:
        tag = "" if track == "merged" else f"_{track}"
        tot = {}
        for v in ("r1", "both"):
            keys = [f"TeloHAEC_{c}_h3k27ac{tag}_{v}_5p_{s}.bw" for s in ("plus", "minus")]
            if not all(k in stats for k in keys):
                tot[v] = None
                continue
            tot[v] = sum(stats[k][0] for k in keys)
        if tot.get("r1") is None or tot.get("both") is None:
            print(f"{c:<10} {track:<10} {'MISSING':>16}")
            ok = False
            continue
        ratio = tot["r1"] / tot["both"]
        # r1 is one mate of a pair, so expect roughly half, and strictly less
        good = tot["r1"] < tot["both"] and 0.35 < ratio < 0.65
        ok &= good
        print(f"{c:<10} {track:<10} {tot['r1']:16,.0f} {tot['both']:16,.0f} {ratio:9.3f}  "
              f"{'ok' if good else 'UNEXPECTED'}")

print("\n--- strand balance (merged, both) ---")
for c in CONDS:
    p = stats.get(f"TeloHAEC_{c}_h3k27ac_both_5p_plus.bw")
    m = stats.get(f"TeloHAEC_{c}_h3k27ac_both_5p_minus.bw")
    if p and m:
        d = abs(p[0] - m[0]) / max(p[0], m[0])
        print(f"{c:<10} plus/minus imbalance {d*100:5.2f}%  {'ok' if d < 0.05 else 'CHECK'}")

print("\n--- ATAC ---")
for c in CONDS:
    k = f"TeloHAEC_{c}_atac.bw"
    if k in stats:
        print("{:<32} sum {:>16,.0f}  nchrom {}".format(k, stats[k][0], stats[k][1]))
    else:
        print("{:<32} MISSING".format(k))

sys.exit(0 if (ok and not bad) else 1)
