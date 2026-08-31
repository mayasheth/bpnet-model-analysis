#!/usr/bin/env python
"""Validate the new 5'-end ATAC tracks, and verify the read-length fix actually worked.

Integrity is necessary but not sufficient: the job ran on `owners` and its build function
skips any output where [[ -s "$out" ]], so a preempted mid-write file would be well-formed
and silently reused. Beyond readability, check the arithmetic that the fix implies:

  sum(atac.bw) / sum(atac_5p.bw) should equal the MEAN READ LENGTH.
  Full-interval coverage counts read_length bases per read; 5'-end counting counts exactly
  one. So the ratio is a direct readout of the smear we removed, and it should come out at
  ~94.5 for K562/GM12878 and ~35.5 for TeloHAEC.

  And the point of the exercise: in the 5' tracks, total signal should now be proportional
  to READ COUNT alone, so the K562-vs-TeloHAEC gap should shrink to the genuine depth
  difference with the read-length component gone.
"""
import os
import pyBigWig

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
P = f"{D}/2026_0824_H3K27ac_model"
PAIRS = [
    ("K562",          f"{D}/2026_0529_multimodal_p300_model/data/atac.bw",
                      f"{D}/2026_0529_multimodal_p300_model/data/atac_5p.bw", 94.5),
    ("GM12878",       f"{D}/2026_0606_GM12878_transferability/data/atac.bw",
                      f"{D}/2026_0606_GM12878_transferability/data/atac_5p.bw", 94.5),
    ("TeloHAEC_ctrl", f"{P}/data/panel/TeloHAEC_ctrl_atac.bw",
                      f"{P}/data/panel/TeloHAEC_ctrl_atac_5p.bw", 35.5),
    ("TeloHAEC_IL1b", f"{P}/data/panel/TeloHAEC_IL1b_atac.bw",
                      f"{P}/data/panel/TeloHAEC_IL1b_atac_5p.bw", 35.5),
    ("TeloHAEC_TNFa", f"{P}/data/panel/TeloHAEC_TNFa_atac.bw",
                      f"{P}/data/panel/TeloHAEC_TNFa_atac_5p.bw", 35.5),
    ("TeloHAEC_VEGF", f"{P}/data/panel/TeloHAEC_VEGF_atac.bw",
                      f"{P}/data/panel/TeloHAEC_VEGF_atac_5p.bw", 35.5),
]

def hdr(p):
    bw = pyBigWig.open(p)
    h, ch = bw.header(), bw.chroms()
    bw.stats("chr1", 1_000_000, 2_000_000)   # force an index read
    bw.close()
    return h, len(ch)

bad = []
rows = []
print("{:<16} {:>12} {:>12} {:>9} {:>9} {:>7}  {}".format(
    "cell type", "old sum", "5p sum (reads)", "ratio", "expected", "nchrom", "status"))
for name, old, new, exp_len in PAIRS:
    try:
        ho, no = hdr(old)
        hn, nn = hdr(new)
    except Exception as e:
        bad.append((name, f"UNREADABLE {e}"))
        print(f"{name:<16} UNREADABLE: {e}")
        continue
    ratio = ho["sumData"] / hn["sumData"]
    ok = nn >= 20 and hn["sumData"] > 0 and abs(ratio - exp_len) / exp_len < 0.15
    if not ok:
        bad.append((name, f"ratio {ratio:.1f} vs expected {exp_len}"))
    rows.append((name, hn["sumData"]))
    print("{:<16} {:>12.3e} {:>12.3e} {:>9.1f} {:>9.1f} {:>7}  {}".format(
        name, ho["sumData"], hn["sumData"], ratio, exp_len, nn, "OK" if ok else "CHECK"))

print("\n--- the point of the fix: 5' totals are read COUNTS, read length removed ---")
k562 = dict(rows).get("K562")
for name, s in rows:
    print("{:<16} {:>12.3e} reads   {:>6.1f}x vs K562".format(name, s, k562 / s if s else 0))

print("\nBAD: " + (", ".join(f"{n} ({w})" for n, w in bad) if bad else "none"))
raise SystemExit(1 if bad else 0)
