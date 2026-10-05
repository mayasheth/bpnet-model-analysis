#!/usr/bin/env python3
"""Pseudobulk 10x fragments -> a coordinate-sorted BAM of 1 bp reads, one per Tn5 insertion, for
TOBIAS ATACorrect (which reads BAMs, not fragments).

10x fragment coordinates are already Tn5-adjusted, so TOBIAS must run with --read_shift 0 0.
Each fragment gives two reads that keep the strand of the cut, because ATACorrect fits its bias
model per strand: the fragment start as a + strand read at `start`, the fragment end as a -
strand read whose 5' end is `end - 1`.

Usage: 2.1.fragments_to_insertion_bam.py --frag F --out insertions.bam [--chroms chr1,...]
"""
import argparse
import subprocess

import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--frag", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--genome-fai", default="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa.fai")
ap.add_argument("--chroms", default=",".join([f"chr{i}" for i in range(1, 23)] + ["chrX"]))
a = ap.parse_args()

chroms = a.chroms.split(",")
sizes = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.genome_fai)}
header = {"HD": {"VN": "1.6", "SO": "coordinate"},
          "SQ": [{"SN": c, "LN": sizes[c]} for c in chroms]}
unsorted = a.out + ".unsorted.bam"
n = 0
with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
    tid = {c: i for i, c in enumerate(chroms)}
    for c in chroms:
        p = subprocess.Popen(["tabix", a.frag, c], stdout=subprocess.PIPE, text=True)
        for line in p.stdout:
            f = line.split("\t")
            s, e = int(f[1]), int(f[2])
            for pos, rev in ((s, False), (e - 1, True)):
                r = pysam.AlignedSegment()
                r.query_name = f"i{n}"
                r.flag = 16 if rev else 0
                r.reference_id = tid[c]
                r.reference_start = pos
                r.mapping_quality = 60
                r.cigartuples = [(0, 1)]
                r.query_sequence = "N"
                r.query_qualities = pysam.qualitystring_to_array("I")
                out.write(r)
                n += 1
        p.wait()
        print(f"{c}: {n:,} insertions so far", flush=True)
pysam.sort("-@", "4", "-m", "2G", "-o", a.out, unsorted)
pysam.index(a.out)
subprocess.run(["rm", "-f", unsorted], check=True)
print(f"wrote {a.out} ({n:,} reads)")
