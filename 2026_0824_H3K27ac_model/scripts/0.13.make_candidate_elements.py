#!/usr/bin/env python
"""
Convert an ENCODE_rE2G `Neighborhoods/EnhancerList.bed` into the narrowPeak form this
project's extraction code expects.

EnhancerList.bed is BED4: chr, start, end, "class|chr:start-end". The trainer and all
evaluators centre windows on `start + summit` (column 10 of a narrowPeak), so the element
midpoint must be written into a summit column — that is the convention the existing K562
and GM12878 element files follow (`summit == width // 2`).

Records the accessibility assay the elements were derived from in the output filename and
in a sidecar `.provenance.txt`, because it is not uniform across the panel: the
ENCODE_rE2G megamap runs are DNase/DHS-derived (`dhs_*` model dirs) while the TeloHAEC
runs are ATAC-derived (`atac_h3k27ac_powerlaw`). Mixing the two silently would confound
element definition with cell type.

Usage:
  python 0.13.make_candidate_elements.py \
      --enhancer-list <.../Neighborhoods/EnhancerList.bed> \
      --label Jurkat --derived-from DNase --outdir reference/celltype_elements
"""

import argparse
import os

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--enhancer-list", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--derived-from", required=True, choices=["DNase", "ATAC"],
                   help="Accessibility assay the element set was called from")
    p.add_argument("--model-dir", default="",
                   help="rE2G model directory the list came from, recorded for provenance")
    p.add_argument("--outdir", required=True)
    p.add_argument("--min-width", type=int, default=50,
                   help="Drop elements narrower than this (default: 50)")
    return p.parse_args()


CANONICAL = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.enhancer_list, sep="\t", header=None,
                     names=["chr", "start", "end", "name"], usecols=[0, 1, 2, 3])
    n0 = len(df)

    df = df[df["chr"].isin(CANONICAL)]
    n_chr = n0 - len(df)
    df = df[(df["end"] - df["start"]) >= args.min_width]
    n_narrow = n0 - n_chr - len(df)
    df = df.drop_duplicates(subset=["chr", "start", "end"]).reset_index(drop=True)
    n_dup = n0 - n_chr - n_narrow - len(df)

    width = df["end"] - df["start"]
    out = pd.DataFrame({
        "chr": df["chr"],
        "start": df["start"],
        "end": df["end"],
        "name": df["chr"].astype(str) + ":" + df["start"].astype(str) + "-"
                + df["end"].astype(str),
        "score": 0,
        "strand": ".",
        "signalValue": 0,
        "pValue": -1,
        "qValue": -1,
        # the load-bearing column: window centring is start + summit
        "summit": (width // 2).astype(int),
    }).sort_values(["chr", "start"]).reset_index(drop=True)

    safe = args.label.replace(" ", "_").replace("+", "plus")
    stem = f"{safe}_{args.derived_from}_candidate_elements"
    path = os.path.join(args.outdir, stem + ".narrowPeak")
    out.to_csv(path, sep="\t", header=False, index=False)

    with open(os.path.join(args.outdir, stem + ".provenance.txt"), "w") as f:
        f.write(f"label:          {args.label}\n")
        f.write(f"derived_from:   {args.derived_from}\n")
        f.write(f"source:         {args.enhancer_list}\n")
        f.write(f"model_dir:      {args.model_dir}\n")
        f.write(f"n_input:        {n0}\n")
        f.write(f"n_output:       {len(out)}\n")
        f.write(f"dropped_non_canonical_chr: {n_chr}\n")
        f.write(f"dropped_width_lt_{args.min_width}: {n_narrow}\n")
        f.write(f"dropped_duplicate_intervals: {n_dup}\n")
        f.write(f"mean_width:     {float(width.mean()):.1f}\n")
        f.write("summit_convention: width // 2 (element midpoint), matching "
                "K562_DNase_candidate_elements.narrowPeak\n")

    print(f"{args.label:<20} {args.derived_from:<6} {n0:>7} -> {len(out):>7} elements "
          f"(dropped {n_chr} non-canonical, {n_narrow} narrow, {n_dup} dup) "
          f"mean_width={width.mean():.0f}")
    print(f"  wrote {path}")


if __name__ == "__main__":
    main()
