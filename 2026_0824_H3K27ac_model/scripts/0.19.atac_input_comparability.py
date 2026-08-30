#!/usr/bin/env python
"""Is the ATAC input comparable across cell types, or does read length distort it?

K562/GM12878 ATAC tagAlign entries are 94-95 bp; TeloHAEC's are 35-36 bp. `genomecov -bg`
makes coverage proportional to interval width, so K562's accessibility is smeared over
~2.6x more bases per insertion. z-normalisation (what the trainer does) fixes mean and sd
but cannot fix a spatial difference -- structurally the same defect as the 250 bp fragment
extension rejected for the H3K27ac target.

Whether it MATTERS is an empirical question, so measure it instead of assuming:

  width_fwhm     full width at half maximum of the normalised ATAC meta-profile around
                 element centres. If read length is smearing the input, K562's peak should
                 be measurably broader than TeloHAEC's.
  shape_r        correlation between normalised meta-profiles, after rescaling each to unit
                 mass. Near 1.0 => the difference is amplitude only and normalisation
                 handles it; clearly below => genuine input-domain shift.

Each cell type uses ITS OWN elements, since element sets differ; the comparison is of
profile shape, not of per-element values.
"""
import numpy as np
import pandas as pd
import pyBigWig

D = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
P = f"{D}/2026_0824_H3K27ac_model"
NP = ["chr", "start", "end", "name", "score", "strand",
      "signalValue", "pValue", "qValue", "summit"]

TRACKS = [
    ("K562",          f"{D}/2026_0529_multimodal_p300_model/data/atac.bw",
     f"{D}/reference/K562_DNase_candidate_elements.narrowPeak"),
    ("GM12878",       f"{D}/2026_0606_GM12878_transferability/data/atac.bw",
     f"{D}/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak"),
    ("TeloHAEC_ctrl", f"{P}/data/panel/TeloHAEC_ctrl_atac.bw",
     f"{P}/reference/celltype_elements/TeloHAEC_ctrl_ATAC_candidate_elements.narrowPeak"),
    ("TeloHAEC_TNFa", f"{P}/data/panel/TeloHAEC_TNFa_atac.bw",
     f"{P}/reference/celltype_elements/TeloHAEC_TNFa_ATAC_candidate_elements.narrowPeak"),
]
FLANK, BIN, NS = 2000, 25, 8000


def load_els(path):
    d = pd.read_csv(path, sep="\t", header=None, names=NP)
    mid = (d["end"] - d["start"]) // 2
    d["center"] = d["start"] + d["summit"].where(d["summit"] >= 0, mid)
    return d


def profile(bw_path, els, rng):
    bw = pyBigWig.open(bw_path)
    ch = bw.chroms()
    nb = (2 * FLANK) // BIN
    acc, n = np.zeros(nb), 0
    idx = rng.choice(len(els), size=min(NS, len(els)), replace=False)
    for i in idx:
        r = els.iloc[i]
        s, e = int(r["center"]) - FLANK, int(r["center"]) + FLANK
        if r["chr"] not in ch or s < 0 or e > ch[r["chr"]]:
            continue
        v = np.nan_to_num(np.array(bw.values(r["chr"], s, e), dtype=np.float64))
        acc += v.reshape(nb, BIN).mean(axis=1)
        n += 1
    bw.close()
    return acc / max(n, 1), n


def fwhm(x, y):
    """Full width at half maximum, above the distal baseline."""
    base = np.median(np.r_[y[:8], y[-8:]])
    pk = y.max()
    half = base + (pk - base) / 2
    above = np.where(y >= half)[0]
    if len(above) < 2:
        return np.nan
    return float(x[above[-1]] - x[above[0]])


nb = (2 * FLANK) // BIN
x = (np.arange(nb) - nb / 2) * BIN + BIN / 2
profs = {}
print("{:<16} {:>8} {:>12} {:>10} {:>12}".format("cell type", "n", "peak_height", "FWHM_bp", "baseline"))
for name, bwp, elp in TRACKS:
    rng = np.random.default_rng(0)
    els = load_els(elp)
    p, n = profile(bwp, els, rng)
    profs[name] = p
    base = np.median(np.r_[p[:8], p[-8:]])
    print("{:<16} {:>8,} {:>12.3f} {:>10.0f} {:>12.3f}".format(name, n, p.max(), fwhm(x, p), base))

print("\n--- shape agreement after normalising each to unit mass ---")
names = list(profs)
norm = {k: v / v.sum() for k, v in profs.items()}
print("{:<16}".format("") + "".join("{:>16}".format(n) for n in names))
for a in names:
    row = "{:<16}".format(a)
    for b in names:
        row += "{:>16.5f}".format(np.corrcoef(norm[a], norm[b])[0, 1])
    print(row)

out = pd.DataFrame({"distance": x, **{k: norm[k] for k in names}})
p1 = f"{P}/results/atac_input_comparability.tsv"
out.to_csv(p1, sep="\t", index=False)
print("\nWrote", p1)
