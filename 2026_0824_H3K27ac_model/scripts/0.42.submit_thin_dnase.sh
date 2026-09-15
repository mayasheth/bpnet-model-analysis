#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 6:00:00
#SBATCH --mem=32G
#SBATCH -o log/thin.%j.txt
#SBATCH -e log/thin.%j.txt
#SBATCH --job-name=thin_dnase
#SBATCH -n 1
#
# Depth-match K562 and THP-1 DNase to GM12878's library depth, the one untested mechanism
# for F-016. GM12878 is already the shallowest, so it is the target and is not thinned; its
# existing models ARE the depth-matched arms for that cell type and are reused.
#
# Main-chromosome 5' totals: K562 293,265,722, THP-1 82,174,264, GM12878 50,299,042.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
TARGET=50299042
cd "$P"

for cell in k562 thp1; do
    OUT=data/${cell}_dnase_5p_depthmatched.bw
    if [[ -s "$OUT" ]]; then echo "exists: $OUT"; continue; fi
    echo "=== thinning $cell to $TARGET ==="
    $PY scripts/0.41.thin_5prime_track.py \
        --in-bw data/${cell}_dnase_5p.bw --out-bw "$OUT" --target-total $TARGET
done

echo "=== per-window means after thinning, each over its OWN elements ==="
$PY - <<'PYIN'
import numpy as np, pandas as pd, pyBigWig
R = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet"
P = f"{R}/2026_0824_H3K27ac_model"
S = {
 "k562 (thinned)":    (f"{P}/data/k562_dnase_5p_depthmatched.bw", f"{R}/reference/K562_DNase_candidate_elements.narrowPeak"),
 "gm12878 (as-is)":   (f"{P}/data/gm12878_dnase_5p.bw",           f"{R}/2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak"),
 "thp1 (thinned)":    (f"{P}/data/thp1_dnase_5p_depthmatched.bw", f"{R}/reference/THP1_DNase_candidate_elements.narrowPeak"),
}
for name, (bw, el) in S.items():
    els = pd.read_csv(el, sep="\t", header=None, usecols=[0,1,2,9],
                      names=["chr","start","end","summit"])
    rng = np.random.default_rng(0)
    els = els.iloc[rng.choice(len(els), min(3000, len(els)), replace=False)]
    b = pyBigWig.open(bw); sizes = b.chroms(); t = []
    for _, r in els.iterrows():
        c = int(r["start"]) + int(r["summit"]); s, e = c - 500, c + 500
        if r["chr"] not in sizes or s < 0 or e > sizes[r["chr"]]: continue
        v = b.values(r["chr"], s, e, numpy=True)
        if v is not None: t.append(float(np.nan_to_num(v).sum()))
    b.close(); t = np.array(t)
    print(f"  {name:<18} mean {t.mean():>8.1f}  median {np.median(t):>7.1f}")
PYIN
echo THIN_DONE
