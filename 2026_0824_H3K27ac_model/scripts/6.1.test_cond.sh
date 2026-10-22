#!/bin/bash
#SBATCH -p gpu,owners
#SBATCH -G 1
#SBATCH -t 1:00:00
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH -o log/test_cond61.%j.txt
#SBATCH -e log/test_cond61.%j.txt
#SBATCH --job-name=test_cond61
#
# Tests for 6.1's FiLM conditioning (--model-dir / --cond-cell), pass and fail directions:
#   1. Corgi, no conditioning: 6.2 block 0 still passes its gate against F-068's 6.1 partial
#   2. Corgi with --cond-cell must be refused
#   3. FiLM, a cell not in its table must be refused
#   4. FiLM per-context rows == fixed-row runs where the rows coincide (bit-identical sums)
#   5. FiLM with a different row in the same context gives different predictions
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
FILM=$P/models/p300_multicell_panelmean_film_chromvar13pc_drop0.1_panel13_holdK562_hw500_clw10
T=$SCRATCH/test_cond61
cd "$P"; rm -f "$T"/*.npz "$T"/*.log; mkdir -p "$T"
SMALL=(--stage run --chunk-start 0 --chunk-size 2 --contexts K562,HepG2
       --anchors GATA1_WGATAA --partners EBOX_TAL1_TCF3 --dmin -10 --dmax 10)
R61=(scripts/6.1.motif_pair_spacing_multimodal.py)

echo "== 1. Corgi regression via the 6.2 gate"
$PY scripts/6.2.motif_pair_seam_control.py --stage run --part-dir "$T/seam" --start 0 --stop 25

expect_fail() {
    local name=$1 pat=$2; shift 2
    if $PY "$@" > "$T/$name.log" 2>&1; then echo "FAIL: $name ran"; exit 1; fi
    grep -q "$pat" "$T/$name.log" || { echo "FAIL: $name wrong error"; cat "$T/$name.log"; exit 1; }
    echo "ok: $name refused ($(grep -o "error: .*" "$T/$name.log" | head -1))"
}
echo "== 2. Corgi + --cond-cell"
expect_fail corgi_cond "not cell-state conditioned" "${R61[@]}" "${SMALL[@]}" --cond-cell K562 --part "$T/x.npz"
echo "== 3. FiLM + unknown cell"
expect_fail film_bogus "not in fold" "${R61[@]}" "${SMALL[@]}" --model-dir "$FILM" --cond-cell NOPE --part "$T/x.npz"

echo "== 4/5. FiLM rows"
$PY "${R61[@]}" "${SMALL[@]}" --model-dir "$FILM" --part "$T/ctx.npz"
$PY "${R61[@]}" "${SMALL[@]}" --model-dir "$FILM" --cond-cell K562 --part "$T/k562.npz"
$PY "${R61[@]}" "${SMALL[@]}" --model-dir "$FILM" --cond-cell HepG2 --part "$T/hepg2.npz"
$PY - "$T" <<'EOF'
import sys, numpy as np
t = sys.argv[1]
z = {k: np.load(f"{t}/{k}.npz") for k in ("ctx", "k562", "hepg2")}
print("cond recorded:", {k: list(v["cond"]) for k, v in z.items()})
K, H = 0, 1                                   # contexts K562, HepG2
ok = True
for key in ("S", "MA", "MP"):
    same_k = np.array_equal(z["ctx"][key][K], z["k562"][key][K])
    same_h = np.array_equal(z["ctx"][key][H], z["hepg2"][key][H])
    print(f"{key}: ctx==fixed K562 row in K562 {same_k}; ctx==fixed HepG2 row in HepG2 {same_h}")
    ok &= same_k and same_h
d = max(np.abs(z["ctx"][k][H] - z["k562"][k][H]).max() for k in ("MA", "MP"))
print(f"HepG2 context, own row vs K562 row: max |marginal difference| {d:.3e}")
ok &= d > 0
print("ALL OK" if ok else "FAILED")
sys.exit(0 if ok else 1)
EOF
echo "Done"
