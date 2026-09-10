#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 2:00:00
#SBATCH --mem=48G
#SBATCH -c 2
#SBATCH -o log/regen_n_figs.%j.txt
#SBATCH -e log/regen_n_figs.%j.txt
#SBATCH --job-name=regen_n_figs
#
# Redraw every figure whose n annotation was corrected.
#
# 3.2 does not cover these: it regenerates the model-free data figures (2, 3, 7, 8's source
# tables). The figures below are drawn straight from the per-fold result tables by 3.1 and
# 3.3-3.10, so they are re-run individually here.
#
# WHAT WAS WRONG. Seven of the eight report figures carrying an n annotation were
# misreporting it:
#   fig4, fig6, fig13  drawn entirely on the top signal quintile, annotated with the
#                      all-element count -- roughly fivefold too large
#   fig1, fig8, fig9, fig12  one panel all elements, one the top quintile, a single
#                      unqualified count covering both
#   fig10              `n_label(df)` read the loop variable after both loops had finished,
#                      so a four-panel figure spanning two cell types with different element
#                      sets was labelled with whichever table was read last
# fig11, fig14, fig15 and fig16 were checked and are correct.
#
# 2.30 must have run first: fig10, fig12 and fig13 now read `n_topq` from the per-fold
# tables, and the scripts will raise rather than fall back if the column is absent.
set -euo pipefail
export PYTHONUNBUFFERED=1
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
PY=$D/.pixi/envs/multimodal/bin/python
cd "$P"

for t in prof_residual_grid fragchan_k562 wide_k562 wide_gm12878; do
    head -1 "results/${t}_per_fold.tsv" | grep -qw n_topq \
        || { echo "ERROR: results/${t}_per_fold.tsv has no n_topq column; run 2.30 first" >&2; exit 1; }
done

echo "### 3.1 -- fig1, fig4, fig6 ###"
$PY scripts/3.1.plot_report_figures.py --results results --figdir figures
echo
echo "### 3.4 -- fig5_residual_grid_both ###"
$PY scripts/3.4.plot_residual_grid_both.py
echo
echo "### 3.5 -- fig9 ###"
$PY scripts/3.5.plot_accs5p_figure.py
echo
echo "### 3.6 -- fig8 ###"
$PY scripts/3.6.plot_coupling_comparison.py
echo
echo "### 3.8 -- fig10, fig11 ###"
$PY scripts/3.8.plot_wide_and_ceiling.py
echo
echo "### 3.9 -- fig12 ###"
$PY scripts/3.9.plot_fragment_channels.py
echo
echo "### 3.10 -- fig13 ###"
$PY scripts/3.10.plot_rc_averaging.py
echo
echo "### the corrected annotations, read back off the rendered figures ###"
$PY - <<'PY'
import re, subprocess
# Pull the text layer out of each PDF so the assertion is against what was actually drawn,
# not against what the plotting code was asked to draw.
FIGS = ["fig1_three_mode_comparison", "fig4_p300_vs_h3k27ac", "fig6_transfer",
        "fig8_coupling_across_celltypes", "fig9_atac_5prime_input",
        "fig10_wide_receptive_field", "fig12_fragment_channels", "fig13_rc_averaging"]
P = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model"
for f in FIGS:
    try:
        txt = subprocess.run(["pdftotext", f"{P}/figures/{f}.pdf", "-"],
                             capture_output=True, text=True, timeout=60).stdout
    except Exception as e:
        print(f"  {f:<34} (could not read pdf: {e})"); continue
    hits = [l.strip() for l in txt.splitlines() if l.strip().startswith("n = ")
            or "top q" in l or "all elements" in l]
    print(f"  {f:<34} {hits if hits else 'NO n ANNOTATION FOUND'}")
PY
echo ALL_DONE
