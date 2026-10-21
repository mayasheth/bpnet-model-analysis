#!/bin/bash
# panel_14: six ENCODE biosamples added to the p300 panel (reference/p300_panel14_manifest.tsv).
# 0.43 fetch to $SCRATCH -> 0.45 tracks -> 1.33 --with-panel14 json; 0.44 CELLS= audit rows;
# panel14_sanity.sh last. engreitz,normal throughout: the downloads and the track builds are
# long and uncheckpointed. BAMs stay on $SCRATCH; delete only after TRACKS_VERIFIED.
set -euo pipefail
P=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model
PY=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/.pixi/envs/multimodal/bin/python
ROOT=$SCRATCH/p300_panel14
CELLS=SK-N-SH_RA,esophagus_muscularis_mucosa,gastroesophageal_sphincter,sigmoid_colon,tibial_nerve,upper_lobe_of_left_lung
cd "$P"
N=$(grep -v '^#' reference/p300_panel14_manifest.tsv | tail -n +2 | grep -c .)
[[ "$N" -eq 27 ]] || { echo "ERROR: expected 27 manifest rows, got $N" >&2; exit 1; }
F=$(DATA_ROOT=$ROOT MANIFEST=$P/reference/p300_panel14_manifest.tsv \
    sbatch --parsable --export=ALL -p engreitz,normal --array=1-$N%6 -J fetch_p14 scripts/0.43.fetch_p300_panel.sh)
S=$(sbatch --parsable -p engreitz,normal -c 1 --mem 2G -t 1:00:00 -J fetch_p14_small \
    -o log/panel14_fetch_small.%j.txt log/panel14_fetch_small.sh)
T=$(DATA_ROOT=$ROOT CELLS=$CELLS sbatch --parsable --export=ALL -p engreitz,normal --array=1-6 \
    -J p300_tracks14 --dependency=afterok:$F scripts/0.45.build_p300_panel_tracks.sh)
A=$(DATA_ROOT=$ROOT CELLS=$CELLS sbatch --parsable --export=ALL -p engreitz,normal \
    --dependency=afterok:$F scripts/0.44.panel_depth_audit.sh)
J=$(sbatch --parsable -p engreitz,normal -c 1 --mem 2G -t 0:20:00 -J p14_json -o log/panel14_json.%j.txt \
    --dependency=afterok:$T --wrap "set -e; cd $P; $PY scripts/1.33.make_panel_json.py --with-panel14 --out config/p300_panel_14.json")
C=$(sbatch --parsable -p engreitz,normal -c 2 --mem 8G -t 1:00:00 -J p14_sanity -o log/panel14_sanity.%j.txt \
    --dependency=afterok:$T:$A log/panel14_sanity.sh)
echo "fetch=$F small=$S tracks=$T audit=$A json=$J sanity=$C"
