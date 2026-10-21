#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 4:00:00
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH -o log/bulk_chromvar.%j.txt
#SBATCH --job-name=bulk_chromvar
#
# Bulk chromVAR TF-activity vector per cell from ATAC only (4.63). Counts/annotation are cached in
# $SCRATCH/bulk_chromvar, so an owners preemption restart only repeats unfinished work.
set -euo pipefail
export PYTHONUNBUFFERED=1
P=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model
PY=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/.pixi/envs/multimodal/bin/python
cd "$P"; mkdir -p log

# guard test, both directions: real inputs pass, a p300 track and a non-ATAC manifest row fail
$PY - <<'EOF'
import importlib.util, pandas as pd
spec = importlib.util.spec_from_file_location("m", "scripts/4.63.bulk_chromvar_panel.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
full = pd.read_csv("config/panel_candidates.tsv", sep="\t")
man = full[full["cell"].isin(m.MANIFEST_CELL)]
m.check_inputs(m.BIGWIGS, man)
try:
    m.check_inputs(m.BIGWIGS, full)
except SystemExit as e:
    print("guard rejects unfiltered manifest (has DNase rows):", e)
else:
    raise SystemExit("guard failed to reject DNase rows")
for bw, mm in [({"x": m.A + "/data/p300_panel/A549_ep300_5p_plus.bw"}, man),
               (m.BIGWIGS, man.assign(assay=["ATAC-seq"] * (len(man) - 1) + ["ChIP-seq"]))]:
    try:
        m.check_inputs(bw, mm)
    except SystemExit as e:
        print("guard rejects as expected:", e)
    else:
        raise SystemExit("guard failed to reject")
print("guard ok")
EOF

$PY scripts/4.63.bulk_chromvar_panel.py --procs 8 "$@"
