#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 6:00:00
#SBATCH --mem=16G
#SBATCH -c 8
#SBATCH -o log/depthaudit.%j.txt
#SBATCH -e log/depthaudit.%j.txt
#SBATCH --job-name=depth_audit
#
# Library depth for every cell type in the multi-cell-type p300 panel, for both assays.
#
# WHY THIS RUNS BEFORE ANY TRAINING. The project has measured the depth question twice and
# got OPPOSITE answers for the two roles, so neither can be assumed here:
#   TARGET depth does NOT matter. The 2026-09-08 subsampling test gave a depth-matched
#     K562 p300 model +0.202 on transfer against full-depth's +0.207, paired -0.005 (p=0.40).
#   INPUT depth DOES matter. F-017 depth-matched the accessibility libraries and turned the
#     GM12878->K562 collapse from -0.133 into +0.003.
# So the number that decides anything here is ATAC depth spread across the panel. Wide
# spread means the model can read library instead of biology; narrow spread means the
# per-cell-type z-scoring already in the trainer is enough.
#
# AND THE DECISION IS NOT AUTOMATIC EITHER WAY, which is why this reports rather than acts.
# F-017's own caveat: "depth-matching is the right control for a transfer question and the
# wrong choice for a production model, which should use all the reads it has." A
# multi-cell-type model is meant to be a production model, so thinning to the shallowest
# cell type spends real reads to remove a confound that may be small.
#
# Counts are alignments in the file as downloaded: ENCODE's filtered 'alignments' for ATAC,
# and for K562/GM12878 the tagAlign the project already uses. Those are the objects the
# tracks are built from, so they are the right unit even though the two formats differ.
set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
# CELLS=a,b (with DATA_ROOT for a $SCRATCH download) APPENDS rows for just those cells to the
# existing table, replacing any earlier rows for them, instead of rebuilding it; the 2026-10-20
# panel_14 biosamples were added this way. AUDIT_OUT redirects the table, for testing.
DATA=${DATA_ROOT:-/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE}
export PATH="$D/.pixi/envs/multimodal/bin:$PATH"
cd "$P"; mkdir -p log results
OUT=${AUDIT_OUT:-results/panel_depth_audit.tsv}

count_bam () { samtools view -c -@ 4 "$1" 2>/dev/null || echo 0; }
count_ta ()  { zcat "$1" 2>/dev/null | wc -l || echo 0; }

if [[ -n "${CELLS:-}" ]]; then
  [[ -s "$OUT" ]] || { echo "ERROR: CELLS appends to $OUT, which is missing" >&2; exit 1; }
  IFS=',' read -ra ADD <<< "$CELLS"
  NEW="$OUT.new.$$"; : > "$NEW"
  for cell in "${ADD[@]}"; do
    for assay in ATAC EP300; do
      nf=0
      for f in "$DATA/$cell/$assay"/*.bam; do
        [[ -e "$f" ]] || continue
        n=$(count_bam "$f"); nf=$((nf + 1))
        # count_bam reports 0 for an unreadable BAM; an appended 0 would read as a real depth
        [[ "$n" -gt 0 ]] || { echo "ERROR: no reads counted in $f" >&2; rm -f "$NEW"; exit 1; }
        echo -e "${cell}\t${assay}\tbam\t$(basename "$f")\t${n}" >> "$NEW"
        echo "  $cell $assay $(basename "$f") = $n"
      done
      [[ $nf -gt 0 ]] || { echo "ERROR: no $assay bams for $cell in $DATA" >&2; rm -f "$NEW"; exit 1; }
    done
  done
  awk -F'\t' -v cs="$CELLS" 'BEGIN{n=split(cs,a,","); for(i=1;i<=n;i++) d[a[i]]=1}
                             !($1 in d)' "$OUT" > "$OUT.keep.$$"
  cat "$OUT.keep.$$" "$NEW" > "$OUT"; rm -f "$OUT.keep.$$" "$NEW"
else
echo -e "cell\tassay\tsource\tfile\treads" > "$OUT"

# --- new cell types, downloaded 2026-09-18 -------------------------------------------
for cell in A549 HepG2 MCF-7; do
  for assay in ATAC EP300; do
    for f in "$DATA/$cell/$assay"/*.bam; do
      [[ -e "$f" ]] || continue
      n=$(count_bam "$f")
      echo -e "${cell}\t${assay}\tbam\t$(basename "$f")\t${n}" >> "$OUT"
      echo "  $cell $assay $(basename "$f") = $n"
    done
  done
done

# --- existing cell types, for the comparison that matters ----------------------------
for f in "$DATA/K562"/*.tn5.sorted.tagAlign.gz; do
  [[ -e "$f" ]] || continue
  n=$(count_ta "$f"); echo -e "K562\tATAC\ttagAlign\t$(basename "$f")\t${n}" >> "$OUT"
  echo "  K562 ATAC $(basename "$f") = $n"
done
for f in "$DATA/GM12878/ATAC"/*.tagAlign.gz; do
  [[ -e "$f" ]] || continue
  n=$(count_ta "$f"); echo -e "GM12878\tATAC\ttagAlign\t$(basename "$f")\t${n}" >> "$OUT"
  echo "  GM12878 ATAC $(basename "$f") = $n"
done
for f in "$DATA/K562/ENCFF466WKF.filtered.sorted.bam" "$DATA/K562/ENCFF163FSR.filtered.sorted.bam" \
         "$DATA/GM12878/EP300"/*.filtered.sorted.bam; do
  [[ -e "$f" ]] || continue
  n=$(count_bam "$f"); cell=K562; [[ "$f" == *GM12878* ]] && cell=GM12878
  echo -e "${cell}\tEP300\tbam\t$(basename "$f")\t${n}" >> "$OUT"
  echo "  $cell EP300 $(basename "$f") = $n"
done
fi

python3 - "$OUT" <<'PYIN'
import sys, collections
rows = [l.rstrip('\n').split('\t') for l in open(sys.argv[1])][1:]
tot = collections.defaultdict(int); nrep = collections.defaultdict(int)
for cell, assay, src, f, n in rows:
    tot[(cell, assay)] += int(n); nrep[(cell, assay)] += 1
print(f"\n{'cell':<9} {'assay':<6} {'reps':>5} {'total reads':>16}")
for assay in ('ATAC', 'EP300'):
    vals = {c: tot[(c, a)] for (c, a) in tot if a == assay}
    if not vals: continue
    for c in sorted(vals, key=lambda k: -vals[k]):
        print(f'{c:<9} {assay:<6} {nrep[(c, assay)]:>5} {vals[c]:>16,}')
    lo, hi = min(vals.values()), max(vals.values())
    print('  %s spread: %.2fx  (%s deepest, %s shallowest)' % (
        assay, hi / max(lo, 1),
        max(vals, key=vals.get), min(vals, key=vals.get)))
    if assay == 'ATAC':
        print('  ATAC is the INPUT, so this is the spread that matters (F-017).')
        print('  Under 2x: per-cell-type z-scoring in the trainer is enough.')
        print('  Well over 2x: thinning becomes a real trade-off, not an obvious win,')
        print('    because it costs the deep cell types reads a production model wants.')
PYIN
echo "wrote $OUT"
