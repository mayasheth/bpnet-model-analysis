#!/bin/bash
#SBATCH -p engreitz,normal,owners
#SBATCH -t 12:00:00
#SBATCH --mem=8G
#SBATCH -c 2
#SBATCH -o log/fetchpanel.%A_%a.txt
#SBATCH -e log/fetchpanel.%A_%a.txt
#SBATCH --job-name=fetch_panel
#SBATCH --array=1-21%6
#
# Download the ENCODE files for the multi-cell-type p300 panel (A549, HepG2, MCF-7),
# one array task per file, six at a time so the portal is not hammered.
#
# WHY AN ARRAY AND NOT A LOOP. The ATAC BAMs are 8-13 GB each and the whole set is 72 GB.
# A single job that dies at 60 GB restarts from zero; an array retries only the file that
# failed, and `owners` preemption is likely over a transfer this long.
#
# IDEMPOTENT BY MD5, NOT BY EXISTENCE. A half-downloaded BAM left by a preempted task looks
# complete to `-f`, and a truncated BAM fails much later and confusingly, in whatever counts
# reads from it. Every task re-checks the md5 recorded by the portal and re-downloads on
# mismatch. This is the same lesson as 4.6's gzip -t gate, applied upstream.
#
# Layout follows the existing convention, Data/ENCODE/<cell>/<assay>/, which GM12878 already
# uses.
set -euo pipefail
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
P=$D/2026_0824_H3K27ac_model
DATA=/oak/stanford/groups/engreitz/Users/sheth/Data/ENCODE
MAN=$P/reference/p300_panel_manifest.tsv
cd "$P"; mkdir -p log

i=${SLURM_ARRAY_TASK_ID:-1}
row=$(grep -v '^#' "$MAN" | tail -n +2 | sed -n "${i}p")
[[ -n "$row" ]] || { echo "no manifest row $i"; exit 0; }
cell=$(cut -f1 <<<"$row"); assay=$(cut -f2 <<<"$row")
kind=$(cut -f4 <<<"$row"); acc=$(cut -f5 <<<"$row")
case "$kind" in bam) ext=bam ;; peak) ext=bed.gz ;; *) echo "bad kind $kind"; exit 1 ;; esac

out="$DATA/$cell/$assay"
mkdir -p "$out"
dest="$out/$acc.$ext"
echo "task $i: $cell $assay $kind $acc -> $dest"

md5=$(curl -s --max-time 60 "https://www.encodeproject.org/files/$acc/?format=json" \
      -H 'Accept: application/json' | python3 -c 'import json,sys; print(json.load(sys.stdin)["md5sum"])')
[[ -n "$md5" ]] || { echo "ERROR: no md5 for $acc" >&2; exit 1; }

if [[ -f "$dest" ]]; then
    have=$(md5sum "$dest" | cut -d' ' -f1)
    if [[ "$have" == "$md5" ]]; then echo "already present and md5 OK, skipping"; exit 0; fi
    echo "present but md5 mismatch ($have != $md5), re-downloading"
fi

curl -L --fail --retry 5 --retry-delay 20 --max-time 36000 \
     -o "$dest.part" "https://www.encodeproject.org/files/$acc/@@download/$acc.$ext"
have=$(md5sum "$dest.part" | cut -d' ' -f1)
if [[ "$have" != "$md5" ]]; then
    echo "ERROR: md5 mismatch after download: got $have want $md5" >&2
    rm -f "$dest.part"; exit 1
fi
mv "$dest.part" "$dest"
echo "OK $dest ($(du -h "$dest" | cut -f1))"
