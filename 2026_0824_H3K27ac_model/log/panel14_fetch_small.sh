#!/bin/bash
# Small ENCODE files for panel_14, md5-checked as in 0.43: SK-N-SH poly(A) RNA-seq gene TPM for
# 4.54 ($SCRATCH/rnaseq/<acc>.tsv) and the released SK-N-SH_RA ATAC IDR peaks for
# $SCRATCH/panel_candidates. Untreated (ENCSR000CTT) and RA-treated (ENCSR000CPN) RNA are both
# fetched, Gingeras whole-cell, GENCODE V29 RSEM, two replicates each.
set -euo pipefail
fetch () {  # fetch <acc> <ext> <dir>
    local acc=$1 ext=$2 dir=$3 dest md5 have
    mkdir -p "$dir"; dest="$dir/$acc.$ext"
    md5=$(curl -s --max-time 60 -H 'Accept: application/json' \
          "https://www.encodeproject.org/files/$acc/?format=json" \
          | python3 -c 'import json,sys; print(json.load(sys.stdin)["md5sum"])')
    [[ -n "$md5" ]] || { echo "ERROR: no md5 for $acc" >&2; return 1; }
    if [[ -f "$dest" ]] && [[ "$(md5sum "$dest" | cut -d' ' -f1)" == "$md5" ]]; then
        echo "present, md5 OK: $dest"; return 0
    fi
    curl -L --fail --retry 5 --retry-delay 20 -o "$dest.part" \
         "https://www.encodeproject.org/files/$acc/@@download/$acc.$ext"
    have=$(md5sum "$dest.part" | cut -d' ' -f1)
    [[ "$have" == "$md5" ]] || { echo "ERROR: md5 $acc got $have want $md5" >&2; rm -f "$dest.part"; return 1; }
    mv "$dest.part" "$dest"; echo "OK $dest"
}
for acc in ENCFF378TYT ENCFF645FNQ ENCFF098ROF ENCFF425UDJ; do
    fetch "$acc" tsv "$SCRATCH/rnaseq"
done
fetch ENCFF855CQR bed.gz "$SCRATCH/panel_candidates"
# the five tissue candidates already listed in config/panel_candidates.tsv: re-check md5 only
for acc in ENCFF272UWG ENCFF580UXE ENCFF697MUB ENCFF623GOH ENCFF656MIU ENCFF195CXM; do
    fetch "$acc" bed.gz "$SCRATCH/panel_candidates"
done
for f in "$SCRATCH"/rnaseq/ENCFF{378TYT,645FNQ,098ROF,425UDJ}.tsv; do
    head -1 "$f" | tr '\t' '\n' | grep -qx TPM || { echo "ERROR: no TPM column in $f" >&2; exit 1; }
    echo "$(basename "$f"): $(($(wc -l < "$f") - 1)) genes"
done
echo "small fetch done"
