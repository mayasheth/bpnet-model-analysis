#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 3:00:00
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH -o log/setup_extra.%j.txt
#SBATCH --job-name=setup_extra
# 4.8 Generic setup for a scE2G run that adds external features to the multiome v3 model.
#   TAG    run tag: run dir $SCRATCH/sce2g_runs/${RUN_DATE}_k562_xu_$TAG, configs config/*_$TAG*
#   PAIRS  E-G pairs file with the feature columns (made beforehand, e.g. by 4.7)
#   SPEC   TSV, no header: model <tab> feature name <tab> input column <tab> aggregate <tab> description
#          [<tab> v3 feature to replace]; one row per added feature; a model listed several times gets
#          several features. The optional 6th column drops that v3 feature from the model (e.g.
#          ARC.E2G.Score, to swap in a recomputed ARC); added 2026-10-04 after the corr run.
# Feature names must map one-to-one to input columns across the run: scE2G keeps one feature table
# per run keyed by feature name (learnings 2026-10-04). The baseline multiome_powerlaw_n6 is included.
# Same seeding as 4.0. Run from the A2 dir.
set -euo pipefail
export PYTHONUNBUFFERED=1
A2=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
V3=/oak/stanford/groups/engreitz/Users/sheth/scE2G_temp/scE2G_main/models/multiome_powerlaw_v3/feature_table.tsv
BASE=$SCRATCH/sce2g_runs/2026_1002_k562_xu
: "${TAG:?}" "${PAIRS:?}" "${SPEC:?}"
RUN=$SCRATCH/sce2g_runs/${RUN_DATE:-2026_1004}_k562_xu_$TAG
MODELS=config/models_$TAG
cd $A2
[[ -e $RUN ]] && { echo "$RUN exists; not overwriting" >&2; exit 1; }
[[ -s $PAIRS && -s $SPEC ]] || { echo "PAIRS or SPEC missing" >&2; exit 1; }

# spec checks: feature <-> input one-to-one, every input column present in PAIRS
[[ $(cut -f2,3 $SPEC | sort -u | wc -l) -eq $(cut -f2 $SPEC | sort -u | wc -l) &&
   $(cut -f2,3 $SPEC | sort -u | wc -l) -eq $(cut -f3 $SPEC | sort -u | wc -l) ]] \
    || { echo "feature names and input columns are not one-to-one in $SPEC" >&2; exit 1; }
HDR=$( (zcat $PAIRS || true) | head -n 1 | tr '\t' '\n')
for c in $(cut -f3 $SPEC | sort -u); do
    grep -qx "$c" <<< "$HDR" || { echo "$c missing from $PAIRS" >&2; exit 1; }
    ! cut -f1 $V3 | grep -qx "$c" || { echo "$c clashes with a v3 feature name" >&2; exit 1; }
done

# results dir
mkdir -p $RUN
cp -a $BASE/tmp $RUN/
mkdir $RUN/K562_Xu_pl
for f in $BASE/K562_Xu_pl/*; do
    b=$(basename $f)
    case $b in
        multiome_powerlaw_n6|scATAC_powerlaw|feature_table.tsv|to_generate.txt|external_features_config.tsv|ActivityOnly_plus_external_features.tsv.gz|genomewide_features.tsv.gz) ;;
        *) cp -a $f $RUN/K562_Xu_pl/ ;;
    esac
done

# configs
EFC=config/external_features_$TAG.tsv
{
    printf 'input_col\tsource_col\taggregate_function\tjoin_by\tsource_file\n'
    cut -f3,4 $SPEC | sort -u | awk -F'\t' -v OFS='\t' -v p=$PAIRS '{print $1, $1, $2, "overlap", p}'
} > $EFC
awk -F'\t' -v OFS='\t' -v efc="$A2/$EFC" 'NR==1 {print $0, "external_features_config"; next} {print $0, efc}' \
    config/cell_clusters_K562_Xu.tsv > config/cell_clusters_K562_Xu_$TAG.tsv
{
    head -n 1 config/models_baseline.tsv
    grep -P '^multiome_powerlaw_n6\t' config/models_baseline.tsv
} > config/models_$TAG.tsv
for m in $(cut -f1 $SPEC | awk '!s[$0]++'); do
    mkdir -p $MODELS/$m
    REPL=$(awk -F'\t' -v m=$m '$1 == m && $6 != "" {print $6}' $SPEC | sort -u)
    for r in $REPL; do
        cut -f1 $V3 | grep -qxF "$r" || { echo "$m: $r to replace is not a v3 feature" >&2; exit 1; }
    done
    awk -F'\t' -v r="$REPL" 'BEGIN {n = split(r, x, "\n"); for (i = 1; i <= n; i++) d[x[i]] = 1} !($1 in d)' $V3 \
        > $MODELS/$m/feature_table.tsv
    sed -i -e '$a\' $MODELS/$m/feature_table.tsv
    awk -F'\t' -v OFS='\t' -v m=$m '$1 == m {print $2, $3, "NA", $4, 0, $5}' $SPEC >> $MODELS/$m/feature_table.tsv
    NREPL=$(grep -c . <<< "$REPL" || true)
    NR_=$(( $(grep -c . $V3) - NREPL + $(awk -F'\t' -v m=$m '$1 == m' $SPEC | wc -l) ))
    [[ $(grep -c . $MODELS/$m/feature_table.tsv) -eq $NR_ ]] || { echo "row count wrong in $m" >&2; exit 1; }
    for r in $REPL; do
        ! cut -f1 $MODELS/$m/feature_table.tsv | grep -qxF "$r" || { echo "$m: $r not removed" >&2; exit 1; }
    done
    printf '%s\tK562_Xu_pl\t\t%s\tFalse\t\n' $m $A2/$MODELS/$m/feature_table.tsv >> config/models_$TAG.tsv
done
sed -e "s#^results_dir: .*#results_dir: \"$RUN\"#" \
    -e "s#^cell_clusters: .*#cell_clusters: \"$A2/config/cell_clusters_K562_Xu_$TAG.tsv\"#" \
    -e "s#^model_config: .*#model_config: \"$A2/config/models_$TAG.tsv\"#" \
    config/config_training.yaml > config/config_training_$TAG.yaml
for k in "results_dir: \"$RUN\"" "cell_clusters_K562_Xu_$TAG.tsv" "models_$TAG.tsv"; do
    grep -qF "$k" config/config_training_$TAG.yaml || { echo "config patch missed: $k" >&2; exit 1; }
done
echo "== models"; cat config/models_$TAG.tsv
echo "== external features"; cat $EFC
echo "setup complete"
