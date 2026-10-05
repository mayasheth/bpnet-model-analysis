#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 3:00:00
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH -o log/setup_featB.%j.txt
#SBATCH --job-name=setup_featB
# 4.0 Set up a scE2G feature-B run (footprints as external features). SRC=tobias|print picks the
# 2.4 element table (TOBIAS scores, or PRINT scores at the same sites); RUN_DATE names the run.
#  1. E-G pairs file of element features (4.1)
#  2. new results dir seeded from the baseline run (cp -a keeps mtimes, so with
#     --rerun-triggers mtime nothing upstream reruns), minus everything that depends on the
#     model set or on external features
#  3. configs: cluster table with external_features_config, model table, per-model feature tables
# Run from the A2 dir.
set -euo pipefail
export PYTHONUNBUFFERED=1
export A2=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
V3=/oak/stanford/groups/engreitz/Users/sheth/scE2G_temp/scE2G_main/models/multiome_powerlaw_v3/feature_table.tsv
BASE=$SCRATCH/sce2g_runs/2026_1002_k562_xu
SRC=${SRC:-tobias}
case $SRC in
    tobias) ELEMENTS=results/2026_1003_footprint_p300_gate.element_features.tsv.gz ;;
    print) ELEMENTS=results/2026_1003_footprint_p300_gate.print.element_features.tsv.gz ;;
    *) echo "SRC must be tobias or print" >&2; exit 1 ;;
esac
RUN=$SCRATCH/sce2g_runs/${RUN_DATE:-2026_1003}_k562_xu_featB_${SRC}
MODELS=config/models_${SRC}  # per-source model dirs (TOBIAS run of 2026-10-03 used config/models/)
EXT=$SCRATCH/sce2g_external
cd $A2
mkdir -p $EXT

# 1. pairs file
PAIRS=$EXT/featB_${SRC}_pairs.tsv.gz
$A2/.pixi/envs/tobias/bin/python scripts/4.1.make_element_feature_pairs.py \
    --elements $ELEMENTS --out $PAIRS

# 2. results dir
[[ -e $RUN ]] && { echo "$RUN exists; not overwriting" >&2; exit 1; }
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
ls $RUN/K562_Xu_pl

# 3. configs
cat > config/external_features_featB_${SRC}.tsv <<EOF
input_col	source_col	aggregate_function	join_by	source_file
erythroid_n	erythroid_n	max	overlap	$PAIRS
erythroid_fp_sum	erythroid_fp_sum	max	overlap	$PAIRS
EOF
awk -F'\t' -v OFS='\t' -v efc="$A2/config/external_features_featB_${SRC}.tsv" 'NR==1 {print $0, "external_features_config"; next} {print $0, efc}' \
    config/cell_clusters_K562_Xu.tsv > config/cell_clusters_K562_Xu_featB_${SRC}.tsv
for m in multiome_n6_ery_motif multiome_n6_ery_fp; do
    mkdir -p $MODELS/$m
    cp $V3 $MODELS/$m/feature_table.tsv
    sed -i -e '$a\' $MODELS/$m/feature_table.tsv  # v3 table lacks a final newline
    printf 'eryMotifCount\terythroid_n\tNA\tmax\t0\terythroid motif sites in E\n' >> $MODELS/$m/feature_table.tsv
done
printf 'eryFootprintSum\terythroid_fp_sum\tNA\tmax\t0\terythroid %s footprint sum in E\n' "$SRC" >> $MODELS/multiome_n6_ery_fp/feature_table.tsv
{
    head -n 1 config/models_baseline.tsv
    grep -P '^multiome_powerlaw_n6\t' config/models_baseline.tsv
    for m in multiome_n6_ery_motif multiome_n6_ery_fp; do
        printf '%s\tK562_Xu_pl\t\t%s\tFalse\t\n' $m $A2/$MODELS/$m/feature_table.tsv
    done
} > config/models_featB_${SRC}.tsv
sed -e "s#^results_dir: .*#results_dir: \"$RUN\"#" \
    -e "s#^cell_clusters: .*#cell_clusters: \"$A2/config/cell_clusters_K562_Xu_featB_${SRC}.tsv\"#" \
    -e "s#^model_config: .*#model_config: \"$A2/config/models_featB_${SRC}.tsv\"#" \
    config/config_training.yaml > config/config_training_featB_${SRC}.yaml
for k in "results_dir: \"$RUN\"" "cell_clusters_K562_Xu_featB_${SRC}.tsv" "models_featB_${SRC}.tsv"; do
    grep -qF "$k" config/config_training_featB_${SRC}.yaml || { echo "config patch missed: $k" >&2; exit 1; }
done
if grep -lF '${SRC}' config/*featB_${SRC}* $MODELS/*/feature_table.tsv; then echo "unexpanded \${SRC} in configs" >&2; exit 1; fi
grep -qF "external_features_featB_${SRC}.tsv" config/cell_clusters_K562_Xu_featB_${SRC}.tsv
echo "== models"; cat config/models_featB_${SRC}.tsv
echo "== feature table (fp)"; cat $MODELS/multiome_n6_ery_fp/feature_table.tsv
echo "setup complete"
