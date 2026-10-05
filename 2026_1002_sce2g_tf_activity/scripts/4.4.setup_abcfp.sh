#!/bin/bash
#SBATCH -p engreitz,normal
#SBATCH -t 3:00:00
#SBATCH --mem=64G
#SBATCH -c 4
#SBATCH -o log/setup_abcfp.%j.txt
#SBATCH --job-name=setup_abcfp
# 4.4 Set up a scE2G run with footprints inside the activity term (4.3): each model is the
# multiome v3 feature set with ARC.E2G.Score swapped for a footprint-weighted ARC-E2G score
# (TOBIAS or PRINT erythroid footprints, lam 1 or 4). Same seeding as 4.0. Run from the A2 dir.
set -euo pipefail
export PYTHONUNBUFFERED=1
export A2=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_1002_sce2g_tf_activity
V3=/oak/stanford/groups/engreitz/Users/sheth/scE2G_temp/scE2G_main/models/multiome_powerlaw_v3/feature_table.tsv
BASE=$SCRATCH/sce2g_runs/2026_1002_k562_xu
RUN=$SCRATCH/sce2g_runs/${RUN_DATE:-2026_1004}_k562_xu_abcfp
MODELS=config/models_abcfp
TAGS=(tobias_l1 tobias_l4 print_l1 print_l4)
EXT=$SCRATCH/sce2g_external
cd $A2
mkdir -p $EXT
[[ -e $RUN ]] && { echo "$RUN exists; not overwriting" >&2; exit 1; }

# 1. pairs file
PAIRS=$EXT/abcfp_pairs.tsv.gz
$A2/.pixi/envs/tobias/bin/python scripts/4.3.make_abc_fp_pairs.py --lams 1,4 --out $PAIRS
HDR=$( (zcat $PAIRS || true) | head -n 1 | tr '\t' '\n')
for t in "${TAGS[@]}"; do
    grep -qx "ARCfp_$t" <<< "$HDR" || { echo "ARCfp_$t missing from $PAIRS" >&2; exit 1; }
done

# 2. results dir
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

# 3. configs
EFC=config/external_features_abcfp.tsv
printf 'input_col\tsource_col\taggregate_function\tjoin_by\tsource_file\n' > $EFC
for t in "${TAGS[@]}"; do printf 'ARCfp_%s\tARCfp_%s\tmean\toverlap\t%s\n' $t $t $PAIRS >> $EFC; done
awk -F'\t' -v OFS='\t' -v efc="$A2/$EFC" 'NR==1 {print $0, "external_features_config"; next} {print $0, efc}' \
    config/cell_clusters_K562_Xu.tsv > config/cell_clusters_K562_Xu_abcfp.tsv
{
    head -n 1 config/models_baseline.tsv
    grep -P '^multiome_powerlaw_n6\t' config/models_baseline.tsv
} > config/models_abcfp.tsv
for t in "${TAGS[@]}"; do
    m=multiome_n6_arcfp_$t
    mkdir -p $MODELS/$m
    awk -F'\t' -v OFS='\t' -v t=$t '$1 == "ARC.E2G.Score" {$1 = "ARCfp_" t; $2 = "ARCfp_" t; $6 = "ARC-E2G score, footprint-weighted activity (" t ")"} {print}' \
        $V3 > $MODELS/$m/feature_table.tsv
    grep -qP "^ARCfp_$t\tARCfp_$t\t" $MODELS/$m/feature_table.tsv || { echo "ARC swap missed in $m" >&2; exit 1; }
    ! grep -q '^ARC.E2G.Score' $MODELS/$m/feature_table.tsv || { echo "ARC.E2G.Score still in $m" >&2; exit 1; }
    [[ $(grep -c . $MODELS/$m/feature_table.tsv) -eq $(grep -c . $V3) ]] || { echo "row count changed in $m" >&2; exit 1; }
    printf '%s\tK562_Xu_pl\t\t%s\tFalse\t\n' $m $A2/$MODELS/$m/feature_table.tsv >> config/models_abcfp.tsv
done
# scE2G keeps one feature table per run keyed by feature name, so a shared name would give every
# model the same column (2026-10-04 run: four identical models)
[[ $(cut -f1 $MODELS/multiome_n6_arcfp_*/feature_table.tsv | grep -c '^ARCfp') -eq ${#TAGS[@]} ]] &&
    [[ $(cut -f1 $MODELS/multiome_n6_arcfp_*/feature_table.tsv | grep '^ARCfp' | sort -u | wc -l) -eq ${#TAGS[@]} ]] \
    || { echo "ARCfp feature names are not unique across models" >&2; exit 1; }
sed -e "s#^results_dir: .*#results_dir: \"$RUN\"#" \
    -e "s#^cell_clusters: .*#cell_clusters: \"$A2/config/cell_clusters_K562_Xu_abcfp.tsv\"#" \
    -e "s#^model_config: .*#model_config: \"$A2/config/models_abcfp.tsv\"#" \
    config/config_training.yaml > config/config_training_abcfp.yaml
for k in "results_dir: \"$RUN\"" "cell_clusters_K562_Xu_abcfp.tsv" "models_abcfp.tsv"; do
    grep -qF "$k" config/config_training_abcfp.yaml || { echo "config patch missed: $k" >&2; exit 1; }
done
grep -qF "external_features_abcfp.tsv" config/cell_clusters_K562_Xu_abcfp.tsv
echo "== models"; cat config/models_abcfp.tsv
echo "== external features"; cat $EFC
echo "== feature table (tobias_l1)"; cat $MODELS/multiome_n6_arcfp_tobias_l1/feature_table.tsv
echo "setup complete"
