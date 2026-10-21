#!/bin/bash
# Score Bolero idea 2 exactly as pre-registered (decisions [2026-10-21]): FILM (1.43 FiLM on the
# panel-13 Corgi, K562 row), its matched BASE (no conditioning) and CLOSEST (FILM with GM12878's
# row), 3 runs each, on the K562 CRISPR benchmark through ABC and rE2G. Copied from
# panelmean_score_launch.sh (F-067). Corgi 5-cell (F-067) is not re-predicted: its ABC outputs in
# 2026_1013_p300_panelmean and its rE2G models in 2026_0923_p300_features are reused, so
# BASE - Corgi5 is paired on the same pairs.
#
# One file, several stages, so no command is nested inside a --wrap string:
#   main      (default) gate on the 30 folds, 4.26 stats, 9 x 4.4, then submit the stages below
#   crispr    after ABC: check the region-calling gate ran, 4.6 --arms film13, 4.7
#   abcboot   after crispr: 4.17, pre-registered pairs and floor context, pair and gene units
#   re2g      after the predictions: wait for no re2g_driver, 4.32 --stage-new, 4.33
#   re2gboot  after the rE2G driver: 4.34, pair and gene units
# Submit: sbatch -p engreitz,normal -c 1 --mem 4G -t 2:00:00 -J film13_score_launch \
#   -o log/film13_score_launch.%j.txt --dependency=afterok:<6 training arrays> log/film13_score_launch.sh
set -euo pipefail
export PYTHONUNBUFFERED=1
P=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2026_0824_H3K27ac_model
D=/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet
PY=$D/.pixi/envs/multimodal/bin/python
CC=/oak/stanford/groups/engreitz/Users/sheth/CRISPR_comparison_v3/CRISPR_comparison
E2G=/oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G
ABCD=/oak/stanford/groups/engreitz/Users/sheth/ABC_working/ABC-Enhancer-Gene-Prediction
SELF=$P/log/film13_score_launch.sh
cd "$P"
# 4.4 is submitted with --export=ALL; nothing from the submitting shell may leak into it
unset COND_CELL PANEL_MEAN ACC_MEAN ACC_STD ACC_BW OUT_PREFIX OFFSET_MODEL MOTIF_BW SIGNAL_BW MODEL_ROOT

FILMM=models/p300_multicell_panelmean_film_chromvar13pc_drop0.1_panel13_holdK562_hw500_clw10
BASEM=models/p300_multicell_panelmean_panel13_holdK562_hw500_clw10
FILMA=multicell_panelmean_film_chromvar13pc_panel13_holdK562
BASEA=multicell_panelmean_panel13_holdK562
ABCRUN=2026_1021_p300_film13
OUT=results/2026_1021_film13
SEEDS=("" _seed1 _seed2)
ABCARMS=()
for SF in "${SEEDS[@]}"; do
    ABCARMS+=(p300pred_${FILMA}_rowK562${SF} p300pred_${FILMA}_rowGM12878${SF} p300pred_${BASEA}${SF})
done
STAGE=${1:-main}
echo "stage $STAGE, job ${SLURM_JOB_ID:-none}, $(date)"

case "$STAGE" in
main)
    for SF in "${SEEDS[@]}"; do
        for f in 0 1 2 3 4; do
            for m in $FILMM$SF/fold$f $BASEM$SF/fold$f; do
                [[ -f $m/training_complete.json && -f $m/panel_mean_stats.json ]] || { echo "ERROR: incomplete $m" >&2; exit 1; }
            done
            [[ -f $FILMM$SF/fold$f/cell_cond.json ]] || { echo "ERROR: no cell_cond.json in $FILMM$SF/fold$f" >&2; exit 1; }
        done
        bw=data/abc_predicted/predp300mc_multicell_panelmean_holdK562${SF}.bw
        [[ -f $bw ]] || { echo "ERROR: Corgi 5-cell bigwig $bw missing (reused, not re-predicted)" >&2; exit 1; }
    done
    echo "all 30 folds complete; Corgi 5-cell bigwigs present"

    REG=$ABCD/results/2026_0721_h3k27ac_counting_comparison/K562_ATAC_only/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed
    K5=$($PY -c "import json; print(json.load(open('config/p300_panel_13.json'))['K562']['accessibility_bw'])")
    K5OLD=$($PY -c "import json; print(json.load(open('config/p300_panel.json'))['K562']['accessibility_bw'])")
    [[ "$K5" == "$K5OLD" ]] || { echo "ERROR: panel 13 K562 ATAC $K5 is not F-067's $K5OLD" >&2; exit 1; }
    ST=$($PY scripts/4.26.acc_norm_stats.py --regions "$REG" --accessibility-bw "$K5" | tail -1)
    M=$(echo "$ST" | awk '{print $2}'); S=$(echo "$ST" | awk '{print $4}')
    [[ -n "$M" && -n "$S" ]] || { echo "ERROR: could not parse 4.26 output: $ST" >&2; exit 1; }
    echo "K562 input stats: $ST (F-067 used 3.7232/1.1428)"

    IDS=
    for SF in "${SEEDS[@]}"; do
        for ROW in K562 GM12878; do
            j=$(COND_CELL=$ROW PANEL_MEAN=model ACC_BW="$K5" ACC_MEAN=$M ACC_STD=$S OUT_PREFIX=predp300mc \
                sbatch --parsable --export=ALL scripts/4.4.submit_predict_h3k27ac_cpu.sh \
                ${FILMA}_row${ROW}${SF} $FILMM$SF multimodal)
            IDS=$IDS:$j; echo "pred FILM row $ROW${SF} $j"
        done
        j=$(PANEL_MEAN=model ACC_BW="$K5" ACC_MEAN=$M ACC_STD=$S OUT_PREFIX=predp300mc \
            sbatch --parsable --export=ALL scripts/4.4.submit_predict_h3k27ac_cpu.sh \
            ${BASEA}${SF} $BASEM$SF multimodal)
        IDS=$IDS:$j; echo "pred BASE${SF} $j"
    done

    AB=$(ABC_ARMS=film13 ABC_RESULTS_DIR=results/$ABCRUN ABC_CFG=config/mine/config_p300_film13.yaml \
         sbatch --parsable --export=ALL --dependency=afterok${IDS} scripts/4.5.submit_abc.sh)
    C1=$(ABCJOB=$AB sbatch --parsable -p engreitz,normal -c 1 --mem 4G -t 24:00:00 -J crispr_film13_launch \
         -o log/crispr_film13_launch.%j.txt --dependency=afterok:$AB --export=ALL "$SELF" crispr)
    B=$(sbatch --parsable -p engreitz,normal -c 2 --mem 16G -t 2:00:00 -J film13_boot \
        -o log/film13_boot.%j.txt --dependency=afterok:$C1 --export=ALL "$SELF" abcboot)
    RE=$(sbatch --parsable -p engreitz,normal -c 1 --mem 4G -t 24:00:00 -J re2g_film13_wait \
         -o log/re2g_film13_wait.%j.txt --dependency=afterok${IDS} --export=ALL "$SELF" re2g)
    echo "preds=${IDS#:} abc=$AB crispr=$C1 boot=$B re2g_wait=$RE"
    ;;

crispr)
    LOG=log/abcrun.${ABCJOB:?}.txt
    grep -q "^region calling is correctly skipped" "$LOG" || { echo "ERROR: the 4.5 region-calling gate did not pass in $LOG" >&2; exit 1; }
    JULYREG=$ABCD/results/2026_0721_h3k27ac_counting_comparison/K562_ATAC_only/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed
    for arm in "${ABCARMS[@]}"; do
        cmp -s "$JULYREG" "$ABCD/results/$ABCRUN/$arm/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed" \
            || { echo "ERROR: $arm candidate regions differ from the July set" >&2; exit 1; }
    done
    echo "all ${#ABCARMS[@]} new arms score the July candidate-region set"
    $PY scripts/4.6.setup_crispr_comparison.py --arms film13
    CC_CFG=$CC/config/config_p300_film13.yml sbatch --wait --export=ALL scripts/4.7.submit_crispr_comparison.sh
    ;;

abcboot)
    MERGED=$CC/workflow/results/$ABCRUN/expt_pred_merged_annot.txt.gz
    g() { echo "--group $1=$2.ABC.Score,$2_seed1.ABC.Score,$2_seed2.ABC.Score"; }
    GR="$(g film p300pred_${FILMA}_rowK562) $(g closest p300pred_${FILMA}_rowGM12878) $(g base p300pred_${BASEA}) $(g corgi5 p300pred_multicell_panelmean_holdK562)"
    # pre-registered: primary FILM - BASE; descriptive FILM - CLOSEST and BASE - Corgi5
    PRE="--pair film base --pair film closest --pair base corgi5"
    # context, not named in the pre-registration: each arm against the ATAC-only floor
    FLR="--pair film K562_ATAC_only.ABC.Score --pair base K562_ATAC_only.ABC.Score"
    $PY scripts/4.17.paired_auprc_bootstrap.py "$MERGED" $GR $PRE --out-tsv ${OUT}_auprc_pairs.tsv
    $PY scripts/4.17.paired_auprc_bootstrap.py "$MERGED" $GR $PRE --cluster gene --out-tsv ${OUT}_auprc_pairs.gene.tsv
    $PY scripts/4.17.paired_auprc_bootstrap.py "$MERGED" $GR $FLR --out-tsv ${OUT}_auprc_floor.tsv
    $PY scripts/4.17.paired_auprc_bootstrap.py "$MERGED" $GR $FLR --cluster gene --out-tsv ${OUT}_auprc_floor.gene.tsv
    ;;

re2g)
    # rE2G drivers share one workdir: only one at a time
    while squeue -u "$USER" -h -o %j | grep -q '^re2g_driver'; do sleep 300; done
    cd "$E2G"
    $PY $P/scripts/4.32.setup_re2g_p300.py --stage-new
    n=$(grep -cE '^ATAC_p300mc(film13|pm13)(seed[12])?_powerlaw	' config/mine/config_datasets_K562_p300.tsv || true)
    [[ $n -eq 6 ]] || { echo "ERROR: 4.32 wrote $n film13 datasets, expected 6" >&2; exit 1; }
    cd "$P"
    DR=$(sbatch --parsable scripts/4.33.submit_re2g_training.sh)
    echo "driver $DR"
    RB=$(sbatch --parsable -p engreitz,normal -c 2 --mem 16G -t 2:00:00 -J film13_re2g_boot \
         -o log/film13_re2g_boot.%j.txt --dependency=afterok:$DR --export=ALL "$SELF" re2gboot)
    echo "re2g boot $RB"
    ;;

re2gboot)
    RD=results/2026_0923_p300_features
    g() { echo "--group $1=ATAC_p300mc$2_powerlaw_all,ATAC_p300mc$2seed1_powerlaw_all,ATAC_p300mc$2seed2_powerlaw_all"; }
    RB="$(g film film13) $(g base pm13) $(g corgi5 panelmean) --pair film base --pair base corgi5"
    $PY scripts/4.34.re2g_paired_bootstrap.py $E2G/$RD $RB --out-tsv ${OUT}_re2g_pairs.tsv
    $PY scripts/4.34.re2g_paired_bootstrap.py $E2G/$RD $RB --cluster gene --out-tsv ${OUT}_re2g_pairs.gene.tsv
    ;;

*) echo "unknown stage $STAGE" >&2; exit 1 ;;
esac
