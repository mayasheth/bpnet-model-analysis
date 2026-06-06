qsjob -j motif_prepare --cpu -m 30G -t 6:00:00 --env bpnet_37 \
  "ml cuda/11.1.1 cudnn/8.1.1.33; \
   python $OAK/Users/sheth/EP300_BPNet/scripts/6.2.motif_pairs.py \
      --mode prepare \
      --analysis-id motif_pairs_v1 \
      --narrow-peak-type all \
      --model-type p300_v1 \
      --out-dir $OAK/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/motif_spacing \
      --spacing-min 0 --spacing-max 40 --spacing-step 9"

NCHUNKS=20
for idx in $(seq 0 $((NCHUNKS-1))); do
  qsjob -j motif_chunk${idx} --cpu -m 80G -t 12:00:00 --env bpnet_37 \
    "ml cuda/11.1.1 cudnn/8.1.1.33; \
     python $OAK/Users/sheth/EP300_BPNet/scripts/6.2.motif_pairs.py \
        --mode chunk \
        --chunk-index $idx \
        --chunk-n-chunks $NCHUNKS \
        --analysis-id motif_pairs_v1 \
        --narrow-peak-type all \
        --model-type p300_v1 \
        --out-dir $OAK/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/motif_spacing \
        --chunk-size 5 \
        --spacing-min 0 --spacing-max 40 --spacing-step 9 --use-existing-singles"
done

qsjob -j motif_collate --cpu -m 20G -t 2:00:00 --env bpnet_37 \
  "ml cuda/11.1.1 cudnn/8.1.1.33; \
   python $OAK/Users/sheth/EP300_BPNet/scripts/6.2.motif_pairs.py \
      --mode collate \
      --analysis-id motif_pairs_v1 \
      --narrow-peak-type all \
      --model-type p300_v1 \
      --out-dir $OAK/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/motif_spacing"