# --- thread/env caps must be FIRST ---
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("TF_NUM_INTRAOP_THREADS", "2")
os.environ.setdefault("TF_NUM_INTEROP_THREADS", "2")
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")

import sys, time, socket, platform, logging, argparse
from logging.handlers import RotatingFileHandler
from functools import partial
import numpy as np
import pandas as pd
import itertools
import multiprocessing as mp
import pickle

from motif_exp_utils import (
    get_model, get_shuffled_peak_sequences, one_hot_encode, pattern_to_string,
    make_model_prediction, get_all_patterns,
    insert_motifs_with_orientation_general, generate_motif_pairs
)

from plot_motif_pairs import plot_motif_pairs

# ---------- argument parsing ----------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze motif spacing effects on model predictions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Core analysis parameters
    parser.add_argument("--analysis-id",
                        help="Unique identifier for this analysis run")
    
    # Spacing parameters
    parser.add_argument("--spacing-min", type=int, default=0,
                        help="Minimum spacing between motifs (bp)")
    parser.add_argument("--spacing-max", type=int, default=0,
                        help="Maximum spacing between motifs (bp)")
    parser.add_argument("--spacing-step", type=int, default=1,
                        help="Step size for spacing values")
    
    # Sequence parameters
    parser.add_argument("--num-bg", type=int, default=1000,
                        help="Number of background sequences to generate")
    parser.add_argument("--seq-len", type=int, default=2114,
                        help="Length of sequences for model input")
    
    # File paths
    parser.add_argument("--reference-genome", type = str, default = '/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/hg38.fa',
                        help="Path to reference genome FASTA file")
    parser.add_argument("--narrow-peak-path",
                        help="Path to narrowPeak or BED file with peak regions")
    parser.add_argument("--narrow-peak-type",
                        choices=["all", "p300_peaks"])
    parser.add_argument("--model-base-path",
                        help="Base path to model directories (should contain fold0, fold1, etc.)")
    parser.add_argument("--model-subpath",
                        help="Subpath within each fold directory to the model")
    parser.add_argument("--model-type",
                        choices=["p300_v1", "p300_v2", "DNase", "GATA1", "GATA2"])
    parser.add_argument("--num-folds", type=int, default=5,
                        help="Number of model folds to use")
    parser.add_argument("--out-dir",
                        help="Output directory for results and plots")
    
    # Processing parameters
    parser.add_argument("--n-proc", type=int, default=20,
                        help="Number of parallel processes to use")
    parser.add_argument("--chunk-size", type=int, default=10,
                        help="Number of spacing values per chunk for parallel processing")
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="Logging level")
    parser.add_argument("--random-seed", type=int, default=17,
                        help="Random seed for reproducibility")
    parser.add_argument("--mode", choices=["prepare", "chunk", "collate"], required=True,
        help="Which pipeline stage: 'prepare' for singles+background, 'chunk' for motif pair chunk, 'collate' to merge results and plot.")
    parser.add_argument("--chunk-index", type=int, default=None,
        help="Chunk index for motif pairs to process (required if --mode chunk)")
    parser.add_argument("--chunk-n-chunks", type=int, default=None,
        help="Total number of chunks (required if --mode chunk)")

    # Plot parameters
    parser.add_argument("--heatmap-cmap", default="RdBu_r",
                        help="Colormap for heatmap plots")
    parser.add_argument("--no-plots", action="store_true",
                        help="Skip generating plots (only save data)")
    parser.add_argument("--use-existing-singles", action="store_true",
                        help="If single motif results exist at the expected path, don't regenerate them")
    
    args = parser.parse_args()
    
    # Validate arguments
    if (args.narrow_peak_path is None) + (args.narrow_peak_type is None) != 1:
        parser.error("Must provide exactly one peak specification")
    if (args.model_base_path is None or args.model_subpath is None)  + (args.model_type is None) != 1:
        parser.error("Must provide exactly one model specification")
    if args.spacing_min < 0:
        parser.error("--spacing-min must be non-negative")
    if args.spacing_max < args.spacing_min:
        parser.error("--spacing-max must be >= --spacing-min")
    if args.num_bg <= 0:
        parser.error("--num-bg must be positive")
    if args.seq_len <= 0:
        parser.error("--seq-len must be positive")
    if args.n_proc <= 0:
        parser.error("--n-proc must be positive")

    # Set default paths
    if args.narrow_peak_type:
        if args.narrow_peak_type == "all":
            args.narrow_peak_path = '/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/reference/K562_DNase_candidate_elements.narrowPeak'
        elif args.narrow_peak_type == "p300_peaks":
            args.narrow_peak_path = '/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/reference/ENCSR000EGE_peaks_inliers.narrowPeak'
        else:
            parser.error("Invalid value provided for --narrow-peak-type")

    if args.model_type:
        if args.model_type == "p300_v1":
            args.model_base_path = '/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/models/release_run_1'
            args.model_subpath = 'ENCSR000EGE/ENCSR000EGE_split000'
        elif args.model_type == "p300_v2":
            args.model_base_path = '/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0703_retrain_p300_model/models'
            args.model_subpath = 'model_split000'
        elif args.model_type == "DNase":
            args.model_base_path = '/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/K562_DNase_ChromBPNet/models'
            args.model_subpath = 'model.chrombpnet_nobias.h5' # 'chrombpnet_wo_bias'
        elif args.model_type in ["GATA1", "GATA2"]:
            args.model_base_path = f'/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/K562_{args.model_type}_BPNet/models'
            args.model_subpath = 'model_split000'
        else:
            parser.error("Invalid value provided for --model-type")

    args.motif_dict = {
        "GATA": "GATAA",
        "AP1_2": "TGACTCA",
        "CREB_ATF_3": "TGAGTCA", # actually AP1 (FOS-JUN)
        "STAT": "TTCCNGGAA",
        "ETS": "CGGAAG",
        "CREB_ATF_1": "TGACGTCA",
        "Ebox_CACGTG": "CACGTG",       
        "CTCF": "CCNNNAGGGGGCG",
        "Grepeats": "GGGGGNGGGGG",
        "NFY": "CCAAT",
        "NRF1": "GCGCATGCGC",
        "NFI": "TGGCNNNNNGCCA",
        "REST": "CAGCACCNNGGACAG",
        "RFX_1": "GTTGCCATGGCAAC",
        "E2F": "TTTCCCGCCAAA",
        "TATAbox": "TATAAAA"
    }
    
    return args


# ---------- chunking utilities ----------
def chunk_motif_pairs(motif_pairs, chunk_size=10):
    """
    Split list of motif pairs into chunks of specified size.
    
    Args:
        motif_pairs: List of (motif1_name, motif2_name) tuples
        chunk_size: Number of pairs per chunk
    
    Returns:
        List of chunks, where each chunk is a list of motif pairs
    """
    chunks = []
    for i in range(0, len(motif_pairs), chunk_size):
        chunk = motif_pairs[i:i + chunk_size]
        chunks.append(chunk)
    return chunks

def estimate_chunk_size(n_pairs, 
                       n_workers,
                       target_chunks_per_worker=2,
                       min_chunk_size=5,
                       max_chunk_size=20):
    """
    Estimate optimal chunk size based on number of pairs and workers.
    
    Args:
        n_pairs: Total number of motif pairs
        n_workers: Number of parallel workers
        target_chunks_per_worker: Target number of chunks per worker
        min_chunk_size: Minimum pairs per chunk
        max_chunk_size: Maximum pairs per chunk
    
    Returns:
        Optimal chunk size
    """
    target_n_chunks = n_workers * target_chunks_per_worker
    chunk_size = max(min_chunk_size, n_pairs // target_n_chunks)
    chunk_size = min(chunk_size, max_chunk_size)
    return chunk_size

# ---------- logging helpers ----------
def sizeof_fmt(num, suffix="B"):
    for unit in ["","K","M","G","T","P","E","Z"]:
        if abs(num) < 1024.0: return f"{num:3.1f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f}Y{suffix}"

def setup_logger(log_dir, analysis_id, level=None):
    os.makedirs(log_dir, exist_ok=True)
    level_name = level or "INFO"
    level_val = getattr(logging, level_name.upper(), logging.INFO)
    logger = logging.getLogger(analysis_id); logger.setLevel(level_val); logger.propagate = False
    fmt = logging.Formatter("%(asctime)s | %(levelname)-8s | pid=%(process)d | %(name)s | %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")
    sh = logging.StreamHandler(sys.stdout); sh.setLevel(level_val); sh.setFormatter(fmt); logger.addHandler(sh)
    fh = RotatingFileHandler(os.path.join(log_dir, f"{analysis_id}.log"), maxBytes=20_000_000, backupCount=3)
    fh.setLevel(level_val); fh.setFormatter(fmt); logger.addHandler(fh)
    return logger

def worker_init(analysis_id, level="INFO"):
    log = logging.getLogger(analysis_id)
    if not log.handlers:
        log.setLevel(getattr(logging, level.upper(), logging.INFO))
        fmt = logging.Formatter("%(asctime)s | %(levelname)-8s | pid=%(process)d | %(name)s | %(message)s",
                                datefmt="%Y-%m-%d %H:%M:%S")
        sh = logging.StreamHandler(sys.stdout); sh.setFormatter(fmt); log.addHandler(sh)

# ---------- process chunk ----------
def process_motif_pair_chunk(pair_info, shared_data):
    """
    Process a chunk of motif pairs across all spacings and orientations.
    
    Args:
        pair_info: list of tuples [(name1, name2, seq1, seq2), ...] for chunk
                   OR single tuple (name1, name2, seq1, seq2) for backward compatibility
        shared_data: dict with bg_encoded, model_paths, baseline_preds, etc.
                     Should also contain 'chunk_idx' for logging
    """
    t0 = time.perf_counter()
    pid = os.getpid()
    
    analysis_id = shared_data['analysis_id']
    log = logging.getLogger(analysis_id)
    
    # Handle both single pair (backward compatible) and chunk of pairs
    if isinstance(pair_info, tuple) and len(pair_info) == 4:
        # Single pair: convert to list for uniform processing
        pair_list = [pair_info]
        chunk_idx = 0
    else:
        # Chunk of pairs
        pair_list = pair_info
        chunk_idx = shared_data.get('chunk_idx', 0)
    
    bg_encoded = shared_data['bg_encoded']
    baseline_preds = shared_data['baseline_preds']
    individual_motif_preds = shared_data['individual_motif_preds']
    model_paths = shared_data['model_paths']
    seq_len = shared_data['seq_len']
    spacing_vals = shared_data['spacing_vals']
    
    # Log chunk info
    pair_names = [f"{n1}+{n2}" for n1, n2, _, _ in pair_list]
    log.info(
        f"[worker pid={pid}] Starting chunk {chunk_idx} | "
        f"pairs={pair_names} | n_pairs={len(pair_list)}"
    )
    
    # Load models ONCE for entire chunk (this is the key optimization)
    t_load = time.perf_counter()
    model_info_list = [get_model(mp) for mp in model_paths]
    load_time = time.perf_counter() - t_load
    
    n_conditions_per_pair = len(list(get_all_patterns(2))) * len(spacing_vals)
    total_conditions = len(pair_list) * n_conditions_per_pair
    
    log.info(
        f"[worker pid={pid}] Loaded {len(model_info_list)} models in {load_time:.2f}s | "
        f"chunk={chunk_idx} | n_pairs={len(pair_list)} | "
        f"total_conds={total_conditions}"
    )
    
    # Process all pairs in chunk
    all_results = []
    
    for pair_idx, (name1, name2, seq1, seq2) in enumerate(pair_list):
        pair_start = time.perf_counter()
        
        # Build all conditions for this pair
        conditions = []
        enc_list = []
        sizes = []
        
        # Loop over all orientations
        for pattern in get_all_patterns(2):
            pstr = pattern_to_string(pattern)
            
            # Loop over all spacings
            for s in spacing_vals:
                t_ins = time.perf_counter()
                enc = insert_motifs_with_orientation_general(
                    bg_encoded, 
                    n=2, 
                    spacing=s, 
                    motif_seq=[seq1, seq2],
                    orientation_pattern=pattern, 
                    seq_len=seq_len
                )
                ins_s = time.perf_counter() - t_ins
                
                conditions.append((name1, name2, seq1, seq2, pstr, s, enc.shape[0], ins_s))
                enc_list.append(enc)
                sizes.append(enc.shape[0])
        
        if not enc_list:
            log.warning(f"[worker pid={pid}] No conditions for pair {name1}+{name2}")
            continue
        
        # Concatenate all encodings for this pair
        big_encoded = np.concatenate(enc_list, axis=0)
        offsets = np.cumsum([0] + sizes)
        
        log.info(
            f"[worker pid={pid}] Processing pair {pair_idx+1}/{len(pair_list)}: "
            f"{name1}+{name2} | {len(conditions)} conditions | "
            f"big_encoded shape={big_encoded.shape}"
        )
        
        # Predict across all models (models already loaded)
        model_full_preds = []
        for mi, model_info in enumerate(model_info_list):
            t_pm = time.perf_counter()
            preds = make_model_prediction(big_encoded, model_info)
            model_full_preds.append(preds)
            pred_time = time.perf_counter() - t_pm
            log.info(f"[worker pid={pid}] model{mi} predict len={sizes[0]} took {pred_time:.2f}s")
        
        # Slice predictions and compute metrics
        for idx, (n1, n2, s1, s2, pstr, s, enc_len, ins_s) in enumerate(conditions):
            i0, i1 = offsets[idx], offsets[idx+1]
            
            for mi in range(len(model_info_list)):
                preds_cond = model_full_preds[mi][i0:i1]
                base = baseline_preds[mi]
                
                # Size check
                if preds_cond.shape[0] != base.shape[0]:
                    log.warning(f"[pid={pid}] Size mismatch for {n1}+{n2} spacing={s} pattern={pstr}")
                
                # Log2 fold changes
                log2_vs_base = float(np.mean(preds_cond - base[:preds_cond.shape[0]]) / np.log(2))
                
                # Compare to individual motifs
                individual1 = individual_motif_preds[n1][mi]
                individual2 = individual_motif_preds[n2][mi]
                log2_vs_ind1 = float(np.mean(preds_cond - individual1[:preds_cond.shape[0]]) / np.log(2))
                log2_vs_ind2 = float(np.mean(preds_cond - individual2[:preds_cond.shape[0]]) / np.log(2))
                
                # Compute synergy metric
                expected_additive = (individual1[:preds_cond.shape[0]] + 
                                   individual2[:preds_cond.shape[0]] - 
                                   base[:preds_cond.shape[0]])
                log2_synergy = float(np.mean(preds_cond - expected_additive) / np.log(2))
                
                result = {
                    "motif1_name": n1,
                    "motif2_name": n2,
                    "motif1_seq": s1,
                    "motif2_seq": s2,
                    "orientation_pattern": pstr,
                    "spacing": s,
                    "model_fold": mi,
                    "log2_fc_vs_baseline": log2_vs_base,
                    f"log2_fc_vs_{n1}_alone": log2_vs_ind1,
                    f"log2_fc_vs_{n2}_alone": log2_vs_ind2,
                    "log2_synergy": log2_synergy,
                    "mean_prediction": float(np.mean(preds_cond)),
                    "n_sequences": int(preds_cond.shape[0])
                }
                
                all_results.append(result)
        
        pair_time = time.perf_counter() - pair_start
        log.info(
            f"[worker pid={pid}] Finished pair {pair_idx+1}/{len(pair_list)}: "
            f"{name1}+{name2} in {pair_time:.2f}s | "
            f"{len(conditions) * len(model_info_list)} results"
        )
    
    chunk_time = time.perf_counter() - t0
    log.info(
        f"[worker pid={pid}] Finished chunk {chunk_idx} | "
        f"pairs={pair_names} | "
        f"time={chunk_time:.2f}s | "
        f"n_results={len(all_results)}"
    )
    
    return all_results

def save_individual_motif_predictions(pred_array, pred_dict, baseline_preds, args, logger):
    # convert array to df, compute log2fc
    all_individual_preds_array = np.stack(pred_array)  # Shape: (n_motifs, n_models, n_sequences)
    pseudocount = 1e-6
    log2fc_array = np.log2(
        (all_individual_preds_array + pseudocount) / (baseline_preds[np.newaxis, :, :] + pseudocount)
    )

    rows = []
    for motif_idx, (motif_name, motif_seq) in enumerate(args.motif_dict.items()):
        for model_idx in range(log2fc_array.shape[1]):
            for seq_idx in range(log2fc_array.shape[2]):
                rows.append({
                    'motif_name': motif_name,
                    'motif_seq': motif_seq,
                    'model_num': model_idx,
                    'seq_num': seq_idx,
                    'log2_fc_vs_baseline': log2fc_array[motif_idx, model_idx, seq_idx]
                })
    df = pd.DataFrame(rows)

    print("\nResults summary:")
    print(df.head(10))
    print(f"Shape: {df.shape}")

    # save raw results
    output_file = os.path.join(args.out_dir, f'individual_motifs.raw_results.tsv.gz')
    df.to_csv(output_file, sep='\t', index=False, compression="gzip")
    print(f"Raw results saved to: {output_file}")

    # summarize per motif
    summary_df = df.groupby(['motif_name', 'motif_seq']).agg(
        mean_log2fc=('log2_fc_vs_baseline', 'mean'),
        std_log2fc=('log2_fc_vs_baseline', 'std'),
        sem_log2fc=('log2_fc_vs_baseline', lambda x: x.std() / np.sqrt(len(x))),
        n=('log2_fc_vs_baseline', 'count')
    ).reset_index()
    summary_df = summary_df.sort_values('mean_log2fc', ascending=False)
    
    print("\nMotif summary:")
    print(summary_df)
    
    # save summary
    summary_file = os.path.join(args.out_dir, f'individual_motifs.summary.tsv')
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"Summary saved to: {summary_file}")

    # save dict 
    dict_file = os.path.join(args.out_dir, 'individual_motif_predictions.pkl')
    with open(dict_file, 'wb') as f:
        pickle.dump(pred_dict, f)
    print(f"Dictionary saved to: {dict_file}")
    
    return df, summary_df

def save_motif_pair_results(df, args, logger):
    """Save and summarize motif pair results"""
    
    logger.info(f"\nMotif pair results summary:")
    logger.info(f"Total combinations: {len(df)}")
    logger.info(f"Unique pairs: {df[['motif1_name', 'motif2_name']].drop_duplicates().shape[0]}")
    logger.info(f"Orientations tested: {df['orientation_pattern'].unique()}")
    logger.info(f"Spacing range: {df['spacing'].min()}-{df['spacing'].max()} bp")
    
    # Save raw results
    output_file = os.path.join(args.out_dir, 'motif_pairs.raw_results.tsv.gz')
    df.to_csv(output_file, sep='\t', index=False, compression='gzip')
    logger.info(f"Raw results saved to: {output_file}")
    
    # Summarize by pair/orientation/spacing (across model folds)
    summary_df = df.groupby(['motif1_name', 'motif2_name', 'orientation_pattern', 'spacing']).agg({
        'log2_fc_vs_baseline': ['mean', 'std', 'sem'],
        'mean_prediction': ['mean', 'std'],
        'n_sequences': 'first'
    }).reset_index()
    
    # Flatten column names
    summary_df.columns = ['_'.join(col).strip('_') for col in summary_df.columns.values]
    
    summary_file = os.path.join(args.out_dir, 'motif_pairs.summary.tsv')
    summary_df.to_csv(summary_file, sep='\t', index=False)
    logger.info(f"Summary saved to: {summary_file}")

    return summary_df


# ---------- main ----------
def main():
    args = parse_args()
    
    # Create output directory
    args.out_dir = os.path.join(args.out_dir, args.analysis_id)
    os.makedirs(args.out_dir, exist_ok=True)
    logger = setup_logger(args.out_dir, args.analysis_id, args.log_level)

    # Uniform model paths and spacing
    model_paths = [os.path.join(args.model_base_path, f'fold{f}', args.model_subpath) 
                   for f in range(args.num_folds)]
    spacing_vals = list(range(args.spacing_min, args.spacing_max + 1, args.spacing_step))

    # Compose filenames for sharing background/singles/baselines
    singles_pickle_file = os.path.join(args.out_dir, 'individual_motif_predictions.pkl')
    singles_raw_file = os.path.join(args.out_dir, 'individual_motifs.raw_results.tsv.gz')
    singles_summary_file = os.path.join(args.out_dir, 'individual_motifs.summary.tsv')
    background_file = os.path.join(args.out_dir, 'background_encoded.npy')
    baseline_preds_file = os.path.join(args.out_dir, 'baseline_predictions.npy')

    ### =============== STEP 1: SHARED PREP (prepare mode) ===============

    if args.mode == "prepare":
        # Save background and singleton motif data to disk for sharing

        # Background sequences/encoding
        logger.info("Generating background sequences...")
        shuffled_sequences = get_shuffled_peak_sequences(
            args.narrow_peak_path, args.reference_genome, args.seq_len, args.num_bg)
        bg_encoded = one_hot_encode(shuffled_sequences, seq_length=args.seq_len)
        logger.info(f"Background shape: {bg_encoded.shape}")
        np.save(background_file, bg_encoded)

        # Compute baseline preds for each model fold
        logger.info("Computing baseline predictions (no motif insertion)...")
        model_info_list = [get_model(p) for p in model_paths]
        baseline_preds = [make_model_prediction(bg_encoded, m) for m in model_info_list]
        baseline_preds = np.stack(baseline_preds)
        np.save(baseline_preds_file, baseline_preds)
        logger.info(f"Saved baseline predictions to: {baseline_preds_file}")

        # Compute single motif predictions
        individual_motif_preds = {}
        all_individual_preds = []
        for motif_name, motif_seq in args.motif_dict.items():
            logger.info(f"Computing single motif for {motif_name}...")
            individual_inserted = insert_motifs_with_orientation_general(
                bg_encoded, n=1, spacing=0, motif_seq=[motif_seq],
                orientation_pattern=('+',), seq_len=args.seq_len
            )
            individual_preds = [make_model_prediction(individual_inserted, m) for m in model_info_list]
            individual_preds_stacked = np.stack(individual_preds)
            individual_motif_preds[motif_name] = individual_preds_stacked
            all_individual_preds.append(individual_preds_stacked)

        # Save picks
        with open(singles_pickle_file, 'wb') as f:
            pickle.dump(individual_motif_preds, f)
        logger.info(f"Saved single-motif predictions to: {singles_pickle_file}")
        
        # Save raw and summary tables
        df_individual, df_individual_summary = save_individual_motif_predictions(
            all_individual_preds, individual_motif_preds, baseline_preds, args, logger
        )
        logger.info(f"Finished prepare mode. Exiting.")
        return

    ### =============== STEP 2: CHUNKED PAIR ANALYSIS ("chunk" mode) ===============
    elif args.mode == "chunk":
        assert args.chunk_index is not None and args.chunk_n_chunks is not None, \
            "Must give --chunk-index and --chunk-n-chunks in chunk mode!"

        # Load background
        background_file = os.path.join(args.out_dir, 'background_encoded.npy')
        logger.info(f"Loading background sequence encodings from {background_file}")
        bg_encoded = np.load(background_file)
        
        # Load model paths
        model_paths = [os.path.join(args.model_base_path, f'fold{f}', args.model_subpath)
                       for f in range(args.num_folds)]
        spacing_vals = list(range(args.spacing_min, args.spacing_max + 1, args.spacing_step))

        # Load baseline predictions from prepare step
        logger.info(f"Loading baseline predictions from {baseline_preds_file}")
        baseline_preds = np.load(baseline_preds_file)

        # Load singles
        with open(singles_pickle_file, 'rb') as f:
            individual_motif_preds = pickle.load(f)
        logger.info(f"Loaded single motif predictions for {len(individual_motif_preds)} motifs")

        # Generate motif pair list and chunk as before
        motif_names = sorted(args.motif_dict.keys())
        all_pair_info = [
            (n1, n2, args.motif_dict[n1], args.motif_dict[n2])
            for n1, n2 in itertools.combinations_with_replacement(motif_names, 2)
        ]

        # Chunking
        total_chunks = args.chunk_n_chunks
        chunk_size = (len(all_pair_info) + total_chunks - 1) // total_chunks  # ceil division
        pair_chunks = chunk_motif_pairs(all_pair_info, chunk_size=chunk_size)

        my_pairs = pair_chunks[args.chunk_index]
        logger.info(f"[Chunk {args.chunk_index}/{total_chunks}] Processing {len(my_pairs)} motif pairs")

        shared_data = {
            'analysis_id': args.analysis_id,
            'bg_encoded': bg_encoded,
            'model_paths': model_paths,
            'baseline_preds': baseline_preds,
            'individual_motif_preds': individual_motif_preds,
            'seq_len': args.seq_len,
            'spacing_vals': spacing_vals,
        }

        # Process my chunk serially (no multiprocessing here, because each job is a chunk)
        results = process_motif_pair_chunk(my_pairs, shared_data)
        chunk_outfile = os.path.join(
            args.out_dir, f"motif_pairs.chunk_{args.chunk_index:03d}.results.tsv.gz")
        pd.DataFrame(results).to_csv(chunk_outfile, index=False, sep='\t', compression="gzip")
        logger.info(f"Saved results for chunk {args.chunk_index} to {chunk_outfile}")
        return

    elif args.mode == "collate":
        import glob
        pattern = os.path.join(args.out_dir, "motif_pairs.chunk_*.results.tsv.gz")
        files = sorted(glob.glob(pattern))
        dfs = [pd.read_csv(f, sep='\t') for f in files]
        df = pd.concat(dfs, ignore_index=True)
        logger.info(f"Found {len(files)} result chunk files: {files}")
        logger.info(f"Loaded and concatenated {len(files)} chunk files. Total results: {len(df)}.")

        # Optionally, also load singles for plot
        with open(os.path.join(args.out_dir, 'individual_motif_predictions.pkl'), 'rb') as f:
            individual_motif_preds = pickle.load(f)
        # If you saved `df_individual`/summary, load with pandas

        # Save & summarize motif pair results
        df_summary = save_motif_pair_results(df, args, logger)

        # Plot results (unless disabled)
        plot_our_dir = os.path.join(args.out_dir, "plots")
        if not args.no_plots:
            # If you saved/need `df_individual`, load it here from file
            df_individual = pd.read_csv(singles_raw_file, sep='\t')
            plot_motif_pairs(df, df_individual, args.motif_dict, plot_our_dir, heatmap_cmap=args.heatmap_cmap, custom_colors=None)
            logger.info("All plots saved successfully!")

        logger.info("Collate complete! All summary tables and plots should be present.")
        return
    else:
        raise ValueError(f"Unrecognized mode: {args.mode!r}")
    
# def main():
#     args = parse_args()
    
#     # Create output directory
#     args.out_dir = os.path.join(args.out_dir, args.analysis_id)
#     os.makedirs(args.out_dir, exist_ok=True)
    
#     # Build model paths
#     model_paths = [os.path.join(args.model_base_path, f'fold{f}', args.model_subpath) 
#                    for f in range(args.num_folds)]
    
#     # Generate spacing values
#     spacing_vals = list(range(args.spacing_min, args.spacing_max + 1, args.spacing_step))
    
#     logger = setup_logger(args.out_dir, args.analysis_id, args.log_level)

#     # start method/context
#     try:
#         if mp.get_start_method(allow_none=True) != "spawn":
#             mp.set_start_method("spawn", force=True)
#             logger.info("Multiprocessing start method set to 'spawn'")
#         else:
#             logger.info("Multiprocessing start method already 'spawn'")
#     except RuntimeError as e:
#         logger.warning(f"Could not set start method: {e}")
#     ctx = mp.get_context("spawn")

#     # run start snapshot
#     logger.info("===== RUN START =====")
#     logger.info(f"Host: {socket.gethostname()} | PID: {os.getpid()} | Python: {platform.python_version()} | System: {platform.platform()}")
#     logger.info(f"CUDA_VISIBLE_DEVICES={os.getenv('CUDA_VISIBLE_DEVICES')}")
#     logger.info(f"Analysis: {args.analysis_id}")
#     logger.info(f"Spacing: {args.spacing_min}–{args.spacing_max} bp (step={args.spacing_step}) | Total: {len(spacing_vals)} values")
#     logger.info(f"Sequences: {args.num_bg} background, length={args.seq_len}")
#     logger.info(f"reference_genome: {args.reference_genome} | exists: {os.path.exists(args.reference_genome)}")
#     logger.info(f"narrow_peak_path: {args.narrow_peak_path} | exists: {os.path.exists(args.narrow_peak_path)}")
#     logger.info(f"Output directory: {args.out_dir}")
#     for i, mpth in enumerate(model_paths):
#         logger.info(f"  Model [{i}] {mpth} | exists={os.path.exists(mpth)}")
#     logger.info(f"Requested processes: {args.n_proc} | chunk size: {args.chunk_size}")

#     # data prep
#     np.random.seed(args.random_seed)
#     t_all = time.perf_counter()
#     logger.info("Generating background sequences...")
#     t0 = time.perf_counter()
#     shuffled_sequences = get_shuffled_peak_sequences(args.narrow_peak_path, args.reference_genome, args.seq_len, args.num_bg)
#     bg_encoded = one_hot_encode(shuffled_sequences, seq_length=args.seq_len)
#     logger.info(f"Background sequences shape: {bg_encoded.shape} | bytes ~ {sizeof_fmt(bg_encoded.nbytes)} | took {time.perf_counter()-t0:.2f}s")

#     # baselines
#     logger.info("Computing baseline predictions (no motif insertion)...")
#     t0 = time.perf_counter()
#     #models = [get_model(p) for p in model_paths]
#     model_info_list = [get_model(mp_) for mp_ in model_paths]
#     baseline_preds = [make_model_prediction(bg_encoded, m) for m in model_info_list]
#     baseline_preds = np.stack(baseline_preds)
#     logger.info(f"Baseline predictions shape: {baseline_preds.shape} | took {time.perf_counter()-t0:.2f}s")

#     # Generate all motif pairs
#     motif_pairs = generate_motif_pairs(args.motif_dict, include_self_pairs=True)
#     logger.info(f"Generated {len(motif_pairs)} motif pairs to test")

#     # Individual motif predictions
#     dict_file = os.path.join(args.out_dir, 'individual_motif_predictions.pkl')
#     if args.use_existing_singles and os.path.exists(dict_file):
#         logger.info('Loading existing individual motif predictions....')
#         with open(dict_file, 'rb') as f:
#             individual_motif_preds = pickle.load(f)
#     else:
#         logger.info('Computing individual motif predictions...')
#         individual_motif_preds = {}
#         all_individual_preds = []

#         for motif_name, motif_seq in args.motif_dict.items():
#             logger.info(f"Computing single motif reference for {motif_name}...")
#             individual_inserted = insert_motifs_with_orientation_general(
#                 bg_encoded, n=1, spacing=0, motif_seq=[motif_seq], 
#                 orientation_pattern=('+',), seq_len=args.seq_len
#             )
#             individual_preds = [make_model_prediction(individual_inserted, m) for m in model_info_list]
#             individual_preds_stacked = np.stack(individual_preds)
#             individual_motif_preds[motif_name] = individual_preds_stacked
#             all_individual_preds.append(individual_preds_stacked)

#         # Save individual motif results
#         df_individual, df_individual_summary = save_individual_motif_predictions(
#             all_individual_preds, individual_motif_preds, baseline_preds, args, logger
#         )

#     # parallel analysis
#     logger.info("Starting parallel analysis...")
#     # Prepare shared data
#     shared_data = {
#         'analysis_id': args.analysis_id,
#         'bg_encoded': bg_encoded,
#         'model_paths': model_paths,
#         'baseline_preds': baseline_preds,
#         'individual_motif_preds': individual_motif_preds,
#         'seq_len': args.seq_len,
#         'spacing_vals': spacing_vals,  # full list of spacings
#     }

#     # Convert motif_pairs to list with sequences
#     motif_names = sorted(args.motif_dict.keys())
#     all_pair_info = [
#         (n1, n2, args.motif_dict[n1], args.motif_dict[n2])
#         for n1, n2 in itertools.combinations_with_replacement(motif_names, 2)
#     ]

#     logger.info(f"Total motif pairs to process: {len(all_pair_info)}")

#     # Determine chunk size
#     if args.chunk_size is None:
#         chunk_size = estimate_chunk_size(
#             n_pairs=len(all_pair_info),
#             n_workers=args.n_proc,
#             target_chunks_per_worker=2,
#             min_chunk_size=5,
#             max_chunk_size=15
#         )
#     else:
#         chunk_size = args.chunk_size
    
#     logger.info(f"Using chunk size: {chunk_size} pairs per chunk")

#     # Create chunks
#     pair_chunks = chunk_motif_pairs(all_pair_info, chunk_size=chunk_size)
#     logger.info(f"Created {len(pair_chunks)} chunks from {len(all_pair_info)} pairs")
    
#     # Prepare arguments for each chunk
#     chunk_args = []
#     for chunk_idx, chunk in enumerate(pair_chunks):
#         # Add chunk_idx to shared_data for this chunk
#         chunk_shared_data = shared_data.copy()
#         chunk_shared_data['chunk_idx'] = chunk_idx
#         chunk_args.append((chunk, chunk_shared_data))
    
#     # Run parallel processing
#     logger.info(f"Starting parallel processing with {args.n_proc} workers...")
    
#     with ctx.Pool(processes=args.n_proc, initializer=worker_init, initargs=(args.analysis_id, args.log_level)) as pool:
#         results = pool.starmap(process_motif_pair_chunk, chunk_args)
    
#     # Flatten results (each worker returns a list)
#     all_results = []
#     for chunk_results in results:
#         all_results.extend(chunk_results)
    
#     logger.info(f"Analysis complete. Total results: {len(all_results)}")

#     # save + plot
#     df = pd.DataFrame(all_results)
#     df_summary = save_motif_pair_results(df, args, logger)

#     # Plot results (unless disabled)
#     plot_our_dir = os.path.join(args.out_dir, "plots")
#     if not args.no_plots:
#         plot_motif_pairs(df, df_individual, args.motif_dict, plot_our_dir,
#                            heatmap_cmap=args.heatmap_cmap, custom_colors=None)
#         logger.info("All plots saved successfully!")

#     logger.info(f"Total runtime: {time.perf_counter()-t_all:.2f}s ({(time.perf_counter()-t_all)/3600:.2f}h)")
#     logger.info("===== RUN END =====")

if __name__ == "__main__":
    mp.freeze_support()
    main()