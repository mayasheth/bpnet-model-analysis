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
import multiprocessing as mp

from motif_exp_utils import (
    get_model, get_shuffled_peak_sequences, one_hot_encode, pattern_to_string,
    make_model_prediction, get_representative_patterns,
    insert_motifs_with_orientation_general
)

from plot_motif_spacing import plot_motif_spacing

# ---------- argument parsing ----------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze motif spacing effects on model predictions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Core analysis parameters
    parser.add_argument("--analysis-id",
                        help="Unique identifier for this analysis run")
    parser.add_argument("--motif-name",
                        help="Name of the motif for labeling")
    parser.add_argument("--motif-seq",
                        help="Motif sequence to insert")
    parser.add_argument("--motif-counts", nargs="+", type=int, default=[2],
                        help="Number of motifs to insert (can specify multiple)")
    
    # Spacing parameters
    parser.add_argument("--spacing-min", type=int, default=0,
                        help="Minimum spacing between motifs (bp)")
    parser.add_argument("--spacing-max", type=int, default=50,
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
    
    # Plot parameters
    parser.add_argument("--heatmap-cmap", default="RdBu_r",
                        help="Colormap for heatmap plots")
    parser.add_argument("--no-plots", action="store_true",
                        help="Skip generating plots (only save data)")
    
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
    
    return args

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

# ---------- worker ----------
def process_spacing_chunk(spacing_chunk, shared_data):
    t0 = time.perf_counter()
    pid = os.getpid()
    analysis_id = shared_data['analysis_id']
    log = logging.getLogger(analysis_id)

    bg_encoded       = shared_data['bg_encoded']
    baseline_preds   = shared_data['baseline_preds']      # (n_models, n_bg)
    single_motif_preds = shared_data['single_motif_preds']
    model_paths      = shared_data['model_paths']
    motif_seq        = shared_data['motif_seq']
    motif_counts     = shared_data['motif_counts']
    seq_len          = shared_data['seq_len']

    # build all conditions for the chunk
    conditions, enc_list, sizes = [], [], []
    for n in motif_counts:
        for pattern in get_representative_patterns(n):
            pstr = pattern_to_string(pattern)
            for s in spacing_chunk:
                t_ins = time.perf_counter()
                enc = insert_motifs_with_orientation_general(bg_encoded, n, s, motif_seq, pattern, seq_len)
                ins_s = time.perf_counter() - t_ins
                conditions.append((n, pstr, s, enc.shape[0], ins_s))
                enc_list.append(enc); sizes.append(enc.shape[0])

    if not enc_list:
        log.info(f"[worker pid={pid}] No conditions in {spacing_chunk}")
        return []

    big_encoded = np.concatenate(enc_list, axis=0)
    offsets = np.cumsum([0] + sizes)

    # load models once
    t_load = time.perf_counter()
    model_info_list = [get_model(mp_) for mp_ in model_paths]
    #models = [get_model(mp_) for mp_ in model_paths]
    log.info(f"[worker pid={pid}] Loaded {len(model_info_list)} models in {time.perf_counter()-t_load:.2f}s | "
             f"chunk={spacing_chunk} | conds={len(conditions)} | big={big_encoded.shape}")

    # one predict per model on the big batch
    model_full_preds = []
    for mi, model_info in enumerate(model_info_list):
        t_pm = time.perf_counter()
        preds = make_model_prediction(big_encoded, model_info)  # (sum sizes,)
        model_full_preds.append(preds)
        log.info(f"[worker pid={pid}] model{mi} predict len={preds.shape[0]} took {time.perf_counter()-t_pm:.2f}s")

    # slice back and compute metrics
    out = []
    for idx, (n, pstr, s, enc_len, ins_s) in enumerate(conditions):
        i0, i1 = offsets[idx], offsets[idx+1]
        for mi in range(len(model_info_list)):
            preds_cond = model_full_preds[mi][i0:i1]
            base = baseline_preds[mi]; single = single_motif_preds[mi]
            # expect enc_len == base.shape[0]
            if preds_cond.shape[0] != base.shape[0]:
                log.warning(f"[pid={pid}] size mismatch spacing={s} pattern={pstr}: {preds_cond.shape[0]} vs {base.shape[0]}")
            log2_vs_base  = float(np.mean(preds_cond - base[:preds_cond.shape[0]]) / np.log(2))
            log2_vs_single= float(np.mean(preds_cond - single[:preds_cond.shape[0]]) / np.log(2))
            out.append({
                "motif_seq": motif_seq, "motif_counts": n, "orientation_pattern": pstr,
                "spacing": s, "model_fold": mi, "log2_fc_vs_baseline": log2_vs_base,
                "log2_fc_vs_single": log2_vs_single, "mean_prediction": float(np.mean(preds_cond)),
                "n_sequences": int(preds_cond.shape[0])
            })
        log.info(f"[worker pid={pid}] done {idx+1}/{len(conditions)} (n={n}, pattern={pstr}, spacing={s}) insert={ins_s:.2f}s")

    log.info(f"[worker pid={pid}] Finished chunk {spacing_chunk} in {time.perf_counter()-t0:.2f}s with {len(out)} rows")
    return out

def chunk_list(lst, n): return [lst[i:i+n] for i in range(0, len(lst), n)]

# ---------- main ----------
def main():
    args = parse_args()
    
    # Create output directory
    args.out_dir = os.path.join(args.out_dir, args.analysis_id)
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Build model paths
    model_paths = [os.path.join(args.model_base_path, f'fold{f}', args.model_subpath) 
                   for f in range(args.num_folds)]
    
    # Generate spacing values
    spacing_vals = list(range(args.spacing_min, args.spacing_max + 1, args.spacing_step))
    
    logger = setup_logger(args.out_dir, args.analysis_id, args.log_level)

    # start method/context
    try:
        if mp.get_start_method(allow_none=True) != "spawn":
            mp.set_start_method("spawn", force=True)
            logger.info("Multiprocessing start method set to 'spawn'")
        else:
            logger.info("Multiprocessing start method already 'spawn'")
    except RuntimeError as e:
        logger.warning(f"Could not set start method: {e}")
    ctx = mp.get_context("spawn")

    # run start snapshot
    logger.info("===== RUN START =====")
    logger.info(f"Host: {socket.gethostname()} | PID: {os.getpid()} | Python: {platform.python_version()} | System: {platform.platform()}")
    logger.info(f"CUDA_VISIBLE_DEVICES={os.getenv('CUDA_VISIBLE_DEVICES')}")
    logger.info(f"Analysis: {args.analysis_id}")
    logger.info(f"Motif: {args.motif_name} ({args.motif_seq}) | Counts: {args.motif_counts}")
    logger.info(f"Spacing: {args.spacing_min}–{args.spacing_max} bp (step={args.spacing_step}) | Total: {len(spacing_vals)} values")
    logger.info(f"Sequences: {args.num_bg} background, length={args.seq_len}")
    logger.info(f"reference_genome: {args.reference_genome} | exists: {os.path.exists(args.reference_genome)}")
    logger.info(f"narrow_peak_path: {args.narrow_peak_path} | exists: {os.path.exists(args.narrow_peak_path)}")
    logger.info(f"Output directory: {args.out_dir}")
    for i, mpth in enumerate(model_paths):
        logger.info(f"  Model [{i}] {mpth} | exists={os.path.exists(mpth)}")
    logger.info(f"Requested processes: {args.n_proc} | chunk size: {args.chunk_size}")

    # data prep
    np.random.seed(args.random_seed)
    t_all = time.perf_counter()
    logger.info("Generating background sequences...")
    t0 = time.perf_counter()
    shuffled_sequences = get_shuffled_peak_sequences(args.narrow_peak_path, args.reference_genome, args.seq_len, args.num_bg)
    bg_encoded = one_hot_encode(shuffled_sequences, seq_length=args.seq_len)
    logger.info(f"Background sequences shape: {bg_encoded.shape} | bytes ~ {sizeof_fmt(bg_encoded.nbytes)} | took {time.perf_counter()-t0:.2f}s")

    # baselines
    logger.info("Computing baseline predictions (no motif insertion)...")
    t0 = time.perf_counter()
    #models = [get_model(p) for p in model_paths]
    model_info_list = [get_model(mp_) for mp_ in model_paths]
    baseline_preds = [make_model_prediction(bg_encoded, m) for m in model_info_list]
    baseline_preds = np.stack(baseline_preds)
    logger.info(f"Baseline predictions shape: {baseline_preds.shape} | took {time.perf_counter()-t0:.2f}s")

    logger.info("Computing single motif reference...")
    t0 = time.perf_counter()
    single_motif_inserted = insert_motifs_with_orientation_general(
        bg_encoded, n=1, spacing=0, motif_seq=args.motif_seq, orientation_pattern=('+',), seq_len=args.seq_len
    )
    single_motif_preds = [make_model_prediction(single_motif_inserted, m) for m in model_info_list]
    single_motif_preds = np.stack(single_motif_preds)
    logger.info(f"Single motif predictions shape: {single_motif_preds.shape} | took {time.perf_counter()-t0:.2f}s")
    del model_info_list  # free parent memory

    # parallel analysis
    logger.info("Starting parallel analysis...")
    shared_data = {
        'analysis_id': args.analysis_id,
        'bg_encoded': bg_encoded,
        'baseline_preds': baseline_preds,
        'single_motif_preds': single_motif_preds,
        'model_paths': model_paths,
        'motif_seq': args.motif_seq,
        'motif_counts': args.motif_counts,
        'seq_len': args.seq_len,
    }
    spacing_chunks = chunk_list(spacing_vals, args.chunk_size)
    logger.info(f"Split {len(spacing_vals)} spacings into {len(spacing_chunks)} chunks")
    all_results = []

    fn = partial(process_spacing_chunk, shared_data=shared_data)
    with ctx.Pool(args.n_proc, initializer=worker_init, initargs=(args.analysis_id, args.log_level)) as pool:
        completed = 0
        for chunk_result in pool.imap_unordered(fn, spacing_chunks):
            all_results.extend(chunk_result)
            completed += 1
            if completed % max(1, len(spacing_chunks)//10) == 0:
                logger.info(f"Progress: {completed}/{len(spacing_chunks)} chunks complete")

    # save + plot
    df = pd.DataFrame(all_results)
    print("\nResults summary:")
    print(df.head(10))
    print(f"Shape: {df.shape}")

    # Add summary statistics
    print(f"\nSummary statistics:")
    print(f"Motif counts tested: {sorted(df['motif_counts'].unique())}")
    print(f"Orientation patterns: {sorted(df['orientation_pattern'].unique())}")
    print(f"Spacing range: {df['spacing'].min()}-{df['spacing'].max()} bp")
    print(f"Log2 FC range (vs baseline): {df['log2_fc_vs_baseline'].min():.3f} to {df['log2_fc_vs_baseline'].max():.3f}")

    # Save raw results
    output_file = os.path.join(args.out_dir, f'raw_results.tsv')
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Raw results saved to: {output_file}")

    # Create summary with confidence intervals
    df_summary = df.groupby(['motif_counts', 'orientation_pattern', 'spacing']).agg({
        'log2_fc_vs_baseline': ['mean', 'std', 'sem'],
        'log2_fc_vs_single': ['mean', 'std', 'sem'],
        'mean_prediction': ['mean', 'std', 'sem']
    }).round(4)

    # Flatten column names
    df_summary.columns = [f"{col[0]}_{col[1]}" for col in df_summary.columns]
    df_summary = df_summary.reset_index()

    # Save summary
    summary_file = os.path.join(args.out_dir, f'results_summary.tsv')
    df_summary.to_csv(summary_file, sep='\t', index=False)
    print(f"Summary saved to: {summary_file}")

    # Plot results (unless disabled)
    if not args.no_plots:
        plot_motif_spacing(df, args.motif_name, args.motif_counts, args.out_dir,
                           heatmap_cmap=args.heatmap_cmap, custom_colors=None)
        logger.info("All plots saved successfully!")

    logger.info(f"Total runtime: {time.perf_counter()-t_all:.2f}s ({(time.perf_counter()-t_all)/3600:.2f}h)")
    logger.info("===== RUN END =====")

if __name__ == "__main__":
    mp.freeze_support()
    main()