import os
import h5py
import numpy as np
import polars as pl

# conda activate finemo_gpu

MODISCO_PATTERN_GROUPS = ["pos_patterns"]

def trim_motif(cwm, trim_threshold):
    """
    Adapted from https://github.com/jmschrei/tfmodisco-lite/blob/570535ee5ccf43d670e898d92d63af43d68c38c5/modiscolite/report.py#L213-L236
    """
    score = np.sum(np.abs(cwm), axis=0) # total importance per position
    trim_thresh = np.max(score) * trim_threshold # trim_threshold is a percent, define threshold relative to max total importance
    pass_inds = np.nonzero(score >= trim_thresh)
    start = max(np.min(pass_inds), 0)
    end = min(np.max(pass_inds) + 1, len(score))

    return start, end

def load_modisco_motifs_one_pattern(modisco_h5_path, modisco_id, pattern_group, pattern_idx, 
	pattern_tag_new, flip_orientation, trim_threshold, start_new_fwd, end_new_fwd, 
	motif_type = "cwm", motif_lambda_default = 0.7, include_rc = True):

	"""
	Adapted from https://github.com/jmschrei/tfmodisco-lite/blob/570535ee5ccf43d670e898d92d63af43d68c38c5/modiscolite/report.py#L252-L272
	"""
	motif_data_lsts = {"motif_name": [], "motif_name_orig": [], "strand": [], "motif_start": [], 
						 "motif_end": [], "motif_scale": [], "lambda": []}
	motif_lst = [] 
	trim_mask_lst = []

	with h5py.File(modisco_h5_path, 'r') as modisco_results:
		pattern_tag_orig = f'{modisco_id}.{pattern_group}.pattern_{pattern_idx}'

		pattern = modisco_results[pattern_group][f'pattern_{pattern_idx}']
		motif_lambda = motif_lambda_default
		pattern_tag = pattern_tag_new

		# get and trim CWM
		cwm_raw = pattern['contrib_scores'][:].T  # assume now 4×L
		if cwm_raw.shape[0] != 4:
			cwm_raw = cwm_raw.T  # ensure 4×L

		L = cwm_raw.shape[1]
		cwm_norm = np.sqrt((cwm_raw**2).sum())

		if not flip_orientation:
			cwm_fwd = cwm_raw / cwm_norm
			cwm_rev = cwm_fwd[::-1,::-1]
		else:
			cwm_rev = cwm_raw / cwm_norm
			cwm_fwd= cwm_rev[::-1,::-1]

		if start_new_fwd is not None and end_new_fwd is not None:
				# Forward slice [start_fwd : end_fwd)
				start_fwd, end_fwd = start_new_fwd, end_new_fwd

				# Reverse slice maps to [L - end_fwd : L - start_fwd)
				start_rev = L - end_fwd
				end_rev   = L - start_fwd
		else:
				start_fwd, end_fwd = trim_motif(cwm_fwd, trim_threshold)
				start_rev, end_rev = trim_motif(cwm_rev, trim_threshold)
		
		trim_mask_fwd = np.zeros(cwm_fwd.shape[1], dtype=np.int8)
		trim_mask_fwd[start_fwd:end_fwd] = 1
		trim_mask_rev = np.zeros(cwm_rev.shape[1], dtype=np.int8)
		trim_mask_rev[start_rev:end_rev] = 1

		if motif_type == "cwm":
			motif_fwd = cwm_fwd
			motif_rev = cwm_rev
			motif_norm = cwm_norm

		else:
			raise ValueError(f"Motif type {motif_type} not supported.")
		

		# add to motif data lists
		motif_data_lsts["motif_name"].append(pattern_tag)
		motif_data_lsts["motif_name_orig"].append(pattern_tag_orig)
		motif_data_lsts["strand"].append('+')
		motif_data_lsts["motif_start"].append(start_fwd)
		motif_data_lsts["motif_end"].append(end_fwd)
		motif_data_lsts["motif_scale"].append(motif_norm)
		motif_data_lsts["lambda"].append(motif_lambda)

		if include_rc:
			motif_data_lsts["motif_name"].append(pattern_tag)
			motif_data_lsts["motif_name_orig"].append(pattern_tag_orig)
			motif_data_lsts["strand"].append('-')
			motif_data_lsts["motif_start"].append(start_rev)
			motif_data_lsts["motif_end"].append(end_rev)
			motif_data_lsts["motif_scale"].append(motif_norm)
			motif_data_lsts["lambda"].append(motif_lambda)

			motif_lst.extend([motif_fwd, motif_rev])
			trim_mask_lst.extend([trim_mask_fwd, trim_mask_rev])

		else:
			motif_lst.append(motif_fwd)
			trim_mask_lst.append(trim_mask_fwd)

	return motif_data_lsts, motif_lst, trim_mask_lst


def load_modisco_motifs_from_table(
    table,                       # polars.DataFrame
    motif_type="cwm",
    trim_threshold_default=0.2,
    motif_lambda_default=0.7,
    include_rc_default=True,
):
    """
    Iterates over a table of motif specs, calls `load_modisco_motifs_one_pattern` per row,
    and returns (motifs_df, cwms, trim_masks, names) in the same shapes as your original function.
    """

    # Polars <-> Pandas compat: get row dicts
    rows = table.to_dicts()

    agg_motif_data = {
        "motif_name": [], "motif_name_orig": [], "strand": [],
        "motif_start": [], "motif_end": [], "motif_scale": [], "lambda": []
    }
    agg_cwms = []
    agg_masks = []

    for r in rows:
        modisco_h5_path = r["modisco_h5_path"]
        modisco_id      = r["modisco_id"]
        pattern_group   = r["pattern_group"]
        pattern_idx     = int(r["pattern_idx"])
        pattern_tag_new = r["pattern_tag_new"]
        flip_orientation= bool(r["flip_orientation"])

        # Per-row overrides (use defaults if missing/None)
        trim_threshold  = r.get("trim_threshold", trim_threshold_default)
        row_lambda      = r.get("lambda", None)
        include_rc      = bool(r.get("include_rc", include_rc_default))

        # Call your single-motif loader
        motif_data_lsts, motif_lst, trim_mask_lst = load_modisco_motifs_one_pattern(
            modisco_h5_path=modisco_h5_path,
            modisco_id=modisco_id,
            pattern_group=pattern_group,
            pattern_idx=pattern_idx,
            pattern_tag_new=pattern_tag_new,
            flip_orientation=flip_orientation,
            trim_threshold=trim_threshold,
            start_new_fwd=r.get("start_new_fwd", None),
            end_new_fwd=r.get("end_new_fwd", None),
            motif_type=motif_type,
            motif_lambda_default=(row_lambda if row_lambda is not None else motif_lambda_default),
            include_rc=include_rc,
        )

        # Append
        for k, v in motif_data_lsts.items():
            agg_motif_data[k].extend(v)
        agg_cwms.extend(motif_lst)
        agg_masks.extend(trim_mask_lst)

    if len(agg_cwms) == 0:
        raise ValueError("No motifs were loaded; check your table filters/paths.")

    # Finalize outputs in the same structure you expect
    motifs_df = pl.DataFrame(agg_motif_data).with_row_index(name="motif_id")
    cwms      = np.stack(agg_cwms,  dtype=np.float16, axis=0)  # [N, 4, L] (or your current orientation)
    trim_masks= np.stack(agg_masks, dtype=np.int8,   axis=0)   # [N, L]

    # Names = forward-strand motif names (to match your original function)
    names = motifs_df.filter(pl.col("strand") == "+").get_column("motif_name").to_numpy()

    return motifs_df, cwms, trim_masks, names


def save_motifs_bundle(
    out_prefix: str,
    motifs_df: pl.DataFrame,
    cwms: np.ndarray,
    trim_masks: np.ndarray,
    names: np.ndarray,
):
    """
    Saves:
      - Polars DataFrame -> Parquet
      - Arrays (cwms, trim_masks, names) -> NPZ (compressed)

    Returns the three paths (parquet, npz).
    """
    # Basic validations
    if not isinstance(motifs_df, pl.DataFrame):
        raise TypeError("motifs_df must be a polars.DataFrame")
    if not (isinstance(cwms, np.ndarray) and cwms.ndim == 3):
        raise ValueError("cwms must be a 3D numpy array [N, 4, L]")
    if not (isinstance(trim_masks, np.ndarray) and trim_masks.ndim == 2):
        raise ValueError("trim_masks must be a 2D numpy array [N, L]")
    if cwms.shape[0] != trim_masks.shape[0] or cwms.shape[2] != trim_masks.shape[1]:
        raise ValueError("cwms and trim_masks shapes are incompatible")
    if not isinstance(names, np.ndarray):
        names = np.asarray(names)

    parquet_path = f"{out_prefix}.motifs.parquet"
    npz_path     = f"{out_prefix}.motifs.npz"

    # Save table
    motifs_df.write_parquet(parquet_path)

    # Save arrays (compressed)
    np.savez_compressed(
        npz_path,
        cwms=cwms,
        trim_masks=trim_masks,
        names=names,
    )

    return parquet_path, npz_path




#### RUN ####

PROJECT_DIR = '/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet'
curated_motif_tsv = os.path.join(PROJECT_DIR, 'reference/curated_motif_data_v2.tsv')
out_prefix = os.path.join(PROJECT_DIR, 'reference/curated_motif_data_for_finemo_v2')

# 1) Read the control table
table = pl.read_csv(curated_motif_tsv, separator="\t")

# 2) Build motifs from table
motifs_df, cwms, trim_masks, names = load_modisco_motifs_from_table(
		table=table,
		motif_type="cwm",
		motif_lambda_default=0.7,
		include_rc_default=True,
    )

# 3) Save the bundle
save_motifs_bundle(out_prefix, motifs_df, cwms, trim_masks, names)

# Later (or in a different script), load them back:
# motifs_df2, cwms2, masks2, names2, meta2 = load_motifs_bundle(out_prefix, expect_include_rc=True)

