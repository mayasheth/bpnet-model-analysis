import h5py
import hdf5plugin
import os
import glob
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Merge SHAP .h5 files across chromosomes for one fold")
parser.add_argument("--input-dir", required=True, help="Directory containing per-chromosome SHAP .h5 files")
parser.add_argument("--h5-filename", required=True, help="Name of h5 file to match in each chromosome directory, e.g. counts_scores.h5")
parser.add_argument("--output-file", required=True, help="Path to write the merged SHAP .h5 file")
args = parser.parse_args()

input_dir = args.input_dir
h5_filename = args.h5_filename
output_file = args.output_file

# Gather all per-chromosome h5 files
pattern = os.path.join(input_dir, "chr*/", h5_filename)
h5_files = sorted(glob.glob(pattern))

if not h5_files:
	raise FileNotFoundError(f"No .h5 files found in {input_dir}")

# Initialize empty lists
hyp_scores_lst=[]
chrom_lst=[]
start_lst=[]
end_lst=[]
input_seqs_lst=[]

# Read data from each file
for shap_h5 in h5_files:
	f = h5py.File(shap_h5, 'r')
	hyp_scores_lst.append(f['hyp_scores'][()])
	chrom_lst.append(f['coords_chrom'][()])
	start_lst.append(f['coords_start'][()])
	end_lst.append(f['coords_end'][()])
	input_seqs_lst.append(f['input_seqs'][()])
	f.close()

# Concatenate all data
[print(f"{chr[0]}: {len(chr)}") for chr in chrom_lst]

hyp_scores_all = np.concatenate(hyp_scores_lst)
print(hyp_scores_all.shape)
chrom_all = np.concatenate(chrom_lst)
start_all = np.concatenate(start_lst)
end_all = np.concatenate(end_lst)
input_seqs_all = np.concatenate(input_seqs_lst)

# Write merged output
num_examples = len(hyp_scores_all)
f = h5py.File(output_file, "w")

coords_chrom_dset = f.create_dataset(
	"coords_chrom", (num_examples,),
	dtype=h5py.string_dtype(encoding="ascii"), **hdf5plugin.Blosc()
)
coords_chrom_dset[:] = chrom_all

coords_start_dset = f.create_dataset(
	"coords_start", (num_examples,), dtype=int, **hdf5plugin.Blosc()
)
coords_start_dset[:] = start_all

coords_end_dset = f.create_dataset(
	"coords_end", (num_examples,), dtype=int, **hdf5plugin.Blosc()
)
coords_end_dset[:] = end_all

hyp_scores_dset = f.create_dataset(
	"hyp_scores", (num_examples, hyp_scores_all.shape[1], hyp_scores_all.shape[2]), **hdf5plugin.Blosc()
)
hyp_scores_dset[:, :, :] = hyp_scores_all

input_seqs_dset = f.create_dataset(
	"input_seqs", (num_examples, hyp_scores_all.shape[1], hyp_scores_all.shape[2]), **hdf5plugin.Blosc()
)
input_seqs_dset[:, :, :] = input_seqs_all

f.close()


print(f"Successfully merged {len(h5_files)} files into {output_file}")
