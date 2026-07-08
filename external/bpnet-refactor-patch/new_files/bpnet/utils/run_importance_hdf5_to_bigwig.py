import os
import pandas as pd
import h5py
import hdf5plugin
import numpy as np
import argparse
import gc
from bpnet.utils.importance_hdf5_to_bigwig import importance_hdf5_to_bigwig

parser = argparse.ArgumentParser(description="Create SHAP bigWig from given HDF5.")
parser.add_argument("--hdf5_path", type=str, default=None,
                    help="SHAP .h5 file path")
parser.add_argument("--regions_path", type=str, default=None,
                    help="Path to narrowPeak file of regions corresponding to SHAP h5")
parser.add_argument("--outfile", type=str, default=None,
                    help="Output file path for bigWig")
parser.add_argument("--outstats", type=str, default=None,
                    help="Output file path for stats")
parser.add_argument("--chrom_sizes", type=str, default=None,
                    help="Path to chromosome sizes")
args = parser.parse_args()


importance_hdf5_to_bigwig(args.hdf5_path, args.regions_path, args.outfile, args.outstats, args.chrom_sizes)
