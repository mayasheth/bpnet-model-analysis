import os
import numpy as np
import pandas as pd
from memelite.io import read_meme
from memelite import fimo

# conda activate tfmodisco

PROJECT_DIR = '/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet'
motif_compendium_file = os.path.join(PROJECT_DIR, 'reference/MotifCompendium-Database-Human.meme.txt')
elements_file = os.path.join(PROJECT_DIR, 'reference/K562_DNase_candidate_elements.fa')
peaks_file = os.path.join(PROJECT_DIR, 'reference/ENCSR000EGE_peaks_inliers.fa') 

# EDIT PER RUN
selected_motifs = ['GATA_1', 'NF2L-NFE_0', 'ATF-CEBP_0', 'STAT_0', 'ZNF_3', 'ELF_0', 'FOS-JUN_0', 'TAL_0',
                   'E2F_1', 'CTCF_0', 'ATF-CREB-JUN_0', 'ETV-ELF-NFAT-NFAC-ELK-SP_0', 'HOX_0', 'SNAI-TCF_0',
                   'SREBF-USF_0', 'CREB_0', 'RUNX_0', 'NFI_0', 'RFX_1', 'NR3C_0','ZNF_0', 'NFI_1',
                   'HNF4_0', 'GLI_1', 'HSF_0'] 

sequence_file = elements_file
out_prefix = os.path.join(PROJECT_DIR, 'FIMO/elements_v1')
os.makedirs(out_prefix, exist_ok = True)

# read in motifs
all_motif_dict = read_meme(motif_compendium_file)
motif_dict = {key: all_motif_dict[key] for key in selected_motifs}


# fimo(motifs, sequences, alphabet=['A', 'C', 'G', 'T'], bin_size=0.1, 
# 	eps=0.0001, threshold=0.0001, reverse_complement=True, return_counts=False, 
# 	dim=0)


fimo_hits = fimo(motif_dict, elements_file, threshold=0.001)

# save
# columns in each df: motif_name	motif_idx	sequence_name start	end	strand	score	p-value
# indexing is 0-based, unlike 1-based in fimo

for df, motif in zip(fimo_hits, selected_motifs):
  out_file = os.path.join(out_prefix, f'{motif}.tsv')
  df.to_csv(out_file, sep = "\t")
