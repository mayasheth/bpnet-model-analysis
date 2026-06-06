import pandas as pd
import pickle
import os

# =============================================================================
# CREATE FIMO TO ANNOTATION REGION MAPPING
# =============================================================================

# Output directory
OUTPUT_DIR = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/FIMO/elements_v1/analysis_v1"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# File paths
REGIONS_FILE = "/oak/stanford/groups/engreitz/Users/sheth/TF_analysis/2025_0609_K562_TF_annotations/finemo_peaks_all_chr.chromatin_annotations.tsv"
FIMO_RES_DIR = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/FIMO/elements_v1/"
MOTIF_NAMES = ["GATA_1", "STAT_0", "E2F_1", "ATF-CEBP_0", "FOS-JUN_0", "ELF_0", "CTCF_0", "ATF-CREB-JUN_0", "NF2L-NFE_0"]
FIMO_FILES = [os.path.join(FIMO_RES_DIR, f"{motif}.tsv") for motif in MOTIF_NAMES]
OUTPUT_MAPPING = os.path.join(OUTPUT_DIR, "fimo_to_annotation_mapping.pkl")

print("="*80)
print("CREATING REGION MAPPING")
print(f"OUTPUT DIR: {OUTPUT_DIR}")
print("="*80)

# Load regions file
print("\nLoading annotation regions...")
regions_df = pd.read_csv(REGIONS_FILE, sep='\t')
regions_df['sequence_name'] = (regions_df['chrom'].astype(str) + ':' + 
                                regions_df['start'].astype(str) + '-' + 
                                regions_df['end'].astype(str))
print(f"Loaded {len(regions_df)} annotation regions")

# Load FIMO files
print("\nLoading FIMO results...")
fimo_dfs = []
for fimo_file in FIMO_FILES:
    df = pd.read_csv(fimo_file, sep='\t')
    fimo_dfs.append(df)
    print(f"  {fimo_file}: {len(df)} hits")

motif_hits = pd.concat(fimo_dfs, ignore_index=True)
print(f"Total motif hits: {len(motif_hits)}")

# Parse regions to dataframes
def parse_regions_to_df(region_names):
    """Convert region names to DataFrame"""
    data = []
    for region in region_names:
        chrom, coords = region.split(':')
        start, end = coords.split('-')
        data.append({'region_name': region, 'chrom': chrom, 
                     'start': int(start), 'end': int(end)})
    return pd.DataFrame(data)

print("\nParsing FIMO regions...")
fimo_regions_df = parse_regions_to_df(motif_hits['sequence_name'].unique())
print(f"FIMO regions: {len(fimo_regions_df)}")

print("\nParsing annotation regions...")
annot_regions_df = regions_df[['sequence_name', 'chrom', 'start', 'end']].drop_duplicates()
annot_regions_df.columns = ['region_name', 'chrom', 'start', 'end']
print(f"Annotation regions: {len(annot_regions_df)}")

# Sort both dataframes
fimo_regions_df = fimo_regions_df.sort_values(['chrom', 'start']).reset_index(drop=True)
annot_regions_df = annot_regions_df.sort_values(['chrom', 'start']).reset_index(drop=True)

# Find overlaps
print("\nFinding overlaps (this may take a minute)...")
mapping = []

for chrom in fimo_regions_df['chrom'].unique():
    print(f"  Processing {chrom}...")
    fimo_chr = fimo_regions_df[fimo_regions_df['chrom'] == chrom]
    annot_chr = annot_regions_df[annot_regions_df['chrom'] == chrom]
    
    for _, fimo_row in fimo_chr.iterrows():
        # Find annotation regions that overlap
        overlaps = annot_chr[
            (annot_chr['start'] < fimo_row['end']) & 
            (annot_chr['end'] > fimo_row['start'])
        ]
        
        if len(overlaps) > 0:
            # If multiple overlaps, take the one with maximum overlap
            overlaps = overlaps.copy()
            overlaps['overlap_len'] = overlaps.apply(
                lambda x: min(x['end'], fimo_row['end']) - max(x['start'], fimo_row['start']), 
                axis=1
            )
            best_match = overlaps.loc[overlaps['overlap_len'].idxmax()]
            
            mapping.append({
                'fimo_region': fimo_row['region_name'],
                'annot_region': best_match['region_name'],
                'overlap': best_match['overlap_len']
            })

mapping_df = pd.DataFrame(mapping)

print(f"\n" + "="*80)
print("MAPPING SUMMARY")
print("="*80)
print(f"FIMO regions: {len(fimo_regions_df)}")
print(f"Annotation regions: {len(annot_regions_df)}")
print(f"FIMO regions mapped: {len(mapping_df)}")
print(f"FIMO regions unmapped: {len(fimo_regions_df) - len(mapping_df)}")

# Create dictionary for easy lookup
mapping_dict = dict(zip(mapping_df['fimo_region'], mapping_df['annot_region']))

# Save both the dataframe and dictionary
print(f"\nSaving mapping to {OUTPUT_MAPPING}...")
with open(OUTPUT_MAPPING, 'wb') as f:
    pickle.dump({
        'mapping_df': mapping_df,
        'mapping_dict': mapping_dict
    }, f)

# Also save as TSV for inspection
mapping_df.to_csv(os.path.join(OUTPUT_DIR, 'fimo_to_annotation_mapping.tsv'), sep='\t', index=False)

print("\n" + "="*80)
print("DONE!")
print("="*80)
print(f"Output directory: {OUTPUT_DIR}")
print(f"Mapping saved to:")
print(f"  - fimo_to_annotation_mapping.pkl (pickle format for Python)")
print(f"  - fimo_to_annotation_mapping.tsv (human-readable)")