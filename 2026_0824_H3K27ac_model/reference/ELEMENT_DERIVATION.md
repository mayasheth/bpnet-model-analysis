# Candidate element provenance across the panel

**Decision (Maya, 2026-08-29): proceed with the element sets we already have and record how
each was derived, rather than re-deriving to a single definition.** K562 and GM12878
ATAC-based element sets are not readily available, and re-deriving would block the panel.
The consequence is a confound that must be stated wherever cell types are compared, not
silently carried.

## What each set is

| Cell type | Element set in use | Derived from | Notes |
|---|---|---|---|
| K562 | `reference/K562_DNase_candidate_elements.narrowPeak` | **DNase** | Separate derivation, NOT its own rE2G EnhancerList. 96% of loci overlap that list, 95% at 50% reciprocal overlap, but boundaries differ, so windows are centred slightly differently than for GM12878. |
| GM12878 | `2026_0606_GM12878_transferability/reference/GM12878_candidate_elements.narrowPeak` | **DNase** | IS its rE2G EnhancerList exactly (154,224 / 154,224). |
| TeloHAEC x4 | `reference/celltype_elements/TeloHAEC_*_ATAC_candidate_elements.narrowPeak` | **ATAC** (`atac_h3k27ac_powerlaw` rE2G model dir) | Converted by `scripts/0.13.make_candidate_elements.py` from `Neighborhoods/EnhancerList.bed`, which is BED4 with no summit column, so `summit = width // 2`. |
| H1, H9, HCT116, Jurkat, THP-1 | `reference/celltype_elements/*_DNase_candidate_elements.narrowPeak` | **DNase** | Built but unused: these cell types have no ATAC (see the ATAC-only rule). H9's elements are much wider (mean 798 bp vs ~570-600), which changes its neighbour-contamination profile at any given window. |

## The confound, stated plainly

The accessibility assay used to CALL elements is not the same across the panel: DNase for
K562/GM12878, ATAC for TeloHAEC. Under the ATAC-only rule for model *inputs*, TeloHAEC is
the one whose element derivation matches its input assay, and K562/GM12878 are the
mismatched ones — the opposite of how this was originally framed.

So any K562 <-> TeloHAEC or GM12878 <-> TeloHAEC comparison differs in **two** ways at once:
cell type, and how its elements were called. Element boundaries determine where every
counting window is centred, so this is not a nuisance detail — it shifts the target.

**Required when reporting:** state the derivation assay alongside any cross-cell-type
number, and do not attribute a K562-vs-TeloHAEC difference to cell type alone. Where the
difference is small relative to the effect being claimed, say so; where it is not,
the comparison is not clean and should be labelled as such.

**Cheapest future fix if this ever blocks a claim:** re-derive TeloHAEC elements from a
DNase-based rE2G run (making all three DNase-derived), or obtain ATAC-based element sets for
K562/GM12878 (making all three ATAC-derived). Either direction is internally consistent; the
current state is not.

## Addendum 2026-08-31: an ATAC-derived K562 element set does exist

Found while auditing disk usage, not while looking for it. `K562_ATAC_ChromBPNet/log.sh`
builds its regions from
`ENCODE_rE2G/results/2025_0226_ATAC_powerlaw_models/ATAC_H3K27ac_powerlaw/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed`
and the converted result sits at `K562_ATAC_ChromBPNet/data/K562_candidate_elements.narrowPeak`
(~154k regions). That is the **same rE2G model type** TeloHAEC's elements come from
(`atac_h3k27ac_powerlaw`), so K562 could be made derivation-consistent with TeloHAEC exactly.

**GM12878 has no equivalent** — that results directory contains K562 only.

So the available options are:
- Leave as is. K562/GM12878 DNase-derived, TeloHAEC ATAC-derived; 1 of 3 matches its input
  assay; confound stated wherever cell types are compared. **Current course.**
- Switch K562 to the ATAC-derived set. 2 of 3 consistent, GM12878 still mismatched — a
  partial fix that would require re-deriving every K562 number in the report and in F-001
  through F-003. Poor trade unless a specific claim turns out to hinge on it.

Recorded so the option is known rather than rediscovered. The earlier statement that
K562/GM12878 ATAC element sets were not readily available was correct for GM12878 and wrong
for K562.
