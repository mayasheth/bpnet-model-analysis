
#!/usr/bin/env Rscript
# Usage: Rscript 5.4.plot_hits_vs_chip_data.R 2025_0703_retrain_p300_model/finemo/pkw_500_curated_motifs/annotated_motifs/hits_renamed.tsv 2025_0703_retrain_p300_model/finemo/pkw_500_curated_motifs/annotated_motifs/hits_per_peak.with_predictions.tsv.gz 2025_0703_retrain_p300_model/finemo/pkw_500_curated_motifs/plot_hits_vs_chip
# Rscript scripts/5.4.plot_hits_vs_chip_data.R 2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs/annotated_motifs/hits_renamed.tsv 2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs/annotated_motifs/hits_per_peak.with_predictions.tsv.gz 2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs/plot_hits_vs_chip

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(viridis)
library(scales)

# Parse command line arguments manually
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript motif_analysis.R hits_file peaks_file output_dir\n")
  cat("  hits_file: Path to hits_renamed.tsv\n")
  cat("  peaks_file: Path to peak_summary.tsv\n") 
  cat("  output_dir: Output directory\n")
  quit(status = 1)
}

hits_file <- args[1]
peaks_file <- args[2] 
outdir <- args[3]

chrom_annot_file <- '/oak/stanford/groups/engreitz/Users/sheth/TF_analysis/2025_0609_K562_TF_annotations/finemo_peaks_all_chr.chromatin_annotations.tsv'

# Create output directory
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load data and merge chromatin data to peaks
hits <- read.table(hits_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
chrom_annot <- read.table(chrom_annot_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
  select(-c("EP300_peak_overlap", "EP300.RPM", "peak_id", "n_hits", "motif_list",	"motif_strand_list", "mean_pred_logcounts", "true_logcounts", "chrom", "start", "end"))
peaks <- read.table(peaks_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% 
  left_join(chrom_annot, by = "peak_name")

## chrom annot columns
# EP300_peak_overlap	CREBBP_peak_overlap	CTCF_peak_overlap	H3K27me3_peak_overlap	H3K4me1_peak_overlap	GATA1_peak_overlap	GATA2_peak_overlap
# ATF4_peak_overlap	CEBPG_peak_overlap	NRF1_peak_overlap	MAFK_peak_overlap	MAFF_peak_overlap	MAFG_peak_overlap	STAT5A_peak_overlap
# STAT5B_peak_overlap	JUN_peak_overlap	JUNB_peak_overlap	JUND_peak_overlap	ZNF281_peak_overlap
# EP300.RPM	CREBBP.RPM	CTCF.RPM	DHS.RPM	H3K27ac.RPM	H3K27me3.RPM	H3K4me1.RPM	GATA1.RPM	GATA2.RPM	ATF4.RPM	CEBPG.RPM	NRF1.RPM
# MAFK.RPM	MAFF.RPM	MAFG.RPM	STAT5A.RPM	STAT5B.RPM	TBL1XR1.RPM	JUN.RPM	JUNB.RPM	JUND.RPM	ZNF281.RPM
# H3K27ac.RPM.expandedRegion	H3K27me3.RPM.expandedRegion	H3K4me1.RPM.expandedRegion

## original peak summary columns
# peak_id	peak_region_start	n_hits	motif_list	motif_strand_list
# sum_hit_coefficient_global	sum_hit_coefficient	sum_hit_importance
# chrom	start	end	mean_pred_logcounts	true_logcounts peak_name	EP300_peak_overlap	EP300.RPM


plot_max_similarity_by_overlap <- function(peaks, hits, motif_names, rpm_col, overlap_col, color, out_file) {
  hits_smry <- hits %>% filter(motif_name %in% motif_names) %>%
    group_by(peak_id) %>% 
    summarize(max_similarity = max(hit_similarity), .groups = "drop")
  
  peaks <- peaks %>% left_join(hits_smry, by = "peak_id") %>% 
    mutate(max_similarity = ifelse(is.na(max_similarity), 0, max_similarity)) %>% 
    mutate(rpm_log1p = log10(!!sym(rpm_col) + 1)) %>%
    mutate(overlap_status = factor(!!sym(overlap_col), levels = c(0, 1), labels = c("No overlap", "Overlap")))

  peaks_nonzero_similarity <- peaks %>%
    filter(max_similarity > 0)

  # Plot 1: RPM as a function of max similarity (2D density plot)
  # Get reasonable y-axis limit
  y_max <- quantile(peaks$rpm_log1p, 0.95)

  # Calculate correlations for finite values
  spearman_r <- cor(peaks$max_similarity, peaks$rpm_log1p, method = "spearman")
  pearson_r <- cor(peaks$max_similarity, peaks$rpm_log1p, method = "pearson")
  
  # Create custom colormap from white to specified color
  n_levels <- 10
  color_palette <- colorRampPalette(c("white", color))(n_levels)
  
  p1 <- ggplot(peaks, aes(x = max_similarity, y = rpm_log1p)) +
    stat_density_2d_filled(contour_var = "ndensity") +
    ylim(c(0, y_max)) +
   scale_fill_manual(values = color_palette) + scale_x_log10(limits=c(1e-2,1)) +
    labs(
      x = "Max hit similarity",
      y = paste0("log1p ", gsub("\\.", " ", rpm_col )),
      title = "All elements"
    ) +
    annotate("text", x = Inf, y = Inf, 
             label = paste0("Spearman r = ", round(spearman_r, 2), "\nPearson r = ", round(pearson_r, 2), "\nN = ", nrow(peaks)),
             hjust = 1.1, vjust = 1.2, size = 3.5) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12),
      plot.title = element_text(color = "black", size = 12),
      axis.ticks = element_line(color = "black")
    )

  # nonzero scatter
  spearman_r <- cor(peaks_nonzero_similarity$max_similarity, peaks_nonzero_similarity$rpm_log1p, method = "spearman")
  pearson_r <- cor(peaks_nonzero_similarity$max_similarity, peaks_nonzero_similarity$rpm_log1p, method = "pearson")

  p1_nonzero <- ggplot(peaks_nonzero_similarity, aes(x = max_similarity, y = rpm_log1p)) +
    stat_density_2d_filled(contour_var = "ndensity") +
    ylim(c(0, y_max)) +
   scale_fill_manual(values = color_palette) + scale_x_log10() +
    labs(
      x = "Max hit similarity",
      y = paste0("log1p ", gsub("\\.", " ", rpm_col )),
      title = "Elements with 1+ motif"
    ) +
    annotate("text", x = Inf, y = Inf, 
             label = paste0("Spearman r = ", round(spearman_r, 2), "\nPearson r = ", round(pearson_r, 2), "\nN = ", nrow(peaks_nonzero_similarity)),
             hjust = 1.1, vjust = 1.2, size = 3.5) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12),
      plot.title = element_text(color = "black", size = 12),
      axis.ticks = element_line(color = "black")
    )

  # Plot 2: Distribution of max similarity by peak overlap (violins with boxplots)  
  p2 <- ggplot(peaks_nonzero_similarity, aes(x = overlap_status, y = max_similarity)) +
    geom_violin(aes(fill = overlap_status), alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", alpha = 1, outliers = FALSE) +
    scale_fill_manual(values = c("No overlap" = "#c5cad7", "Overlap" = color)) +
    labs(
      x = "Element overlaps ChIP peak?",
      y = "Max hit similarity (nonzero)"
    ) +
    scale_y_log10() +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12),
      axis.ticks = element_line(color = "black")
    )

  # Plot 3: Distribution of RPM by peak overlap (violins with boxplots)
  p3 <- ggplot(peaks, aes(x = overlap_status, y = rpm_log1p)) +
    geom_violin(aes(fill = overlap_status), alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.1, color = "#000000", fill = "#ffffff", alpha = 1, outliers = FALSE) +
    scale_fill_manual(values = c("No overlap" = "#c5cad7", "Overlap" = color)) +
    labs(
      x = "Element overlaps ChIP peak?",
      y = paste0("log1p ", gsub("\\.", " ", rpm_col))
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12),
      axis.ticks = element_line(color = "black")
    )

  # Combine plots
  combined_plot <- plot_grid(p1, p1_nonzero, p2, p3, nrow = 1, align = "hv")
  
  # Save plot
  ggsave(out_file, combined_plot, width = 15, height = 4)
  
  return(combined_plot)
}

# Main execution
main <- function() {

  cp <- c(GATA = "#0096a0", AP1_1 = "#429130", REPEAT_G = "#6e788d",
    GATA_TAL1 = "#006eae", ETS_1 = "#ca9b23", STAT_2 = "#e96a00",
    CREB_ATF_3 = "#a64791", CTCF = "#96a008", CREB_ATF_1 = "#c5373d"
  )

  peaks <- peaks %>% 
    mutate(GATA.RPM = sqrt(GATA1.RPM + GATA2.RPM), GATA_peak_overlap = ifelse(GATA1_peak_overlap + GATA2_peak_overlap > 0, 1, 0),
      ATF4_CEBPG.RPM = sqrt(ATF4.RPM + CEBPG.RPM), ATF4_CEBPG_peak_overlap = ifelse(ATF4_peak_overlap + CEBPG_peak_overlap > 0, 1, 0),
      MAF.RPM = (MAFG.RPM + MAFF.RPM + MAFK.RPM) ^ (1/3),
        MAF_peak_overlap = ifelse(MAFF_peak_overlap + MAFG_peak_overlap + MAFK_peak_overlap > 0, 1, 0))

  plot_max_similarity_by_overlap(
    peaks = peaks,
    hits = hits,
    motif_names = c("GATA"),
    rpm_col = "GATA1.RPM",
    overlap_col = "GATA1_peak_overlap",
    color = "#0096a0",
    out_file = file.path(outdir, "GATA_motif.GATA1_ChIP.pdf")
  )

  plot_max_similarity_by_overlap(
    peaks = peaks,
    hits = hits,
    motif_names = c("GATA"),
    rpm_col = "GATA2.RPM",
    overlap_col = "GATA2_peak_overlap",
    color = "#0096a0",
    out_file = file.path(outdir, "GATA_motif.GATA2_ChIP.pdf")
  )

  plot_max_similarity_by_overlap(
    peaks = peaks,
    hits = hits,
    motif_names = c("GATA"),
    rpm_col = "GATA.RPM",
    overlap_col = "GATA_peak_overlap",
    color = "#006479",
    out_file = file.path(outdir, "GATA_motif.GATA1_GATA2_ChIP.pdf")
  )

  plot_max_similarity_by_overlap(
    peaks = peaks,
    hits = hits,
    motif_names = c("CREB_ATF_3"),
    rpm_col = "ATF4.RPM",
    overlap_col = "ATF4_peak_overlap",
    color = "#006eae",
    out_file = file.path(outdir, "CARE_motif.ATF4_ChIP.pdf")
  )

  plot_max_similarity_by_overlap(
    peaks = peaks,
    hits = hits,
    motif_names = c("CREB_ATF_3"),
    rpm_col = "CEBPG.RPM",
    overlap_col = "CEBPG_peak_overlap",
    color = "#006eae",
    out_file = file.path(outdir, "CARE_motif.CEBPG_ChIP.pdf")
  )

  plot_max_similarity_by_overlap(
    peaks = peaks,
    hits = hits,
    motif_names = c("CREB_ATF_3"),
    rpm_col = "ATF4_CEBPG.RPM",
    overlap_col = "ATF4_CEBPG_peak_overlap",
    color = "#00488d",
    out_file = file.path(outdir, "CARE_motif.ATF4_CEBPG_ChIP.pdf")
  )

  plot_max_similarity_by_overlap(
    peaks = peaks,
    hits = hits,
    motif_names = c("STAT_2"),
    rpm_col = "STAT5A.RPM",
    overlap_col = "STAT5A_peak_overlap",
    color = "#e96a00",
    out_file = file.path(outdir, "STAT_motif.STAT5A_ChIP.pdf")
  )

  plot_max_similarity_by_overlap(
    peaks = peaks,
    hits = hits,
    motif_names = c("CREB_ATF_3"),
    rpm_col = "JUN.RPM",
    overlap_col = "JUN_peak_overlap",
    color = "#429130",
    out_file = file.path(outdir, "CREB_ATF_3_motif.JUN_ChIP.pdf")
  )

  plot_max_similarity_by_overlap(
    peaks = peaks,
    hits = hits,
    motif_names = c("CREB_ATF_3"),
    rpm_col = "JUND.RPM",
    overlap_col = "JUND_peak_overlap",
    color = "#429130",
    out_file = file.path(outdir, "CREB_ATF_3_motif.JUND_ChIP.pdf")
  )

  plot_max_similarity_by_overlap(
    peaks = peaks,
    hits = hits,
    motif_names = c("CREB_ATF_3"),
    rpm_col = "MAF.RPM",
    overlap_col = "MAF_peak_overlap",
    color = "#429130",
    out_file = file.path(outdir, "CREB_ATF_3_motif.MAF_ChIP.pdf")
  )

  plot_max_similarity_by_overlap(
    peaks = peaks,
    hits = hits,
    motif_names = c("AP1_1"),
    rpm_col = "MAFG.RPM",
    overlap_col = "MAFG_peak_overlap",
    color = "#435369",
    out_file = file.path(outdir, "AP1_1_motif.MAFG_ChIP.pdf")
  )

  plot_max_similarity_by_overlap(
    peaks = peaks,
    hits = hits,
    motif_names = c("AP1_1"),
    rpm_col = "MAF.RPM",
    overlap_col = "MAF_peak_overlap",
    color = "#435369",
    out_file = file.path(outdir, "AP1_1_motif.MAF_ChIP.pdf")
  )

  plot_max_similarity_by_overlap(
    peaks = peaks,
    hits = hits,
    motif_names = c("AP1_1"),
    rpm_col = "JUND.RPM",
    overlap_col = "JUND_peak_overlap",
    color = "#435369",
    out_file = file.path(outdir, "AP1_1_motif.JUND_ChIP.pdf")
  )

}





# Run main function
main()

