
# Rscript scripts/5.3.plot_finemo_hits.R 2025_0703_retrain_p300_model/finemo/pkw_500_curated_motifs/annotated_motifs/hits_renamed.tsv 2025_0703_retrain_p300_model/finemo/pkw_500_curated_motifs/annotated_motifs/hits_per_peak.with_predictions.tsv.gz 2025_0703_retrain_p300_model/finemo/pkw_500_curated_motifs/plot_annotated_motifs
# Rscript scripts/5.3.plot_finemo_hits.R 2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs/annotated_motifs/hits_renamed.tsv 2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs/annotated_motifs/hits_per_peak.with_predictions.tsv.gz 2025_0517_official_EP300_K562_model/finemo/pkw_500_curated_motifs/plot_annotated_motifs

#!/usr/bin/env Rscript

# R script for motif frequency and importance analysis

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

# Create output directory
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load data
hits <- read.table(hits_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
peaks <- read.table(peaks_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Function 1: Bar plot of motif frequency and importance
plot_motif_frequency_importance <- function(hits, cp, output_path) {
  
  # Calculate frequency and sum of importance for each motif
  motif_stats <- hits %>%
    group_by(motif_name) %>%
    summarise(
      frequency = n(),
      total_importance = sum(hit_importance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(motif_name = factor(motif_name, levels = rev(names(cp)), ordered = TRUE))
  
  # Create dual-axis plot
  p1 <- ggplot(motif_stats, aes(x = frequency, y = motif_name)) +
    geom_col(fill = "#006eae", alpha = 0.7) +
    labs(
      x = "Frequency (count)",
      y = "Motif name",
      title = "Motif frequency"
    ) +
    scale_x_continuous(labels=comma_format()) +
    theme_classic() +
    theme(
      axis.text = element_text(color = "#000000", size = 12),
      axis.title = element_text(color = "#000000", size = 14),
      plot.title = element_text(color = "#000000", size = 14, hjust = 0.5)
    )
  
  p2 <- ggplot(motif_stats, aes(x = total_importance, y = motif_name)) +
    geom_col(fill = "#792374", alpha = 0.7) +
    labs(
      x = "Total hit importance",
      y = "",
      title = "Total importance"
    ) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text = element_text(color = "#000000", size = 12),
      axis.title = element_text(color = "#000000", size = 14),
      plot.title = element_text(color = "#000000", size = 14, hjust = 0.5)
    )
  
  # Combine plots
  combined_plot <- plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 0.8), align = "h")
  
  # Save plot
  ggsave(output_path, combined_plot, width = 8, height = max(6, nrow(motif_stats) * 0.2))
  
  return(motif_stats)
}

# Function 2: Line chart of motif frequency vs P300 signal percentiles
plot_motif_frequency_by_percentile <- function(hits, peaks, cp, output_path) {
  
  # Define percentile thresholds (top 100%, 90%, 80%, ..., 10%)
  percentiles <- seq(10, 100, by = 10)
  
  # Get motif frequencies for each percentile threshold
  frequency_data <- data.frame()
  
  for (pct in percentiles) {
    # Calculate threshold for this percentile
    threshold <- quantile(peaks$true_logcounts, (100 - pct) / 100, na.rm = TRUE)
    
    # Filter peaks above threshold
    peaks_filtered <- peaks[peaks$true_logcounts >= threshold, ]
    
    # Get hits in these peaks
    hits_filtered <- hits[hits$peak_id %in% peaks_filtered$peak_id, ]
    
    # Calculate motif frequencies
    motif_freq <- hits_filtered %>%
      count(motif_name) %>%
      mutate(
        percentile = pct,
        relative_frequency = n / sum(n)
      )
    
    frequency_data <- rbind(frequency_data, motif_freq)
  }
  
  # Filter to only include motifs that appear in at least 3 percentile bins
  motifs_to_keep <- frequency_data %>%
    group_by(motif_name) %>%
    summarise(n_percentiles = n(), .groups = "drop") %>%
    filter(n_percentiles >= 3) %>%
    pull(motif_name)
  
  frequency_data <- frequency_data %>%
    filter(motif_name %in% motifs_to_keep)
  
  # Create line plot
  p <- ggplot(frequency_data, aes(x = percentile, y = relative_frequency, color = motif_name)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 3, alpha = 1, shape = 16) +
    scale_x_continuous(
      breaks = seq(10, 100, by = 10),
      labels = paste0("Top ", seq(10, 100, by = 10), "%")
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_color_manual(values = cp, name = "Motif") +
    labs(
      x = "Percentile of p300 signal in element",
      y = "Relative frequency of motifs",
      title = "Motif frequency by p300 signal percentile"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(color = "#000000", angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(color = "#000000", size = 12),
      axis.title = element_text(color = "#000000", size = 12),
      plot.title = element_text(color = "#000000", size = 14, hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(color = "#000000", size = 10),
      legend.title = element_text(color = "#000000", size = 12)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  # Save plot
  ggsave(output_path, p, width = 6, height = 4)
  
  return(frequency_data)
}

# Function 2.1: Line chart of motif scores vs P300 signal percentiles
plot_motif_score_by_percentile <- function(hits, peaks, cp, output_path) {
  
  # Define percentile thresholds (top 100%, 90%, 80%, ..., 10%)
  percentiles <- seq(10, 100, by = 10)
  
  # Get motif frequencies for each percentile threshold
  avg_data <- data.frame()
  
  for (pct in percentiles) {
    # Calculate threshold for this percentile
    threshold <- quantile(peaks$true_logcounts, (100 - pct) / 100, na.rm = TRUE)
    
    # Filter peaks above threshold
    peaks_filtered <- peaks[peaks$true_logcounts >= threshold, ]
    
    # Get hits in these peaks
    hits_filtered <- hits[hits$peak_id %in% peaks_filtered$peak_id, ]
    
    # Calculate motif frequencies
    motif_avg <- hits_filtered %>%
      group_by(motif_name) %>% 
      summarize(
        avg_importance = mean(hit_importance),
        avg_similarity = mean(hit_similarity),
        avg_coef_global = mean(hit_coefficient_global)
      ) %>%
      mutate(percentile = pct)
    
    avg_data <- rbind(motif_avg, avg_data)
  }
  
  # Create line plot
  p1 <- ggplot(avg_data, aes(x = percentile, y = avg_coef_global, color = motif_name)) +
    geom_line(size = 1, alpha = 0.8) +
    geom_point(size = 3, alpha = 1, shape = 16) +
    scale_x_continuous(
      breaks = seq(10, 100, by = 10),
      labels = paste0("Top ", seq(10, 100, by = 10), "%")
    ) +
    scale_color_manual(values = cp, name = "Motif") +
    labs(
      x = "Percentile of p300 signal in element",
      y = "Average hit importance",
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(color = "#000000", angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(color = "#000000", size = 12),
      axis.title = element_text(color = "#000000", size = 12),
      legend.position = "none",
    )

  p2 <- ggplot(avg_data, aes(x = percentile, y = avg_importance, color = motif_name)) +
    geom_line(size = 1, alpha = 0.8) +
    geom_point(size = 3, alpha = 1, shape = 16) +
    scale_x_continuous(
      breaks = seq(10, 100, by = 10),
      labels = paste0("Top ", seq(10, 100, by = 10), "%")
    ) +
    scale_color_manual(values = cp, name = "Motif") +
    labs(
      x = "Percentile of p300 signal in element",
      y = "Average hit importance",
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(color = "#000000", angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(color = "#000000", size = 12),
      axis.title = element_text(color = "#000000", size = 12),
      legend.position = "none",
    )

  p3 <- ggplot(avg_data, aes(x = percentile, y = avg_similarity, color = motif_name)) +
    geom_line(size = 1, alpha = 0.8) +
    geom_point(size = 3, alpha = 1, shape = 16) +
    scale_x_continuous(
      breaks = seq(10, 100, by = 10),
      labels = paste0("Top ", seq(10, 100, by = 10), "%")
    ) +
    scale_color_manual(values = cp, name = "Motif") +
    labs(
      x = "Percentile of p300 signal in element",
      y = "Average similarity of motifs",
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(color = "#000000", angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(color = "#000000", size = 12),
      axis.title = element_text(color = "#000000", size = 12),
      legend.position = "right",
      legend.text = element_text(color = "#000000", size = 10),
      legend.title = element_text(color = "#000000", size = 12)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))

  # Save plot
  grid <- cowplot::plot_grid(p1, p2, p3, ncol = 1, nrow = 3, align = "hv")
  ggsave(output_path, grid, width = 6, height = 10)
  
  return(avg_data)
}

# Function 3: Bubble plot combining frequency and importance
plot_motif_bubble <- function(hits, cp, output_path) {
  
  # Calculate motif statistics
  motif_stats <- hits %>%
    group_by(motif_name) %>%
    summarise(
      frequency = n(),
      avg_importance = mean(hit_importance, na.rm = TRUE),
      total_importance = sum(hit_importance, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create bubble plot
  p <- ggplot(motif_stats, aes(x = frequency, y = avg_importance)) +
    geom_point(aes(size = total_importance, color = motif_name), alpha = 0.7) +
    scale_size_continuous(
      name = "Total\nimportance",
      range = c(2, 15),
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    scale_color_manual(values = cp, name = "Motif", guide = "none") + 
    scale_x_log10(labels = comma_format()) +
    scale_y_log10(labels = comma_format()) +
    labs(
      x = "Frequency (log10)",
      y = "Average hit importance (log10)",
      title = "Motif frequency vs. importance"
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(color = "#000000", size = 12),
      axis.title = element_text(color = "#000000", size = 14),
      plot.title = element_text(color = "#000000", size = 14, hjust = 0.5),
      legend.position = "right"
    )
  
  # Add text labels for top motifs
  top_motifs <- motif_stats %>%
    arrange(desc(total_importance)) %>%
    head(8)
  
  p <- p + 
    ggrepel::geom_text_repel(
      data = top_motifs,
      aes(label = motif_name),
      size = 3,
      max.overlaps = 15,
      box.padding = 0.3,
      point.padding = 0.3
    )
  
  # Save plot
  ggsave(output_path, p, width = 5, height = 4)
  
  return(motif_stats)
}

# Function 4: Heatmap of motif importance by P300 signal bins
plot_motif_heatmap_by_signal <- function(hits, peaks, cp, colname, collabel, output_path) {
  
  # Create P300 signal bins
  peaks$signal_bin <- cut(
    peaks$true_logcounts,
    breaks = quantile(peaks$true_logcounts, probs = seq(0, 1, by = 0.2), na.rm = TRUE),
    labels = c("Bottom 20%", "20-40%", "40-60%", "60-80%", "Top 20%"),
    include.lowest = TRUE
  )
  
  # Merge hits with peak signal bins
  hits_with_bins <- hits %>%
    left_join(peaks[, c("peak_id", "signal_bin")], by = "peak_id") %>%
    filter(!is.na(signal_bin))
  
  # Calculate average importance for each motif-bin combination
  heatmap_data <- hits_with_bins %>%
    group_by(motif_name, signal_bin) %>%
    summarise(
      avg_value = mean(!!sym(colname), na.rm = TRUE),
      count = n(),
      .groups = "drop"
    ) %>%
    filter(count >= 5)  # Only include combinations with at least 5 hits
  
  # Filter to top motifs by overall frequency
  top_motifs <- hits %>%
    count(motif_name, sort = TRUE) %>%
    head(20) %>%
    pull(motif_name)
  
  heatmap_data <- heatmap_data %>%
    filter(motif_name %in% top_motifs) %>% 
    mutate(motif_name = factor(motif_name, levels = rev(names(cp)), ordered = TRUE))
  
  # Create heatmap
  cp = ifelse(colname == "hit_importance", "mako", "rocket")
  p <- ggplot(heatmap_data, aes(x = signal_bin, y = motif_name)) +
    geom_tile(aes(fill = avg_value), color = "white", size = 0.5) +
    scale_fill_viridis_c(
      option = cp,
      name = collabel,
      trans = "log10",
      labels = comma_format()
    ) +
    labs(
      x = "p300 signal bin",
      y = "Motif name",
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      axis.text.x = element_text(color = "#000000", angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(color = "#000000", size = 12),
      axis.title = element_text(color = "#000000", size = 14),
      plot.title = element_text(color = "#000000", size = 14, hjust = 0.5),
      legend.position = "right"
    )
  
  # Save plot
  ggsave(output_path, p, width = 6, height = 6)
  
  return(heatmap_data)
}

# Function 4: Per-peak analysis with hit distributions and boxplots
plot_per_peak_analysis <- function(hits, peaks, output_path) {
  
  # Define quantile thresholds (all hits, top 80%, 60%, 40%, 20%)
  quantiles <- c(1.0, 0.8, 0.6, 0.4, 0.2)
  quantile_metric <- "hit_coefficient_global"
  max_n <- 9
  
  # Check if predictions are available
  has_predictions <- all(c("true_logcounts", "mean_pred_logcounts") %in% colnames(peaks))
  
  # Create list to store plots
  plot_list <- list()
  
  for (i in seq_along(quantiles)) {
    top_pct <- quantiles[i]
    
    # Calculate threshold for this percentile
    threshold <- quantile(hits[[quantile_metric]], 1.0 - top_pct, na.rm = TRUE)
    
    # Filter hits above threshold
    hits_filtered <- hits[hits[[quantile_metric]] >= threshold, ]
    
    # Recompute hits per peak for this subset
    peak_hits_summary <- hits_filtered %>%
      group_by(peak_id) %>%
      summarise(n_hits = n(), .groups = "drop")
    
    # Merge with peaks data
    peaks_subset <- peaks %>% select(-n_hits) %>%
      left_join(peak_hits_summary, by = "peak_id") %>%
      mutate(n_hits = ifelse(is.na(n_hits), 0, n_hits)) %>% 
      filter(n_hits <= max_n)
    
    # Create label
    subset_label <- ifelse(top_pct == 1.0, 
                          paste0("All hits (n=", nrow(hits_filtered), ")"),
                          paste0("Top ", round(top_pct * 100), "% (n=", nrow(hits_filtered), ")"))
    
    # Plot 1: Distribution of hits per peak
    peak_freq <- peaks_subset %>% 
      group_by(n_hits) %>% 
      summarize(freq = n())

    p1 <- ggplot(peak_freq, aes(x = n_hits, y = freq)) +
      geom_col(fill = "#b778b3", color = "#430b4e") +
      scale_x_continuous(n.breaks=10) +
      labs(x = "Motif hits per peak", y = "Number of peaks", title = subset_label) +
      theme_classic() +
      theme(
        axis.text = element_text(color = "#000000", size = 9),
        axis.title = element_text(color = "#000000", size = 10),
        plot.title = element_text(color = "#000000", size = 10)
      )
    
    # Plot 2: True counts vs hits
    p2 <- ggplot(peaks_subset %>% filter(n_hits <= 10), 
                 aes(x = factor(n_hits), y = true_logcounts)) +
      geom_boxplot(fill = "#b778b3", alpha = 0.7, outlier.size = 0.5) +
      scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 3)) +
      labs(x = "Motif hits per peak", y = "Observed log counts", title = "") +
      theme_classic() +
      theme(
        axis.text = element_text(color = "#000000", size = 9),
        axis.title = element_text(color = "#000000", size = 10)
      )
    
    if (has_predictions) {
      # Plot 3: Predicted counts vs hits
      p3 <- ggplot(peaks_subset %>% filter(n_hits <= 10), 
                   aes(x = factor(n_hits), y = mean_pred_logcounts)) +
        geom_boxplot(fill = "#b778b3", alpha = 0.7, outlier.size = 0.2, outlier.alpha = 0.5) +
        scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 3)) +
        labs(x = "Motif hits per peak", y = "Predicted log counts", title = "") +
        theme_classic() +
        theme(
          axis.text = element_text(color = "#000000", size = 9),
          axis.title = element_text(color = "#000000", size = 10)
        )
      
      # Calculate delta if not present
      if (!"delta_logcounts" %in% colnames(peaks_subset)) {
        peaks_subset$delta_logcounts <- peaks_subset$mean_pred_logcounts - peaks_subset$true_logcounts
      }
      
      # Plot 4: Delta vs hits
      p4 <- ggplot(peaks_subset %>% filter(n_hits <= 10), 
                   aes(x = factor(n_hits), y = delta_logcounts)) +
        geom_boxplot(fill = "#b778b3", alpha = 0.7, outlier.size = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "#9b241c") +
        scale_y_continuous(limits = c(-8, 8), breaks = seq(-9, 9, 3)) +
        labs(x = "Motif hits per peak", y = "Predicted - observed log counts", title = "") +
        theme_classic() +
        theme(
          axis.text = element_text(color = "#000000", size = 9),
          axis.title = element_text(color = "#000000", size = 10)
        )
      
      # Combine 4 plots per row
      row_plot <- plot_grid(p1, p2, p3, p4, nrow = 1, rel_widths = c(1, 1, 1, 1))
      
    } else {
      # Only 2 plots if no predictions
      row_plot <- plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1))
    }
    
    plot_list[[i]] <- row_plot
  }
  
  # Combine all rows
  final_plot <- plot_grid(plotlist = plot_list, ncol = 1)
  
  # Save plot
  plot_height <- ifelse(has_predictions, 2.2 * length(quantiles), 4)
  plot_width <- ifelse(has_predictions, 8, 6)
  
  ggsave(output_path, final_plot, width = plot_width, height = plot_height)
  
  return(invisible())
}

# Main execution
main <- function() {
  cat("Loading data...\n")
  cat("Hits data dimensions:", nrow(hits), "x", ncol(hits), "\n")
  cat("Peaks data dimensions:", nrow(peaks), "x", ncol(peaks), "\n")

  cp <- c(GATA = "#0096a0", AP1_1 = "#429130", REPEAT_G = "#6e788d",
    GATA_TAL1 = "#006eae", ETS_1 = "#ca9b23", STAT_2 = "#e96a00",
    CREB_ATF_3 = "#a64791", CTCF = "#96a008", CREB_ATF_1 = "#c5373d"
  )
  
  # Generate plots
  cat("Creating motif frequency and importance bar plot...\n")
  plot_motif_frequency_importance(
    hits = hits,
    cp = cp,
    output_path = file.path(outdir, "motif_frequency_importance.pdf")
  )
  
  cat("Creating motif frequency by percentile line plot...\n")
  plot_motif_frequency_by_percentile(
    hits = hits,
    peaks = peaks,
    cp = cp,
    output_path = file.path(outdir, "motif_frequency_by_percentile.pdf")
  )

  plot_motif_score_by_percentile(
    hits = hits,
    peaks = peaks,
    cp = cp,
    output_path = file.path(outdir, "motif_scores_by_percentile.pdf")
  )
  
  cat("Creating motif bubble plot...\n")
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    plot_motif_bubble(
      hits = hits,
      cp = cp,
      output_path = file.path(outdir, "motif_bubble_plot.pdf")
    )
  } else {
    cat("Note: ggrepel package not available, skipping bubble plot with labels\n")
  }
  
  cat("Creating motif importance heatmap by signal...\n")
  plot_motif_heatmap_by_signal(
    hits = hits,
    peaks = peaks,
    cp = cp,
    colname = "hit_importance",
    collabel = "Average\nimportance",
    output_path = file.path(outdir, "motif_importance_heatmap.pdf")
  )

  cat("Creating motif similarity heatmap by signal...\n")
  plot_motif_heatmap_by_signal(
    hits = hits,
    peaks = peaks,
    cp = cp,
    colname = "hit_similarity",
    collabel = "Average\nsimilarity",
    output_path = file.path(outdir, "motif_similarity_heatmap.pdf")
  )

  cat("Creating motif counts per peak plot...\n")
  plot_per_peak_analysis(
    hits = hits,
    peaks = peaks,
    output_path = file.path(outdir, "counts_per_peak.pdf"))
  
  cat("All plots completed successfully!\n")
}



# Run main function
main()