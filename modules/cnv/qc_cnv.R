# ==============================================================================
# qc_cnv.R — Quality Control & Visualization for CNV
# OmicsFlow — CNV Module
# ==============================================================================

source("modules/cnv/utils_cnv.R")

generate_cnv_qc_plots <- function(cnv_matrix, outdir) {
  cnv_step("QC", "Generating QC Visualizations")
  
  plot_dir <- file.path(outdir, "plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  # 1. Density Plot
  set.seed(42)
  sample_cols <- sample(colnames(cnv_matrix), min(20, ncol(cnv_matrix)))
  cnv_melt <- reshape2::melt(cnv_matrix[, sample_cols])
  colnames(cnv_melt) <- c("Gene", "Sample", "CNV_Value")
  
  fig_cnv_dist <- ggplot(cnv_melt, aes(x = CNV_Value, group = Sample)) +
    geom_density(color = "darkblue", alpha = 0.2) +
    theme_bw() +
    labs(title = "CNV Log2 Copy Ratio Distribution (Subset of 20 Samples)",
         subtitle = "Values centered at 0 (Neutral Diploid state), clipped at \u00b15",
         x = "Log2 Copy Ratio", y = "Density") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1)
  
  ggsave(file.path(plot_dir, "cnv_density.pdf"), fig_cnv_dist, width = 8, height = 5)
  ggsave(file.path(plot_dir, "cnv_density.png"), fig_cnv_dist, width = 8, height = 5, dpi = 300)
  cnv_msg("Saved: plots/cnv_density.pdf and .png")
  
  # 2. Correlation Histogram
  sample_cor <- cor(cnv_matrix, method = "pearson")
  cor_values <- sample_cor[lower.tri(sample_cor)]
  
  fig_cor <- ggplot(data.frame(Cor = cor_values), aes(x = Cor)) +
    geom_histogram(binwidth = 0.02, fill = "#8e44ad", color = "white") +
    theme_classic() +
    labs(title = "Inter-Sample Pearson Correlation (CNV Layer)",
         subtitle = paste0("Median Correlation: ", round(median(cor_values, na.rm=TRUE), 3)),
         x = "Pearson R", y = "Count")
         
  ggsave(file.path(plot_dir, "cnv_correlation_histogram.pdf"), fig_cor, width = 8, height = 5)
  ggsave(file.path(plot_dir, "cnv_correlation_histogram.png"), fig_cor, width = 8, height = 5, dpi = 300)
  cnv_msg("Saved: plots/cnv_correlation_histogram.pdf and .png")
  
  # 3. PCA Plot
  pca_res <- prcomp(t(cnv_matrix), scale. = FALSE)
  pca_var <- summary(pca_res)$importance
  
  df_pca <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2])
  fig_pca <- ggplot(df_pca, aes(x = PC1, y = PC2)) +
    geom_point(color = "#2980b9", size = 2, alpha = 0.6) +
    theme_minimal() +
    labs(title = "PCA of Copy Number Variation Profiles",
         x = sprintf("PC1 (%.1f%%)", pca_var[2,1]*100),
         y = sprintf("PC2 (%.1f%%)", pca_var[2,2]*100))
         
  ggsave(file.path(plot_dir, "cnv_pca.pdf"), fig_pca, width = 6, height = 5)
  ggsave(file.path(plot_dir, "cnv_pca.png"), fig_pca, width = 6, height = 5, dpi = 300)
  cnv_msg("Saved: plots/cnv_pca.pdf and .png")
  
  # 4. Boxplot Distribution
  pdf(file.path(plot_dir, "cnv_distribution_per_sample.pdf"), width = 14, height = 6)
  n_plot <- min(50, ncol(cnv_matrix))
  boxplot(
    cnv_matrix[, 1:n_plot, drop = FALSE],
    las = 2, cex.axis = 0.5,
    col = "steelblue", border = "gray30",
    main = sprintf("CNV Distribution (log2 Copy Ratio) - %d Samples", n_plot),
    ylab = "log2 Copy Ratio",
    outline = FALSE
  )
  abline(h = 0, col = "red", lty = 2, lwd = 1.5)
  legend("topright", legend = "Neutral (0)", col = "red", lty = 2, cex = 0.8)
  dev.off()
  png(file.path(plot_dir, "cnv_distribution_per_sample.png"), width = 14, height = 6, units = "in", res = 300)
  boxplot(
    cnv_matrix[, 1:n_plot, drop = FALSE],
    las = 2, cex.axis = 0.5,
    col = "steelblue", border = "gray30",
    main = sprintf("CNV Distribution (log2 Copy Ratio) - %d Samples", n_plot),
    ylab = "log2 Copy Ratio",
    outline = FALSE
  )
  abline(h = 0, col = "red", lty = 2, lwd = 1.5)
  legend("topright", legend = "Neutral (0)", col = "red", lty = 2, cex = 0.8)
  dev.off()
  cnv_msg("Saved: plots/cnv_distribution_per_sample.pdf and .png")
  
  # 5. Combined Validation Figure
  final_cnv_plot <- grid.arrange(fig_cnv_dist, fig_cor, fig_pca, ncol = 1)
  ggsave(file.path(plot_dir, "CNV_Research_Validation_Figures.png"),
         plot = final_cnv_plot, width = 8, height = 12, dpi = 300)
  cnv_msg("Saved: plots/CNV_Research_Validation_Figures.png")
}
