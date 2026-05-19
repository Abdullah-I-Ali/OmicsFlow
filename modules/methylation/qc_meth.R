# ==============================================================================
# qc_meth.R — QC Plots & Statistical Validation Figures for Methylation Module
# OmicsFlow | Phase 1: Methylation Module
# Scientific source: data/full_scripts.R (lines 1124-1147 & 1187-1270)
# ==============================================================================

#' Generate technical QC plots for the methylation preprocessing pipeline
#'
#' Produces:
#'   plots/pca_after_combat.pdf    — PCA colored by Batch
#'   plots/beta_distribution.pdf   — Beta value distribution
generate_meth_qc_plots <- function(pca_final, pca_var, batch_vec, met_beta_top, outdir) {
  plots_dir <- file.path(outdir, "plots")
  meth_msg("Generating QC plots...")

  # PCA plot
  batch_for_plot <- batch_vec[rownames(pca_final$x)]
  pca_df <- data.frame(PC1 = pca_final$x[,1], PC2 = pca_final$x[,2], Batch = batch_for_plot)
  p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Batch)) + 
    geom_point(alpha = 0.8, size = 1.5) +
    theme_minimal() + theme(legend.position = "none") +
    labs(title = "PCA after ComBat", 
         x = sprintf("PC1 (%.1f%%)", pca_var[2,1]*100),
         y = sprintf("PC2 (%.1f%%)", pca_var[2,2]*100))
  
  ggsave(file.path(plots_dir, "pca_after_combat.pdf"), p_pca, width = 6, height = 5)
  ggsave(file.path(plots_dir, "pca_after_combat.png"), p_pca, width = 6, height = 5, dpi = 300)
  meth_msg("Saved: plots/pca_after_combat.pdf and .png")

  # Beta density
  set.seed(42)
  n_probes_plot <- min(500, nrow(met_beta_top))
  sample_probes <- sample(rownames(met_beta_top), n_probes_plot)
  beta_long <- data.frame(beta = as.vector(met_beta_top[sample_probes, ]))
  p_beta <- ggplot(beta_long, aes(x = beta)) + 
    geom_density(fill = "steelblue", alpha = 0.5) +
    theme_minimal() + 
    labs(title = "Beta Value Distribution", x = "Beta", y = "Density")
  
  ggsave(file.path(plots_dir, "beta_distribution.pdf"), p_beta, width = 6, height = 4)
  ggsave(file.path(plots_dir, "beta_distribution.png"), p_beta, width = 6, height = 4, dpi = 300)
  meth_msg("Saved: plots/beta_distribution.pdf and .png")
  
  invisible(NULL)
}

#' Generate publication-ready statistical validation figures
#'
#' Produces:
#'   plots/Methylation_Research_Validation_Figures.png (3 panels)
generate_meth_validation_figures <- function(meth_m, meth_beta, outdir) {
  plots_dir <- file.path(outdir, "plots")
  meth_msg("Generating publication-ready validation figures...")

  set.seed(42)
  sample_probes <- sample(rownames(meth_beta), min(10000, nrow(meth_beta)))

  # Figure 1: Beta Value Distribution (Bimodal Check)
  beta_melt <- melt(meth_beta[sample_probes, 1:min(20, ncol(meth_beta))])
  colnames(beta_melt) <- c("Probe", "Sample", "Beta")

  fig_beta <- ggplot(beta_melt, aes(x = Beta, group = Sample)) +
    geom_density(color = "steelblue", alpha = 0.1) +
    theme_bw() +
    labs(title = "DNA Methylation Beta Value Distribution (Top 20 Samples)",
         subtitle = "Showing the characteristic biological bimodal shape (unmethylated vs methylated)",
         x = "Beta Value (0 to 1)", y = "Density")

  # Figure 2: M-Value Distribution (Gaussian Check for ML)
  m_melt <- melt(meth_m[sample_probes, 1:min(20, ncol(meth_m))])
  colnames(m_melt) <- c("Probe", "Sample", "M_Value")

  fig_m <- ggplot(m_melt, aes(x = M_Value, group = Sample)) +
    geom_density(color = "#e74c3c", alpha = 0.1) +
    theme_bw() +
    labs(title = "Transformed M-Value Distribution (Top 20 Samples)",
         subtitle = "Gaussian-like distribution optimized for Machine Learning assumptions",
         x = "M-Value [ log2(Beta / (1-Beta)) ]", y = "Density")

  # Figure 3: Inter-Sample Correlation Histogram
  var_probes <- rownames(meth_m)[order(rowVars(meth_m), decreasing = TRUE)[1:min(5000, nrow(meth_m))]]
  sample_cor <- cor(meth_m[var_probes, ], method = "pearson")
  cor_values <- sample_cor[lower.tri(sample_cor)]

  fig_cor <- ggplot(data.frame(Cor = cor_values), aes(x = Cor)) +
    geom_histogram(binwidth = 0.01, fill = "#2ecc71", color = "white") +
    theme_classic() +
    labs(title = "Inter-Sample Pearson Correlation (Methylation Layer)",
         subtitle = paste0("Based on Top ", length(var_probes), " Variable Probes | Median Correlation: ", round(median(cor_values), 3)),
         x = "Pearson R", y = "Count")

  # Combine and Save
  final_meth_plot <- grid.arrange(fig_beta, fig_m, fig_cor, ncol = 1)

  ggsave(file.path(plots_dir, "Methylation_Research_Validation_Figures.png"), 
         plot = final_meth_plot, width = 8, height = 12, dpi = 300)
  meth_msg("Saved: plots/Methylation_Research_Validation_Figures.png")

  # Statistical Summary
  pca_res <- prcomp(t(meth_m[var_probes, ]), scale. = TRUE)
  pca_var <- summary(pca_res)$importance

  meth_msg("--- Statistical Report ---", level = "DETAIL")
  cat(sprintf("  1. Total High-Quality Patients Matched : %d\n", ncol(meth_m)))
  cat(sprintf("  2. Total Clean Probes Retained         : %s\n", format(nrow(meth_m), big.mark = ",")))
  cat(sprintf("  3. Median Inter-sample Correlation (R) : %.4f\n", median(cor_values)))
  cat(sprintf("  4. Mean M-Value Range                  : [%.2f to %.2f]\n", min(rowMeans(meth_m)), max(rowMeans(meth_m))))
  cat(sprintf("  5. Variance Explained by PC1           : %.2f%%\n", pca_var[2,1]*100))

  invisible(NULL)
}
