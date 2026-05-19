# ==============================================================================
# qc_rna.R — QC Plots & Statistical Validation Figures for the RNA Module
# OmicsFlow | Phase 1: RNA Module
# Scientific source: data/full_scripts.R (lines 462-628)
# ==============================================================================
# PURPOSE:
#   Two exported functions, sourced by preprocess_rna.R:
#
#   generate_rna_qc_plots()          — Step 13: limma::plotDensities,
#                                      pheatmap sample correlation heatmap,
#                                      PCA before/after batch correction,
#                                      zero-proportion check.
#
#   generate_rna_validation_figures() — Publication-ready ggplot2 figures:
#                                       Figure 1: Expression density (top 20 samples)
#                                       Figure 2: Mean-SD relationship (CV analysis)
#                                       Figure 3: Inter-sample correlation histogram
#
# SCIENTIFIC RULES:
#   - plotDensities uses rainbow colors, max 50 samples, no legend.
#   - pheatmap uses blue-white-red palette, no row/col names.
#   - PCA uses scale. = FALSE (data already on log2 CPM scale).
#   - Validation figures use exactly 20 samples for density, top-20 subset.
#   - Correlation histogram uses lower triangle of full Pearson matrix.
#   - CV figure uses GAM smoother (geom_smooth method = "gam").
#   - All figures saved to outdir/plots/.
# ==============================================================================


# ==============================================================================
# FUNCTION: generate_rna_qc_plots
# ==============================================================================

#' Generate technical QC plots for the RNA preprocessing pipeline
#'
#' Produces:
#'   plots/pca_before_after_batch.pdf  — PCA before and after correction
#'   plots/density_log2cpm.pdf         — limma::plotDensities
#'   plots/sample_correlation.pdf      — pheatmap correlation heatmap
#'   Prints zero-proportion statistic to console.
#'
#' @param rna_top    Numeric matrix — log2 CPM, top variable genes x samples
#' @param df_before  data.frame with PC1, PC2, Batch — before correction
#' @param df_after   data.frame with PC1, PC2, Batch — after correction
#' @param pca_before prcomp object — before correction (for variance labels)
#' @param pca_after  prcomp object — after correction (for variance labels)
#' @param batch_info factor — plate batch labels aligned to rna_top columns
#' @param outdir     Character — base output directory
generate_rna_qc_plots <- function(rna_top, df_before, df_after,
                                   pca_before, pca_after, batch_info,
                                   outdir) {

  plots_dir <- file.path(outdir, "plots")

  # --- PCA Before / After Batch Correction ---
  rna_msg("Generating PCA before/after batch correction...")

  var_before <- summary(pca_before)$importance
  var_after  <- summary(pca_after)$importance

  p1 <- ggplot(df_before, aes(PC1, PC2, color = Batch)) +
    geom_point(size = 2, alpha = 0.7) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      title    = "PCA Before Batch Correction",
      subtitle = "Colored by Plate ID",
      x = sprintf("PC1 (%.1f%%)", var_before[2, 1] * 100),
      y = sprintf("PC2 (%.1f%%)", var_before[2, 2] * 100)
    )

  p2 <- ggplot(df_after, aes(PC1, PC2, color = Batch)) +
    geom_point(size = 2, alpha = 0.7) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      title    = "PCA After Batch Correction",
      subtitle = "Batch effects should be reduced",
      x = sprintf("PC1 (%.1f%%)", var_after[2, 1] * 100),
      y = sprintf("PC2 (%.1f%%)", var_after[2, 2] * 100)
    )

  pca_combined <- grid.arrange(p1, p2, ncol = 2)
  ggsave(
    file.path(plots_dir, "pca_before_after_batch.pdf"),
    plot   = pca_combined,
    width  = 12,
    height = 5
  )
  ggsave(
    file.path(plots_dir, "pca_before_after_batch.png"),
    plot   = pca_combined,
    width  = 12,
    height = 5,
    dpi    = 300
  )
  rna_msg("Saved: plots/pca_before_after_batch.pdf and .png")

  # --- Density Plot (limma::plotDensities) ---
  rna_msg("Generating density plot (limma::plotDensities)...")

  pdf(file.path(plots_dir, "density_log2cpm.pdf"), width = 9, height = 6)
  plotDensities(
    rna_top,
    legend = FALSE,
    main   = "Density Plot: log2 CPM (After Batch Correction)",
    col    = rainbow(min(ncol(rna_top), 50), alpha = 0.3)
  )
  dev.off()
  png(file.path(plots_dir, "density_log2cpm.png"), width = 9, height = 6, units = "in", res = 300)
  plotDensities(
    rna_top,
    legend = FALSE,
    main   = "Density Plot: log2 CPM (After Batch Correction)",
    col    = rainbow(min(ncol(rna_top), 50), alpha = 0.3)
  )
  dev.off()
  rna_msg("Saved: plots/density_log2cpm.pdf and .png")

  # --- Sample Correlation Heatmap (pheatmap) ---
  rna_msg("Generating sample correlation heatmap...")

  sample_cor_final <- cor(rna_top, method = "pearson")

  pdf(file.path(plots_dir, "sample_correlation_heatmap.pdf"), width = 8, height = 7)
  pheatmap(
    sample_cor_final,
    main          = "Sample Correlation Heatmap (Final)",
    show_colnames = FALSE,
    show_rownames = FALSE,
    color         = colorRampPalette(c("blue", "white", "red"))(100)
  )
  dev.off()
  png(file.path(plots_dir, "sample_correlation_heatmap.png"), width = 8, height = 7, units = "in", res = 300)
  pheatmap(
    sample_cor_final,
    main          = "Sample Correlation Heatmap (Final)",
    show_colnames = FALSE,
    show_rownames = FALSE,
    color         = colorRampPalette(c("blue", "white", "red"))(100)
  )
  dev.off()
  rna_msg("Saved: plots/sample_correlation_heatmap.pdf and .png")

  # --- Zero Proportion Check ---
  zero_prop <- sum(rna_top == 0) / length(rna_top)
  rna_msg(sprintf("Zero proportion in final matrix: %.2f%% (expected: very low)",
                  zero_prop * 100))

  invisible(NULL)
}


# ==============================================================================
# FUNCTION: generate_rna_validation_figures
# ==============================================================================

#' Generate publication-ready statistical validation figures
#'
#' Produces:
#'   plots/Research_Validation_Figures.png — 3-panel figure (8x12 in, 300 DPI)
#'     Panel 1: Gene expression density (top 20 samples)
#'     Panel 2: Mean vs SD scatter with GAM smoother (CV analysis)
#'     Panel 3: Inter-sample Pearson correlation histogram
#'
#' Prints statistical summary to console (for thesis/paper writing).
#'
#' @param rna_ml   Numeric matrix — log2 CPM ML version (genes x samples)
#' @param outdir   Character — base output directory
generate_rna_validation_figures <- function(rna_ml, outdir) {

  plots_dir <- file.path(outdir, "plots")
  rna_msg("Generating publication-ready validation figures...")

  # --------------------------------------------------------------------------
  # Figure 1: Global Expression Distribution (top 20 samples)
  # --------------------------------------------------------------------------
  n_samples_plot <- min(20, ncol(rna_ml))
  melted_rna     <- melt(rna_ml[, 1:n_samples_plot])
  colnames(melted_rna) <- c("Gene", "Sample", "Log2CPM")

  fig_density <- ggplot(melted_rna, aes(x = Log2CPM, group = Sample)) +
    geom_density(color = "steelblue", alpha = 0.1) +
    theme_bw() +
    labs(
      title    = "Gene Expression Distribution (Top 20 Samples)",
      subtitle = "Post-TMM Normalization & Batch Correction",
      x        = "Log2(CPM + 2)",
      y        = "Density"
    )

  # --------------------------------------------------------------------------
  # Figure 2: Coefficient of Variation — Mean vs SD (GAM smoother)
  # --------------------------------------------------------------------------
  gene_means <- rowMeans(rna_ml)
  gene_sds   <- apply(rna_ml, 1, sd)
  cv_values  <- (gene_sds / gene_means) * 100

  df_cv <- data.frame(Mean = gene_means, SD = gene_sds, CV = cv_values)

  fig_cv <- ggplot(df_cv, aes(x = Mean, y = SD)) +
    geom_point(alpha = 0.4, color = "#e67e22") +
    geom_smooth(method = "gam", color = "black") +
    theme_minimal() +
    labs(
      title    = "Mean-Standard Deviation Relationship",
      subtitle = "Validating Feature Selection (Top 2000 Genes)",
      x        = "Mean Expression (Log2 CPM)",
      y        = "Standard Deviation"
    )

  # --------------------------------------------------------------------------
  # Figure 3: Inter-Sample Pearson Correlation Histogram
  # --------------------------------------------------------------------------
  sample_cor <- cor(rna_ml, method = "pearson")
  cor_values <- sample_cor[lower.tri(sample_cor)]

  fig_cor_dist <- ggplot(data.frame(Cor = cor_values), aes(x = Cor)) +
    geom_histogram(binwidth = 0.01, fill = "#27ae60", color = "white") +
    theme_classic() +
    labs(
      title    = "Inter-Sample Pearson Correlation",
      subtitle = paste0("Median Correlation: ", round(median(cor_values), 3)),
      x        = "Pearson R",
      y        = "Count"
    )

  # --------------------------------------------------------------------------
  # Combine and Save (8x12 inches, 300 DPI)
  # --------------------------------------------------------------------------
  final_report_plot <- grid.arrange(fig_density, fig_cv, fig_cor_dist, ncol = 1)

  ggsave(
    file.path(plots_dir, "Research_Validation_Figures.png"),
    plot   = final_report_plot,
    width  = 8,
    height = 12,
    dpi    = 300
  )
  rna_msg("Saved: plots/Research_Validation_Figures.png")

  # --------------------------------------------------------------------------
  # Statistical Summary (printed for thesis/paper writing)
  # --------------------------------------------------------------------------
  pca_summary <- summary(prcomp(t(rna_ml)))$importance

  rna_msg("--- Statistical Report ---", level = "DETAIL")
  cat(sprintf("  1. Total High-Quality Samples  : %d\n",    ncol(rna_ml)))
  cat(sprintf("  2. Final Gene Feature Set      : %d genes\n", nrow(rna_ml)))
  cat(sprintf("  3. Median Inter-sample R       : %.4f\n",  median(cor_values)))
  cat(sprintf("  4. Mean Expression Range       : [%.2f to %.2f]\n",
              min(gene_means), max(gene_means)))
  cat(sprintf("  5. PC1 Variance Explained      : %.2f%%\n",
              pca_summary[2, 1] * 100))

  invisible(NULL)
}
