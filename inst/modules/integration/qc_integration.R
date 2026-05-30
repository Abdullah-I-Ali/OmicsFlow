# ==============================================================================
# qc_integration.R — Quality Control & Visualization for Integration
# OmicsFlow — Integration Module
# ==============================================================================

source("modules/integration/utils_integration.R")

generate_int_qc_plots <- function(mofa, outdir) {
  int_step("QC", "Generating QC Visualizations")
  
  plot_dir <- file.path(outdir, "plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  # 1. ELBO Convergence Plot
  if (!is.null(mofa@training_stats$elbo)) {
    elbo <- mofa@training_stats$elbo
    pdf(file.path(plot_dir, "mofa_elbo_convergence.pdf"), width = 6, height = 4)
    plot(elbo, type = "l", col = "steelblue", lwd = 2,
         xlab = "Iteration", ylab = "ELBO",
         main = "MOFA Training Convergence")
    dev.off()
    png(file.path(plot_dir, "mofa_elbo_convergence.png"), width = 6, height = 4, units = "in", res = 300)
    plot(elbo, type = "l", col = "steelblue", lwd = 2,
         xlab = "Iteration", ylab = "ELBO",
         main = "MOFA Training Convergence")
    dev.off()
    int_msg("Saved: plots/mofa_elbo_convergence.pdf and .png")
  }
  
  # 2. Total Variance Explained
  p_var_total <- plot_variance_explained(mofa, plot_total = TRUE)[[2]] +
    theme_bw(base_size = 12) +
    labs(title = "Total Variance Explained", 
         subtitle = "Contribution of all factors per omics layer") +
    theme(plot.title = element_text(face = "bold"))
  
  # 3. Variance per Factor
  p_var_factor <- plot_variance_explained(mofa, x="view", y="factor") +
    theme_bw(base_size = 12) +
    labs(title = "Variance Explained per Factor", 
         subtitle = "Identifying layer-specific vs. shared biological drivers") +
    theme(plot.title = element_text(face = "bold"))
  
  fig_variance <- grid.arrange(p_var_total, p_var_factor, ncol = 2, widths = c(1, 2.5))
  ggsave(file.path(plot_dir, "MOFA_Fig1_Variance_Landscape.png"), 
         plot = fig_variance, width = 14, height = 6, dpi = 300)
  int_msg("Saved: plots/MOFA_Fig1_Variance_Landscape.png")
  
  # 4. Factor Correlation Matrix
  pdf(file.path(plot_dir, "MOFA_Fig2_Factor_Correlation.pdf"), width = 6, height = 5)
  plot_factor_cor(mofa)
  dev.off()
  png(file.path(plot_dir, "MOFA_Fig2_Factor_Correlation.png"), width = 6, height = 5, units = "in", res = 300)
  plot_factor_cor(mofa)
  dev.off()
  int_msg("Saved: plots/MOFA_Fig2_Factor_Correlation.pdf and .png")
  
  # 5. Top Molecular Drivers for Factor 1
  # Find which view is most strongly driven by Factor 1
  var_exp <- get_variance_explained(mofa)
  f1_var <- var_exp$r2_per_factor[[1]]["Factor1", ]
  best_view_f1 <- names(which.max(f1_var))
  
  p_weights <- plot_weights(mofa,
                            view = best_view_f1,
                            factor = 1,
                            nfeatures = 15,
                            scale = TRUE,
                            abs = FALSE) +
    theme_minimal(base_size = 12) +
    labs(title = paste("Top 15 Molecular Drivers in Factor 1\n(Layer:", best_view_f1, ")"),
         subtitle = "Features with highest absolute weights",
         x = "Molecular Feature", y = "Scaled MOFA Weight") +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(plot_dir, "MOFA_Fig3_Factor1_Drivers.png"), 
         plot = p_weights, width = 7, height = 6, dpi = 300)
  int_msg("Saved: plots/MOFA_Fig3_Factor1_Drivers.png")
}
