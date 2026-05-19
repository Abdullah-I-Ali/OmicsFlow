# ==============================================================================
# qc_ml.R — Quality Control & Visualization for ML
# OmicsFlow — Machine Learning Module
# ==============================================================================

source("modules/ml/utils_ml.R")

generate_ml_qc_plots <- function(results_summary, top_features, test_data, xgb_risk, os_time_col = "os_time", os_event_col = "os_event", outdir) {
  ml_step("QC", "Generating ML Visualizations")
  
  plot_dir <- file.path(outdir, "plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  # 1. Model Performance Comparison (C-Index)
  results_summary$model <- factor(results_summary$model, 
                                  levels = results_summary$model[order(results_summary$c_index)])
  
  fig_cindex <- ggplot(results_summary, aes(x = model, y = c_index, fill = model)) +
    geom_bar(stat = "identity", width = 0.6, color = "black", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.3f", c_index)), hjust = -0.2, size = 5, fontface = "bold") +
    coord_flip() +
    scale_fill_manual(values = c("LASSO Cox" = "#bdc3c7", 
                                 "Random Forest" = "#3498db", 
                                 "XGBoost Cox" = "#e74c3c")) +
    theme_classic(base_size = 14) +
    scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, by = 0.2)) +
    theme(legend.position = "none") +
    labs(title = "Prognostic Model Performance",
         subtitle = "Concordance Index (C-index) on hold-out test set",
         x = "Survival Algorithm", y = "C-index (Higher is better)")
  
  # 2. Top Predictive Features (Biomarkers)
  top15_genes <- top_features %>% 
    arrange(desc(train_variance)) %>% 
    head(15)
  
  top15_genes$feature <- factor(top15_genes$feature, levels = rev(top15_genes$feature))
  
  fig_features <- ggplot(top15_genes, aes(x = train_variance, y = feature)) +
    geom_segment(aes(x = 0, xend = train_variance, y = feature, yend = feature), color = "gray50", linewidth = 1) +
    geom_point(color = "#27ae60", size = 4) +
    theme_minimal(base_size = 14) +
    labs(title = "Top 15 Prognostic Biomarkers",
         subtitle = "Genes with highest predictive variance derived from MOFA latent factors",
         x = "Feature Importance (Training Variance)", y = "Gene Symbol") +
    theme(panel.grid.major.y = element_blank())
  
  # Combine
  final_ml_plot <- grid.arrange(fig_cindex, fig_features, ncol = 2, widths = c(1, 1.2))
  ggsave(file.path(plot_dir, "ML_Research_Validation_Figures.png"), 
         plot = final_ml_plot, width = 14, height = 6, dpi = 300)
  ml_msg("Saved: plots/ML_Research_Validation_Figures.png")
  
  # 3. Kaplan-Meier by Risk Group (XGBoost risk score)
  test_data$risk_group <- factor(
    ifelse(xgb_risk < median(xgb_risk), "Low Risk", "High Risk"),
    levels = c("Low Risk", "High Risk")
  )
  
  surv_fit <- survfit(Surv(os_time, os_event) ~ risk_group, data = test_data)
  
  p_km <- ggsurvplot(
    surv_fit,
    data       = test_data,
    pval       = TRUE,
    risk.table = TRUE,
    palette    = c("steelblue", "firebrick"),
    title      = "Kaplan-Meier: High vs Low Risk (MOFA-Hybrid -> XGBoost Cox)",
    xlab       = "Time (days)",
    ylab       = "Overall Survival Probability"
  )
  
  pdf(file.path(plot_dir, "kaplan_meier_risk_groups.pdf"), width = 8, height = 6)
  print(p_km, newpage = FALSE)
  dev.off()
  png(file.path(plot_dir, "kaplan_meier_risk_groups.png"), width = 8, height = 6, units = "in", res = 300)
  print(p_km, newpage = FALSE)
  dev.off()
  ml_msg("Saved: plots/kaplan_meier_risk_groups.pdf and .png")
}
