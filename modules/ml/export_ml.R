# ==============================================================================
# export_ml.R — Output Generation for ML Module
# OmicsFlow — Machine Learning Module
# ==============================================================================

source("modules/ml/utils_ml.R")

export_ml_results <- function(results_summary, top_features_export, rf_model, xgb_model, cv_lasso, mofa_top_df, selected_genes, lasso_genes, rna, outdir) {
  ml_step("Export", "Saving Results")
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  write.csv(results_summary, file.path(outdir, "ml_results_summary.csv"), row.names = FALSE)
  write.csv(top_features_export, file.path(outdir, "ml_top_features.csv"), row.names = FALSE)
  ml_msg("Saved: ml_results_summary.csv")
  ml_msg("Saved: ml_top_features.csv")
  
  saveRDS(rf_model,  file.path(outdir, "rf_survival_model.rds"))
  saveRDS(xgb_model, file.path(outdir, "xgb_cox_model.rds"))
  saveRDS(cv_lasso,  file.path(outdir, "lasso_cox_model.rds"))
  ml_msg("Saved: Trained Models (RF, XGBoost, LASSO)")
  
  saveRDS(mofa_top_df, file.path(outdir, "mofa_top_genes.rds"))
  saveRDS(selected_genes, file.path(outdir, "rf_top_genes.rds"))
  saveRDS(lasso_genes, file.path(outdir, "lasso_selected_genes.rds"))
  saveRDS(rna, file.path(outdir, "rna_for_pathway.rds"))
  ml_msg("Saved: Downstream data components")
  
  invisible(list(
    summary = file.path(outdir, "ml_results_summary.csv"),
    rf = file.path(outdir, "rf_survival_model.rds"),
    xgb = file.path(outdir, "xgb_cox_model.rds")
  ))
}
