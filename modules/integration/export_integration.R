# ==============================================================================
# export_integration.R — Output Generation for Integration Module
# OmicsFlow — Integration Module
# ==============================================================================

source("modules/integration/utils_integration.R")

export_int_results <- function(mofa, active_factors, outdir) {
  int_step("Export", "Saving Processed Data")
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # Saved as named list
  path_mofa <- file.path(outdir, "mofa_model.rds")
  saveRDS(
    list(model = mofa, active_factors = active_factors),
    path_mofa
  )
  
  int_msg(sprintf("Saved: mofa_model.rds"))
  int_msg(sprintf("  $model          - MOFA2 trained object"), level = "DETAILS")
  int_msg(sprintf("  $active_factors - Factor indices (mean R\u00b2 > 1%%)"), level = "DETAILS")
  
  invisible(list(
    model = path_mofa
  ))
}
