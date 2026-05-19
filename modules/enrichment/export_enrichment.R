# ==============================================================================
# export_enrichment.R — Output Generation for Enrichment Module
# OmicsFlow — Enrichment Module
# ==============================================================================

source("modules/enrichment/utils_enrichment.R")

export_pathway_results <- function(go_bp, go_mf, go_cc, kegg, outdir) {
  path_step("Export", "Saving Enrichment Results")
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  saveRDS(
    list(go_bp = go_bp, go_mf = go_mf, go_cc = go_cc, kegg = kegg),
    file.path(outdir, "enrichment_results.rds")
  )
  path_msg("Saved: enrichment_results.rds")
  
  if (nrow(go_bp) > 0) {
    write.csv(as.data.frame(go_bp), file.path(outdir, "go_bp_results.csv"), row.names = FALSE)
    path_msg("Saved: go_bp_results.csv")
  }
  
  if (nrow(go_mf) > 0) {
    write.csv(as.data.frame(go_mf), file.path(outdir, "go_mf_results.csv"), row.names = FALSE)
    path_msg("Saved: go_mf_results.csv")
  }
  
  if (nrow(go_cc) > 0) {
    write.csv(as.data.frame(go_cc), file.path(outdir, "go_cc_results.csv"), row.names = FALSE)
    path_msg("Saved: go_cc_results.csv")
  }
  
  if (nrow(kegg) > 0) {
    write.csv(as.data.frame(kegg), file.path(outdir, "kegg_results.csv"), row.names = FALSE)
    path_msg("Saved: kegg_results.csv")
  }
  
  invisible(list(
    results_rds = file.path(outdir, "enrichment_results.rds")
  ))
}
