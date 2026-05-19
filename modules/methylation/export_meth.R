# ==============================================================================
# export_meth.R — Output Export for the Methylation Preprocessing Module
# OmicsFlow | Phase 1: Methylation Module
# Scientific source: data/full_scripts.R (lines 1079-1095 & 1150-1175)
# ==============================================================================

#' Export all Methylation module outputs
#'
#' Writes:
#'   - methylation_processed_matrix.rds     (M-values, MOFA+/DIABLO, renamed to standardized output)
#'   - methylation_beta_top_variable.rds    (Beta values, MOFA+/DIABLO)
#'   - methylation_m_FULL_ML_Ready.rds      (M-values, full matrix)
#'   - methylation_beta_FULL_ML_Ready.rds   (Beta values, full matrix)
#'   - sample_metadata.csv                  (Patient ID + batch)
#'   - sample_metadata.rds
#'   - qc_metrics.json
#'   - probe_filter_log.csv
#'   - final_sample_ids.txt
export_meth_results <- function(meth_m_top, met_beta_top, meth_m_full, met_beta_full, 
                                sample_info, qc_metrics, probe_log, outdir) {

  meth_validate_barcodes(meth_m_top,  "meth_m_top (MOFA)")
  meth_validate_barcodes(meth_m_full, "meth_m_full (ML)")

  # --- MOFA / DIABLO matrices ---
  path_mofa_m <- file.path(outdir, "methylation_processed_matrix.rds")
  saveRDS(meth_m_top, path_mofa_m)
  meth_msg(sprintf("Saved: methylation_processed_matrix.rds (%d probes x %d samples, M-values)", 
                   nrow(meth_m_top), ncol(meth_m_top)))

  path_mofa_beta <- file.path(outdir, "methylation_beta_top_variable.rds")
  saveRDS(met_beta_top, path_mofa_beta)
  meth_msg(sprintf("Saved: methylation_beta_top_variable.rds (%d probes x %d samples, Beta)", 
                   nrow(met_beta_top), ncol(met_beta_top)))

  # --- ML matrices ---
  path_ml_m <- file.path(outdir, "methylation_m_FULL_ML_Ready.rds")
  saveRDS(meth_m_full, path_ml_m, compress = "xz")
  meth_msg(sprintf("Saved: methylation_m_FULL_ML_Ready.rds (%d probes x %d samples, M-values)", 
                   nrow(meth_m_full), ncol(meth_m_full)))

  path_ml_beta <- file.path(outdir, "methylation_beta_FULL_ML_Ready.rds")
  saveRDS(met_beta_full, path_ml_beta, compress = "xz")
  meth_msg(sprintf("Saved: methylation_beta_FULL_ML_Ready.rds (%d probes x %d samples, Beta)", 
                   nrow(met_beta_full), ncol(met_beta_full)))

  # --- Metadata ---
  path_meta_csv <- file.path(outdir, "sample_metadata.csv")
  path_meta_rds <- file.path(outdir, "sample_metadata.rds")
  write.csv(sample_info, path_meta_csv, row.names = FALSE)
  saveRDS(sample_info, path_meta_rds)
  meth_msg(sprintf("Saved: sample_metadata.csv (%d samples)", nrow(sample_info)))

  # --- QC & Audit ---
  meth_save_qc_json(qc_metrics, outdir)

  probe_df <- data.frame(Step = names(probe_log), Probes = unlist(probe_log))
  write.csv(probe_df, file.path(outdir, "probe_filter_log.csv"), row.names = FALSE)
  meth_msg("Saved: probe_filter_log.csv")

  writeLines(colnames(meth_m_top), file.path(outdir, "final_sample_ids.txt"))
  meth_msg("Saved: final_sample_ids.txt")

  invisible(list(
    processed_matrix = path_mofa_m,
    ml_m_matrix      = path_ml_m,
    qc_metrics_json  = file.path(outdir, "qc_metrics.json")
  ))
}
