# ==============================================================================
# export_snv.R — Output Generation for SNV Module
# OmicsFlow — SNV Module
# ==============================================================================

source("modules/snv/utils_snv.R")

export_snv_results <- function(snv_matrix, maf_filtered, outdir, sample_info = NULL) {
  snv_step("Export", "Saving Processed Data")
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # Standardized integration matrix filename
  path_mofa <- file.path(outdir, "snv_processed_matrix.rds")
  saveRDS(snv_matrix, path_mofa)
  snv_msg(sprintf("Saved: snv_processed_matrix.rds  (%d genes x %d samples)",
                  nrow(snv_matrix), ncol(snv_matrix)))
  
  # Filtered MAF object for downstream maftools visualization
  path_maf <- file.path(outdir, "snv_maf_clean.rds")
  saveRDS(maf_filtered, path_maf)
  snv_msg("Saved: snv_maf_clean.rds")
  
  # Sample metadata
  if (is.null(sample_info)) {
    sample_info <- data.frame(
      patient_id = colnames(snv_matrix),
      stringsAsFactors = FALSE
    )
  }
  path_meta_csv <- file.path(outdir, "sample_metadata.csv")
  path_meta_rds <- file.path(outdir, "sample_metadata.rds")
  write.csv(sample_info, path_meta_csv, row.names = FALSE)
  saveRDS(sample_info, path_meta_rds)
  snv_msg("Saved: sample_metadata.csv and .rds")
  
  invisible(list(
    processed_matrix = path_mofa,
    maf_clean        = path_maf,
    sample_metadata  = path_meta_csv
  ))
}
