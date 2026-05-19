# ==============================================================================
# export_cnv.R — Output Generation for CNV Module
# OmicsFlow — CNV Module
# ==============================================================================

source("modules/cnv/utils_cnv.R")

export_cnv_results <- function(cnv_matrix, gene_var_df, outdir) {
  cnv_step("Export", "Saving Processed Data")
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # Standardized integration matrix filename
  path_mofa <- file.path(outdir, "cnv_processed_matrix.rds")
  saveRDS(cnv_matrix, path_mofa)
  cnv_msg(sprintf("Saved: cnv_processed_matrix.rds  (%d genes x %d samples)",
                  nrow(cnv_matrix), ncol(cnv_matrix)))
  
  # Gene variance information
  path_var <- file.path(outdir, "cnv_gene_variance_info.rds")
  saveRDS(gene_var_df, path_var)
  cnv_msg("Saved: cnv_gene_variance_info.rds")
  
  # Sample metadata
  # CNV data doesn't have an explicit batch vector derived within this module,
  # but we provide patient_id for uniform structure.
  sample_info <- data.frame(
    patient_id = colnames(cnv_matrix),
    stringsAsFactors = FALSE
  )
  path_meta_csv <- file.path(outdir, "sample_metadata.csv")
  path_meta_rds <- file.path(outdir, "sample_metadata.rds")
  write.csv(sample_info, path_meta_csv, row.names = FALSE)
  saveRDS(sample_info, path_meta_rds)
  cnv_msg("Saved: sample_metadata.csv and .rds")
  
  invisible(list(
    processed_matrix = path_mofa,
    sample_metadata  = path_meta_csv
  ))
}
