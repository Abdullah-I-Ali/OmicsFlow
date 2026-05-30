# ==============================================================================
# export_rna.R — Output Export for the RNA Preprocessing Module
# OmicsFlow | Phase 1: RNA Module
# Scientific source: data/full_scripts.R (lines 492-537)
# ==============================================================================
# PURPOSE:
#   Single exported function sourced by preprocess_rna.R.
#   Writes all standardized RNA module outputs:
#
#     rna_processed_matrix.rds — Z-scored log2 CPM (MOFA+/DIABLO)
#     rna_ml.rds              — Raw log2 CPM (Machine Learning)
#     sample_metadata.csv     — Patient ID + batch (CSV)
#     sample_metadata.rds     — Patient ID + batch (RDS)
#     qc_metrics.json         — Full QC statistics JSON
#
# NOTE:
#   No scientific logic lives here. All preprocessing has already occurred
#   in preprocess_rna.R. This file only handles serialization.
# ==============================================================================


#' Export all RNA module outputs to the standardized output directory
#'
#' Writes five files to outdir/:
#'   - rna_processed_matrix.rds (Z-scored, for MOFA+/DIABLO)
#'   - rna_ml.rds             (log2 CPM, for Machine Learning)
#'   - sample_metadata.csv
#'   - sample_metadata.rds
#'   - qc_metrics.json
#'
#' @param rna_scaled   Numeric matrix — Z-scored log2 CPM (genes x samples)
#' @param rna_ml       Numeric matrix — Raw log2 CPM (genes x samples)
#' @param sample_info  data.frame — columns: patient_id, batch
#' @param qc_metrics   Named list — accumulated QC statistics
#' @param outdir       Character — base output directory
export_rna_results <- function(rna_scaled, rna_ml, sample_info,
                                qc_metrics, outdir, metadata_supplied = FALSE) {

  # --- Validate barcodes before saving ---
  rna_validate_barcodes(rna_scaled, "rna_scaled (MOFA)", metadata_supplied)
  rna_validate_barcodes(rna_ml,     "rna_ml (ML)", metadata_supplied)

  # --------------------------------------------------------------------------
  # rna_processed_matrix.rds — Z-scored matrix (primary integration output)
  # Named "rna_processed_matrix" for standardized inter-module compatibility.
  # --------------------------------------------------------------------------
  path_mofa <- file.path(outdir, "rna_processed_matrix.rds")
  saveRDS(rna_scaled, path_mofa)
  rna_msg(sprintf("Saved: rna_processed_matrix.rds  (%d genes x %d samples, Z-scored)",
                  nrow(rna_scaled), ncol(rna_scaled)))

  # --------------------------------------------------------------------------
  # rna_ml.rds — log2 CPM matrix (Machine Learning)
  # Separate file so ML pipelines can access without re-running preprocessing.
  # --------------------------------------------------------------------------
  path_ml <- file.path(outdir, "rna_ml.rds")
  saveRDS(rna_ml, path_ml)
  rna_msg(sprintf("Saved: rna_ml.rds  (%d genes x %d samples, log2 CPM)",
                  nrow(rna_ml), ncol(rna_ml)))

  # --------------------------------------------------------------------------
  # sample_metadata.csv + sample_metadata.rds
  # Both formats required: CSV for human inspection, RDS for downstream R use.
  # --------------------------------------------------------------------------
  path_meta_csv <- file.path(outdir, "sample_metadata.csv")
  path_meta_rds <- file.path(outdir, "sample_metadata.rds")

  write.csv(sample_info, path_meta_csv, row.names = FALSE)
  saveRDS(sample_info,   path_meta_rds)

  rna_msg(sprintf("Saved: sample_metadata.csv  (%d samples)", nrow(sample_info)))
  rna_msg(sprintf("Saved: sample_metadata.rds  (%d samples)", nrow(sample_info)))

  # --------------------------------------------------------------------------
  # qc_metrics.json — standardized QC output required by all modules
  # --------------------------------------------------------------------------
  rna_save_qc_json(qc_metrics, outdir)

  # --------------------------------------------------------------------------
  # Print file manifest
  # --------------------------------------------------------------------------
  rna_msg("Output manifest:")
  rna_msg(sprintf("  %s", file.path(outdir, "rna_processed_matrix.rds")),
          level = "DETAIL")
  rna_msg(sprintf("  %s", file.path(outdir, "rna_ml.rds")),
          level = "DETAIL")
  rna_msg(sprintf("  %s", file.path(outdir, "sample_metadata.csv")),
          level = "DETAIL")
  rna_msg(sprintf("  %s", file.path(outdir, "sample_metadata.rds")),
          level = "DETAIL")
  rna_msg(sprintf("  %s", file.path(outdir, "qc_metrics.json")),
          level = "DETAIL")
  rna_msg(sprintf("  %s", file.path(outdir, "plots/")),
          level = "DETAIL")

  invisible(list(
    processed_matrix = path_mofa,
    rna_ml           = path_ml,
    sample_metadata_csv = path_meta_csv,
    sample_metadata_rds = path_meta_rds,
    qc_metrics_json  = file.path(outdir, "qc_metrics.json")
  ))
}
