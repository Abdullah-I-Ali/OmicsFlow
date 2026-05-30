# ==============================================================================
# R/build_metadata_table.R — Construct sample_metadata.csv Data Frame
# ==============================================================================

#' Construct a sample metadata data frame from detected sample IDs
#'
#' Unions all non-NULL sample ID vectors from the detected modalities,
#' resolves patient IDs (via \code{infer_patient_ids} if not pre-computed),
#' and returns a standardized metadata data frame with conservative defaults.
#'
#' @param sample_ids_list Named list of sample ID character vectors
#'   (from \code{detect_sample_ids})
#' @param patient_ids_df Data frame containing \code{sample_id},
#'   \code{patient_id}, \code{method} mapping (optional; inferred if NULL)
#' @return A data frame with columns: \code{sample_id}, \code{patient_id},
#'   \code{sample_class}, \code{batch}, \code{center}
#' @keywords internal
build_metadata_table <- function(sample_ids_list, patient_ids_df = NULL) {
  # Union all non-NULL sample ID vectors
  all_samples <- unique(unlist(sample_ids_list))
  if (is.null(all_samples) || length(all_samples) == 0) {
    return(data.frame(
      sample_id = character(0),
      patient_id = character(0),
      sample_class = character(0),
      batch = character(0),
      center = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Resolve patient IDs
  if (is.null(patient_ids_df)) {
    # Call infer_patient_ids on all sample IDs
    patient_ids_df <- infer_patient_ids(all_samples)
  }
  
  # Ensure we have mapping for all samples
  matched_idx <- match(all_samples, patient_ids_df$sample_id)
  patient_ids <- patient_ids_df$patient_id[matched_idx]
  # Fallback for any unmatched (should not happen if mapped properly)
  patient_ids[is.na(patient_ids)] <- all_samples[is.na(patient_ids)]
  
  # Apply conservative defaults
  df_meta <- data.frame(
    sample_id = all_samples,
    patient_id = patient_ids,
    sample_class = "Tumor",
    batch = "B1",
    center = "CenterA",
    stringsAsFactors = FALSE
  )
  
  # TCGA Barcode Inference
  # A standard TCGA barcode is 28 characters (e.g. TCGA-CC-A7IJ-01A-41D-A915-36)
  is_tcga <- startsWith(as.character(df_meta$sample_id), "TCGA-") & nchar(as.character(df_meta$sample_id)) >= 28
  
  if (any(is_tcga)) {
    tcga_samples <- as.character(df_meta$sample_id[is_tcga])
    
    # Extract tissue type code (positions 14-15)
    tissue_codes <- substr(tcga_samples, 14, 15)
    class_map <- ifelse(tissue_codes == "01", "Primary Tumor",
                 ifelse(tissue_codes == "11", "Solid Tissue Normal",
                 ifelse(tissue_codes == "02", "Recurrent Tumor",
                 ifelse(tissue_codes == "06", "Metastatic", "Tumor"))))
    df_meta$sample_class[is_tcga] <- class_map
    
    # Extract Plate/Batch (positions 22-25)
    plate_codes <- substr(tcga_samples, 22, 25)
    df_meta$batch[is_tcga] <- paste0("Plate_", plate_codes)
    
    # Extract Center (positions 27-28)
    center_codes <- substr(tcga_samples, 27, 28)
    df_meta$center[is_tcga] <- paste0("Center_", center_codes)
  }
  
  return(df_meta)
}
