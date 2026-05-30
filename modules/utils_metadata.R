# ==============================================================================
# utils_metadata.R — Shared Utilities for Universal Metadata Layer
# OmicsFlow | Phase 1: Universal Metadata Layer
# ==============================================================================

#' Validate that a data frame conforms to the standardized metadata schema
#' Standardized schema columns:
#'   - sample_id (required)
#'   - patient_id (required)
#'   - sample_class (required)
#'   - batch (required)
#'   - center (optional)
#'
#' @param df Data frame to validate
#' @return Logical TRUE if valid, otherwise throws an error
validate_metadata_schema <- function(df) {
  if (is.null(df) || !is.data.frame(df)) {
    stop("Metadata object must be a valid data frame")
  }
  
  required_cols <- c("sample_id", "patient_id", "sample_class", "batch")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Metadata is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Ensure no NA or empty values in required columns
  for (col in required_cols) {
    vals <- df[[col]]
    if (any(is.na(vals) | trimws(as.character(vals)) == "")) {
      stop("Metadata column contains missing or empty values: ", col)
    }
  }
  
  return(TRUE)
}

#' Load and validate a standardized metadata schema
#'
#' @param metadata_file Path to metadata CSV file
#' @return A validated data frame, or NULL if metadata_file is NULL or absent
load_metadata <- function(metadata_file) {
  if (is.null(metadata_file) || metadata_file == "") {
    return(NULL)
  }
  if (!file.exists(metadata_file)) {
    stop("Metadata file not found: ", metadata_file)
  }
  
  df <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Validate schema
  validate_metadata_schema(df)
  
  # Standardize types and trim whitespace
  required_cols <- c("sample_id", "patient_id", "sample_class", "batch")
  for (col in required_cols) {
    df[[col]] <- trimws(as.character(df[[col]]))
  }
  if ("center" %in% colnames(df)) {
    df[["center"]] <- trimws(as.character(df[["center"]]))
  }
  
  # Deduplicate metadata rows to be clean
  cols_to_keep <- intersect(colnames(df), c(required_cols, "center"))
  df <- unique(df[, cols_to_keep, drop = FALSE])
  
  return(df)
}

#' Map sample IDs to patient IDs using metadata or TCGA barcode parsing
#'
#' @param sample_ids Vector of sample IDs/barcodes
#' @param metadata Data frame loaded via load_metadata(), or NULL
#' @return Vector of mapped patient IDs
get_patient_id <- function(sample_ids, metadata = NULL) {
  if (is.null(metadata)) {
    # Fallback to TCGA: 12-character prefix
    return(substr(as.character(sample_ids), 1, 12))
  } else {
    matched_idx <- match(as.character(sample_ids), metadata$sample_id)
    if (any(is.na(matched_idx))) {
      missing_samples <- unique(sample_ids[is.na(matched_idx)])
      stop("The following sample IDs in raw data are missing from metadata: ", 
           paste(head(missing_samples, 5), collapse = ", "))
    }
    return(metadata$patient_id[matched_idx])
  }
}

#' Map sample IDs to sample classes/types
#'
#' @param sample_ids Vector of sample IDs/barcodes
#' @param metadata Data frame loaded via load_metadata(), or NULL
#' @return Vector of mapped sample classes/types
get_sample_class <- function(sample_ids, metadata = NULL) {
  if (is.null(metadata)) {
    # Fallback to TCGA: positions 14-15
    return(substr(as.character(sample_ids), 14, 15))
  } else {
    matched_idx <- match(as.character(sample_ids), metadata$sample_id)
    if (any(is.na(matched_idx))) {
      missing_samples <- unique(sample_ids[is.na(matched_idx)])
      stop("The following sample IDs in raw data are missing from metadata: ", 
           paste(head(missing_samples, 5), collapse = ", "))
    }
    return(metadata$sample_class[matched_idx])
  }
}

#' Map sample IDs to batches
#'
#' @param sample_ids Vector of sample IDs/barcodes
#' @param metadata Data frame loaded via load_metadata(), or NULL
#' @param omics_type Character, one of "rna", "meth"
#' @return Vector of mapped batches
get_batch <- function(sample_ids, metadata = NULL, omics_type = "rna") {
  if (is.null(metadata)) {
    # Fallback to TCGA
    if (omics_type == "rna") {
      # RNA Plate ID: positions 22-25
      res <- substr(as.character(sample_ids), 22, 25)
      res[is.na(res) | res == ""] <- "batch1"
      return(res)
    } else if (omics_type == "meth") {
      # Methylation Plate ID: position 6 after splitting by "-"
      res <- sapply(strsplit(as.character(sample_ids), "-"),
                    function(x) if (length(x) >= 6) x[6] else NA_character_)
      res[is.na(res) | res == ""] <- "batch1"
      return(res)
    } else {
      return(rep("batch1", length(sample_ids)))
    }
  } else {
    matched_idx <- match(as.character(sample_ids), metadata$sample_id)
    if (any(is.na(matched_idx))) {
      missing_samples <- unique(sample_ids[is.na(matched_idx)])
      stop("The following sample IDs in raw data are missing from metadata: ", 
           paste(head(missing_samples, 5), collapse = ", "))
    }
    return(metadata$batch[matched_idx])
  }
}

#' Map sample IDs to centers (optional biological covariate)
#'
#' @param sample_ids Vector of sample IDs/barcodes
#' @param metadata Data frame loaded via load_metadata(), or NULL
#' @return Vector of centers, or NULL if not present in metadata
get_center <- function(sample_ids, metadata = NULL) {
  if (is.null(metadata) || !"center" %in% colnames(metadata)) {
    return(NULL)
  } else {
    matched_idx <- match(as.character(sample_ids), metadata$sample_id)
    if (any(is.na(matched_idx))) {
      missing_samples <- unique(sample_ids[is.na(matched_idx)])
      stop("The following sample IDs in raw data are missing from metadata: ", 
           paste(head(missing_samples, 5), collapse = ", "))
    }
    return(metadata$center[matched_idx])
  }
}
