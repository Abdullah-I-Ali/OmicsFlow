# ==============================================================================
# R/infer_patient_ids.R — Conservative Patient ID Inference
# ==============================================================================

#' Conservative regex-based patient ID inference from sample IDs
#'
#' Attempts to extract patient identifiers from sample barcodes using a
#' tiered strategy: TCGA barcode parsing, common suffix stripping, and
#' identity fallback.
#'
#' @param sample_ids Character vector of sample IDs
#' @return A data frame with columns: \code{sample_id}, \code{patient_id},
#'   \code{method}
#' @keywords internal
infer_patient_ids <- function(sample_ids) {
  if (is.null(sample_ids) || length(sample_ids) == 0) {
    return(data.frame(
      sample_id = character(0),
      patient_id = character(0),
      method = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  sample_ids <- unique(as.character(sample_ids))
  patient_ids <- rep(NA_character_, length(sample_ids))
  methods <- rep(NA_character_, length(sample_ids))
  
  # TCGA pattern: e.g. TCGA-02-0001-01A-...
  # We match the first 12 characters: TCGA-XX-XXXX
  tcga_pattern <- "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}"
  
  # Suffix patterns: _T, _N, _01A, -01 at the end of the string
  suffix_pattern <- "(_[TN]$|_[0-9]{2}[A-Z]$|-[0-9]{2}$)"
  
  for (i in seq_along(sample_ids)) {
    sid <- sample_ids[i]
    if (grepl(tcga_pattern, sid, ignore.case = TRUE)) {
      # Extract first 12 chars
      patient_ids[i] <- substr(sid, 1, 12)
      methods[i] <- "tcga"
    } else if (grepl(suffix_pattern, sid, ignore.case = TRUE)) {
      # Strip suffix
      patient_ids[i] <- sub(suffix_pattern, "", sid, ignore.case = TRUE)
      methods[i] <- "suffix_strip"
    } else {
      # Fallback
      patient_ids[i] <- sid
      methods[i] <- "identity"
    }
  }
  
  data.frame(
    sample_id = sample_ids,
    patient_id = patient_ids,
    method = methods,
    stringsAsFactors = FALSE
  )
}
