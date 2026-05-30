# ==============================================================================
# R/validate_inputs.R — Pre-flight Validation of Pipeline Inputs
# ==============================================================================

#' Pre-flight validation of OmicsFlow pipeline inputs
#'
#' Performs comprehensive validation of user-prepared pipeline inputs including
#' file existence checks, modality counting (minimum 2 required), metadata
#' schema validation, sample ID alignment, and clinical column mapping
#' verification.
#'
#' @param rna Path to RNA-seq matrix file (.rds, .csv, .tsv, .txt), or NULL
#' @param meth Path to DNA methylation matrix file (.rds, .csv, .tsv, .txt), or NULL
#' @param cnv Path to CNV data file (.rds, .csv, .tsv, .txt), or NULL
#' @param snv Path to SNV data file (.rds, .csv, .tsv, .txt), or NULL
#' @param metadata Path to sample metadata CSV file (required for execution)
#' @param clinical Path to clinical TSV/CSV file, or NULL
#' @param clinical_map Path to clinical mapping JSON file, or NULL
#' @param cnv_cache Path to Ensembl gene coordinates cache for CNV preprocessing.
#' @param meth_cross_react Path to cross-reactive probes list for Methylation preprocessing.
#'
#' @return An S3 object of class \code{omicsflow_validation} containing:
#'   \describe{
#'     \item{status}{Character, either "success" or "failed"}
#'     \item{valid}{Logical indicating whether all checks passed}
#'     \item{detected_modalities}{Character vector of executed/detected modules}
#'     \item{skipped_modalities}{Character vector of skipped modules}
#'     \item{messages}{Character vector of validation log messages}
#'     \item{warnings}{Character vector of validation warnings}
#'     \item{errors}{Character vector of validation errors}
#'     \item{runtime}{Numeric execution time in seconds}
#'   }
#'
#' @examples
#' \dontrun{
#' result <- validate_inputs(
#'   rna = "data/rna.rds",
#'   meth = "data/meth.rds",
#'   metadata = "data/sample_metadata.csv"
#' )
#' if (result$valid) message("Ready to run!")
#' }
#'
#' @export
validate_inputs <- function(rna = NULL, meth = NULL, cnv = NULL, snv = NULL,
                            metadata = NULL, clinical = NULL, clinical_map = NULL,
                            cnv_cache = NULL, meth_cross_react = NULL) {
  # validate_metadata_schema() is now a package-internal function in

  # R/utils_metadata.R — no source() call needed.
  
  start_time <- Sys.time()
  valid <- TRUE
  messages <- character(0)
  warnings_list <- character(0)
  errors_list <- character(0)
  
  # Helper to log validation checks
  log_ok <- function(msg) {
    messages <<- c(messages, sprintf("[OK] %s", msg))
    msg_ok(msg)
  }
  log_warn <- function(msg) {
    warnings_list <<- c(warnings_list, msg)
    messages <<- c(messages, sprintf("[WARNING] %s", msg))
    msg_warn(msg)
  }
  log_fail <- function(msg) {
    errors_list <<- c(errors_list, msg)
    messages <<- c(messages, sprintf("[FAIL] %s", msg))
    msg_fail(msg)
    valid <<- FALSE
  }
  
  # 1. Modality and file checks
  provided_modalities <- 0
  detected_mod_names <- character(0)
  skipped_mod_names <- character(0)
  
  check_modality <- function(label, path, name) {
    if (!is.null(path)) {
      if (!file.exists(path)) {
        log_fail(sprintf("%s file does not exist: %s", label, path))
      } else {
        log_ok(sprintf("%s detected: %s", name, path))
        provided_modalities <<- provided_modalities + 1
        detected_mod_names <<- c(detected_mod_names, name)
      }
    } else {
      log_warn(sprintf("%s not provided — skipping %s module", name, name))
      skipped_mod_names <<- c(skipped_mod_names, name)
    }
  }
  
  check_exists <- function(label, path) {
    if (!is.null(path)) {
      if (!file.exists(path)) {
        log_fail(sprintf("%s file does not exist: %s", label, path))
      } else {
        log_ok(sprintf("%s file exists: %s", label, path))
      }
    }
  }
  
  check_modality("RNA-seq matrix", rna, "RNA")
  check_modality("Methylation matrix", meth, "Methylation")
  
  if (!is.null(meth)) {
    if (is.null(meth_cross_react)) {
      try({
        candidate_cross <- omicsflow_path("configs", "cross_reactive_probes.csv")
        if (file.exists(candidate_cross)) {
          meth_cross_react <- candidate_cross
        }
      }, silent = TRUE)
    }
    if (is.null(meth_cross_react) || !file.exists(meth_cross_react)) {
      log_warning("Methylation cross-reactive probes list not found locally. Preprocessing will attempt to download it from GitHub.")
    } else {
      log_ok(sprintf("Methylation cross-reactive probes resolved: %s", meth_cross_react))
    }
  }

  check_modality("CNV data", cnv, "CNV")
  
  if (!is.null(cnv)) {
    if (is.null(cnv_cache)) {
      # Try to auto-resolve using project path logic
      try({
        candidate_cache <- omicsflow_path("configs", "gene_coordinates.rds")
        candidate_cache2 <- omicsflow_path("realistic_cache.rds")
        if (file.exists(candidate_cache)) {
          cnv_cache <- candidate_cache
        } else if (file.exists(candidate_cache2)) {
          cnv_cache <- candidate_cache2
        }
      }, silent = TRUE)
    }
    if (is.null(cnv_cache) || !file.exists(cnv_cache)) {
      log_fail("CNV enabled but no valid cache available. Please provide 'cnv_cache' or run 'generate_cnv_cache()'.")
    } else {
      log_ok(sprintf("CNV cache resolved: %s", cnv_cache))
    }
  }
  check_modality("SNV data", snv, "SNV")
  
  check_exists("Metadata CSV", metadata)
  check_exists("Clinical table", clinical)
  check_exists("Clinical mapping", clinical_map)
  
  if (valid && provided_modalities < 2) {
    log_fail(sprintf("Fewer than 2 modalities provided (detected %d). OmicsFlow requires at least 2 modalities.", provided_modalities))
  }
  
  # Stop early if critical files are missing or minimum modality rule is violated
  if (!valid) {
    res <- list(
      status = "failed",
      valid = FALSE,
      detected_modalities = detected_mod_names,
      skipped_modalities = skipped_mod_names,
      messages = messages,
      warnings = warnings_list,
      errors = errors_list,
      runtime = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    )
    class(res) <- "omicsflow_validation"
    return(res)
  }
  
  # 2. Metadata Schema Validation
  df_meta <- NULL
  if (!is.null(metadata)) {
    tryCatch({
      df_meta <- read.csv(metadata, stringsAsFactors = FALSE, check.names = FALSE)
      validate_metadata_schema(df_meta)
      log_ok("Metadata matches required schema (sample_id, patient_id, sample_class, batch)")
    }, error = function(e) {
      log_fail(sprintf("Metadata schema invalid: %s", e$message))
    })
  } else {
    log_fail("Metadata path is required for pipeline execution.")
  }
  
  # 3. Align Omics Sample IDs with Metadata
  if (!is.null(df_meta)) {
    meta_sids <- unique(df_meta$sample_id)
    
    check_sids_alignment <- function(label, path, sids_extractor) {
      if (!is.null(path)) {
        extracted_sids <- sids_extractor(path)
        if (is.null(extracted_sids) || length(extracted_sids) == 0) {
          log_fail(sprintf("No samples detected in %s file: %s", label, path))
        } else {
          # Check duplicates
          if (length(extracted_sids) != length(unique(extracted_sids))) {
            log_warn(sprintf("Duplicate sample IDs found in %s file: %s", label, path))
          }
          
          # Check alignment
          missing_sids <- setdiff(extracted_sids, meta_sids)
          if (length(missing_sids) > 0) {
            log_fail(sprintf("%d sample IDs in %s raw data are missing from metadata: %s", 
                              length(missing_sids), label, paste(head(missing_sids, 5), collapse = ", ")))
          } else {
            log_ok(sprintf("All %d sample IDs in %s align correctly with metadata", length(extracted_sids), label))
          }
          
          # 6. Check for zero-variance or all-NA rows (Lightweight, first 100 rows only)
          # Only check for matrices (RNA, Methylation)
          if (label %in% c("RNA-seq matrix", "Methylation matrix")) {
            mat <- read_matrix_safe(path)
            if (!is.null(mat)) {
              sub_mat <- mat[1:min(100, nrow(mat)), , drop = FALSE]
              na_counts <- rowSums(is.na(sub_mat))
              all_na_rows <- sum(na_counts == ncol(sub_mat))
              if (all_na_rows > 0) {
                log_warn(sprintf("Detected %d completely NA rows in first 100 features of %s", all_na_rows, label))
              }
              
              row_vars <- apply(sub_mat, 1, var, na.rm = TRUE)
              zero_var_rows <- sum(row_vars == 0, na.rm = TRUE)
              if (zero_var_rows > 0) {
                log_warn(sprintf("Detected %d zero-variance rows in first 100 features of %s", zero_var_rows, label))
              }
            }
          }
        }
      }
    }
    
    # RNA
    check_sids_alignment("RNA-seq matrix", rna, function(p) {
      m <- read_matrix_safe(p)
      if (is.null(m)) NULL else colnames(m)
    })
    
    # Methylation
    check_sids_alignment("Methylation matrix", meth, function(p) {
      m <- read_matrix_safe(p)
      if (is.null(m)) NULL else colnames(m)
    })
    
    # CNV
    check_sids_alignment("CNV data", cnv, function(p) {
      d <- read_table_safe(p)
      if (is.null(d)) return(NULL)
      if ("Sample" %in% colnames(d)) unique(as.character(d$Sample)) else unique(as.character(d[[1]]))
    })
    
    # SNV
    check_sids_alignment("SNV data", snv, function(p) {
      d <- read_table_safe(p)
      if (is.null(d)) return(NULL)
      if ("Tumor_Sample_Barcode" %in% colnames(d)) unique(as.character(d$Tumor_Sample_Barcode)) else unique(as.character(d[[1]]))
    })
  }
  
  # 4. Clinical mapping checks
  if (!is.null(clinical)) {
    if (!file.exists(clinical)) {
      log_fail(sprintf("Clinical file does not exist: %s", clinical))
    } else {
      df_cl <- read_table_safe(clinical)
      if (is.null(df_cl)) {
        log_fail("Failed to parse clinical table.")
      } else {
        # Check clinical map
        if (!is.null(clinical_map) && file.exists(clinical_map)) {
          tryCatch({
            cmap <- jsonlite::fromJSON(clinical_map)
            # Verify clinical column mapping keys resolve to real columns
            required_keys <- c("patient_id", "os_time", "os_event")
            for (k in required_keys) {
              if (k %in% names(cmap)) {
                mapped_col <- cmap[[k]]
                if (!mapped_col %in% colnames(df_cl)) {
                  log_fail(sprintf("Clinical column map key '%s' resolves to column '%s' which does not exist in clinical file.", k, mapped_col))
                } else {
                  log_ok(sprintf("Clinical map key '%s' matches column '%s'", k, mapped_col))
                }
              } else {
                log_fail(sprintf("Clinical map is missing key: %s", k))
              }
            }
          }, error = function(e) {
            log_fail(sprintf("Failed to read clinical map JSON: %s", e$message))
          })
        } else {
          log_warn("No clinical mapping file provided. Clinical parser will fallback to auto-detection.")
        }
      }
    }
  }
  
  if (valid) {
    log_ok("Input validation PASSED. Ready for pipeline execution.")
  } else {
    log_fail("Input validation FAILED. Resolve errors before running the pipeline.")
  }
  
  res <- list(
    status = if (valid) "success" else "failed",
    valid = valid,
    detected_modalities = detected_mod_names,
    skipped_modalities = skipped_mod_names,
    messages = messages,
    warnings = warnings_list,
    errors = errors_list,
    runtime = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  )
  class(res) <- "omicsflow_validation"
  return(res)
}
