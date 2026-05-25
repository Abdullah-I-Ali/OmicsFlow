# ==============================================================================
# utils_integration.R — Shared Utilities for the Integration Layer (MOFA+)
# OmicsFlow — Phase 5: Integration Module
# ==============================================================================

# ==============================================================================
# LOGGING
# ==============================================================================

#' Print a top-level pipeline banner
int_banner <- function(msg) {
  width <- 64
  border <- paste0("\u2554", strrep("\u2550", width), "\u2557")
  inner  <- sprintf("\u2551  %-*s\u2551", width - 2, msg)
  footer <- paste0("\u255a", strrep("\u2550", width), "\u255d")
  message(border)
  message(inner)
  message(footer)
}

#' Print a section step header
int_step <- function(step_num, title) {
  width <- 62
  message(sprintf("\n\u250c%s\u2510", strrep("\u2500", width)))
  message(sprintf("\u2502  STEP %s: %-*s\u2502", step_num, width - 10, title))
  message(sprintf("\u2514%s\u2518", strrep("\u2500", width)))
}

#' Print an indented info message
int_msg <- function(..., level = "INFO") {
  prefix <- switch(level,
                   "INFO"    = "  \u2713",
                   "WARN"    = "  \u26a0",
                   "ERROR"   = "  \u2716",
                   "DETAILS" = "    \u2022",
                   "  \u2794")
  message(sprintf("%s %s", prefix, paste0(...)))
}

# ==============================================================================
# PACKAGE MANAGEMENT
# ==============================================================================

#' Load all packages required for Integration preprocessing
load_int_packages <- function() {
  suppressPackageStartupMessages({
    library(MOFA2)
    library(dplyr)
    library(ggplot2)
    library(gridExtra)
    library(jsonlite)
  })
}

# ==============================================================================
# DATA LOADING
# ==============================================================================

#' Read omics matrix
load_omics_matrix <- function(file_path, name, metadata = NULL) {
  if (!file.exists(file_path)) {
    stop(sprintf("Input %s matrix not found: %s", name, file_path))
  }
  mat <- readRDS(file_path)
  
  # Standardize barcodes to 12 chars if needed
  if (is.matrix(mat) || is.data.frame(mat)) {
    if (is.null(metadata)) {
      colnames(mat) <- substr(colnames(mat), 1, 12)
    } else {
      cols <- colnames(mat)
      meta_pids <- unique(metadata$patient_id)
      if (all(cols %in% meta_pids)) {
        # Already standardized to patient_ids
      } else {
        # Try remapping via metadata sample_id -> patient_id
        # Source modules/utils_metadata.R if functions are not available
        if (!exists("get_patient_id")) {
          # Try to find utils_metadata.R
          for (p in c("modules/utils_metadata.R", "../utils_metadata.R", "utils_metadata.R")) {
            if (file.exists(p)) {
              source(p)
              break
            }
          }
        }
        
        if (exists("get_patient_id")) {
          tryCatch({
            colnames(mat) <- get_patient_id(cols, metadata)
          }, error = function(e) {
            # Fallback to TCGA truncation if mapping fails
            colnames(mat) <- substr(cols, 1, 12)
          })
        } else {
          colnames(mat) <- substr(cols, 1, 12)
        }
      }
    }
  }
  
  return(mat)
}

# ==============================================================================
# QC TRACKING SYSTEM
# ==============================================================================

#' Initialize the QC tracking list
init_int_qc <- function() {
  list(
    samples = list(),
    features = list(),
    model = list(),
    metadata = list(
      pipeline_version = "2.0",
      timestamp = as.character(Sys.time())
    )
  )
}

#' Add a key-value pair to a specific QC category
add_int_qc <- function(qc_list, category, key, value) {
  qc_list[[category]][[key]] <- value
  return(qc_list)
}

#' Export the complete QC tracking list to a JSON file
export_int_qc <- function(qc_list, outdir) {
  out_path <- file.path(outdir, "qc_metrics.json")
  jsonlite::write_json(
    qc_list,
    path = out_path,
    auto_unbox = TRUE,
    pretty = TRUE
  )
  int_msg(sprintf("Saved QC metrics: %s", out_path))
}
