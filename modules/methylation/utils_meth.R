# ==============================================================================
# utils_meth.R — Shared Utilities for the Methylation Preprocessing Module
# OmicsFlow — Phase 1: Methylation Module
# ==============================================================================
# PURPOSE:
#   Provides helper functions used across all methylation module scripts:
#     - Formatted console logging (same style as utils_rna.R)
#     - TCGA barcode parsing (sample type, plate ID, TSS, patient ID)
#     - Directory management
#     - Package loading
#     - Probe filter audit log (accumulated and exported as CSV)
#     - QC metrics accumulation and JSON export
#
# SCIENTIFIC NOTE:
#   This file contains NO scientific logic. All preprocessing thresholds
#   and methods live exclusively in preprocess_meth.R.
#
# USAGE:
#   source("modules/methylation/utils_meth.R")
# ==============================================================================


# ==============================================================================
# LOGGING
# ==============================================================================

meth_banner <- function(msg) {
  width <- 64
  message(paste0("\u2554", strrep("\u2550", width), "\u2557"))
  message(sprintf("\u2551  %-*s\u2551", width - 2, msg))
  message(paste0("\u255a", strrep("\u2550", width), "\u255d"))
}

meth_step <- function(step_num, title) {
  width <- 62
  message(sprintf("\n\u250c%s\u2510", strrep("\u2500", width)))
  message(sprintf("\u2502  STEP %s: %-*s\u2502", step_num, width - 10, title))
  message(sprintf("\u2514%s\u2518", strrep("\u2500", width)))
}

meth_msg <- function(..., level = "INFO") {
  prefix <- switch(level,
    "INFO"   = "  \u2713",
    "WARN"   = "  \u26a0",
    "ERROR"  = "  \u274c",
    "DETAIL" = "   ",
    "  \u2022"
  )
  message(paste0(prefix, " ", paste(..., sep = "")))
}

meth_stat <- function(label, value) {
  message(sprintf("    %-40s %s", label, value))
}


# ==============================================================================
# TCGA BARCODE PARSERS
# (Exact logic preserved from master script helper functions, lines 720-737)
# ==============================================================================

#' Extract TCGA sample type code (e.g. "01" = primary tumor)
#' Uses strsplit on "-" and takes positions 4[1:2]
meth_extract_sample_type <- function(barcodes) {
  sapply(strsplit(as.character(barcodes), "-"),
         function(x) if (length(x) >= 4) substr(x[4], 1, 2) else NA_character_)
}

#' Extract Plate ID for batch correction
#' Position 6 after splitting on "-"
meth_extract_plate_id <- function(barcodes) {
  sapply(strsplit(as.character(barcodes), "-"),
         function(x) if (length(x) >= 6) x[6] else NA_character_)
}

#' Extract Tissue Source Site (TSS) — used as biological covariate in ComBat
#' Position 2 after splitting on "-"
meth_extract_tss <- function(barcodes) {
  sapply(strsplit(as.character(barcodes), "-"),
         function(x) if (length(x) >= 2) x[2] else NA_character_)
}

#' Truncate any TCGA barcode to the 12-char patient ID
meth_extract_patient_id <- function(barcodes) {
  substr(as.character(barcodes), 1, 12)
}

#' Safely retrieve a column from the 450k annotation by trying candidate names.
#' Fails hard if none found (preserves probe filtering correctness).
meth_safe_get_col <- function(anno, candidates) {
  hit <- intersect(candidates, names(anno))
  if (length(hit) == 0) {
    stop("Cannot find annotation column. Tried: ", paste(candidates, collapse = ", "))
  }
  anno[[hit[1]]]
}

#' Map clinical stage strings to standardized 4-level factor
#' Exact logic from master script map_stage(), lines 746-754
meth_map_stage <- function(stage_vec) {
  stage_vec <- toupper(as.character(stage_vec))
  out <- rep(NA_character_, length(stage_vec))
  out[grepl("STAGE IV",  stage_vec)] <- "Stage IV"
  out[grepl("STAGE III", stage_vec)] <- "Stage III"
  out[grepl("STAGE II",  stage_vec)] <- "Stage II"
  out[grepl("STAGE I",   stage_vec)] <- "Stage I"
  return(out)
}


# ==============================================================================
# PROBE FILTER AUDIT LOG
# ==============================================================================

#' Initialize an empty probe filter log list (populated via meth_log_probes)
meth_init_probe_log <- function() {
  list()
}

#' Record current probe and sample count for an audit step
#'
#' @param probe_log  Named list accumulated across pipeline steps
#' @param step       Step label string (e.g. "After Probe filter")
#' @param mat        Numeric matrix at this checkpoint
#' @return Updated probe_log list (invisibly)
meth_log_probes <- function(probe_log, step, mat) {
  probe_log[[step]] <- nrow(mat)
  meth_msg(sprintf("[%s] Probes: %d | Samples: %d", step, nrow(mat), ncol(mat)),
           level = "DETAIL")
  invisible(probe_log)
}


# ==============================================================================
# DIRECTORY MANAGEMENT
# ==============================================================================

#' Create all required methylation output subdirectories
meth_ensure_dirs <- function(outdir) {
  dirs <- c(outdir, file.path(outdir, "plots"))
  for (d in dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(dirs)
}


# ==============================================================================
# PACKAGE LOADING
# ==============================================================================

#' Load all packages required by the methylation module
meth_load_packages <- function() {
  required <- c(
    "matrixStats", "limma", "sva", "impute",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "SummarizedExperiment",
    "ggplot2", "gridExtra", "reshape2",
    "optparse", "jsonlite"
  )

  missing_pkgs <- required[!vapply(required, requireNamespace,
                                    quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required packages for methylation module:\n  ",
      paste(missing_pkgs, collapse = ", "),
      "\nInstall with: BiocManager::install(c('",
      paste(missing_pkgs, collapse = "','"), "'))"
    )
  }

  suppressPackageStartupMessages({
    library(matrixStats)
    library(limma)
    library(sva)
    library(impute)
    library(IlluminaHumanMethylation450kmanifest)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(SummarizedExperiment)
    library(ggplot2)
    library(gridExtra)
    library(reshape2)
    library(jsonlite)
  })

  invisible(TRUE)
}


# ==============================================================================
# QC METRICS ACCUMULATOR
# ==============================================================================

#' Initialize an empty QC metrics list for the methylation module
meth_init_qc <- function() {
  list(
    pipeline_version = "2.0.1",
    omics_layer      = "DNA Methylation (450k)",
    timestamp        = as.character(Sys.time()),

    # Sample-level counts
    samples = list(
      raw                      = NA_integer_,
      primary_tumor            = NA_integer_,
      after_deduplication      = NA_integer_,
      after_batch_na_removal   = NA_integer_,
      after_singleton_removal  = NA_integer_,
      after_clinical_filter    = NA_integer_,
      final                    = NA_integer_
    ),

    # Probe-level counts at each filter step
    probes = list(
      raw                      = NA_integer_,
      after_dedup              = NA_integer_,
      after_detection_pval     = NA_integer_,
      after_probe_filter       = NA_integer_,
      after_na_filter          = NA_integer_,
      selected_top_variable    = NA_integer_
    ),

    # Filter breakdown
    filters = list(
      cross_reactive_removed   = NA_integer_,
      snp_probes_removed       = NA_integer_,
      sex_chromosome_removed   = NA_integer_,
      non_cpg_removed          = NA_integer_,
      total_bad_probes         = NA_integer_,
      detection_pval_removed   = NA_integer_,
      na_fraction_removed      = NA_integer_
    ),

    # Outliers (PCA-based)
    outliers = list(
      pc1_outliers             = NA_integer_,
      pc2_outliers             = NA_integer_,
      total_flagged            = NA_integer_,
      outlier_ids              = list()
    ),

    # Batch correction
    batches = list(
      total_detected           = NA_integer_,
      singleton_batches        = NA_integer_,
      combat_applied           = NA,
      covariates_protected     = NA_character_
    ),

    # Imputation
    imputation = list(
      method_used              = NA_character_,
      missing_fraction         = NA_real_,
      probes_imputed           = NA_integer_
    ),

    # M-value range statistics
    m_value_stats = list(
      global_min               = NA_real_,
      global_max               = NA_real_,
      global_mean              = NA_real_
    ),

    # Output dimensions
    output_matrices = list(
      mofa_probes              = NA_integer_,
      mofa_samples             = NA_integer_,
      ml_full_probes           = NA_integer_,
      ml_full_samples          = NA_integer_
    )
  )
}

#' Update a nested field in the QC metrics list
meth_update_qc <- function(qc, section, key, value) {
  qc[[section]][[key]] <- value
  qc
}


# ==============================================================================
# JSON EXPORT
# ==============================================================================

#' Serialize QC metrics to qc_metrics.json
meth_save_qc_json <- function(qc_metrics, outdir) {
  clean <- rapply(qc_metrics, function(x) {
    if (length(x) == 1 && is.na(x)) NULL else x
  }, how = "replace")

  json_path <- file.path(outdir, "qc_metrics.json")
  jsonlite::write_json(clean, json_path, pretty = TRUE, auto_unbox = TRUE,
                       null = "null")
  meth_msg(sprintf("Saved: %s", json_path))
  invisible(json_path)
}


# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

#' Validate that a required file path exists
meth_check_file <- function(path, label = "Input file") {
  if (is.null(path) || !nzchar(path)) {
    stop(label, " path is required but was not provided.")
  }
  if (!file.exists(path)) {
    stop(label, " not found: ", path)
  }
  invisible(TRUE)
}

#' Check that column names look like 12-char TCGA patient barcodes
meth_validate_barcodes <- function(mat, label = "matrix", metadata_supplied = FALSE) {
  if (isTRUE(metadata_supplied)) {
    return(invisible(TRUE))
  }
  ids <- colnames(mat)
  if (!all(nchar(ids) == 12)) {
    warning(sprintf("[%s] %d / %d column names are not exactly 12 characters.",
                    label, sum(nchar(ids) != 12), length(ids)))
  }
  if (!all(grepl("^TCGA-", ids))) {
    warning(sprintf("[%s] Some column names do not start with 'TCGA-'.", label))
  }
  invisible(TRUE)
}
