# ==============================================================================
# utils_rna.R — Shared Utilities for the RNA Preprocessing Module
# OmicsFlow — Phase 1: RNA Module
# ==============================================================================
#
# PURPOSE:
#   Provides helper functions used across all RNA module scripts:
#     - Formatted console logging
#     - Directory management
#     - Package loading
#     - QC metrics accumulation and JSON export
#
# SCIENTIFIC NOTE:
#   This file contains NO scientific logic. All preprocessing thresholds
#   and methods live exclusively in preprocess_rna.R.
#
# USAGE:
#   source("modules/rna/utils_rna.R")
#
# ==============================================================================


# ==============================================================================
# LOGGING
# ==============================================================================

#' Print a top-level pipeline banner
rna_banner <- function(msg) {
  width <- 64
  border <- paste0("\u2554", strrep("\u2550", width), "\u2557")
  inner  <- sprintf("\u2551  %-*s\u2551", width - 2, msg)
  footer <- paste0("\u255a", strrep("\u2550", width), "\u255d")
  message(border)
  message(inner)
  message(footer)
}

#' Print a section step header
rna_step <- function(step_num, title) {
  width <- 62
  message(sprintf("\n\u250c%s\u2510", strrep("\u2500", width)))
  message(sprintf("\u2502  STEP %s: %-*s\u2502", step_num, width - 10, title))
  message(sprintf("\u2514%s\u2518", strrep("\u2500", width)))
}

#' Print an indented info message
rna_msg <- function(..., level = "INFO") {
  prefix <- switch(level,
    "INFO"    = "  \u2713",
    "WARN"    = "  \u26a0",
    "ERROR"   = "  \u274c",
    "DETAIL"  = "   ",
    "  \u2022"
  )
  message(paste0(prefix, " ", paste(..., sep = "")))
}

#' Print a key-value statistic line
rna_stat <- function(label, value) {
  message(sprintf("    %-40s %s", label, value))
}


# ==============================================================================
# DIRECTORY MANAGEMENT
# ==============================================================================

#' Create all required RNA output subdirectories
#'
#' @param outdir  Base output directory (e.g. "results/rna/")
#' @return Invisibly returns the created directory paths
rna_ensure_dirs <- function(outdir) {
  dirs <- c(
    outdir,
    file.path(outdir, "plots")
  )
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
  }
  invisible(dirs)
}


# ==============================================================================
# PACKAGE LOADING
# ==============================================================================

#' Load all packages required by the RNA module
#'
#' Suppresses startup messages for clean pipeline output.
#' Fails fast with a clear message if any package is missing.
rna_load_packages <- function() {
  required <- c(
    "edgeR", "limma", "dplyr", "org.Hs.eg.db", "AnnotationDbi",
    "SummarizedExperiment", "ggplot2", "gridExtra", "pheatmap",
    "RColorBrewer", "ggpubr", "reshape2", "optparse", "jsonlite"
  )

  missing_pkgs <- required[!vapply(required, requireNamespace,
                                    quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required packages for RNA module:\n  ",
      paste(missing_pkgs, collapse = ", "),
      "\nInstall with: BiocManager::install(c('",
      paste(missing_pkgs, collapse = "','"), "'))"
    )
  }

  suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
    library(dplyr)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(SummarizedExperiment)
    library(ggplot2)
    library(gridExtra)
    library(pheatmap)
    library(RColorBrewer)
    library(ggpubr)
    library(reshape2)
    library(jsonlite)
  })

  invisible(TRUE)
}


# ==============================================================================
# QC METRICS ACCUMULATOR
# ==============================================================================

#' Initialize an empty QC metrics list
#'
#' This list is populated incrementally throughout the pipeline and
#' serialized to qc_metrics.json at the end.
#'
#' @return Named list with empty sub-lists for each QC category
rna_init_qc <- function() {
  list(
    pipeline_version = "2.0.1",
    omics_layer      = "RNA-seq",
    timestamp        = as.character(Sys.time()),

    # Sample-level counts (populated during Steps 2, 9, 10)
    samples = list(
      raw                    = NA_integer_,
      primary_tumor          = NA_integer_,
      after_deduplication    = NA_integer_,
      low_depth_flagged      = NA_integer_,   # warned, NOT removed
      after_outlier_removal  = NA_integer_,
      after_singleton_removal= NA_integer_,
      final                  = NA_integer_
    ),

    # Gene/feature counts (populated during Steps 3–8, 11)
    genes = list(
      raw                    = NA_integer_,
      after_zero_removal     = NA_integer_,
      after_cpm_filter       = NA_integer_,
      after_ensembl_mapping  = NA_integer_,
      after_symbol_dedup     = NA_integer_,
      after_na_filter        = NA_integer_,
      after_variance_filter  = NA_integer_,
      selected_top_variable  = NA_integer_,
      final_after_mean_filter= NA_integer_
    ),

    # Outlier detection (Step 9)
    outliers = list(
      low_correlation_threshold  = NA_real_,
      high_correlation_threshold = NA_real_,
      low_outliers_count         = NA_integer_,
      high_outliers_count        = NA_integer_,
      total_removed              = NA_integer_
    ),

    # Batch correction (Step 10)
    batches = list(
      total_batches_detected  = NA_integer_,
      singleton_batches       = NA_integer_,
      batches_corrected       = NA_integer_,
      correction_applied      = NA
    ),

    # Library size statistics (Step 3)
    library_size = list(
      min_reads    = NA_real_,
      max_reads    = NA_real_,
      median_reads = NA_real_
    ),

    # Expression range of final matrices
    output_matrices = list(
      zscore_mean  = NA_real_,
      zscore_sd    = NA_real_,
      ml_range_min = NA_real_,
      ml_range_max = NA_real_,
      genes_mofa   = NA_integer_,
      genes_ml     = NA_integer_,
      samples_mofa = NA_integer_,
      samples_ml   = NA_integer_
    )
  )
}

#' Update a nested field in the QC metrics list
#'
#' @param qc      Existing QC metrics list
#' @param section Top-level section name (e.g. "samples")
#' @param key     Field name within section
#' @param value   Value to set
#' @return Updated QC metrics list
rna_update_qc <- function(qc, section, key, value) {
  qc[[section]][[key]] <- value
  qc
}


# ==============================================================================
# JSON EXPORT
# ==============================================================================

#' Serialize QC metrics to qc_metrics.json
#'
#' @param qc_metrics  Named list produced by rna_init_qc() and populated
#'                    throughout the pipeline
#' @param outdir      Output directory (file written as outdir/qc_metrics.json)
#' @return Invisibly returns the full path of the written file
rna_save_qc_json <- function(qc_metrics, outdir) {
  # Replace NA with null for valid JSON
  clean <- rapply(qc_metrics, function(x) {
    if (length(x) == 1 && is.na(x)) NULL else x
  }, how = "replace")

  json_path <- file.path(outdir, "qc_metrics.json")
  jsonlite::write_json(clean, json_path, pretty = TRUE, auto_unbox = TRUE,
                       null = "null")
  rna_msg(sprintf("Saved: %s", json_path))
  invisible(json_path)
}


# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

#' Validate that a required file path exists
#'
#' @param path  File path to check
#' @param label Human-readable label for error messages
rna_check_file <- function(path, label = "Input file") {
  if (is.null(path) || !nzchar(path)) {
    stop(label, " path is required but was not provided.")
  }
  if (!file.exists(path)) {
    stop(label, " not found: ", path)
  }
  invisible(TRUE)
}

#' Check that column names look like 12-char TCGA patient barcodes
#'
#' @param mat  Matrix whose column names will be checked
#' @param label  Descriptive label for error reporting
rna_validate_barcodes <- function(mat, label = "matrix", metadata_supplied = FALSE) {
  if (isTRUE(metadata_supplied)) {
    return(invisible(TRUE))
  }
  ids <- colnames(mat)
  if (!all(nchar(ids) == 12)) {
    warning(sprintf(
      "[%s] %d / %d column names are not exactly 12 characters.",
      label, sum(nchar(ids) != 12), length(ids)
    ))
  }
  if (!all(grepl("^TCGA-", ids))) {
    warning(sprintf("[%s] Some column names do not start with 'TCGA-'.", label))
  }
  invisible(TRUE)
}
