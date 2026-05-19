# ==============================================================================
# utils_cnv.R — Shared Utilities for the CNV Preprocessing Module
# OmicsFlow — CNV Module
# ==============================================================================

# ==============================================================================
# LOGGING
# ==============================================================================

#' Print a top-level pipeline banner
cnv_banner <- function(msg) {
  width <- 64
  border <- paste0("\u2554", strrep("\u2550", width), "\u2557")
  inner  <- sprintf("\u2551  %-*s\u2551", width - 2, msg)
  footer <- paste0("\u255a", strrep("\u2550", width), "\u255d")
  message(border)
  message(inner)
  message(footer)
}

#' Print a section step header
cnv_step <- function(step_num, title) {
  width <- 62
  message(sprintf("\n\u250c%s\u2510", strrep("\u2500", width)))
  message(sprintf("\u2502  STEP %s: %-*s\u2502", step_num, width - 10, title))
  message(sprintf("\u2514%s\u2518", strrep("\u2500", width)))
}

#' Print an indented info message
cnv_msg <- function(..., level = "INFO") {
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

#' Load all packages required for CNV preprocessing
load_cnv_packages <- function() {
  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(matrixStats)
    library(biomaRt)
    library(data.table)
    library(ggplot2)
    library(gridExtra)
    library(reshape2)
    library(jsonlite)
  })
}

# ==============================================================================
# DATA LOADING
# ==============================================================================

#' Read raw CNV segment data and Ensembl gene coordinates
#' Ensures consistent structures before preprocessing
load_cnv_data <- function(cnv_file, gene_cache_file) {
  if (!file.exists(cnv_file)) {
    stop("Input CNV file not found: ", cnv_file)
  }
  
  cnv_raw <- readRDS(cnv_file)
  cnv_data <- as.data.frame(cnv_raw)
  
  # Ensure necessary columns are present
  required_cols <- c("GDC_Aliquot", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean", "Sample")
  missing_cols <- setdiff(required_cols, colnames(cnv_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in CNV data: ", paste(missing_cols, collapse = ", "))
  }
  
  gene_coords <- NULL
  if (file.exists(gene_cache_file)) {
    cnv_msg("Loading cached gene coordinates...", level = "INFO")
    gene_coords <- readRDS(gene_cache_file)
  }
  
  return(list(
    cnv = cnv_data,
    gene_coords = gene_coords
  ))
}

# ==============================================================================
# QC TRACKING SYSTEM
# ==============================================================================

#' Initialize the QC tracking list
init_cnv_qc <- function() {
  list(
    samples = list(),
    genes = list(),
    filters = list(),
    parameters = list(),
    metadata = list(
      pipeline_version = "2.0",
      timestamp = as.character(Sys.time())
    )
  )
}

#' Add a key-value pair to a specific QC category
add_cnv_qc <- function(qc_list, category, key, value) {
  qc_list[[category]][[key]] <- value
  return(qc_list)
}

#' Export the complete QC tracking list to a JSON file
export_cnv_qc <- function(qc_list, outdir) {
  out_path <- file.path(outdir, "qc_metrics.json")
  jsonlite::write_json(
    qc_list,
    path = out_path,
    auto_unbox = TRUE,
    pretty = TRUE
  )
  cnv_msg(sprintf("Saved QC metrics: %s", out_path))
}
