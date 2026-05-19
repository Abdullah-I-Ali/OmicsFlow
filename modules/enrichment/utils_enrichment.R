# ==============================================================================
# utils_enrichment.R — Shared Utilities for Pathway Enrichment Module
# OmicsFlow — Phase 7: Enrichment Module
# ==============================================================================

# ==============================================================================
# LOGGING
# ==============================================================================

#' Print a top-level pipeline banner
path_banner <- function(msg) {
  width <- 64
  border <- paste0("\u2554", strrep("\u2550", width), "\u2557")
  inner  <- sprintf("\u2551  %-*s\u2551", width - 2, msg)
  footer <- paste0("\u255a", strrep("\u2550", width), "\u255d")
  message(border)
  message(inner)
  message(footer)
}

#' Print a section step header
path_step <- function(step_num, title) {
  width <- 62
  message(sprintf("\n\u250c%s\u2510", strrep("\u2500", width)))
  message(sprintf("\u2502  STEP %s: %-*s\u2502", step_num, width - 10, title))
  message(sprintf("\u2514%s\u2518", strrep("\u2500", width)))
}

#' Print an indented info message
path_msg <- function(..., level = "INFO") {
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

#' Load all packages required for Pathway enrichment
load_pathway_packages <- function() {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(ggplot2)
    library(stringr)
    library(jsonlite)
  })
}

# ==============================================================================
# QC TRACKING SYSTEM
# ==============================================================================

#' Initialize the QC tracking list
init_pathway_qc <- function() {
  list(
    genes = list(),
    enrichment = list(),
    metadata = list(
      pipeline_version = "1.1",
      timestamp = as.character(Sys.time())
    )
  )
}

#' Add a key-value pair to a specific QC category
add_pathway_qc <- function(qc_list, category, key, value) {
  qc_list[[category]][[key]] <- value
  return(qc_list)
}

#' Export the complete QC tracking list to a JSON file
export_pathway_qc <- function(qc_list, outdir) {
  out_path <- file.path(outdir, "qc_metrics.json")
  jsonlite::write_json(
    qc_list,
    path = out_path,
    auto_unbox = TRUE,
    pretty = TRUE
  )
  path_msg(sprintf("Saved QC metrics: %s", out_path))
}
