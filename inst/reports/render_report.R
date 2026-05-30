#!/usr/bin/env Rscript
# =============================================================================
# reports/render_report.R — OmicsFlow Report Rendering Script
# =============================================================================
#
# PURPOSE:
#   Renders the OmicsFlow_Report.qmd into a self-contained HTML file.
#   Works with either Quarto CLI or rmarkdown::render() as fallback.
#
# USAGE:
#   Rscript reports/render_report.R [--results_dir results] [--output_dir results/reports]
#
# =============================================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--results_dir", type = "character", default = "results",
              help = "Path to pipeline results directory [default: results]"),
  make_option("--output_dir", type = "character", default = "results/reports",
              help = "Output directory for the HTML report [default: results/reports]"),
  make_option("--report_qmd", type = "character", default = NULL,
              help = "Path to QMD template [default: auto-detect]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Locate report template ---
if (is.null(opt$report_qmd)) {
  script_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) "reports"
  )
  candidates <- c(
    "reports/OmicsFlow_Report.qmd",
    file.path(script_dir, "OmicsFlow_Report.qmd")
  )
  opt$report_qmd <- Filter(file.exists, candidates)[1]
}

if (is.null(opt$report_qmd) || !file.exists(opt$report_qmd)) {
  stop("Report template not found. Specify --report_qmd path.", call. = FALSE)
}

# --- Create output directory ---
if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  OmicsFlow v1.0.0 — Report Generation                      ║\n")
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Template   : %s\n", opt$report_qmd))
cat(sprintf("║  Results    : %s\n", opt$results_dir))
cat(sprintf("║  Output     : %s\n", opt$output_dir))
cat(sprintf("║  Timestamp  : %s\n", Sys.time()))
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# --- Check for Quarto CLI ---
quarto_available <- tryCatch({
  ret <- system2("quarto", "--version", stdout = TRUE, stderr = TRUE)
  TRUE
}, error = function(e) FALSE,
   warning = function(w) FALSE)

output_file <- file.path(normalizePath(opt$output_dir), "OmicsFlow_Report.html")

if (quarto_available) {
  cat("  ✓ Quarto CLI detected. Rendering with Quarto...\n")

  # Copy files to temp location for rendering
  tmp_dir <- tempdir()
  file.copy(opt$report_qmd, file.path(tmp_dir, "OmicsFlow_Report.qmd"), overwrite = TRUE)

  css_src <- file.path(dirname(opt$report_qmd), "styles.css")
  if (file.exists(css_src)) {
    file.copy(css_src, file.path(tmp_dir, "styles.css"), overwrite = TRUE)
  }

  system2("quarto", c(
    "render", file.path(tmp_dir, "OmicsFlow_Report.qmd"),
    "--to", "html",
    "-P", sprintf("results_dir:%s", normalizePath(opt$results_dir)),
    "--output-dir", normalizePath(opt$output_dir)
  ))

} else {
  cat("  ⚠ Quarto CLI not found. Rendering with rmarkdown/knitr...\n")

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Neither Quarto nor rmarkdown are available. Install one of them.", call. = FALSE)
  }

  # rmarkdown can render .qmd files with knitr
  rmarkdown::render(
    input       = opt$report_qmd,
    output_file = output_file,
    output_format = rmarkdown::html_document(
      theme       = "cosmo",
      toc         = TRUE,
      toc_depth   = 3,
      toc_float   = TRUE,
      number_sections = TRUE,
      self_contained  = TRUE,
      css         = file.path(dirname(opt$report_qmd), "styles.css"),
      code_folding = "hide"
    ),
    params = list(results_dir = normalizePath(opt$results_dir)),
    envir  = new.env(parent = globalenv()),
    quiet  = FALSE
  )
}

if (file.exists(output_file)) {
  cat(sprintf("\n  ✓ Report generated: %s\n", output_file))
  cat(sprintf("  ✓ File size: %s bytes\n", format(file.size(output_file), big.mark = ",")))
  cat("\n  ╔══════════════════════════════════════════════╗\n")
  cat("  ║  REPORT GENERATION COMPLETE ✔                ║\n")
  cat("  ╚══════════════════════════════════════════════╝\n\n")
} else {
  cat("\n  ✖ Report generation failed.\n")
  quit(status = 1)
}
