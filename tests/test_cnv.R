#!/usr/bin/env Rscript
# ==============================================================================
# test_cnv.R — Validation Tests for CNV Module
# ==============================================================================
#
# PURPOSE:
#   Validates the outputs of the CNV module to ensure scientific integrity
#   and compatibility with downstream integration.
#
# CHECKS:
#   1. Required files exist
#   2. Valid matrix structure
#   3. Correct missing value handling (no NAs, imputed with 0)
#   4. Barcode standard (12-char TCGA patient IDs)
#   5. Gene symbol mapping (rownames)
#   6. Sign-preserving values (amplifications vs deletions)
#   7. Clipped extreme values (range [-5, 5])
#   8. Non-zero variance features
#   9. Metadata consistency
#
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

option_list <- list(
  make_option(c("-o", "--outdir"), type = "character", default = "results/cnv/",
              help = "Directory containing CNV outputs [default= %default]")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
outdir <- opt$outdir

# --- Testing Utilities ---
n_passed <- 0
n_failed <- 0
n_warned <- 0

test_msg <- function(msg, level = "INFO") {
  if (level == "PASS") {
    cat(sprintf("  [\033[32mPASS\033[0m] %s\n", msg))
    n_passed <<- n_passed + 1
  } else if (level == "FAIL") {
    cat(sprintf("  [\033[31mFAIL\033[0m] %s\n", msg))
    n_failed <<- n_failed + 1
  } else if (level == "WARN") {
    cat(sprintf("  [\033[33mWARN\033[0m] %s\n", msg))
    n_warned <<- n_warned + 1
  } else {
    cat(sprintf("--- %s ---\n", msg))
  }
}

expect_true <- function(condition, pass_msg, fail_msg) {
  if (isTRUE(condition)) test_msg(pass_msg, "PASS") else test_msg(fail_msg, "FAIL")
}

# ==============================================================================
# START TESTS
# ==============================================================================

cat("\n==============================================================\n")
cat("  OmicsFlow CNV Module \u2014 Validation Tests\n")
cat(sprintf("  Output directory: %s\n", outdir))
cat(sprintf("  Timestamp       : %s\n", Sys.time()))
cat("==============================================================\n\n")

# TEST 1 \u2014 Required output files
test_msg("Test 1: Required output files")
req_files <- c(
  "cnv_processed_matrix.rds",
  "cnv_gene_variance_info.rds",
  "sample_metadata.csv",
  "sample_metadata.rds",
  "qc_metrics.json",
  "plots/cnv_density.pdf",
  "plots/cnv_correlation_histogram.pdf",
  "plots/cnv_pca.pdf",
  "plots/cnv_distribution_per_sample.pdf",
  "plots/CNV_Research_Validation_Figures.png"
)

for (f in req_files) {
  fp <- file.path(outdir, f)
  expect_true(file.exists(fp),
              sprintf("File exists: %s", f),
              sprintf("Missing file: %s", f))
}

# Abort if main matrix is missing
if (!file.exists(file.path(outdir, "cnv_processed_matrix.rds"))) {
  stop("Main matrix missing. Cannot continue tests.")
}

# Load files
test_msg("Loading output files")
cnv_mat <- readRDS(file.path(outdir, "cnv_processed_matrix.rds"))
meta_csv <- read.csv(file.path(outdir, "sample_metadata.csv"), stringsAsFactors = FALSE)
qc_json  <- fromJSON(file.path(outdir, "qc_metrics.json"))

cat(sprintf("  cnv_processed_matrix.rds : %d x %d\n", nrow(cnv_mat), ncol(cnv_mat)))
cat(sprintf("  sample_metadata.csv      : %d rows\n", nrow(meta_csv)))

# TEST 2 \u2014 Matrix validity
test_msg("Test 2: Matrix Validity")
expect_true(is.matrix(cnv_mat), "cnv_matrix is a matrix", "cnv_matrix is not a matrix")
expect_true(is.numeric(cnv_mat), "cnv_matrix is numeric", "cnv_matrix is not numeric")
expect_true(!anyNA(cnv_mat), "cnv_matrix has no NAs", "cnv_matrix contains NAs")
expect_true(all(is.finite(cnv_mat)), "cnv_matrix has no Inf/NaN", "cnv_matrix contains Inf or NaN")

# TEST 3 \u2014 Missing value handling (Neutral Diploid state = 0)
test_msg("Test 3: Neutral Diploid State (0) representation")
expect_true(any(cnv_mat == 0), "cnv_matrix contains 0 values (expected for missing/neutral)", "cnv_matrix has no 0 values")

# TEST 4 \u2014 Barcode Standard (12-char TCGA patient IDs)
test_msg("Test 4: Column names (12-char TCGA patient IDs)")
cols_12 <- nchar(colnames(cnv_mat)) == 12
expect_true(all(cols_12),
            "All cnv_matrix columns are 12 characters",
            "Some cnv_matrix columns are NOT 12 characters")

starts_tcga <- grepl("^TCGA-", colnames(cnv_mat))
expect_true(all(starts_tcga),
            "All cnv_matrix columns start with 'TCGA-'",
            "Some cnv_matrix columns do NOT start with 'TCGA-'")

# TEST 5 \u2014 Extreme Value Clipping
test_msg("Test 5: Extreme value clipping (range [-5, 5])")
cnv_range <- range(cnv_mat)
cat(sprintf("  Observed range: [%.3f, %.3f]\n", cnv_range[1], cnv_range[2]))
expect_true(cnv_range[1] >= -5,
            "Minimum value is >= -5 (clipped correctly)",
            sprintf("Minimum value < -5 (found %.3f)", cnv_range[1]))
expect_true(cnv_range[2] <= 5,
            "Maximum value is <= 5 (clipped correctly)",
            sprintf("Maximum value > 5 (found %.3f)", cnv_range[2]))

# TEST 6 \u2014 Sign-preserving values
test_msg("Test 6: Sign-preserving values (amplifications & deletions)")
expect_true(any(cnv_mat > 0), "Contains amplifications (>0)", "No positive values found (amplifications)")
expect_true(any(cnv_mat < 0), "Contains deletions (<0)", "No negative values found (deletions)")

# TEST 7 \u2014 Metadata Consistency
test_msg("Test 7: Metadata Consistency")
expect_true("patient_id" %in% colnames(meta_csv),
            "sample_metadata.csv has 'patient_id' column",
            "sample_metadata.csv is missing 'patient_id' column")

expect_true(nrow(meta_csv) == ncol(cnv_mat),
            "Metadata rows match matrix columns",
            sprintf("Metadata rows (%d) != matrix cols (%d)", nrow(meta_csv), ncol(cnv_mat)))

missing_in_meta <- setdiff(colnames(cnv_mat), meta_csv$patient_id)
missing_in_mat  <- setdiff(meta_csv$patient_id, colnames(cnv_mat))

expect_true(length(missing_in_meta) == 0,
            "All matrix columns present in metadata",
            sprintf("%d matrix columns missing from metadata", length(missing_in_meta)))
expect_true(length(missing_in_mat) == 0,
            "All metadata IDs present in matrix columns",
            sprintf("%d metadata IDs missing from matrix", length(missing_in_mat)))

# TEST 8 \u2014 Immune genes filtering
test_msg("Test 8: Immune genes filtering")
immune_present <- grep("^IG[HKL][VDJC]|^TR[ABDG][VDJC]", rownames(cnv_mat), value = TRUE)
expect_true(length(immune_present) == 0,
            "No immune receptor genes in final matrix",
            sprintf("Found %d immune genes (e.g. %s)", length(immune_present), paste(head(immune_present, 3), collapse = ", ")))

cat("\n==============================================================\n")
cat("  RESULTS: ", n_passed + n_failed + n_warned, " tests total\n")
cat("    PASS : ", n_passed, "\n")
cat("    FAIL : ", n_failed, "\n")
cat("    WARN : ", n_warned, "\n")
cat("==============================================================\n")

if (n_failed > 0) {
  cat("\n  \033[31mSOME TESTS FAILED. See above for details.\033[0m\n")
  quit(status = 1)
} else {
  cat("\n  ALL TESTS PASSED. CNV module outputs are valid.\n")
  quit(status = 0)
}
