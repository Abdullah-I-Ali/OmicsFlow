#!/usr/bin/env Rscript
# ==============================================================================
# test_snv.R — Validation Tests for SNV Module
# ==============================================================================
#
# PURPOSE:
#   Validates the outputs of the SNV module to ensure scientific integrity
#   and compatibility with downstream integration.
#
# CHECKS:
#   1. Required files exist
#   2. Valid binary matrix structure
#   3. Correct missing value handling (no NAs)
#   4. Barcode standard (12-char TCGA patient IDs)
#   5. Valid binary content (only 0 and 1)
#   6. Dimensionality checks
#   7. Metadata consistency
#
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

option_list <- list(
  make_option(c("-o", "--outdir"), type = "character", default = "results/snv/",
              help = "Directory containing SNV outputs [default= %default]")
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
cat("  OmicsFlow SNV Module \u2014 Validation Tests\n")
cat(sprintf("  Output directory: %s\n", outdir))
cat(sprintf("  Timestamp       : %s\n", Sys.time()))
cat("==============================================================\n\n")

# TEST 1 \u2014 Required output files
test_msg("Test 1: Required output files")
req_files <- c(
  "snv_processed_matrix.rds",
  "snv_maf_clean.rds",
  "sample_metadata.csv",
  "sample_metadata.rds",
  "qc_metrics.json",
  "plots/SNV_Oncoplot.png",
  "plots/SNV_Research_Validation_Figures.png"
)

for (f in req_files) {
  fp <- file.path(outdir, f)
  expect_true(file.exists(fp),
              sprintf("File exists: %s", f),
              sprintf("Missing file: %s", f))
}

if (!file.exists(file.path(outdir, "snv_processed_matrix.rds"))) {
  stop("Main matrix missing. Cannot continue tests.")
}

# Load files
test_msg("Loading output files")
snv_mat <- readRDS(file.path(outdir, "snv_processed_matrix.rds"))
meta_csv <- read.csv(file.path(outdir, "sample_metadata.csv"), stringsAsFactors = FALSE)
qc_json  <- fromJSON(file.path(outdir, "qc_metrics.json"))

cat(sprintf("  snv_processed_matrix.rds : %d x %d\n", nrow(snv_mat), ncol(snv_mat)))
cat(sprintf("  sample_metadata.csv      : %d rows\n", nrow(meta_csv)))

# TEST 2 \u2014 Matrix validity
test_msg("Test 2: Matrix Validity")
expect_true(is.matrix(snv_mat), "snv_matrix is a matrix", "snv_matrix is not a matrix")
expect_true(is.numeric(snv_mat), "snv_matrix is numeric", "snv_matrix is not numeric")
expect_true(!anyNA(snv_mat), "snv_matrix has no NAs", "snv_matrix contains NAs")

# TEST 3 \u2014 Binary values check
test_msg("Test 3: Binary Values Only")
all_binary <- all(snv_mat %in% c(0, 1))
expect_true(all_binary,
            "snv_matrix contains only 0s and 1s",
            "snv_matrix contains non-binary values")

# TEST 4 \u2014 Barcode Standard (12-char TCGA patient IDs)
test_msg("Test 4: Column names (12-char TCGA patient IDs)")
cols_12 <- nchar(colnames(snv_mat)) == 12
expect_true(all(cols_12),
            "All snv_matrix columns are 12 characters",
            "Some snv_matrix columns are NOT 12 characters")

starts_tcga <- grepl("^TCGA-", colnames(snv_mat))
expect_true(all(starts_tcga),
            "All snv_matrix columns start with 'TCGA-'",
            "Some snv_matrix columns do NOT start with 'TCGA-'")

# TEST 5 \u2014 Metadata Consistency
test_msg("Test 5: Metadata Consistency")
expect_true("patient_id" %in% colnames(meta_csv),
            "sample_metadata.csv has 'patient_id' column",
            "sample_metadata.csv is missing 'patient_id' column")

expect_true(nrow(meta_csv) == ncol(snv_mat),
            "Metadata rows match matrix columns",
            sprintf("Metadata rows (%d) != matrix cols (%d)", nrow(meta_csv), ncol(snv_mat)))

missing_in_meta <- setdiff(colnames(snv_mat), meta_csv$patient_id)
missing_in_mat  <- setdiff(meta_csv$patient_id, colnames(snv_mat))

expect_true(length(missing_in_meta) == 0,
            "All matrix columns present in metadata",
            sprintf("%d matrix columns missing from metadata", length(missing_in_meta)))
expect_true(length(missing_in_mat) == 0,
            "All metadata IDs present in matrix columns",
            sprintf("%d metadata IDs missing from matrix", length(missing_in_mat)))

# TEST 6 \u2014 Sparsity
test_msg("Test 6: Mutation Sparsity")
mean_mut_rate <- mean(snv_mat)
cat(sprintf("  Mean mutation rate: %.4f\n", mean_mut_rate))
expect_true(mean_mut_rate < 0.5,
            "Matrix is appropriately sparse for SNV (< 50%)",
            "Matrix is surprisingly dense (> 50%)")

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
  cat("\n  ALL TESTS PASSED. SNV module outputs are valid.\n")
  quit(status = 0)
}
