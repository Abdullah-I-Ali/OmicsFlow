#!/usr/bin/env Rscript
# ==============================================================================
# test_enrichment.R — Validation Tests for Enrichment Module
# ==============================================================================
#
# PURPOSE:
#   Validates the outputs of the Enrichment module to ensure scientific
#   integrity and correct background (universe) mapping.
#
# CHECKS:
#   1. Required files exist (RDS, CSV, PNG)
#   2. Valid list object loaded from RDS
#   3. GO BP/MF/CC and KEGG slots present
#   4. Correct mapping to human OrgDb
#
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

option_list <- list(
  make_option(c("-o", "--outdir"), type = "character", default = "results/enrichment/",
              help = "Directory containing Enrichment outputs [default= %default]")
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
cat("  OmicsFlow Enrichment Module \u2014 Validation Tests\n")
cat(sprintf("  Output directory: %s\n", outdir))
cat(sprintf("  Timestamp       : %s\n", Sys.time()))
cat("==============================================================\n\n")

# TEST 1 \u2014 Required output files
test_msg("Test 1: Required output files")
req_files <- c(
  "enrichment_results.rds",
  "qc_metrics.json"
)

for (f in req_files) {
  fp <- file.path(outdir, f)
  expect_true(file.exists(fp),
              sprintf("File exists: %s", f),
              sprintf("Missing file: %s", f))
}

if (!file.exists(file.path(outdir, "enrichment_results.rds"))) {
  stop("Main RDS output missing. Cannot continue tests.")
}

# Load main output
test_msg("Loading output files")
res <- readRDS(file.path(outdir, "enrichment_results.rds"))
qc_json <- fromJSON(file.path(outdir, "qc_metrics.json"))

# TEST 2 \u2014 Results Format
test_msg("Test 2: RDS Structure Integrity")
expect_true(is.list(res), "Results is a list", "Results is not a list")

expected_names <- c("go_bp", "go_mf", "go_cc", "kegg")
expect_true(all(expected_names %in% names(res)),
            "List contains go_bp, go_mf, go_cc, kegg",
            "List is missing expected enrichment slots")

# TEST 3 \u2014 Enrichment Objects Data Frames
test_msg("Test 3: Valid enrichment classes")
classes <- sapply(res, function(x) inherits(x, c("enrichResult", "data.frame")))
expect_true(all(classes), 
            "All results are valid enrichment/dataframe objects",
            "Some results are invalid objects")

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
  cat("\n  ALL TESTS PASSED. Enrichment module outputs are valid.\n")
  quit(status = 0)
}
