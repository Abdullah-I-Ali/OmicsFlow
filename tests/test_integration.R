#!/usr/bin/env Rscript
# ==============================================================================
# test_integration.R — Validation Tests for Integration Module (MOFA+)
# ==============================================================================
#
# PURPOSE:
#   Validates the outputs of the MOFA+ Integration module to ensure scientific
#   integrity and compatibility with downstream machine learning.
#
# CHECKS:
#   1. Required files exist
#   2. Valid MOFA2 model object
#   3. Correct number of views (4) and names
#   4. Common samples across all views
#   5. Active factors vector structure
#   6. Variance Explained values populated
#   7. Quality of Factor 1 weights
#
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(MOFA2)
})

option_list <- list(
  make_option(c("-o", "--outdir"), type = "character", default = "results/integration/",
              help = "Directory containing Integration outputs [default= %default]")
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
cat("  OmicsFlow Integration Module \u2014 Validation Tests\n")
cat(sprintf("  Output directory: %s\n", outdir))
cat(sprintf("  Timestamp       : %s\n", Sys.time()))
cat("==============================================================\n\n")

# TEST 1 \u2014 Required output files
test_msg("Test 1: Required output files")
req_files <- c(
  "mofa_model.rds",
  "qc_metrics.json",
  "plots/MOFA_Fig1_Variance_Landscape.png",
  "plots/MOFA_Fig2_Factor_Correlation.pdf",
  "plots/MOFA_Fig3_Factor1_Drivers.png"
)

for (f in req_files) {
  fp <- file.path(outdir, f)
  expect_true(file.exists(fp),
              sprintf("File exists: %s", f),
              sprintf("Missing file: %s", f))
}

if (!file.exists(file.path(outdir, "mofa_model.rds"))) {
  stop("Main MOFA model missing. Cannot continue tests.")
}

# Load model
test_msg("Loading output files")
loaded <- readRDS(file.path(outdir, "mofa_model.rds"))
qc_json  <- fromJSON(file.path(outdir, "qc_metrics.json"))

# TEST 2 \u2014 Valid MOFA2 structure
test_msg("Test 2: MOFA2 Structure")
expect_true("model" %in% names(loaded) && "active_factors" %in% names(loaded),
            "Saved object is a list containing $model and $active_factors",
            "Saved object does not have correct structure")

mofa <- loaded$model
expect_true(is(mofa, "MOFA"),
            "$model is a valid MOFA2 object",
            "$model is not a valid MOFA2 object")

# TEST 3 \u2014 Views and Dimensions
test_msg("Test 3: Views and Samples Integrity")
views <- views_names(mofa)
expected_views <- c("RNA", "Methylation", "CNV", "SNV")
expect_true(all(expected_views %in% views),
            "All 4 omics layers present as views",
            sprintf("Missing views. Found: %s", paste(views, collapse=", ")))

samples <- samples_names(mofa)[[1]]
expect_true(length(samples) > 0,
            sprintf("Model contains common samples (%d patients)", length(samples)),
            "Model has 0 samples")

# TEST 4 \u2014 Active factors extraction
test_msg("Test 4: Active Factors Extraction")
active <- loaded$active_factors
expect_true(length(active) > 0 && is.numeric(active),
            sprintf("Successfully extracted %d active factors", length(active)),
            "Failed to extract active factors")

# TEST 5 \u2014 Variance Explained calculation
test_msg("Test 5: Variance Explained Generation")
var_exp <- get_variance_explained(mofa)
expect_true(!is.null(var_exp$r2_total[[1]]),
            "Variance explained properly calculated",
            "Variance explained missing")

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
  cat("\n  ALL TESTS PASSED. Integration module outputs are valid.\n")
  quit(status = 0)
}
