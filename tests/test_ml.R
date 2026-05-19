#!/usr/bin/env Rscript
# ==============================================================================
# test_ml.R — Validation Tests for Machine Learning Module
# ==============================================================================
#
# PURPOSE:
#   Validates the outputs of the ML module to ensure scientific
#   integrity and successful prediction models.
#
# CHECKS:
#   1. Required files exist
#   2. Models successfully created
#   3. Results summary integrity
#   4. Extracted features logic
#
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(randomForestSRC)
  library(xgboost)
  library(glmnet)
})

option_list <- list(
  make_option(c("-o", "--outdir"), type = "character", default = "results/ml/",
              help = "Directory containing ML outputs [default= %default]")
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
cat("  OmicsFlow ML Module \u2014 Validation Tests\n")
cat(sprintf("  Output directory: %s\n", outdir))
cat(sprintf("  Timestamp       : %s\n", Sys.time()))
cat("==============================================================\n\n")

# TEST 1 \u2014 Required output files
test_msg("Test 1: Required output files")
req_files <- c(
  "ml_results_summary.csv",
  "ml_top_features.csv",
  "rf_survival_model.rds",
  "xgb_cox_model.rds",
  "lasso_cox_model.rds",
  "qc_metrics.json",
  "plots/ML_Research_Validation_Figures.png",
  "plots/kaplan_meier_risk_groups.pdf"
)

for (f in req_files) {
  fp <- file.path(outdir, f)
  expect_true(file.exists(fp),
              sprintf("File exists: %s", f),
              sprintf("Missing file: %s", f))
}

# Load main output
test_msg("Loading output files")
results <- read.csv(file.path(outdir, "ml_results_summary.csv"), stringsAsFactors = FALSE)
qc_json <- fromJSON(file.path(outdir, "qc_metrics.json"))

# TEST 2 \u2014 Results Format
test_msg("Test 2: Results Summary Integrity")
expect_true(all(c("model", "c_index", "top_5_genes") %in% colnames(results)),
            "ml_results_summary.csv has correct columns",
            "ml_results_summary.csv has missing/incorrect columns")

expect_true(nrow(results) == 3,
            "Summary contains 3 models (RF, XGBoost, LASSO)",
            sprintf("Summary contains %d rows, expected 3", nrow(results)))

# TEST 3 \u2014 Models loadable
test_msg("Test 3: Saved Models Structure")
rf <- readRDS(file.path(outdir, "rf_survival_model.rds"))
xgb <- readRDS(file.path(outdir, "xgb_cox_model.rds"))
lasso <- readRDS(file.path(outdir, "lasso_cox_model.rds"))

expect_true(inherits(rf, "rfsrc"), "RF model is valid class (rfsrc)", "RF model is invalid")
expect_true(inherits(xgb, "xgb.Booster"), "XGBoost model is valid class (xgb.Booster)", "XGBoost model is invalid")
expect_true(inherits(lasso, "cv.glmnet"), "LASSO model is valid class (cv.glmnet)", "LASSO model is invalid")

# TEST 4 \u2014 Quality Check on Performance
test_msg("Test 4: Performance Validation")
expect_true(all(results$c_index >= 0 & results$c_index <= 1),
            "All C-index values are between 0 and 1",
            "Some C-index values are outside valid [0,1] range")

# TEST 5 \u2014 Top features
test_msg("Test 5: Top Features Metadata")
features <- read.csv(file.path(outdir, "ml_top_features.csv"))
expect_true("train_variance" %in% colnames(features),
            "Features table includes 'train_variance' from Hybrid selection",
            "Features table missing 'train_variance'")
expect_true(nrow(features) > 0,
            "Features table is populated",
            "Features table is empty")

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
  cat("\n  ALL TESTS PASSED. ML module outputs are valid.\n")
  quit(status = 0)
}
