#!/usr/bin/env Rscript
# ==============================================================================
# test_rna.R — Validation Tests for the RNA Preprocessing Module
# OmicsFlow | Phase 1: RNA Module
# ==============================================================================
# PURPOSE:
#   Validates that all RNA module outputs are scientifically correct and
#   structurally complete. Compares outputs against expected properties
#   derived from the master script (data/full_scripts.R).
#
# USAGE:
#   Rscript tests/test_rna.R --outdir results/rna/
#
# TESTS:
#   1.  All required output files exist
#   2.  rna_processed_matrix.rds is a numeric matrix (no NAs)
#   3.  rna_ml.rds is a numeric matrix (no NAs)
#   4.  Column names are exactly 12-character TCGA patient barcodes
#   5.  Row names are gene symbols (no ENSG IDs remain)
#   6.  Z-scored matrix has mean ~ 0 and sd ~ 1 (per-gene)
#   7.  ML matrix is on log2 CPM scale (range check)
#   8.  No zero-variance genes in either output
#   9.  Both matrices have the same column set
#   10. sample_metadata.csv has required columns (patient_id, batch)
#   11. sample_metadata patient IDs match matrix columns
#   12. qc_metrics.json is valid JSON with all required keys
#   13. Dimensions: rna_processed_matrix cols <= rna_ml cols (Z-score may lose rows)
#   14. No duplicate row names in either matrix
#   15. No duplicate column names in either matrix
# ==============================================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--outdir", type = "character", default = "results/rna/",
              help = "RNA module output directory to validate [default: results/rna/]")
)
parser <- OptionParser(option_list = option_list,
                       description = "OmicsFlow RNA Module Validation Tests")
args   <- parse_args(parser)

# Load jsonlite for QC JSON validation
suppressPackageStartupMessages(library(jsonlite))

# ------------------------------------------------------------------------------
# TEST HARNESS
# ------------------------------------------------------------------------------
PASS <- 0L
FAIL <- 0L
WARN <- 0L

test_pass <- function(msg) {
  PASS <<- PASS + 1L
  cat(sprintf("  [PASS] %s\n", msg))
}

test_fail <- function(msg) {
  FAIL <<- FAIL + 1L
  cat(sprintf("  [FAIL] %s\n", msg))
}

test_warn <- function(msg) {
  WARN <<- WARN + 1L
  cat(sprintf("  [WARN] %s\n", msg))
}

expect_true <- function(condition, pass_msg, fail_msg) {
  if (isTRUE(condition)) test_pass(pass_msg) else test_fail(fail_msg)
}

outdir <- args$outdir
cat("\n")
cat("==============================================================\n")
cat("  OmicsFlow RNA Module — Validation Tests\n")
cat(sprintf("  Output directory: %s\n", outdir))
cat(sprintf("  Timestamp       : %s\n", Sys.time()))
cat("==============================================================\n\n")

# ==============================================================================
# TEST 1 — Required output files exist
# ==============================================================================
cat("--- Test 1: Required output files ---\n")

required_files <- c(
  "rna_processed_matrix.rds",
  "rna_ml.rds",
  "sample_metadata.csv",
  "sample_metadata.rds",
  "qc_metrics.json",
  file.path("plots", "pca_before_after_batch.pdf"),
  file.path("plots", "density_log2cpm.pdf"),
  file.path("plots", "sample_correlation_heatmap.pdf"),
  file.path("plots", "Research_Validation_Figures.png")
)

all_exist <- TRUE
for (f in required_files) {
  full_path <- file.path(outdir, f)
  if (file.exists(full_path)) {
    test_pass(sprintf("File exists: %s", f))
  } else {
    test_fail(sprintf("MISSING: %s", f))
    all_exist <- FALSE
  }
}

if (!all_exist) {
  cat("\nFATAL: One or more required output files are missing.\n")
  cat("Cannot proceed with further validation.\n\n")
  quit(status = 1)
}

# ==============================================================================
# LOAD OUTPUTS
# ==============================================================================
cat("\n--- Loading output files ---\n")
mofa_mat    <- readRDS(file.path(outdir, "rna_processed_matrix.rds"))
ml_mat      <- readRDS(file.path(outdir, "rna_ml.rds"))
meta_csv    <- read.csv(file.path(outdir, "sample_metadata.csv"),
                        stringsAsFactors = FALSE)
meta_rds    <- readRDS(file.path(outdir, "sample_metadata.rds"))
qc_json     <- jsonlite::read_json(file.path(outdir, "qc_metrics.json"))

cat(sprintf("  rna_processed_matrix.rds : %d x %d\n", nrow(mofa_mat), ncol(mofa_mat)))
cat(sprintf("  rna_ml.rds           : %d x %d\n", nrow(ml_mat),   ncol(ml_mat)))
cat(sprintf("  sample_metadata.csv  : %d rows\n",  nrow(meta_csv)))

# ==============================================================================
# TEST 2 — rna_processed_matrix.rds is a valid numeric matrix
# ==============================================================================
cat("\n--- Test 2: rna_processed_matrix.rds (Z-scored matrix) ---\n")

expect_true(is.matrix(mofa_mat),
            "rna_processed_matrix is a matrix",
            "rna_processed_matrix is NOT a matrix")
expect_true(is.numeric(mofa_mat),
            "rna_processed_matrix is numeric",
            "rna_processed_matrix is NOT numeric")
expect_true(!anyNA(mofa_mat),
            "rna_processed_matrix has no NAs",
            "rna_processed_matrix contains NAs")
expect_true(all(is.finite(mofa_mat)),
            "rna_processed_matrix has no Inf/NaN",
            "rna_processed_matrix contains Inf or NaN values")
expect_true(nrow(mofa_mat) > 0 && ncol(mofa_mat) > 0,
            "rna_processed_matrix has non-zero dimensions",
            "rna_processed_matrix has zero rows or columns")

# ==============================================================================
# TEST 3 — rna_ml.rds is a valid numeric matrix
# ==============================================================================
cat("\n--- Test 3: rna_ml.rds (log2 CPM matrix) ---\n")

expect_true(is.matrix(ml_mat),
            "rna_ml is a matrix",
            "rna_ml is NOT a matrix")
expect_true(is.numeric(ml_mat),
            "rna_ml is numeric",
            "rna_ml is NOT numeric")
expect_true(!anyNA(ml_mat),
            "rna_ml has no NAs",
            "rna_ml contains NAs")
expect_true(all(is.finite(ml_mat)),
            "rna_ml has no Inf/NaN",
            "rna_ml contains Inf or NaN values")

# ==============================================================================
# TEST 4 — Column names are 12-char TCGA barcodes
# ==============================================================================
cat("\n--- Test 4: Column names (12-char TCGA patient IDs) ---\n")

mofa_col_nchar <- all(nchar(colnames(mofa_mat)) == 12)
ml_col_nchar   <- all(nchar(colnames(ml_mat))   == 12)

expect_true(mofa_col_nchar,
            "All rna_processed_matrix columns are 12 characters",
            "Some rna_processed_matrix columns are NOT 12 characters")
expect_true(ml_col_nchar,
            "All rna_ml columns are 12 characters",
            "Some rna_ml columns are NOT 12 characters")

mofa_tcga <- all(grepl("^TCGA-", colnames(mofa_mat)))
ml_tcga   <- all(grepl("^TCGA-", colnames(ml_mat)))

expect_true(mofa_tcga,
            "All rna_processed_matrix columns start with 'TCGA-'",
            "Some rna_processed_matrix columns do NOT start with 'TCGA-'")
expect_true(ml_tcga,
            "All rna_ml columns start with 'TCGA-'",
            "Some rna_ml columns do NOT start with 'TCGA-'")

# ==============================================================================
# TEST 5 — Row names are gene symbols (no ENSG IDs)
# ==============================================================================
cat("\n--- Test 5: Row names are gene symbols (no Ensembl IDs) ---\n")

mofa_ensg <- sum(grepl("^ENSG", rownames(mofa_mat)))
ml_ensg   <- sum(grepl("^ENSG", rownames(ml_mat)))

expect_true(mofa_ensg == 0,
            "rna_processed_matrix has no ENSG row names",
            sprintf("rna_processed_matrix has %d ENSG row names (mapping failed)", mofa_ensg))
expect_true(ml_ensg == 0,
            "rna_ml has no ENSG row names",
            sprintf("rna_ml has %d ENSG row names (mapping failed)", ml_ensg))

# ==============================================================================
# TEST 6 — Z-scored matrix: per-gene mean ~ 0, sd ~ 1
# ==============================================================================
cat("\n--- Test 6: Z-score validation (mean ~ 0, sd ~ 1) ---\n")

gene_means_z <- rowMeans(mofa_mat)
gene_sds_z   <- apply(mofa_mat, 1, sd)

mean_of_means <- mean(gene_means_z)
mean_of_sds   <- mean(gene_sds_z)

cat(sprintf("  Per-gene mean — grand mean : %.6f (expected ~ 0)\n", mean_of_means))
cat(sprintf("  Per-gene sd   — grand mean : %.6f (expected ~ 1)\n", mean_of_sds))

expect_true(abs(mean_of_means) < 0.01,
            sprintf("Z-score grand mean is near 0 (%.6f)", mean_of_means),
            sprintf("Z-score grand mean is NOT near 0 (%.6f)", mean_of_means))
expect_true(abs(mean_of_sds - 1) < 0.05,
            sprintf("Z-score grand SD is near 1 (%.6f)", mean_of_sds),
            sprintf("Z-score grand SD is NOT near 1 (%.6f)", mean_of_sds))

# ==============================================================================
# TEST 7 — ML matrix is on log2 CPM scale (biologically plausible range)
# ==============================================================================
cat("\n--- Test 7: ML matrix expression range (log2 CPM scale) ---\n")

ml_min <- min(ml_mat)
ml_max <- max(ml_mat)
cat(sprintf("  Expression range: [%.2f, %.2f]\n", ml_min, ml_max))

expect_true(ml_min > -12,
            sprintf("ML matrix min is within plausible range (actual: %.2f, expected > -12)", ml_min),
            sprintf("ML matrix min is extremely low (%.2f), check normalization", ml_min))
expect_true(ml_max < 25,
            sprintf("ML matrix max is within plausible range (actual: %.2f, expected < 25)", ml_max),
            sprintf("ML matrix max is extremely high (%.2f), check normalization", ml_max))
# Since log2 CPM with prior.count=2 can naturally be negative for low-expression genes, we check if it's biologically plausible
expect_true(any(ml_mat < 0),
            "ML matrix contains negative values (biologically expected for log2 CPM with prior.count=2)",
            "ML matrix contains no negative values (unusual for log2 CPM with prior.count=2)")


# ==============================================================================
# TEST 8 — No zero-variance genes in either matrix
# ==============================================================================
cat("\n--- Test 8: Zero-variance genes ---\n")

mofa_zero_var <- sum(apply(mofa_mat, 1, var) < 1e-8)
ml_zero_var   <- sum(apply(ml_mat,   1, var) < 1e-8)

expect_true(mofa_zero_var == 0,
            "rna_processed_matrix has no zero-variance genes",
            sprintf("rna_processed_matrix has %d zero-variance genes", mofa_zero_var))
expect_true(ml_zero_var == 0,
            "rna_ml has no zero-variance genes",
            sprintf("rna_ml has %d zero-variance genes", ml_zero_var))

# ==============================================================================
# TEST 9 — Both matrices share the same sample set
# ==============================================================================
cat("\n--- Test 9: Sample consistency between matrices ---\n")

mofa_samples <- sort(colnames(mofa_mat))
ml_samples   <- sort(colnames(ml_mat))

expect_true(identical(mofa_samples, ml_samples),
            "rna_processed_matrix and rna_ml have identical sample sets",
            sprintf("Sample sets differ: %d in MOFA, %d in ML",
                    length(mofa_samples), length(ml_samples)))

# ==============================================================================
# TEST 10 — sample_metadata.csv has required columns
# ==============================================================================
cat("\n--- Test 10: sample_metadata.csv structure ---\n")

required_cols <- c("patient_id", "batch")
for (col in required_cols) {
  expect_true(col %in% colnames(meta_csv),
              sprintf("sample_metadata.csv has column '%s'", col),
              sprintf("sample_metadata.csv MISSING column '%s'", col))
}

expect_true(nrow(meta_csv) > 0,
            sprintf("sample_metadata.csv has %d rows", nrow(meta_csv)),
            "sample_metadata.csv is empty")

# Check CSV and RDS agree
expect_true(identical(meta_csv, meta_rds),
            "sample_metadata.csv and .rds are identical",
            "sample_metadata.csv and .rds differ")

# ==============================================================================
# TEST 11 — Metadata patient IDs match matrix columns
# ==============================================================================
cat("\n--- Test 11: Metadata patient IDs match matrix columns ---\n")

meta_ids   <- sort(meta_csv$patient_id)
matrix_ids <- sort(colnames(ml_mat))

n_in_meta_not_matrix <- sum(!meta_ids %in% matrix_ids)
n_in_matrix_not_meta <- sum(!matrix_ids %in% meta_ids)

expect_true(n_in_meta_not_matrix == 0,
            "All metadata IDs present in matrix columns",
            sprintf("%d metadata IDs NOT in matrix columns", n_in_meta_not_matrix))
expect_true(n_in_matrix_not_meta == 0,
            "All matrix columns present in metadata",
            sprintf("%d matrix columns NOT in metadata", n_in_matrix_not_meta))

# ==============================================================================
# TEST 12 — qc_metrics.json is valid and contains required keys
# ==============================================================================
cat("\n--- Test 12: qc_metrics.json structure ---\n")

required_qc_sections <- c("samples", "genes", "outliers", "batches",
                           "library_size", "output_matrices")
for (section in required_qc_sections) {
  expect_true(section %in% names(qc_json),
              sprintf("qc_metrics.json has section '%s'", section),
              sprintf("qc_metrics.json MISSING section '%s'", section))
}

# Check key fields within sections
if ("samples" %in% names(qc_json)) {
  expect_true("final" %in% names(qc_json$samples),
              "qc_metrics.json$samples$final is present",
              "qc_metrics.json$samples$final is MISSING")
}
if ("output_matrices" %in% names(qc_json)) {
  expect_true("zscore_mean" %in% names(qc_json$output_matrices),
              "qc_metrics.json$output_matrices$zscore_mean is present",
              "qc_metrics.json$output_matrices$zscore_mean is MISSING")
}

# ==============================================================================
# TEST 13 — Dimension consistency
# ==============================================================================
cat("\n--- Test 13: Matrix dimensions ---\n")

expect_true(nrow(mofa_mat) <= nrow(ml_mat),
            sprintf("rna_processed_matrix rows (%d) <= rna_ml rows (%d)",
                    nrow(mofa_mat), nrow(ml_mat)),
            sprintf("rna_processed_matrix rows (%d) > rna_ml rows (%d) — unexpected",
                    nrow(mofa_mat), nrow(ml_mat)))

if (nrow(mofa_mat) < 100) {
  test_warn(sprintf("rna_processed_matrix has only %d genes — verify feature selection",
                    nrow(mofa_mat)))
}
if (ncol(ml_mat) < 30) {
  test_warn(sprintf("rna_ml has only %d samples — small cohort", ncol(ml_mat)))
}

# ==============================================================================
# TEST 14 — No duplicate row names
# ==============================================================================
cat("\n--- Test 14: Unique gene names ---\n")

mofa_dup_rows <- sum(duplicated(rownames(mofa_mat)))
ml_dup_rows   <- sum(duplicated(rownames(ml_mat)))

expect_true(mofa_dup_rows == 0,
            "rna_processed_matrix has no duplicate gene names",
            sprintf("rna_processed_matrix has %d duplicate gene names", mofa_dup_rows))
expect_true(ml_dup_rows == 0,
            "rna_ml has no duplicate gene names",
            sprintf("rna_ml has %d duplicate gene names", ml_dup_rows))

# ==============================================================================
# TEST 15 — No duplicate column names
# ==============================================================================
cat("\n--- Test 15: Unique sample IDs ---\n")

mofa_dup_cols <- sum(duplicated(colnames(mofa_mat)))
ml_dup_cols   <- sum(duplicated(colnames(ml_mat)))

expect_true(mofa_dup_cols == 0,
            "rna_processed_matrix has no duplicate sample IDs",
            sprintf("rna_processed_matrix has %d duplicate sample IDs", mofa_dup_cols))
expect_true(ml_dup_cols == 0,
            "rna_ml has no duplicate sample IDs",
            sprintf("rna_ml has %d duplicate sample IDs", ml_dup_cols))

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
total <- PASS + FAIL + WARN
cat("\n")
cat("==============================================================\n")
cat(sprintf("  RESULTS: %d tests total\n", total))
cat(sprintf("    PASS : %d\n", PASS))
cat(sprintf("    FAIL : %d\n", FAIL))
cat(sprintf("    WARN : %d\n", WARN))
cat("==============================================================\n\n")

if (FAIL > 0) {
  cat(sprintf("VALIDATION FAILED: %d test(s) failed.\n\n", FAIL))
  quit(status = 1)
} else if (WARN > 0) {
  cat(sprintf("VALIDATION PASSED with %d warning(s).\n\n", WARN))
  quit(status = 0)
} else {
  cat("ALL TESTS PASSED. RNA module outputs are valid.\n\n")
  quit(status = 0)
}
