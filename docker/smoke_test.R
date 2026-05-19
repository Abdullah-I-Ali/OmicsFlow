#!/usr/bin/env Rscript
# =============================================================================
# docker/smoke_test.R — Container Validation Smoke Test
# =============================================================================
#
# PURPOSE:
#   Validates that ALL required R/Bioconductor packages are installed and
#   loadable inside the Docker container. Runs automatically during image build.
#   If any package fails to load, the Docker build will abort.
#
# =============================================================================

cat("\n")
cat("==============================================================\n")
cat("  OmicsFlow Docker Smoke Test\n")
cat(sprintf("  R version    : %s\n", R.version.string))
cat(sprintf("  Timestamp    : %s\n", Sys.time()))
cat("==============================================================\n\n")

n_pass <- 0
n_fail <- 0

check_pkg <- function(pkg) {
  result <- tryCatch({
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    TRUE
  }, error = function(e) FALSE)

  if (result) {
    ver <- as.character(packageVersion(pkg))
    cat(sprintf("  [\033[32mOK\033[0m]   %-40s %s\n", pkg, ver))
    n_pass <<- n_pass + 1
  } else {
    cat(sprintf("  [\033[31mFAIL\033[0m] %-40s MISSING\n", pkg))
    n_fail <<- n_fail + 1
  }
}

# --- CRAN Packages ---
cat("--- CRAN Packages ---\n")
cran_pkgs <- c(
  "optparse", "jsonlite", "data.table", "matrixStats", "stringr",
  "dplyr", "tidyr", "tibble", "scales", "RColorBrewer", "pheatmap",
  "ggplot2", "survival", "caret", "glmnet", "xgboost", "randomForestSRC"
)
for (p in cran_pkgs) check_pkg(p)

# --- Bioconductor: Core Genomics ---
cat("\n--- Bioconductor: Core Genomics ---\n")
bioc_core <- c(
  "SummarizedExperiment", "edgeR", "limma", "sva", "impute",
  "GenomicRanges", "IRanges", "GenomeInfoDb", "biomaRt"
)
for (p in bioc_core) check_pkg(p)

# --- Bioconductor: Methylation ---
cat("\n--- Bioconductor: Methylation ---\n")
bioc_meth <- c(
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"
)
for (p in bioc_meth) check_pkg(p)

# --- Bioconductor: Mutations ---
cat("\n--- Bioconductor: Mutations ---\n")
bioc_mut <- c("maftools")
for (p in bioc_mut) check_pkg(p)

# --- Bioconductor: Integration ---
cat("\n--- Bioconductor: Integration ---\n")
bioc_int <- c("MOFA2", "reticulate", "rhdf5")
for (p in bioc_int) check_pkg(p)

# --- Bioconductor: Enrichment ---
cat("\n--- Bioconductor: Enrichment ---\n")
bioc_enrich <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE", "AnnotationDbi")
for (p in bioc_enrich) check_pkg(p)

# --- Python Backend (mofapy2) ---
cat("\n--- Python Backend ---\n")
mofapy2_ok <- tryCatch({
  ret <- system2("python3", c("-c", "'import mofapy2; print(mofapy2.__version__)'"),
                  stdout = TRUE, stderr = TRUE)
  TRUE
}, error = function(e) FALSE)

if (mofapy2_ok) {
  cat("  [\033[32mOK\033[0m]   mofapy2 (Python)\n")
  n_pass <- n_pass + 1
} else {
  cat("  [\033[31mFAIL\033[0m] mofapy2 (Python)                         MISSING\n")
  n_fail <- n_fail + 1
}

# --- Summary ---
cat("\n==============================================================\n")
cat(sprintf("  RESULTS: %d packages checked\n", n_pass + n_fail))
cat(sprintf("    PASS : %d\n", n_pass))
cat(sprintf("    FAIL : %d\n", n_fail))
cat("==============================================================\n")

if (n_fail > 0) {
  cat("\n  \033[31mSMOKE TEST FAILED. Docker image is incomplete.\033[0m\n")
  quit(status = 1)
} else {
  cat("\n  \033[32mALL PACKAGES VERIFIED. Image is production-ready.\033[0m\n")
  quit(status = 0)
}
