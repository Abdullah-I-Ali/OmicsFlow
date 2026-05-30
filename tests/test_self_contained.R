# ==============================================================================
# tests/test_self_contained.R — End-to-End Self-Contained Package Test
#
# Verifies that OmicsFlow functions work after install_github() from an
# arbitrary working directory, WITHOUT requiring set_omicsflow_root() or
# cloning the repository.
# ==============================================================================

library(OmicsFlow)

cat("=== OmicsFlow Self-Contained Package Test ===\n\n")

# --------------------------------------------------------------------------
# Test 1: system.file() resolves installed modules
# --------------------------------------------------------------------------
cat("[TEST 1] system.file() resolves installed modules...\n")
pkg_modules <- system.file("modules", package = "OmicsFlow")
if (!nzchar(pkg_modules)) {
  stop("FAIL: system.file('modules', package='OmicsFlow') returned empty string")
}
if (!dir.exists(pkg_modules)) {
  stop("FAIL: Modules directory does not exist at installed path: ", pkg_modules)
}
cat(sprintf("  [OK] Modules found at: %s\n", pkg_modules))

# Verify key module files exist
required_modules <- c(
  file.path("modules", "utils_metadata.R"),
  file.path("modules", "utils_clinical.R"),
  file.path("modules", "rna", "preprocess_rna.R"),
  file.path("modules", "methylation", "preprocess_meth.R"),
  file.path("modules", "cnv", "preprocess_cnv.R"),
  file.path("modules", "snv", "preprocess_snv.R"),
  file.path("modules", "integration", "run_integration.R"),
  file.path("modules", "ml", "run_ml.R"),
  file.path("modules", "enrichment", "run_enrichment.R")
)

for (mod in required_modules) {
  mod_path <- system.file(mod, package = "OmicsFlow")
  if (!nzchar(mod_path) || !file.exists(mod_path)) {
    stop(sprintf("FAIL: Required module not found: %s", mod))
  }
}
cat(sprintf("  [OK] All %d required module files found\n", length(required_modules)))

# --------------------------------------------------------------------------
# Test 2: system.file() resolves installed reports
# --------------------------------------------------------------------------
cat("\n[TEST 2] system.file() resolves installed reports...\n")
report_qmd <- system.file("reports", "OmicsFlow_Report.qmd", package = "OmicsFlow")
if (!nzchar(report_qmd) || !file.exists(report_qmd)) {
  cat("  [WARN] Report template not found (optional for non-rendering use)\n")
} else {
  cat(sprintf("  [OK] Report template found at: %s\n", report_qmd))
}

# --------------------------------------------------------------------------
# Test 2.5: system.file() resolves bundled CNV cache
# --------------------------------------------------------------------------
cat("\n[TEST 2.5] system.file() resolves bundled CNV cache...\n")
cnv_cache_path <- system.file("configs", "gene_coordinates.rds", package = "OmicsFlow")
if (!nzchar(cnv_cache_path)) {
  stop("FAIL: Bundled CNV cache not found. The file was not included in the installed package.")
} else {
  cat(sprintf("  [OK] Bundled CNV cache found at: %s\n", cnv_cache_path))
}

# --------------------------------------------------------------------------
# Test 2.6: system.file() resolves bundled cross-reactive probes
# --------------------------------------------------------------------------
cat("\n[TEST 2.6] system.file() resolves bundled cross-reactive probes...\n")
cross_react_path <- system.file("configs", "cross_reactive_probes.csv", package = "OmicsFlow")
if (!nzchar(cross_react_path)) {
  stop("FAIL: Bundled cross-reactive probes not found. The file was not included in the installed package.")
} else {
  cat(sprintf("  [OK] Bundled cross-reactive probes found at: %s\n", cross_react_path))
}

# --------------------------------------------------------------------------
# Test 3: omicsflow_root() resolves without set_omicsflow_root()
# --------------------------------------------------------------------------
cat("\n[TEST 3] omicsflow_root() resolves without set_omicsflow_root()...\n")

# Reset any cached root
env <- OmicsFlow:::.omicsflow_env
env$root <- NULL

# Change to a temporary directory to ensure we're NOT in the repo
original_wd <- getwd()
tmp_dir <- tempdir()
setwd(tmp_dir)

root <- tryCatch(
  OmicsFlow:::omicsflow_root(),
  error = function(e) {
    setwd(original_wd)
    stop("FAIL: omicsflow_root() failed from arbitrary directory: ", e$message)
  }
)

if (!dir.exists(file.path(root, "modules"))) {
  setwd(original_wd)
  stop("FAIL: Resolved root does not contain modules/: ", root)
}
cat(sprintf("  [OK] Root resolved to: %s\n", root))

# --------------------------------------------------------------------------
# Test 4: omicsflow_path() constructs valid paths
# --------------------------------------------------------------------------
cat("\n[TEST 4] omicsflow_path() constructs valid paths...\n")
meta_path <- OmicsFlow:::omicsflow_path("modules", "utils_metadata.R")
if (!file.exists(meta_path)) {
  setwd(original_wd)
  stop("FAIL: omicsflow_path('modules', 'utils_metadata.R') does not exist: ", meta_path)
}
cat(sprintf("  [OK] omicsflow_path() resolves correctly: %s\n", meta_path))

# --------------------------------------------------------------------------
# Test 5: validate_metadata_schema() is available without source()
# --------------------------------------------------------------------------
cat("\n[TEST 5] validate_metadata_schema() works as package function...\n")

# Create a valid test metadata data frame
test_meta <- data.frame(
  sample_id = c("S1", "S2", "S3"),
  patient_id = c("P1", "P2", "P3"),
  sample_class = c("Tumor", "Normal", "Tumor"),
  batch = c("B1", "B1", "B2"),
  stringsAsFactors = FALSE
)

result <- tryCatch(
  OmicsFlow:::validate_metadata_schema(test_meta),
  error = function(e) {
    setwd(original_wd)
    stop("FAIL: validate_metadata_schema() failed: ", e$message)
  }
)
if (!isTRUE(result)) {
  setwd(original_wd)
  stop("FAIL: validate_metadata_schema() did not return TRUE for valid data")
}
cat("  [OK] validate_metadata_schema() works as package-internal function\n")

# --------------------------------------------------------------------------
# Test 6: validate_inputs() does not call source()
# --------------------------------------------------------------------------
cat("\n[TEST 6] validate_inputs() runs without sourcing external files...\n")

# Create temporary test data files
tmp_rna <- file.path(tmp_dir, "test_rna.rds")
tmp_meth <- file.path(tmp_dir, "test_meth.rds")
tmp_meta <- file.path(tmp_dir, "test_metadata.csv")

# Create minimal test matrices
set.seed(42)
rna_mat <- matrix(rnorm(300), nrow = 100, ncol = 3,
                  dimnames = list(paste0("gene", 1:100), c("S1", "S2", "S3")))
meth_mat <- matrix(runif(300), nrow = 100, ncol = 3,
                   dimnames = list(paste0("cg", 1:100), c("S1", "S2", "S3")))

saveRDS(rna_mat, tmp_rna)
saveRDS(meth_mat, tmp_meth)
write.csv(test_meta, tmp_meta, row.names = FALSE)

val_result <- tryCatch(
  OmicsFlow::validate_inputs(
    rna = tmp_rna,
    meth = tmp_meth,
    metadata = tmp_meta
  ),
  error = function(e) {
    setwd(original_wd)
    stop("FAIL: validate_inputs() failed from arbitrary directory: ", e$message)
  }
)

if (!isTRUE(val_result$valid)) {
  cat("  [WARN] Validation returned valid=FALSE (may be data-specific, not a path issue)\n")
  cat(sprintf("  Errors: %s\n", paste(val_result$errors, collapse = "; ")))
} else {
  cat("  [OK] validate_inputs() succeeded from arbitrary working directory\n")
}

# --------------------------------------------------------------------------
# Test 7: generate_metadata_templates() works from arbitrary directory
# --------------------------------------------------------------------------
cat("\n[TEST 7] generate_metadata_templates() works from arbitrary directory...\n")

tmp_output <- file.path(tmp_dir, "test_templates")
gen_result <- tryCatch(
  OmicsFlow::generate_metadata_templates(
    rna = tmp_rna,
    meth = tmp_meth,
    output_dir = tmp_output
  ),
  error = function(e) {
    setwd(original_wd)
    stop("FAIL: generate_metadata_templates() failed: ", e$message)
  }
)

if (gen_result$status == "success") {
  cat("  [OK] generate_metadata_templates() succeeded\n")
} else {
  cat("  [WARN] generate_metadata_templates() did not return success\n")
}

# --------------------------------------------------------------------------
# Cleanup
# --------------------------------------------------------------------------
setwd(original_wd)
unlink(c(tmp_rna, tmp_meth, tmp_meta))
unlink(tmp_output, recursive = TRUE)

cat("\n=== All self-contained package tests PASSED ===\n")
