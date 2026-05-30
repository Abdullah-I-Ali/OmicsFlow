#!/usr/bin/env Rscript
# ==============================================================================
# tests/test_clinical_ingestion.R — Unit Tests for Dataset-Agnostic Clinical Ingestion
# ==============================================================================

library(OmicsFlow)

cat("=== OmicsFlow Clinical Ingestion Unit Tests ===\n\n")

# Create a temporary directory for test files
tmp_dir <- tempdir()
tmp_clinical_tcga <- file.path(tmp_dir, "clinical_tcga.tsv")
tmp_clinical_geo <- file.path(tmp_dir, "clinical_geo.csv")
tmp_clinical_custom <- file.path(tmp_dir, "clinical_custom.tsv")
tmp_metadata <- file.path(tmp_dir, "metadata.csv")

# ------------------------------------------------------------------------------
# Mock Data Preparation
# ------------------------------------------------------------------------------

# 1. TCGA-like Clinical File (Dual-column survival time, age in days, barcode)
df_tcga <- data.frame(
  bcr_patient_barcode = c("TCGA-01", "TCGA-02", "TCGA-03"),
  days_to_death = c(NA, 450, NA),
  days_to_last_follow_up = c(1200, NA, 850),
  vital_status = c("Alive", "Dead", "Alive"),
  age_at_diagnosis = c(23450, 18900, 21200),
  gender = c("male", "female", "female"),
  stringsAsFactors = FALSE
)
write.table(df_tcga, tmp_clinical_tcga, sep = "\t", row.names = FALSE, quote = FALSE)

# 2. GEO/METABRIC-like Clinical File (Single survival column, age in years)
df_geo <- data.frame(
  patient = c("GEO-01", "GEO-02", "GEO-03"),
  survival_months = c(24.5, 12.1, 45.0),
  status = c("alive", "deceased", "alive"),
  age_years = c(64, 51, 72),
  sex = c("M", "F", "F"),
  stringsAsFactors = FALSE
)
write.csv(df_geo, tmp_clinical_geo, row.names = FALSE)

# 3. Custom Institutional Clinical File (Custom headers)
df_custom <- data.frame(
  Subject_ID = c("PAT-1", "PAT-2", "PAT-3"),
  Time_Days = c(300, 150, 800),
  Event_Occurred = c(0, 1, 0),
  Age_Years = c(45, 62, 55),
  Sex = c("Male", "Female", "Male"),
  stringsAsFactors = FALSE
)
write.table(df_custom, tmp_clinical_custom, sep = "\t", row.names = FALSE, quote = FALSE)

# 4. Sample Metadata File
df_meta <- data.frame(
  sample_id = c("S1", "S2", "S3"),
  patient_id = c("TCGA-01", "TCGA-02", "TCGA-03"),
  sample_class = c("primary_tumor", "primary_tumor", "primary_tumor"),
  batch = c("B1", "B1", "B1"),
  stringsAsFactors = FALSE
)
write.csv(df_meta, tmp_metadata, row.names = FALSE)

# ------------------------------------------------------------------------------
# Test 1: Auto-detection on TCGA-like cohort (Dual-column Survival)
# ------------------------------------------------------------------------------
cat("[TEST 1] Auto-detecting TCGA dual-column survival columns...\n")
det_tcga <- OmicsFlow:::detect_clinical_columns(tmp_clinical_tcga)
cmap_tcga <- det_tcga$mapping

if (cmap_tcga$patient_id != "bcr_patient_barcode") stop("Fail: patient_id not mapped to bcr_patient_barcode")
if (cmap_tcga$os_time != "days_to_death,days_to_last_follow_up") stop("Fail: os_time dual column detection failed. Got: ", cmap_tcga$os_time)
if (cmap_tcga$os_event != "vital_status") stop("Fail: os_event not mapped to vital_status")
if (cmap_tcga$age != "age_at_diagnosis") stop("Fail: age not mapped to age_at_diagnosis")
if (cmap_tcga$gender != "gender") stop("Fail: gender not mapped to gender")
cat("  [OK] Successfully auto-detected TCGA schema, including dual-column 'days_to_death,days_to_last_follow_up'\n")

# ------------------------------------------------------------------------------
# Test 2: Auto-detection on GEO/METABRIC-like cohort
# ------------------------------------------------------------------------------
cat("\n[TEST 2] Auto-detecting GEO/METABRIC-like columns...\n")
det_geo <- OmicsFlow:::detect_clinical_columns(tmp_clinical_geo)
cmap_geo <- det_geo$mapping

if (cmap_geo$patient_id != "patient") stop("Fail: patient_id not mapped to patient")
if (cmap_geo$os_time != "survival_months") stop("Fail: os_time not mapped to survival_months")
if (cmap_geo$os_event != "status") stop("Fail: os_event not mapped to status")
if (cmap_geo$age != "age_years") stop("Fail: age not mapped to age_years")
if (cmap_geo$gender != "sex") stop("Fail: gender not mapped to sex")
cat("  [OK] Successfully auto-detected GEO/METABRIC-like schema\n")

# ------------------------------------------------------------------------------
# Test 3: Standarized dual-column survival ingestion
# ------------------------------------------------------------------------------
cat("\n[TEST 3] Loading and standardizing dual-column survival clinical data...\n")
# Simulating load_clinical_data call using our auto-detected mapping
std_tcga <- OmicsFlow:::load_clinical_data(tmp_clinical_tcga, cmap_tcga, df_meta)

if (nrow(std_tcga) != 3) stop("Fail: incorrect number of parsed rows")
if (!all(std_tcga$patient_id == c("TCGA-01", "TCGA-02", "TCGA-03"))) stop("Fail: patient_id mismatch")
if (!all(is.na(std_tcga$os_time) == FALSE)) stop("Fail: os_time should not be NA")
# dead has 450, alive has 1200 and 850
if (!all(std_tcga$os_time == c(1200, 450, 850))) stop("Fail: dual-column os_time values resolved incorrectly")
if (!all(std_tcga$os_event == c(0, 1, 0))) stop("Fail: vital_status event conversion failed")
# age in days should be divided by 365.25
if (abs(std_tcga$age[1] - (23450/365.25)) > 0.01) stop("Fail: age conversion from days to years failed")
cat("  [OK] Successfully loaded and standardized TCGA clinical data with dual columns and age-in-days conversion\n")

# ------------------------------------------------------------------------------
# Test 4: Pre-flight validation with comma-separated mappings
# ------------------------------------------------------------------------------
cat("\n[TEST 4] Validating pre-flight input validation with comma-separated mappings...\n")
# Create mock RNA and methylation files
tmp_rna <- file.path(tmp_dir, "test_rna.rds")
tmp_meth <- file.path(tmp_dir, "test_meth.rds")
saveRDS(matrix(rnorm(300), nrow = 100, ncol = 3, dimnames = list(paste0("g", 1:100), c("S1", "S2", "S3"))), tmp_rna)
saveRDS(matrix(runif(300), nrow = 100, ncol = 3, dimnames = list(paste0("cg", 1:100), c("S1", "S2", "S3"))), tmp_meth)

# Write a clinical map JSON referencing dual columns
tmp_map_json <- file.path(tmp_dir, "clinical_map.json")
write(jsonlite::toJSON(cmap_tcga, auto_unbox = TRUE, pretty = TRUE), tmp_map_json)

val_res <- validate_inputs(
  rna = tmp_rna,
  meth = tmp_meth,
  metadata = tmp_metadata,
  clinical = tmp_clinical_tcga,
  clinical_map = tmp_map_json
)

if (!isTRUE(val_res$valid)) {
  stop("Fail: validation failed for valid dual-column clinical mapping: ", paste(val_res$errors, collapse = "; "))
}
cat("  [OK] Pre-flight validation succeeded for multi-column mapping and patient overlap verification\n")

# ------------------------------------------------------------------------------
# Test 5: Patient overlap verification failure trigger
# ------------------------------------------------------------------------------
cat("\n[TEST 5] Verifying that validation fails if patient overlap is zero...\n")
# Write metadata with non-overlapping patient IDs (e.g. PAT-X)
df_meta_bad <- data.frame(
  sample_id = c("S1", "S2", "S3"),
  patient_id = c("PAT-99", "PAT-98", "PAT-97"),
  sample_class = c("primary_tumor", "primary_tumor", "primary_tumor"),
  batch = c("B1", "B1", "B1"),
  stringsAsFactors = FALSE
)
tmp_metadata_bad <- file.path(tmp_dir, "metadata_bad.csv")
write.csv(df_meta_bad, tmp_metadata_bad, row.names = FALSE)

val_res_bad <- validate_inputs(
  rna = tmp_rna,
  meth = tmp_meth,
  metadata = tmp_metadata_bad,
  clinical = tmp_clinical_tcga,
  clinical_map = tmp_map_json
)

if (isTRUE(val_res_bad$valid)) {
  stop("Fail: validation should have failed due to zero patient overlap between metadata and clinical files")
}

# Verify that the expected error message is present
has_error <- any(grepl("No overlapping patient IDs found between sample metadata and clinical data", val_res_bad$errors))
if (!has_error) {
  stop("Fail: could not find expected patient overlap error message in validation errors")
}
cat("  [OK] Successfully failed validation with clear warning when patient overlap is 0\n")

# ------------------------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------------------------
unlink(c(tmp_clinical_tcga, tmp_clinical_geo, tmp_clinical_custom, tmp_metadata, tmp_metadata_bad, tmp_rna, tmp_meth, tmp_map_json))
cat("\n=== All clinical ingestion unit tests PASSED successfully! ===\n")
