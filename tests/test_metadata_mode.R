#!/usr/bin/env Rscript
# ==============================================================================
# test_metadata_mode.R — Integration Tests for Universal Metadata Layer
# OmicsFlow | Phase 1: Universal Metadata Layer
# ==============================================================================
# PURPOSE:
#   Validates that all preprocessing modules run successfully in metadata-driven
#   mode using custom, non-TCGA patient IDs and batch mappings, and verifies
#   scientific and structural correctness of outputs.
# ==============================================================================

set.seed(42)
cat("\n==============================================================\n")
cat("  OmicsFlow — Metadata-Driven Integration Tests\n")
cat("==============================================================\n\n")

# Create output directories
test_dir <- "results/test_metadata"
if (!dir.exists(test_dir)) dir.create(test_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. SYNTHESIZE MOCK METADATA FROM RAW DATA BARCODES
# ------------------------------------------------------------------------------
cat("--- Step 1: Synthesizing Mock Metadata CSV ---\n")

# Read some raw columns from RNA to get real barcodes
rna_raw <- readRDS("data/rna_expression_raw.rds")
if (inherits(rna_raw, "SummarizedExperiment")) {
  rna_barcodes <- colnames(assay(rna_raw))
} else {
  rna_barcodes <- colnames(rna_raw)
}

# Keep only primary tumors for synthesis to ensure enough matches
sample_types <- substr(rna_barcodes, 14, 15)
tumor_barcodes <- rna_barcodes[sample_types == "01"]

# Standardize to a subset of 30 samples for fast test runs
subset_barcodes <- head(tumor_barcodes, 30)
tcga_patients <- substr(subset_barcodes, 1, 12)

# Generate custom patient IDs and batches
unique_patients <- unique(tcga_patients)
custom_patient_map <- setNames(sprintf("PAT_%03d", seq_along(unique_patients)), unique_patients)

# Plate ID batches for RNA
tcga_plates <- substr(subset_barcodes, 22, 25)
unique_plates <- unique(tcga_plates)
custom_batch_map <- setNames(sprintf("BATCH_%d", seq_along(unique_plates) %% 2 + 1), unique_plates)

# 1. RNA metadata dataframe
rna_meta <- data.frame(
  sample_id    = subset_barcodes,
  patient_id   = custom_patient_map[tcga_patients],
  sample_class = "primary_tumor",
  batch        = custom_batch_map[tcga_plates],
  center       = "CENTER_A",
  stringsAsFactors = FALSE
)

# 2. DNA Methylation metadata dataframe
meth_raw <- readRDS("data/methylation_beta_raw.rds")
if (inherits(meth_raw, "SummarizedExperiment")) {
  meth_barcodes <- colnames(assay(meth_raw))
} else {
  meth_barcodes <- colnames(meth_raw)
}
meth_patients <- substr(meth_barcodes, 1, 12)
meth_keep_idx <- which(meth_patients %in% unique_patients)
meth_keep_barcodes <- meth_barcodes[meth_keep_idx]
meth_keep_patients <- meth_patients[meth_keep_idx]

meth_plates <- substr(meth_keep_barcodes, 22, 25)
unique_meth_plates <- unique(meth_plates)
meth_batch_map <- setNames(sprintf("BATCH_%d", seq_along(unique_meth_plates) %% 2 + 1), unique_meth_plates)

meth_meta <- data.frame(
  sample_id    = meth_keep_barcodes,
  patient_id   = custom_patient_map[meth_keep_patients],
  sample_class = "primary_tumor",
  batch        = meth_batch_map[meth_plates],
  center       = "CENTER_A",
  stringsAsFactors = FALSE
)

# 3. CNV metadata dataframe
cnv_raw <- readRDS("data/cnv_segment_raw.rds")
cnv_barcodes <- unique(unlist(strsplit(as.character(cnv_raw$Sample), ";")))
cnv_patients <- substr(cnv_barcodes, 1, 12)
cnv_keep_idx <- which(cnv_patients %in% unique_patients)
cnv_keep_barcodes <- cnv_barcodes[cnv_keep_idx]
cnv_keep_patients <- cnv_patients[cnv_keep_idx]

cnv_meta <- data.frame(
  sample_id    = cnv_keep_barcodes,
  patient_id   = custom_patient_map[cnv_keep_patients],
  sample_class = "primary_tumor",
  batch        = "BATCH_1",
  center       = "CENTER_A",
  stringsAsFactors = FALSE
)

# 4. SNV metadata dataframe
snv_raw <- readRDS("data/snv_mutation_raw.rds")
snv_barcodes <- unique(as.character(snv_raw$Tumor_Sample_Barcode))
snv_patients <- substr(snv_barcodes, 1, 12)
snv_keep_idx <- which(snv_patients %in% unique_patients)
snv_keep_barcodes <- snv_barcodes[snv_keep_idx]
snv_keep_patients <- snv_patients[snv_keep_idx]

snv_meta <- data.frame(
  sample_id    = snv_keep_barcodes,
  patient_id   = custom_patient_map[snv_keep_patients],
  sample_class = "primary_tumor",
  batch        = "BATCH_1",
  center       = "CENTER_A",
  stringsAsFactors = FALSE
)

# Combine and unique
meta_df <- rbind(rna_meta, meth_meta, cnv_meta, snv_meta)
meta_df <- unique(meta_df)

metadata_file <- file.path(test_dir, "sample_metadata.csv")
write.csv(meta_df, metadata_file, row.names = FALSE)
cat(sprintf("  Successfully synthesized metadata for %d samples across %d patients\n", 
            nrow(meta_df), length(unique(meta_df$patient_id))))
cat(sprintf("  Saved to: %s\n\n", metadata_file))


# ------------------------------------------------------------------------------
# 2. RUN RNA MODULE IN METADATA-DRIVEN MODE
# ------------------------------------------------------------------------------
cat("--- Step 2: Running RNA Preprocessing ---\n")
rna_out <- file.path(test_dir, "rna")
if (!dir.exists(rna_out)) dir.create(rna_out, recursive = TRUE)

rna_cmd <- sprintf(
  "Rscript modules/rna/preprocess_rna.R --input data/rna_expression_raw.rds --metadata %s --outdir %s --n-top 100 --cor-low 0.10 --cor-high 0.99",
  metadata_file, rna_out
)
cat(sprintf("  Executing: %s\n", rna_cmd))
system(rna_cmd)

# ------------------------------------------------------------------------------
# 3. RUN METHYLATION MODULE IN METADATA-DRIVEN MODE
# ------------------------------------------------------------------------------
cat("\n--- Step 3: Running DNA Methylation Preprocessing ---\n")
meth_out <- file.path(test_dir, "methylation")
if (!dir.exists(meth_out)) dir.create(meth_out, recursive = TRUE)

# Use lower knn-k and n-top to keep speed high
meth_cmd <- sprintf(
  "Rscript modules/methylation/preprocess_meth.R --input data/methylation_beta_raw.rds --metadata %s --outdir %s --n-top 100 --knn-k 2 --cross-react data/Chen_2013_cross_reactive_probes.csv",
  metadata_file, meth_out
)
cat(sprintf("  Executing: %s\n", meth_cmd))
system(meth_cmd)

# ------------------------------------------------------------------------------
# 4. RUN CNV MODULE IN METADATA-DRIVEN MODE
# ------------------------------------------------------------------------------
cat("\n--- Step 4: Running CNV Preprocessing ---\n")
cnv_out <- file.path(test_dir, "cnv")
if (!dir.exists(cnv_out)) dir.create(cnv_out, recursive = TRUE)

cnv_cmd <- sprintf(
  "Rscript modules/cnv/preprocess_cnv.R --input data/cnv_segment_raw.rds --cache data/gene_coords_hg38.rds --metadata %s --outdir %s --ntop 100",
  metadata_file, cnv_out
)
cat(sprintf("  Executing: %s\n", cnv_cmd))
system(cnv_cmd)

# ------------------------------------------------------------------------------
# 5. RUN SNV MODULE IN METADATA-DRIVEN MODE
# ------------------------------------------------------------------------------
cat("\n--- Step 5: Running SNV Preprocessing ---\n")
snv_out <- file.path(test_dir, "snv")
if (!dir.exists(snv_out)) dir.create(snv_out, recursive = TRUE)

snv_cmd <- sprintf(
  "Rscript modules/snv/preprocess_snv.R --input data/snv_mutation_raw.rds --metadata %s --outdir %s --freq 0.0 --topn 100",
  metadata_file, snv_out
)
cat(sprintf("  Executing: %s\n", snv_cmd))
system(snv_cmd)

# ------------------------------------------------------------------------------
# 6. RUN AND AGGREGATE VALIDATION TESTS
# ------------------------------------------------------------------------------
cat("\n--- Step 6: Running Validation Verification ---\n")

test_errors <- FALSE

# Check RNA outputs
cat("\n  [Verification] Validating RNA Metadata outputs...\n")
rna_scaled <- readRDS(file.path(rna_out, "rna_processed_matrix.rds"))
rna_meta <- read.csv(file.path(rna_out, "sample_metadata.csv"))

cat(sprintf("    RNA Scaled Matrix dimensions : %d genes x %d samples\n", nrow(rna_scaled), ncol(rna_scaled)))
cat(sprintf("    RNA metadata rows            : %d\n", nrow(rna_meta)))

if (!all(colnames(rna_scaled) %in% meta_df$patient_id)) {
  cat("    [FAIL] Column names in rna_processed_matrix are not mapped to custom patient IDs!\n")
  test_errors <- TRUE
} else {
  cat("    [PASS] Column names mapped successfully to custom patient IDs (e.g. PAT_001)\n")
}

if (!all(c("sample_id", "patient_id", "sample_class", "batch", "center") %in% colnames(rna_meta))) {
  cat("    [FAIL] RNA sample_metadata.csv is missing universal schema columns!\n")
  test_errors <- TRUE
} else {
  cat("    [PASS] RNA sample_metadata.csv has all standardized columns\n")
}

# Check Methylation outputs
cat("\n  [Verification] Validating Methylation Metadata outputs...\n")
meth_scaled <- readRDS(file.path(meth_out, "methylation_processed_matrix.rds"))
meth_meta <- read.csv(file.path(meth_out, "sample_metadata.csv"))

cat(sprintf("    Methylation Matrix dimensions: %d probes x %d samples\n", nrow(meth_scaled), ncol(meth_scaled)))
cat(sprintf("    Methylation metadata rows    : %d\n", nrow(meth_meta)))

if (!all(colnames(meth_scaled) %in% meta_df$patient_id)) {
  cat("    [FAIL] Column names in methylation_processed_matrix are not mapped to custom patient IDs!\n")
  test_errors <- TRUE
} else {
  cat("    [PASS] Column names mapped successfully to custom patient IDs\n")
}

if (!all(c("sample_id", "patient_id", "sample_class", "batch", "center") %in% colnames(meth_meta))) {
  cat("    [FAIL] Methylation sample_metadata.csv is missing universal schema columns!\n")
  test_errors <- TRUE
} else {
  cat("    [PASS] Methylation sample_metadata.csv has all standardized columns\n")
}

# Check CNV outputs
cat("\n  [Verification] Validating CNV Metadata outputs...\n")
cnv_scaled <- readRDS(file.path(cnv_out, "cnv_processed_matrix.rds"))
cnv_meta <- read.csv(file.path(cnv_out, "sample_metadata.csv"))

cat(sprintf("    CNV Matrix dimensions: %d genes x %d samples\n", nrow(cnv_scaled), ncol(cnv_scaled)))
if (!all(colnames(cnv_scaled) %in% meta_df$patient_id)) {
  cat("    [FAIL] Column names in cnv_processed_matrix are not mapped to custom patient IDs!\n")
  test_errors <- TRUE
} else {
  cat("    [PASS] Column names mapped successfully to custom patient IDs\n")
}

# Check SNV outputs
cat("\n  [Verification] Validating SNV Metadata outputs...\n")
snv_scaled <- readRDS(file.path(snv_out, "snv_processed_matrix.rds"))
snv_meta <- read.csv(file.path(snv_out, "sample_metadata.csv"))

cat(sprintf("    SNV Matrix dimensions: %d genes x %d samples\n", nrow(snv_scaled), ncol(snv_scaled)))
if (!all(colnames(snv_scaled) %in% meta_df$patient_id)) {
  cat("    [FAIL] Column names in snv_processed_matrix are not mapped to custom patient IDs!\n")
  test_errors <- TRUE
} else {
  cat("    [PASS] Column names mapped successfully to custom patient IDs\n")
}

cat("\n==============================================================\n")
if (test_errors) {
  cat("  METADATA INTEGRATION TESTS FAILED!\n")
  cat("==============================================================\n\n")
  quit(status = 1)
} else {
  cat("  ALL METADATA INTEGRATION TESTS PASSED SUCCESSFULLY!\n")
  cat("==============================================================\n\n")
  quit(status = 0)
}
