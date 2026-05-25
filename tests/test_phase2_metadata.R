#!/usr/bin/env Rscript
# ==============================================================================
# test_phase2_metadata.R — End-to-End Validation Suite for Phase 2
# OmicsFlow | Phase 2: Clinical Metadata Generalization
# ==============================================================================
# PURPOSE:
#   Validates the entire OmicsFlow pipeline in metadata-driven and generalized
#   clinical mode. Runs preprocessing, multi-omics integration, ML survival
#   modeling (with custom mapping), pathway enrichment (with custom keywords),
#   and report generation.
# ==============================================================================

set.seed(42)
cat("\n==============================================================\n")
cat("  OmicsFlow — Phase 2 Clinical & Metadata End-to-End Tests\n")
cat("==============================================================\n\n")

# Create output directories
test_dir <- "results/test_phase2"
if (!dir.exists(test_dir)) dir.create(test_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. SYNTHESIZE MOCK METADATA CSV
# ------------------------------------------------------------------------------
cat("--- Step 1: Synthesizing Mock Metadata CSV ---\n")

# Read some raw columns from RNA to get real barcodes
rna_raw <- readRDS("data/rna_expression_raw.rds")
if (inherits(rna_raw, "SummarizedExperiment")) {
  rna_barcodes <- colnames(assay(rna_raw))
} else {
  rna_barcodes <- colnames(rna_raw)
}

# Keep only primary tumors for synthesis
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
cat(sprintf("  Successfully synthesized metadata for %d samples\n", nrow(meta_df)))
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
# 6. RUN MULTI-OMICS INTEGRATION (MOFA+) IN METADATA MODE
# ------------------------------------------------------------------------------
cat("\n--- Step 6: Running Multi-Omics Integration (MOFA+) ---\n")
int_out <- file.path(test_dir, "integration")
if (!dir.exists(int_out)) dir.create(int_out, recursive = TRUE)

int_cmd <- sprintf(
  "Rscript modules/integration/run_integration.R --rna %s --meth %s --cnv %s --snv %s --metadata %s --outdir %s --factors 3 --iter 10",
  file.path(rna_out, "rna_processed_matrix.rds"),
  file.path(meth_out, "methylation_processed_matrix.rds"),
  file.path(cnv_out, "cnv_processed_matrix.rds"),
  file.path(snv_out, "snv_processed_matrix.rds"),
  metadata_file, int_out
)
cat(sprintf("  Executing: %s\n", int_cmd))
system(int_cmd)


# ------------------------------------------------------------------------------
# 7. SYNTHESIZE CUSTOM CLINICAL DATA WITH MAPPED PATIENT IDS
# ------------------------------------------------------------------------------
cat("\n--- Step 7: Synthesizing Custom Clinical Data ---\n")
clin_raw <- read.delim("data/clinical_data.tsv", stringsAsFactors = FALSE)

clin_patients <- clin_raw$bcr_patient_barcode
clin_standard_patients <- substr(clin_patients, 1, 12)
clin_keep_idx <- which(clin_standard_patients %in% unique_patients)

custom_clin_raw <- clin_raw[clin_keep_idx, , drop = FALSE]
custom_clin_pids <- custom_patient_map[clin_standard_patients[clin_keep_idx]]

# Map days to single survival time column
vital_status <- custom_clin_raw$vital_status
days_to_death <- as.numeric(custom_clin_raw$days_to_death)
days_to_follow <- as.numeric(custom_clin_raw$days_to_last_follow_up)
os_time <- ifelse(vital_status == "Dead", days_to_death, days_to_follow)
os_time[is.na(os_time)] <- 100 # fallback

# Rename/synthesize custom clinical columns
custom_clin <- data.frame(
  Subject_ID = custom_clin_pids,
  Time_Days = os_time,
  Event_Occurred = ifelse(vital_status == "Dead", 1, 0),
  Age_Years = as.numeric(custom_clin_raw$age_at_diagnosis) / 365.25,
  Sex = custom_clin_raw$gender,
  stringsAsFactors = FALSE
)

custom_clinical_file <- file.path(test_dir, "custom_clinical.csv")
write.csv(custom_clin, custom_clinical_file, row.names = FALSE)

# Write custom clinical mapping JSON
clinical_map_file <- file.path(test_dir, "clinical_map.json")
clinical_map <- list(
  patient_id = "Subject_ID",
  os_time    = "Time_Days",
  os_event   = "Event_Occurred",
  age        = "Age_Years",
  gender     = "Sex"
)
write(jsonlite::toJSON(clinical_map, auto_unbox = TRUE, pretty = TRUE), clinical_map_file)
cat(sprintf("  Successfully synthesized custom clinical and mapping JSON\n"))


# ------------------------------------------------------------------------------
# 8. RUN MACHINE LEARNING MODULE IN METADATA MODE
# ------------------------------------------------------------------------------
cat("\n--- Step 8: Running Machine Learning ---\n")
ml_out <- file.path(test_dir, "ml")
if (!dir.exists(ml_out)) dir.create(ml_out, recursive = TRUE)

ml_cmd <- sprintf(
  "Rscript modules/ml/run_ml.R --mofa %s --rna %s --clinical %s --clinical_map %s --metadata %s --outdir %s --mofa_prefilter 10 --final_features 5 --train_ratio 0.8",
  file.path(int_out, "mofa_model.rds"),
  file.path(rna_out, "rna_ml.rds"),
  custom_clinical_file, clinical_map_file, metadata_file, ml_out
)
cat(sprintf("  Executing: %s\n", ml_cmd))
system(ml_cmd)


# ------------------------------------------------------------------------------
# 9. RUN ENRICHMENT MODULE IN METADATA MODE
# ------------------------------------------------------------------------------
cat("\n--- Step 9: Running Pathway Enrichment ---\n")
enr_out <- file.path(test_dir, "pathway")
if (!dir.exists(enr_out)) dir.create(enr_out, recursive = TRUE)

keywords_file <- file.path(test_dir, "keywords.json")
write(jsonlite::toJSON(c("extracellular", "collagen", "matrix"), auto_unbox = TRUE), keywords_file)

enr_cmd <- sprintf(
  "Rscript modules/enrichment/run_enrichment.R --mofa %s --rf %s --lasso %s --rna %s --validation_keywords %s --outdir %s",
  file.path(ml_out, "mofa_top_genes.rds"),
  file.path(ml_out, "rf_top_genes.rds"),
  file.path(ml_out, "lasso_selected_genes.rds"),
  file.path(ml_out, "rna_for_pathway.rds"),
  keywords_file, enr_out
)
cat(sprintf("  Executing: %s\n", enr_cmd))
system(enr_cmd)


# ------------------------------------------------------------------------------
# 10. GENERATE AUTOMATED CLINICAL REPORT
# ------------------------------------------------------------------------------
cat("\n--- Step 10: Generating Dynamic HTML Report ---\n")
report_out <- file.path(test_dir, "reports")
if (!dir.exists(report_out)) dir.create(report_out, recursive = TRUE)

report_cmd <- sprintf(
  "Rscript reports/render_report.R --results_dir %s --output_dir %s",
  test_dir, report_out
)
cat(sprintf("  Executing: %s\n", report_cmd))
system(report_cmd)


# ------------------------------------------------------------------------------
# 11. VERIFY ALL PHASES OF PIPELINE RUN SUCCESSFULLY
# ------------------------------------------------------------------------------
cat("\n--- Step 11: Running Validation Verification ---\n")

test_errors <- FALSE

# 11.1 Check RNA outputs
rna_scaled <- readRDS(file.path(rna_out, "rna_processed_matrix.rds"))
if (!all(colnames(rna_scaled) %in% meta_df$patient_id)) {
  cat("    [FAIL] RNA matrix columns are not remapped to patient IDs!\n")
  test_errors <- TRUE
} else {
  cat("    [PASS] RNA matrix columns successfully remapped to patient IDs\n")
}

# 11.2 Check Integration outputs
if (!file.exists(file.path(int_out, "mofa_model.rds"))) {
  cat("    [FAIL] MOFA model RDS does not exist!\n")
  test_errors <- TRUE
} else {
  cat("    [PASS] MOFA model successfully generated in metadata mode\n")
  mofa_loaded <- readRDS(file.path(int_out, "mofa_model.rds"))
  mofa_pids <- colnames(mofa_loaded$model@data$RNA[[1]])
  if (!all(mofa_pids %in% meta_df$patient_id)) {
    cat("    [FAIL] MOFA model contains untruncated/unmapped patient IDs!\n")
    test_errors <- TRUE
  } else {
    cat("    [PASS] MOFA model contains correct custom patient IDs (no substr 1-12 corruption)\n")
  }
}

# 11.3 Check Machine Learning outputs
if (!file.exists(file.path(ml_out, "ml_results_summary.csv"))) {
  cat("    [FAIL] ML results summary does not exist!\n")
  test_errors <- TRUE
} else {
  cat("    [PASS] ML survival analysis executed successfully using custom clinical mapping\n")
  ml_summary <- read.csv(file.path(ml_out, "ml_results_summary.csv"))
  if (nrow(ml_summary) != 3) {
    cat("    [FAIL] ML summary is incomplete!\n")
    test_errors <- TRUE
  } else {
    cat("    [PASS] All 3 survival models trained and evaluated successfully\n")
  }
}

# 11.4 Check Enrichment outputs
if (!file.exists(file.path(enr_out, "enrichment_results.rds"))) {
  cat("    [FAIL] Enrichment results RDS does not exist!\n")
  test_errors <- TRUE
} else {
  cat("    [PASS] Enrichment module executed successfully using custom keywords validation\n")
}

# 11.5 Check Report output
if (!file.exists(file.path(report_out, "OmicsFlow_Report.html"))) {
  cat("    [FAIL] Dynamic HTML report was not rendered!\n")
  test_errors <- TRUE
} else {
  cat("    [PASS] Dynamic cohort-agnostic HTML report rendered successfully\n")
  
  # Read rendered HTML file and verify no hardcoded TCGA mentions
  report_lines <- readLines(file.path(report_out, "OmicsFlow_Report.html"), warn = FALSE)
  report_text <- paste(report_lines, collapse = " ")
  if (grepl("TCGA omics data", report_text)) {
    cat("    [FAIL] Found hardcoded 'TCGA omics data' inside figure captions!\n")
    test_errors <- TRUE
  } else {
    cat("    [PASS] Verified no hardcoded 'TCGA omics data' in figure captions\n")
  }
}

cat("\n==============================================================\n")
if (test_errors) {
  cat("  PHASE 2 END-TO-END VALIDATION TESTS FAILED!\n")
  cat("==============================================================\n\n")
  quit(status = 1)
} else {
  cat("  ALL PHASE 2 END-TO-END VALIDATION TESTS PASSED SUCCESSFULLY!\n")
  cat("==============================================================\n\n")
  quit(status = 0)
}
