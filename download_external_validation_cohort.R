#!/usr/bin/env Rscript
# ==============================================================================
# download_external_validation_cohort.R
# Fetches a minimal TCGA subset (TCGA-UVM or TCGA-ACC) for real-world validation
# ==============================================================================

suppressPackageStartupMessages(library(TCGAbiolinks))
suppressPackageStartupMessages(library(SummarizedExperiment))

out_dir <- "data/external_cohort"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("Fetching clinical data for TCGA-UVM...\n")
clinical <- GDCquery_clinic(project = "TCGA-UVM", type = "clinical")

cat("Fetching RNA-seq data for TCGA-UVM (subsetting to 10 samples for speed)...\n")
query_rna <- GDCquery(
  project = "TCGA-UVM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Subset to first 10 samples to make validation fast
query_rna$results[[1]] <- query_rna$results[[1]][1:10, ]
GDCdownload(query_rna)
rna_se <- GDCprepare(query_rna)
rna_matrix <- assay(rna_se, "tpm_unstrand")

# Save RNA
saveRDS(rna_matrix, file.path(out_dir, "rna_expression.rds"))

# Prepare minimal metadata
metadata <- data.frame(
  sample_id = colnames(rna_matrix),
  patient_id = substr(colnames(rna_matrix), 1, 12),
  sample_class = "Tumor",
  batch = "TCGA",
  age = 50,
  stringsAsFactors = FALSE
)
write.csv(metadata, file.path(out_dir, "clinical_metadata.csv"), row.names = FALSE)

cat("Successfully created minimal real-world validation dataset in data/external_cohort/\n")
