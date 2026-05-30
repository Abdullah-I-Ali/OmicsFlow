#!/usr/bin/env Rscript
# ==============================================================================
# preprocess_snv.R — SNV Mutation Preprocessing Pipeline
# OmicsFlow | Phase 4: SNV Module
# ==============================================================================
#
# PURPOSE:
#   Executes the scientific steps for filtering and transforming SNV MAF data
#   into a binary gene x patient matrix.
#   
# METHODOLOGY:
#   1. Load raw mutation data (MAF format).
#   2. Filter functional mutations (remove silent/intron/IGR).
#   3. Standardize barcodes to 12-char patient level.
#   4. Collapse duplicate patient records (any mutation = mutated).
#   5. Build Gene x Patient binary matrix (1=mutated, 0=WT).
#   6. Remove hypermutated patients (Top 1% by mutation count).
#   7. Frequency-based feature selection (>= 2% prevalence).
#
# ==============================================================================
# USAGE:
#   Rscript preprocess_snv.R \
#     --input      data/snv_mutation_raw.rds \
#     --outdir     results/snv/ \
#     --freq       0.02
# ==============================================================================

suppressPackageStartupMessages(library(optparse))

# ------------------------------------------------------------------------------
# COMMAND-LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "data/snv_mutation_raw.rds",
              help = "Path to raw SNV mutation data (.rds) [default= %default]"),
  make_option(c("-o", "--outdir"), type = "character", default = "results/snv/",
              help = "Output directory for SNV results [default= %default]"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Path to sample metadata CSV file (optional)"),
  make_option(c("-f", "--freq"), type = "numeric", default = 0.02,
              help = "Minimum mutation frequency threshold [default= %default]"),
  make_option(c("-t", "--topn"), type = "integer", default = 3000,
              help = "Max number of genes to retain [default= %default]"),
  make_option(c("-s", "--seed"), type = "integer", default = 42,
              help = "Random seed for reproducibility [default= %default]")
)

opt_parser <- OptionParser(
  usage = "Usage: %prog [options]",
  option_list = option_list,
  description = "OmicsFlow Phase 4: SNV Mutation Preprocessing Module"
)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("--input is required.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# SOURCE MODULE FILES
# ------------------------------------------------------------------------------
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) "modules/snv"
)
source(file.path(script_dir, "utils_snv.R"))
source(file.path(script_dir, "qc_snv.R"))
source(file.path(script_dir, "export_snv.R"))
source(file.path(dirname(script_dir), "utils_metadata.R"))

# ------------------------------------------------------------------------------
# INITIALISE
# ------------------------------------------------------------------------------
tryCatch({
  load_snv_packages()
  set.seed(opt$seed)
  
  if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
  }
  
  qc_metrics <- init_snv_qc()
  
  metadata <- load_metadata(opt$metadata)
  
  snv_banner("OmicsFlow v1.0.0 | SNV Preprocessing")
  snv_msg(sprintf("Start time : %s", Sys.time()))
  snv_msg(sprintf("Input      : %s", opt$input))
  snv_msg(sprintf("Output dir : %s", opt$outdir))
  if (!is.null(metadata)) {
    snv_msg(sprintf("Metadata   : %s (%d samples)", opt$metadata, nrow(metadata)))
  }
  
  # ==============================================================================
  # 1. LOAD DATA
  # ==============================================================================
  snv_step(1, "Loading Raw SNV Data")
  snv_raw <- load_snv_data(opt$input)
  
  maf_obj <- read.maf(maf = snv_raw)
  n_initial_mutations <- nrow(maf_obj@data)
  n_initial_samples <- length(unique(maf_obj@data$Tumor_Sample_Barcode))
  
  snv_msg(sprintf("Loaded MAF: %s mutations across %d samples",
                  format(n_initial_mutations, big.mark = ","),
                  n_initial_samples))
  qc_metrics <- add_snv_qc(qc_metrics, "samples", "initial_samples", n_initial_samples)
  qc_metrics <- add_snv_qc(qc_metrics, "filters", "initial_mutations", n_initial_mutations)
  
  # ==============================================================================
  # 2. FILTER MUTATIONS (FUNCTIONAL VARIANTS ONLY)
  # ==============================================================================
  snv_step(2, "Filtering Functional Mutations")
  
  valid_variants <- c(
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Splice_Site",
    "Nonstop_Mutation"
  )
  
  maf_filtered <- subsetMaf(
    maf   = maf_obj,
    query = paste0("Variant_Classification %in% c('",
                   paste(valid_variants, collapse = "','"), "')")
  )
  
  n_retained_variants <- nrow(maf_filtered@data)
  n_removed_variants <- n_initial_mutations - n_retained_variants
  
  snv_msg(sprintf("Retained %s functional mutations", format(n_retained_variants, big.mark = ",")))
  snv_msg(sprintf("Removed %s silent/non-coding mutations", format(n_removed_variants, big.mark = ",")), level="DETAILS")
  qc_metrics <- add_snv_qc(qc_metrics, "filters", "retained_functional_mutations", n_retained_variants)
  qc_metrics <- add_snv_qc(qc_metrics, "filters", "removed_nonfunctional_mutations", n_removed_variants)
  
  # ==============================================================================
  # 3. EXTRACT MUTATION DATA
  # ==============================================================================
  snv_step(3, "Extracting Mutation Table")
  
  mut_df <- as.data.frame(maf_filtered@data)
  mut_df <- mut_df[, c("Hugo_Symbol", "Tumor_Sample_Barcode")]
  colnames(mut_df) <- c("gene", "sample_id")
  
  snv_msg(sprintf("Extracted %s mutation records", format(nrow(mut_df), big.mark = ",")))
  
  # ==============================================================================
  # 4. STANDARDIZE BARCODES TO 12-CHAR (PATIENT LEVEL) & COLLAPSE
  # ==============================================================================
  snv_step(4, "Standardizing Barcodes to Patient Level")
  
  if (!is.null(metadata)) {
    # Filter to samples present in metadata
    mut_df <- mut_df %>% filter(sample_id %in% metadata$sample_id)
    n_samples_before <- length(unique(mut_df$sample_id))
    mut_df$sample_id <- get_patient_id(mut_df$sample_id, metadata)
    n_patients_after <- length(unique(mut_df$sample_id))
    snv_msg(sprintf("Filtered & mapped to metadata: %d unique samples -> %d unique patients", 
                    n_samples_before, n_patients_after))
  } else {
    n_samples_before <- length(unique(mut_df$sample_id))
    mut_df$sample_id <- substr(mut_df$sample_id, 1, 12)
    n_patients_after <- length(unique(mut_df$sample_id))
    snv_msg(sprintf("Standardized to 12-char patient barcodes: %d unique samples -> %d unique patients", 
                    n_samples_before, n_patients_after))
  }
  
  n_before_collapse <- nrow(mut_df)
  mut_df <- mut_df %>% distinct(gene, sample_id)
  n_after_collapse <- nrow(mut_df)
  
  if (n_after_collapse < n_before_collapse) {
    snv_msg(sprintf("Collapsed duplicate records: %s -> %s (any mutation = mutated)",
                    format(n_before_collapse, big.mark = ","),
                    format(n_after_collapse, big.mark = ",")))
  }
  qc_metrics <- add_snv_qc(qc_metrics, "samples", "patient_level_standardized", n_patients_after)
  
  # ==============================================================================
  # 5. BUILD GENE-LEVEL BINARY MUTATION MATRIX
  # ==============================================================================
  snv_step(5, "Building Gene x Patient Binary Matrix")
  
  mut_df$mutated <- 1
  
  snv_matrix <- mut_df %>%
    distinct() %>%
    pivot_wider(
      names_from  = sample_id,
      values_from = mutated,
      values_fill = 0
    ) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  rownames(snv_matrix) <- make.unique(rownames(snv_matrix))
  
  snv_msg(sprintf("Binary matrix: %d genes x %d patients", nrow(snv_matrix), ncol(snv_matrix)))
  
  top_mutated <- sort(rowSums(snv_matrix), decreasing = TRUE)[1:5]
  snv_msg("Top 5 most mutated genes:", level = "DETAILS")
  for (i in seq_along(top_mutated)) {
    snv_msg(sprintf("%d. %-12s - %d patients (%.1f%%)",
                    i, names(top_mutated)[i], top_mutated[i],
                    100 * top_mutated[i] / ncol(snv_matrix)), level = "DETAILS")
  }
  
  # ==============================================================================
  # 6. REMOVE HYPERMUTATED PATIENTS
  # ==============================================================================
  snv_step(6, "Removing Hypermutated Patients")
  
  patient_mut_count <- colSums(snv_matrix)
  threshold         <- quantile(patient_mut_count, 0.99)
  n_removed         <- sum(patient_mut_count >= threshold)
  
  if (n_removed > 0) {
    snv_matrix <- snv_matrix[, patient_mut_count < threshold, drop = FALSE]
    snv_msg(sprintf("Mutation burden threshold (99th percentile): %.0f mutations", threshold))
    snv_msg(sprintf("Removed %d hypermutated patient(s)", n_removed))
    snv_msg(sprintf("Patients retained: %d", ncol(snv_matrix)))
  } else {
    snv_msg("No hypermutated patients detected (all below 99th percentile)")
  }
  qc_metrics <- add_snv_qc(qc_metrics, "filters", "hypermutated_patients_removed", n_removed)
  qc_metrics <- add_snv_qc(qc_metrics, "samples", "final", ncol(snv_matrix))
  
  # ==============================================================================
  # 7. FREQUENCY-BASED FEATURE SELECTION
  # ==============================================================================
  snv_step(7, sprintf("Selecting Recurrently Mutated Genes (>= %.1f%%)", opt$freq * 100))
  
  mutation_freq <- rowMeans(snv_matrix)
  snv_filtered  <- snv_matrix[mutation_freq >= opt$freq, , drop = FALSE]
  
  if (nrow(snv_filtered) == 0) {
    current_max_freq <- max(mutation_freq)
    stop(paste(
      sprintf("\n[ERROR] No genes passed the %.1f%% mutation frequency threshold.", opt$freq * 100),
      "\n   Current max mutation frequency:", round(current_max_freq, 3),
      "\n   Consider lowering the cutoff.",
      "\n   Current patient count:", ncol(snv_matrix)
    ))
  }
  
  filtered_freq <- mutation_freq[mutation_freq >= opt$freq]
  snv_msg(sprintf("Genes mutated in >= %.1f%% of patients: %d", opt$freq * 100, nrow(snv_filtered)))
  
  top_n <- min(opt$topn, nrow(snv_filtered))
  
  if (nrow(snv_filtered) <= 500) {
    top_genes <- names(sort(filtered_freq, decreasing = TRUE))
    snv_msg(sprintf("Keeping ALL %d genes (discovery mode)", nrow(snv_filtered)))
  } else {
    top_genes <- names(sort(filtered_freq, decreasing = TRUE))[1:top_n]
    snv_msg(sprintf("Selected top %d genes by mutation frequency", top_n))
  }
  
  snv_final <- snv_filtered[top_genes, , drop = FALSE]
  snv_final <- snv_final[order(rownames(snv_final)), , drop = FALSE]
  
  snv_msg(sprintf("Final gene set: %d genes", nrow(snv_final)))
  snv_msg(sprintf("Frequency range: [%.4f, %.4f]", min(filtered_freq[top_genes]), max(filtered_freq[top_genes])))
  
  qc_metrics <- add_snv_qc(qc_metrics, "genes", "final_genes", nrow(snv_final))
  
  # ==============================================================================
  # 8. QC CHECKS
  # ==============================================================================
  snv_step(8, "Quality Control Checks")
  
  stopifnot("NAs should not exist in binary matrix" = !any(is.na(snv_final)))
  stopifnot("Gene names must be unique" = !any(duplicated(rownames(snv_final))))
  
  mean_mut_rate <- mean(snv_final)
  snv_msg("No missing values")
  snv_msg("No duplicate gene names")
  snv_msg(sprintf("Matrix dimensions: %d genes x %d patients", nrow(snv_final), ncol(snv_final)))
  snv_msg(sprintf("Matrix sparsity (mean mutation rate): %.4f", mean_mut_rate))
  snv_msg(sprintf("Total mutation events: %s", format(sum(snv_final), big.mark = ",")))
  
  qc_metrics <- add_snv_qc(qc_metrics, "parameters", "mean_mutation_rate", mean_mut_rate)
  
  if (mean_mut_rate < 0.10) {
    snv_msg("Sparse matrix is expected for SNV data", level = "DETAILS")
  }
  
  zero_mut_patients <- sum(colSums(snv_final) == 0)
  if (zero_mut_patients > 0) {
    snv_msg(sprintf("%d patients have zero mutations in selected genes", zero_mut_patients), level = "WARN")
  }
  qc_metrics <- add_snv_qc(qc_metrics, "samples", "zero_mutation_patients", zero_mut_patients)
  
  # ==============================================================================
  # 9. VISUALIZATION & OUTPUT
  # ==============================================================================
  generate_snv_qc_plots(snv_final, maf_filtered, opt$outdir)
  if (!is.null(metadata)) {
    sample_info <- metadata[match(colnames(snv_final), metadata$patient_id), ]
  } else {
    sample_info <- data.frame(
      sample_id        = colnames(snv_final),
      patient_id       = colnames(snv_final),
      sample_class     = "primary_tumor",
      batch            = "batch1",
      stringsAsFactors = FALSE
    )
  }
  export_snv_results(snv_final, maf_filtered, opt$outdir, sample_info)
  export_snv_qc(qc_metrics, opt$outdir)
  
  snv_banner("SNV MODULE COMPLETE \u2714")
  
}, error = function(e) {
  snv_msg(sprintf("FATAL ERROR: %s", e$message), level = "ERROR")
  quit(status = 1)
})
