#!/usr/bin/env Rscript
# ==============================================================================
# preprocess_rna.R — RNA-seq Preprocessing Pipeline
# OmicsFlow | Phase 1: RNA Module
# Scientific source: data/full_scripts.R (lines 1-628)
# ==============================================================================
# USAGE:
#   Rscript preprocess_rna.R \
#     --input      path/to/rna_expression_raw.rds \
#     --outdir     results/rna/ \
#     --n-top      2000 \
#     --cor-low    0.75 \
#     --cor-high   0.895
# ==============================================================================

suppressPackageStartupMessages(library(optparse))

# ------------------------------------------------------------------------------
# CLI ARGUMENT PARSING
# ------------------------------------------------------------------------------
option_list <- list(
  make_option("--input",        type="character", default=NULL,
              help="Path to input RNA-seq RDS (SummarizedExperiment or matrix)"),
  make_option("--outdir",       type="character", default="results/rna/",
              help="Output directory [default: results/rna/]"),
  make_option("--metadata",     type="character", default=NULL,
              help="Path to sample metadata CSV file (optional)"),
  make_option("--cpm-threshold",type="double",    default=1.0,
              help="CPM threshold for low-expression filtering [default: 1]"),
  make_option("--cpm-min-pct",  type="double",    default=0.20,
              help="Min fraction of samples above CPM threshold [default: 0.20]"),
  make_option("--prior-count",  type="double",    default=2.0,
              help="Prior count for log2 CPM stabilization [default: 2]"),
  make_option("--n-top",        type="integer",   default=2000L,
              help="Number of top variable genes to select [default: 2000]"),
  make_option("--mean-expr-min",type="double",    default=1.0,
              help="Min mean log2 CPM after feature selection [default: 1]"),
  make_option("--cor-low",      type="double",    default=0.75,
              help="Low Pearson correlation outlier threshold [default: 0.75]"),
  make_option("--cor-high",     type="double",    default=0.895,
              help="High Pearson correlation outlier threshold [default: 0.895]"),
  make_option("--sd-threshold", type="double",    default=1e-8,
              help="Min SD for zero-variance removal [default: 1e-8]"),
  make_option("--na-threshold", type="double",    default=0.20,
              help="Max NA fraction per gene before removal [default: 0.20]")
)

parser <- OptionParser(option_list = option_list,
                       description = "OmicsFlow RNA-seq Preprocessing Pipeline")
args   <- parse_args(parser)

if (is.null(args$input)) {
  print_help(parser)
  stop("--input is required.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# SOURCE MODULE FILES
# ------------------------------------------------------------------------------
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) "modules/rna"
)
source(file.path(script_dir, "utils_rna.R"))
source(file.path(script_dir, "qc_rna.R"))
source(file.path(script_dir, "export_rna.R"))
source(file.path(dirname(script_dir), "utils_metadata.R"))

# ------------------------------------------------------------------------------
# INITIALISE
# ------------------------------------------------------------------------------
rna_load_packages()
set.seed(42)
rna_ensure_dirs(args$outdir)
qc <- rna_init_qc()

metadata <- load_metadata(args$metadata)

rna_banner("OmicsFlow v1.0.0 | RNA-seq Preprocessing")
rna_msg(sprintf("Start time : %s", Sys.time()))
rna_msg(sprintf("Input      : %s", args$input))
rna_msg(sprintf("Output dir : %s", args$outdir))
if (!is.null(metadata)) {
  rna_msg(sprintf("Metadata   : %s (%d samples)", args$metadata, nrow(metadata)))
}

# ==============================================================================
# STEP 1 — DATA LOADING
# ==============================================================================
rna_step(1, "Loading RNA-seq Data")
rna_check_file(args$input, "RNA input")

rna_raw <- readRDS(args$input)
if (inherits(rna_raw, "SummarizedExperiment")) {
  rna_counts <- assay(rna_raw)
} else {
  rna_counts <- as.matrix(rna_raw)
}
storage.mode(rna_counts) <- "double"

qc <- rna_update_qc(qc, "genes", "raw", nrow(rna_counts))
qc <- rna_update_qc(qc, "samples", "raw", ncol(rna_counts))
rna_msg(sprintf("Raw matrix : %s genes x %s samples",
                format(nrow(rna_counts), big.mark=","),
                format(ncol(rna_counts), big.mark=",")))

# ==============================================================================
# STEP 2 — SAMPLE SELECTION: PRIMARY TUMOR ONLY + DEDUPLICATION
# ==============================================================================
rna_step(2, "Sample Selection & Deduplication")

if (!is.null(metadata)) {
  # --- Filter columns to those in metadata ---
  common_samples <- intersect(colnames(rna_counts), metadata$sample_id)
  if (length(common_samples) == 0) {
    stop("None of the raw count matrix column names match the sample_id column in metadata.")
  }
  rna_counts <- rna_counts[, common_samples, drop = FALSE]
  
  qc <- rna_update_qc(qc, "samples", "primary_tumor", ncol(rna_counts))
  rna_msg(sprintf("Samples matched to metadata : %d", ncol(rna_counts)))
  
  # --- Save FULL barcodes/sample IDs before any truncation ---
  full_barcodes <- colnames(rna_counts)
  
  # --- Deduplicate: per patient keep aliquot with highest library size ---
  patient_ids   <- get_patient_id(colnames(rna_counts), metadata)
  lib_sizes_all <- colSums(rna_counts)
  
  keep_samples <- tapply(
    seq_len(ncol(rna_counts)),
    patient_ids,
    function(idx) idx[which.max(lib_sizes_all[idx])]
  )
  keep_idx <- sort(unlist(keep_samples))
  
  rna_counts <- rna_counts[, keep_idx, drop = FALSE]
  
  # Preserve original barcodes/sample IDs for kept samples
  original_barcodes        <- full_barcodes[keep_idx]
  names(original_barcodes) <- get_patient_id(original_barcodes, metadata)
  
  # Rename columns to patient IDs
  colnames(rna_counts) <- get_patient_id(colnames(rna_counts), metadata)
  
  qc <- rna_update_qc(qc, "samples", "after_deduplication", ncol(rna_counts))
  rna_msg(sprintf("After deduplication   : %d unique patients", ncol(rna_counts)))
} else {
  # --- Keep only primary tumor samples (TCGA barcode positions 14-15 = "01") ---
  sample_type_code <- substr(colnames(rna_counts), 14, 15)
  rna_counts       <- rna_counts[, sample_type_code == "01", drop = FALSE]
  qc <- rna_update_qc(qc, "samples", "primary_tumor", ncol(rna_counts))
  rna_msg(sprintf("Primary tumor samples : %d", ncol(rna_counts)))
  
  # --- Save FULL barcodes before any truncation (needed for batch extraction) ---
  full_barcodes <- colnames(rna_counts)
  
  # --- Deduplicate: per patient keep aliquot with highest library size ---
  patient_ids   <- substr(colnames(rna_counts), 1, 12)
  lib_sizes_all <- colSums(rna_counts)
  
  keep_samples <- tapply(
    seq_len(ncol(rna_counts)),
    patient_ids,
    function(idx) idx[which.max(lib_sizes_all[idx])]
  )
  keep_idx <- sort(unlist(keep_samples))
  
  rna_counts <- rna_counts[, keep_idx, drop = FALSE]
  
  # Preserve full barcodes for the kept samples (needed for Plate ID extraction)
  original_barcodes        <- full_barcodes[keep_idx]
  names(original_barcodes) <- substr(original_barcodes, 1, 12)
  
  # Rename columns to 12-char patient IDs for downstream integration
  colnames(rna_counts) <- substr(colnames(rna_counts), 1, 12)
  
  qc <- rna_update_qc(qc, "samples", "after_deduplication", ncol(rna_counts))
  rna_msg(sprintf("After deduplication   : %d unique patients", ncol(rna_counts)))
  rna_msg("Full barcodes saved for batch extraction")
}

# ==============================================================================
# STEP 3 — QUALITY CONTROL: LIBRARY SIZE + ZERO-GENE REMOVAL
# ==============================================================================
rna_step(3, "Quality Control")

# Library size statistics
lib_sizes <- colSums(rna_counts)
qc <- rna_update_qc(qc, "library_size", "min_reads",    min(lib_sizes))
qc <- rna_update_qc(qc, "library_size", "max_reads",    max(lib_sizes))
qc <- rna_update_qc(qc, "library_size", "median_reads", median(lib_sizes))
rna_msg(sprintf("Library size : %.1fM - %.1fM reads (median %.1fM)",
                min(lib_sizes)/1e6, max(lib_sizes)/1e6, median(lib_sizes)/1e6))

# Flag low-depth samples (warn only — do NOT remove)
low_depth <- names(lib_sizes[lib_sizes < 1e6])
qc <- rna_update_qc(qc, "samples", "low_depth_flagged", length(low_depth))
if (length(low_depth) > 0) {
  warning(sprintf("  %d low-depth sample(s) (<1M reads) flagged (not removed)",
                  length(low_depth)))
}

# Remove all-zero genes
n_before   <- nrow(rna_counts)
rna_counts <- rna_counts[rowSums(rna_counts) > 0, , drop = FALSE]
qc <- rna_update_qc(qc, "genes", "after_zero_removal", nrow(rna_counts))
rna_msg(sprintf("Removed %d all-zero genes", n_before - nrow(rna_counts)))

# ==============================================================================
# STEP 4 — LOW-EXPRESSION FILTERING (CPM > 1 in >= 20% of samples)
# ==============================================================================
rna_step(4, "Low-Expression Gene Filtering")
rna_msg(sprintf("Criterion: CPM > %.0f in >= %.0f%% of samples",
                args[["cpm-threshold"]], args[["cpm-min-pct"]] * 100))

cpm_matrix <- cpm(rna_counts)
keep_genes <- rowSums(cpm_matrix > args[["cpm-threshold"]]) >=
              args[["cpm-min-pct"]] * ncol(rna_counts)
rna_counts <- rna_counts[keep_genes, , drop = FALSE]

qc <- rna_update_qc(qc, "genes", "after_cpm_filter", nrow(rna_counts))
rna_msg(sprintf("Genes retained : %s", format(nrow(rna_counts), big.mark=",")))

# ==============================================================================
# STEP 5 — NORMALIZATION: TMM + log2 CPM (prior.count = 2)
# ==============================================================================
rna_step(5, "TMM Normalization + log2 CPM")
rna_msg(sprintf("Method: TMM + cpm(log=TRUE, prior.count=%.0f)",
                args[["prior-count"]]))

dge      <- DGEList(counts = rna_counts)
dge      <- calcNormFactors(dge, method = "TMM")
rna_lcpm <- cpm(dge, log = TRUE, prior.count = args[["prior-count"]])

rna_msg(sprintf("Normalized matrix : %s genes x %d samples",
                format(nrow(rna_lcpm), big.mark=","), ncol(rna_lcpm)))
rna_msg(sprintf("Expression range  : [%.2f, %.2f]",
                min(rna_lcpm), max(rna_lcpm)))

# ==============================================================================
# STEP 6 — GENE ID MAPPING: ENSEMBL -> SYMBOL (multiVals = "asNA")
# ==============================================================================
rna_step(6, "Gene ID Mapping (Ensembl -> Symbol)")

# Strip Ensembl version suffix (e.g. ENSG00000001234.5 -> ENSG00000001234)
ensembl_ids <- gsub("\\..*", "", rownames(rna_lcpm))

# Ambiguous 1:many mappings -> NA (not collapsed) — multiVals = "asNA" is critical
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys      = ensembl_ids,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "asNA"
)

valid_idx    <- !is.na(gene_symbols)
rna_mapped   <- rna_lcpm[valid_idx, , drop = FALSE]
gene_symbols <- gene_symbols[valid_idx]

qc <- rna_update_qc(qc, "genes", "after_ensembl_mapping", nrow(rna_mapped))
rna_msg(sprintf("Valid symbol mappings : %s",
                format(nrow(rna_mapped), big.mark=",")))

# ==============================================================================
# STEP 7 — RESOLVE DUPLICATE GENE SYMBOLS (keep highest mean expression)
# ==============================================================================
rna_step(7, "Resolving Duplicate Gene Symbols")
rna_msg("Strategy: keep isoform with highest mean expression per symbol")

rna_df           <- as.data.frame(rna_mapped)
rna_df$symbol    <- gene_symbols
rna_df$mean_expr <- rowMeans(rna_df[, colnames(rna_df) != "symbol", drop=FALSE],
                              na.rm = TRUE)

rna_unique <- rna_df %>%
  group_by(symbol) %>%
  slice_max(order_by = mean_expr, n = 1, with_ties = FALSE) %>%
  ungroup()

rna_unique$mean_expr <- NULL

symbols_vec             <- rna_unique$symbol
expr_data               <- as.matrix(rna_unique[, colnames(rna_unique) != "symbol"])
storage.mode(expr_data) <- "double"
rownames(expr_data)     <- symbols_vec
rna_mat                 <- expr_data

qc <- rna_update_qc(qc, "genes", "after_symbol_dedup", nrow(rna_mat))
rna_msg(sprintf("Unique genes after deduplication : %s",
                format(nrow(rna_mat), big.mark=",")))

# ==============================================================================
# STEP 8 — MATRIX QUALITY FILTERS (NA, imputation, zero-variance)
# ==============================================================================
rna_step(8, "Matrix Quality Filters")

# Remove genes with > na-threshold fraction of NA values
na_fraction <- rowMeans(is.na(rna_mat))
n_removed_na <- sum(na_fraction >= args[["na-threshold"]])
rna_mat     <- rna_mat[na_fraction < args[["na-threshold"]], , drop = FALSE]
rna_msg(sprintf("Removed %d genes with > %.0f%% NA",
                n_removed_na, args[["na-threshold"]] * 100))

# Impute remaining NAs with per-gene (row) mean
n_na_total <- sum(is.na(rna_mat))
if (n_na_total > 0) {
  for (i in seq_len(nrow(rna_mat))) {
    na_idx <- is.na(rna_mat[i, ])
    if (any(na_idx)) {
      rna_mat[i, na_idx] <- mean(rna_mat[i, ], na.rm = TRUE)
    }
  }
  rna_msg(sprintf("Imputed %d remaining NAs with row means", n_na_total))
}

# Remove non-finite rows
finite_rows <- is.finite(rowSums(rna_mat))
rna_mat     <- rna_mat[finite_rows, , drop = FALSE]

# Remove zero-variance genes (sd <= sd-threshold)
gene_sd <- apply(rna_mat, 1, sd)
rna_mat <- rna_mat[gene_sd > args[["sd-threshold"]], , drop = FALSE]

qc <- rna_update_qc(qc, "genes", "after_variance_filter", nrow(rna_mat))
rna_msg(sprintf("After NA filter        : %s genes",
                format(nrow(rna_mat), big.mark=",")))
rna_msg(sprintf("Final clean matrix     : %s genes x %d samples",
                format(nrow(rna_mat), big.mark=","), ncol(rna_mat)))

# ==============================================================================
# STEP 9 — OUTLIER SAMPLE DETECTION (Pearson correlation thresholds)
# ==============================================================================
rna_step(9, "Outlier Sample Detection")
rna_msg(sprintf("Method     : Mean inter-sample Pearson correlation"))
rna_msg(sprintf("Thresholds : low < %.3f | high > %.3f",
                args[["cor-low"]], args[["cor-high"]]))

sample_cor <- cor(rna_mat, method = "pearson")
mean_cor   <- rowMeans(sample_cor)

low_outliers      <- names(mean_cor[mean_cor < args[["cor-low"]]])
high_outliers     <- names(mean_cor[mean_cor > args[["cor-high"]]])
samples_to_remove <- unique(c(low_outliers, high_outliers))

qc <- rna_update_qc(qc, "outliers", "low_correlation_threshold",  args[["cor-low"]])
qc <- rna_update_qc(qc, "outliers", "high_correlation_threshold", args[["cor-high"]])
qc <- rna_update_qc(qc, "outliers", "low_outliers_count",         length(low_outliers))
qc <- rna_update_qc(qc, "outliers", "high_outliers_count",        length(high_outliers))
qc <- rna_update_qc(qc, "outliers", "total_removed",              length(samples_to_remove))

rna_msg(sprintf("Low-correlation outliers  (<%.3f) : %d", args[["cor-low"]],
                length(low_outliers)))
rna_msg(sprintf("High-correlation outliers (>%.3f) : %d", args[["cor-high"]],
                length(high_outliers)))
rna_msg(sprintf("Total samples to remove           : %d",
                length(samples_to_remove)))

if (length(samples_to_remove) > 0) {
  rna_mat <- rna_mat[, !colnames(rna_mat) %in% samples_to_remove, drop = FALSE]
}
qc <- rna_update_qc(qc, "samples", "after_outlier_removal", ncol(rna_mat))
rna_msg(sprintf("Samples remaining : %d", ncol(rna_mat)))

# ==============================================================================
# STEP 10 — BATCH EFFECT CORRECTION (Plate ID, positions 22-25)
# ==============================================================================
rna_step(10, "Batch Effect Correction")
if (!is.null(metadata)) {
  rna_msg("Batch variable: From Metadata")
} else {
  rna_msg("Batch variable: Plate ID (TCGA barcode positions 22-25)")
}

# Recover full barcodes for surviving samples using the saved named vector
matched_barcodes <- original_barcodes[colnames(rna_mat)]
plate_id         <- get_batch(matched_barcodes, metadata, omics_type = "rna")
batch_info       <- as.factor(plate_id)

n_batches_total <- nlevels(batch_info)
qc <- rna_update_qc(qc, "batches", "total_batches_detected", n_batches_total)

# Remove singleton batches (batches with only 1 sample cannot be corrected)
plate_table       <- table(batch_info)
singleton_batches <- names(plate_table[plate_table == 1])
qc <- rna_update_qc(qc, "batches", "singleton_batches", length(singleton_batches))

if (length(singleton_batches) > 0) {
  rna_msg(sprintf("Removing %d singleton batch(es): %s",
                  length(singleton_batches),
                  paste(singleton_batches, collapse=", ")))
  keep_samp  <- !batch_info %in% singleton_batches
  rna_mat    <- rna_mat[, keep_samp, drop = FALSE]
  batch_info <- droplevels(batch_info[keep_samp])
}
qc <- rna_update_qc(qc, "samples", "after_singleton_removal", ncol(rna_mat))

# PCA BEFORE batch correction (scale. = FALSE — data is already log2 CPM)
pca_before <- prcomp(t(rna_mat), scale. = FALSE)
df_before  <- data.frame(
  PC1   = pca_before$x[, 1],
  PC2   = pca_before$x[, 2],
  Batch = batch_info
)

# Apply removeBatchEffect (limma) — only when >= 2 batches remain
if (nlevels(batch_info) > 1) {
  rna_msg(sprintf("Running removeBatchEffect on %d batches...",
                  nlevels(batch_info)))
  rna_mat <- removeBatchEffect(rna_mat, batch = batch_info)
  qc <- rna_update_qc(qc, "batches", "correction_applied",  TRUE)
  qc <- rna_update_qc(qc, "batches", "batches_corrected",   nlevels(batch_info))
  rna_msg("Batch correction complete.")
} else {
  rna_msg("Only one batch — correction skipped.")
  qc <- rna_update_qc(qc, "batches", "correction_applied", FALSE)
}

# PCA AFTER batch correction
pca_after <- prcomp(t(rna_mat), scale. = FALSE)
df_after  <- data.frame(
  PC1   = pca_after$x[, 1],
  PC2   = pca_after$x[, 2],
  Batch = batch_info
)

# ==============================================================================
# STEP 11 — FEATURE SELECTION: TOP n VARIABLE GENES + MEAN EXPR FILTER
# ==============================================================================
rna_step(11, sprintf("Feature Selection — Top %d Variable Genes", args[["n-top"]]))
rna_msg("CRITICAL: Feature selection BEFORE Z-scoring (post Z-score all vars ~ 1)")

# Variance calculated on log2 CPM (not Z-scored) — correct order is critical
gene_var <- apply(rna_mat, 1, var)

n_top   <- min(args[["n-top"]], nrow(rna_mat))
top_idx <- order(gene_var, decreasing = TRUE)[1:n_top]
rna_top <- rna_mat[top_idx, , drop = FALSE]

qc <- rna_update_qc(qc, "genes", "selected_top_variable", nrow(rna_top))

# Apply mean expression filter: keep genes with mean log2 CPM >= mean-expr-min
mean_expr <- rowMeans(rna_top)
rna_top   <- rna_top[mean_expr >= args[["mean-expr-min"]], , drop = FALSE]

qc <- rna_update_qc(qc, "genes", "final_after_mean_filter", nrow(rna_top))
rna_msg(sprintf("Selected %d top variable genes", nrow(rna_top)))
rna_msg(sprintf("Variance range     : [%.4f, %.4f]",
                min(gene_var[top_idx]), max(gene_var[top_idx])))
rna_msg(sprintf("Mean expr range    : [%.2f, %.2f]",
                min(mean_expr), max(mean_expr)))

# Report top 10 most variable genes
top10 <- names(sort(gene_var, decreasing = TRUE))[1:10]
rna_msg("Top 10 most variable genes:")
for (i in seq_along(top10)) {
  rna_msg(sprintf("  %2d. %-12s (variance = %.4f)", i, top10[i],
                  gene_var[top10[i]]), level = "DETAIL")
}

# ==============================================================================
# STEP 12 — CREATE TWO OUTPUT VERSIONS (Z-scored + log2 CPM)
# ==============================================================================
rna_step(12, "Creating Output Versions")

# --- Version A: Z-scored per gene (mean=0, sd=1) — for MOFA+/DIABLO ---
# Z-scoring is applied ROW-WISE (per gene) using t(scale(t(...)))
rna_msg("12a. Z-scored version for MOFA+/DIABLO...")
rna_scaled <- t(scale(t(rna_top)))
rna_scaled <- rna_scaled[is.finite(rowSums(rna_scaled)), , drop = FALSE]
rna_scaled <- rna_scaled[complete.cases(rna_scaled), , drop = FALSE]

z_mean <- mean(rna_scaled)
z_sd   <- sd(rna_scaled)
qc <- rna_update_qc(qc, "output_matrices", "zscore_mean",  z_mean)
qc <- rna_update_qc(qc, "output_matrices", "zscore_sd",    z_sd)
qc <- rna_update_qc(qc, "output_matrices", "genes_mofa",   nrow(rna_scaled))
qc <- rna_update_qc(qc, "output_matrices", "samples_mofa", ncol(rna_scaled))
rna_msg(sprintf("Z-score validation : mean = %.5f, sd = %.5f", z_mean, z_sd))
rna_msg(sprintf("MOFA/DIABLO matrix : %d genes x %d samples",
                nrow(rna_scaled), ncol(rna_scaled)))

# --- Version B: log2 CPM — for Machine Learning (preserves biological scale) ---
rna_msg("12b. log2 CPM version for Machine Learning...")
rna_ml <- rna_top   # NO transformation — ML needs absolute expression values

qc <- rna_update_qc(qc, "output_matrices", "ml_range_min", min(rna_ml))
qc <- rna_update_qc(qc, "output_matrices", "ml_range_max", max(rna_ml))
qc <- rna_update_qc(qc, "output_matrices", "genes_ml",     nrow(rna_ml))
qc <- rna_update_qc(qc, "output_matrices", "samples_ml",   ncol(rna_ml))
rna_msg(sprintf("ML matrix          : %d genes x %d samples",
                nrow(rna_ml), ncol(rna_ml)))
rna_msg(sprintf("Expression range   : [%.2f, %.2f]", min(rna_ml), max(rna_ml)))

# Final sample count
qc <- rna_update_qc(qc, "samples", "final", ncol(rna_top))

# ==============================================================================
# STEP 13 — QC PLOTS + STATISTICAL VALIDATION FIGURES
# ==============================================================================
rna_step(13, "Technical Validation Plots")

generate_rna_qc_plots(
  rna_top   = rna_top,
  df_before = df_before,
  df_after  = df_after,
  pca_before= pca_before,
  pca_after = pca_after,
  batch_info= batch_info,
  outdir    = args$outdir
)

generate_rna_validation_figures(
  rna_ml  = rna_ml,
  outdir  = args$outdir
)

# ==============================================================================
# STEP 14 — EXPORT ALL OUTPUTS
# ==============================================================================
rna_step(14, "Saving Outputs")

if (!is.null(metadata)) {
  sample_info <- metadata[match(original_barcodes[colnames(rna_top)], metadata$sample_id), ]
} else {
  sample_info <- data.frame(
    sample_id        = original_barcodes[colnames(rna_top)],
    patient_id       = colnames(rna_top),
    sample_class     = "primary_tumor",
    batch            = as.character(batch_info),
    stringsAsFactors = FALSE
  )
}

export_rna_results(
  rna_scaled  = rna_scaled,
  rna_ml      = rna_ml,
  sample_info = sample_info,
  qc_metrics  = qc,
  outdir      = args$outdir,
  metadata_supplied = !is.null(metadata)
)

# ==============================================================================
# FINAL REPORT
# ==============================================================================
rna_banner("RNA-seq PREPROCESSING COMPLETE")
rna_msg(sprintf("Completion time  : %s", Sys.time()))
rna_msg(sprintf("Samples processed: %d", ncol(rna_top)))
rna_msg(sprintf("Genes selected   : %d (top variable)", nrow(rna_top)))
rna_msg(sprintf("MOFA/DIABLO      : %d genes x %d samples (Z-scored)",
                nrow(rna_scaled), ncol(rna_scaled)))
rna_msg(sprintf("ML               : %d genes x %d samples (log2 CPM)",
                nrow(rna_ml), ncol(rna_ml)))
rna_msg(sprintf("Batches corrected: %d", nlevels(batch_info)))
rna_msg(sprintf("Outliers removed : %d", length(samples_to_remove)))
rna_msg(sprintf("Outputs in       : %s", args$outdir))
