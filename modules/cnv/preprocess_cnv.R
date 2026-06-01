#!/usr/bin/env Rscript
# ==============================================================================
# preprocess_cnv.R — CNV Segment Preprocessing Pipeline
# OmicsFlow | Phase 3: CNV Module
# ==============================================================================
#
# PURPOSE:
#   Executes the scientific steps for filtering, normalizing, and transforming
#   CNV segment data into a gene-level matrix.
#   
# METHODOLOGY:
#   1. Barcode parsing (16-char -> 12-char), separate merged barcodes.
#   2. Standardize chromosomes and check hg38 genome build.
#   3. Build CNV GRanges object.
#   4. Fetch/Load Ensembl gene coordinates (hg38).
#   5. Map segments to genes (Genomic overlap).
#   6. Tumor-only filtering (code 01) and deduplication (1:1 mapping max segments).
#   7. Aggregate gene-level matrix preserving CNV sign (using data.table).
#   8. Missing value imputation (0).
#   9. Feature selection (Top 5000 variable genes) + Post-processing.
#
# ==============================================================================
# USAGE:
#   Rscript preprocess_cnv.R \
#     --input      data/cnv_segment_raw.rds \
#     --cache      data/gene_coords_hg38.rds \
#     --outdir     results/cnv/ \
#     --ntop       5000
# ==============================================================================

suppressPackageStartupMessages(library(optparse))

# ------------------------------------------------------------------------------
# COMMAND-LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "data/cnv_segment_raw.rds",
              help = "Path to raw CNV segment data (.rds) [default= %default]"),
  make_option(c("-c", "--cache"), type = "character", default = "data/gene_coords_hg38.rds",
              help = "Path to Ensembl gene coordinates cache (.rds) [default= %default]"),
  make_option(c("-o", "--outdir"), type = "character", default = "results/cnv/",
              help = "Output directory for CNV results [default= %default]"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Path to sample metadata CSV file (optional)"),
  make_option(c("-n", "--ntop"), type = "integer", default = 5000,
              help = "Number of top variable genes to retain [default= %default]"),
  make_option(c("-s", "--seed"), type = "integer", default = 42,
              help = "Random seed for reproducibility [default= %default]")
)

opt_parser <- OptionParser(
  usage = "Usage: %prog [options]",
  option_list = option_list,
  description = "OmicsFlow Phase 3: CNV Segment Preprocessing Module"
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
  error = function(e) "modules/cnv"
)
source(file.path(script_dir, "utils_cnv.R"))
source(file.path(script_dir, "qc_cnv.R"))
source(file.path(script_dir, "export_cnv.R"))
source(file.path(dirname(script_dir), "utils_metadata.R"))

# ------------------------------------------------------------------------------
# INITIALISE
# ------------------------------------------------------------------------------
tryCatch({
  load_cnv_packages()
  set.seed(opt$seed)
  
  if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
  }
  
  qc_metrics <- init_cnv_qc()
  
  metadata <- load_metadata(opt$metadata)
  
  cnv_banner("OmicsFlow v2.0.1 | CNV Preprocessing")
  cnv_msg(sprintf("Start time : %s", Sys.time()))
  cnv_msg(sprintf("Input      : %s", opt$input))
  cnv_msg(sprintf("Output dir : %s", opt$outdir))
  if (!is.null(metadata)) {
    cnv_msg(sprintf("Metadata   : %s (%d samples)", opt$metadata, nrow(metadata)))
  }
  
  # ==============================================================================
  # 1. DATA LOADING & BARCODE PARSING
  # ==============================================================================
  cnv_step(1, "Data Loading & Barcode Parsing")
  cnv_msg("Loading inputs...", level = "DETAILS")
  
  data_list <- load_cnv_data(opt$input, opt$cache)
  cnv_data <- data_list$cnv
  gene_coords_cache <- data_list$gene_coords
  
  n_initial_segments <- nrow(cnv_data)
  qc_metrics <- add_cnv_qc(qc_metrics, "samples", "initial_segments", n_initial_segments)
  
  n_before_sep <- nrow(cnv_data)
  cnv_data <- cnv_data %>%
    tidyr::separate_rows(Sample, sep = ";") %>%
    as.data.frame()
  cnv_data$Sample <- trimws(cnv_data$Sample)
  n_after_sep <- nrow(cnv_data)
  
  if (n_after_sep > n_before_sep) {
    cnv_msg(sprintf("Separated merged barcodes: %s -> %s segments",
                    format(n_before_sep, big.mark = ","),
                    format(n_after_sep, big.mark = ",")))
  } else {
    cnv_msg("No merged barcodes found")
  }
  
  if (!is.null(metadata)) {
    cnv_data <- cnv_data[cnv_data$Sample %in% metadata$sample_id, , drop = FALSE]
    if (nrow(cnv_data) == 0) {
      stop("None of the raw CNV segments match the sample_id column in metadata.")
    }
  }
  
  cnv_data$patient_id  <- get_patient_id(cnv_data$Sample, metadata)
  cnv_data$sample_type <- get_sample_class(cnv_data$Sample, metadata)
  
  n_patients <- length(unique(cnv_data$patient_id))
  n_samples  <- length(unique(cnv_data$Sample))
  cnv_msg(sprintf("Parsed barcodes: %d unique patients, %d unique samples", n_patients, n_samples))
  qc_metrics <- add_cnv_qc(qc_metrics, "samples", "initial_patients", n_patients)
  qc_metrics <- add_cnv_qc(qc_metrics, "samples", "initial_samples", n_samples)
  
  # Genome build validation
  max_coord <- max(cnv_data$End, na.rm = TRUE)
  if (max_coord < 240e6) {
    cnv_msg(sprintf("GENOME BUILD: Max coordinate = %s bp - verify hg38 compatibility", format(max_coord, big.mark = ",")), level = "WARN")
  } else {
    cnv_msg(sprintf("Genome build: max coordinate = %s bp (consistent with hg38)", format(max_coord, big.mark = ",")))
  }
  
  # ==============================================================================
  # 2. CHROMOSOME STANDARDIZATION
  # ==============================================================================
  cnv_step(2, "Chromosome Name Standardization")
  
  has_chr_prefix <- grepl("^chr", cnv_data$Chromosome)
  cnv_data$Chromosome <- ifelse(
    has_chr_prefix,
    cnv_data$Chromosome,
    paste0("chr", cnv_data$Chromosome)
  )
  cnv_data$Chromosome <- gsub("^chrchr", "chr", cnv_data$Chromosome)
  
  valid_chrs <- c(paste0("chr", 1:22), "chrX", "chrY")
  n_before_chr_filter <- nrow(cnv_data)
  cnv_data <- cnv_data[cnv_data$Chromosome %in% valid_chrs, ]
  n_removed_chr <- n_before_chr_filter - nrow(cnv_data)
  
  chr_final <- sort(unique(cnv_data$Chromosome))
  cnv_msg(sprintf("Standardized chromosomes (%d): %s", length(chr_final), paste(chr_final, collapse = ", ")))
  if (n_removed_chr > 0) {
    cnv_msg(sprintf("Removed %s segments from non-standard contigs", format(n_removed_chr, big.mark = ",")), level = "WARN")
  }
  
  # ==============================================================================
  # 3. BUILD CNV GRANGES OBJECT
  # ==============================================================================
  cnv_step(3, "Building CNV GRanges Object")
  
  cnv_gr <- GRanges(
    seqnames     = cnv_data$Chromosome,
    ranges       = IRanges(start = cnv_data$Start, end = cnv_data$End),
    segment_mean = cnv_data$Segment_Mean,
    num_probes   = cnv_data$Num_Probes,
    sample_id    = cnv_data$Sample,
    patient_id   = cnv_data$patient_id
  )
  cnv_msg(sprintf("GRanges built: %s segments", format(length(cnv_gr), big.mark = ",")))
  
  # ==============================================================================
  # 4. ENSEMBL GENE COORDINATES (hg38)
  # ==============================================================================
  cnv_step(4, "Retrieving Gene Coordinates")
  
  cache_file <- opt$cache
  if (is.null(gene_coords_cache) || (!is.data.frame(gene_coords_cache) && !file.exists(gene_coords_cache))) {
    cnv_msg("Downloading from Ensembl (1-3 minutes)...", level = "INFO")
    mart <- useEnsembl(
      biomart = "genes",
      dataset = "hsapiens_gene_ensembl",
      mirror  = "useast"
    )
    gene_coords <- getBM(
      attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
      filters    = "chromosome_name",
      values     = c(1:22, "X", "Y"),
      mart       = mart
    )
    gene_coords <- gene_coords[gene_coords$hgnc_symbol != "", ]
    
    if (is.character(cache_file)) {
      saveRDS(gene_coords, cache_file)
      cnv_msg(sprintf("Cached to: %s", cache_file))
    }
  } else if (is.data.frame(gene_coords_cache)) {
    gene_coords <- gene_coords_cache
    cnv_msg(sprintf("Using provided gene coordinates: %d genes", nrow(gene_coords)))
  } else {
    gene_coords <- readRDS(cache_file)
    cnv_msg(sprintf("Loaded %s genes from cache", format(nrow(gene_coords), big.mark = ",")))
  }
  
  # ==============================================================================
  # 5. BUILD GENE GRANGES OBJECT
  # ==============================================================================
  cnv_step(5, "Building Gene GRanges Object")
  
  genes_gr <- GRanges(
    seqnames    = paste0("chr", gene_coords$chromosome_name),
    ranges      = IRanges(start = gene_coords$start_position, end = gene_coords$end_position),
    gene_symbol = gene_coords$hgnc_symbol
  )
  cnv_msg(sprintf("Gene GRanges built: %s genes", format(length(genes_gr), big.mark = ",")))
  
  # ==============================================================================
  # 6. MAP SEGMENTS TO GENES VIA OVERLAP
  # ==============================================================================
  cnv_step(6, "Mapping Segments -> Genes")
  
  overlaps <- findOverlaps(query = genes_gr, subject = cnv_gr, type = "any")
  
  matched_data <- data.frame(
    gene_symbol  = genes_gr$gene_symbol[queryHits(overlaps)],
    sample_id    = cnv_gr$sample_id[subjectHits(overlaps)],
    patient_id   = cnv_gr$patient_id[subjectHits(overlaps)],
    cnv_value    = cnv_gr$segment_mean[subjectHits(overlaps)],
    num_probes   = cnv_gr$num_probes[subjectHits(overlaps)],
    stringsAsFactors = FALSE
  )
  
  n_pairs <- nrow(matched_data)
  n_genes_mapped <- length(unique(matched_data$gene_symbol))
  cnv_msg(sprintf("Overlaps found: %s gene-segment pairs", format(n_pairs, big.mark = ",")))
  qc_metrics <- add_cnv_qc(qc_metrics, "genes", "genes_mapped", n_genes_mapped)
  
  rm(overlaps, genes_gr, cnv_gr)
  gc()
  
  # ==============================================================================
  # 7. AGGREGATE TO GENE-LEVEL MATRIX (TUMOR ONLY, 1:1)
  # ==============================================================================
  cnv_step(7, "Aggregation - Gene x Tumor Matrix (1:1)")
  
  if (!is.null(metadata)) {
    tumor_sample_ids <- unique(cnv_data$Sample)
  } else {
    tumor_sample_ids <- unique(cnv_data$Sample[cnv_data$sample_type == "01"])
  }
  matched_data_tumor <- matched_data %>% filter(sample_id %in% tumor_sample_ids)
  
  tumor_samples_per_patient <- matched_data_tumor %>%
    distinct(patient_id, sample_id) %>%
    group_by(patient_id) %>%
    summarise(n_tumor_samples = n(), .groups = "drop")
  
  patients_with_multiple <- tumor_samples_per_patient %>% filter(n_tumor_samples > 1)
  
  if (nrow(patients_with_multiple) > 0) {
    cnv_msg(sprintf("%d patients with multiple tumors - resolving...", nrow(patients_with_multiple)), level = "WARN")
    segments_per_sample <- matched_data_tumor %>%
      group_by(sample_id) %>%
      summarise(n_segments = n(), .groups = "drop")
    
    best_tumor_per_patient <- matched_data_tumor %>%
      distinct(patient_id, sample_id) %>%
      left_join(segments_per_sample, by = "sample_id") %>%
      group_by(patient_id) %>%
      slice_max(order_by = n_segments, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      pull(sample_id)
    
    matched_data <- matched_data_tumor %>% filter(sample_id %in% best_tumor_per_patient)
  } else {
    matched_data <- matched_data_tumor
  }
  
  n_final_samples  <- length(unique(matched_data$sample_id))
  n_final_patients <- length(unique(matched_data$patient_id))
  cnv_msg(sprintf("Filtered to %d tumor samples, %d unique patients (1:1 mapping)", n_final_samples, n_final_patients))
  qc_metrics <- add_cnv_qc(qc_metrics, "samples", "final", n_final_patients)
  
  cnv_msg("Aggregating gene-segment pairs (data.table method)...", level = "DETAILS")
  matched_dt <- as.data.table(matched_data)
  gene_sample_cnv <- matched_dt[, {
    best_idx <- which.max(abs(cnv_value))
    list(cnv = cnv_value[best_idx], n_segments = .N)
  }, by = .(gene_symbol, sample_id)]
  gene_sample_cnv <- as.data.frame(gene_sample_cnv)
  
  cnv_msg("Reshaping to matrix...", level = "DETAILS")
  cnv_matrix <- gene_sample_cnv[, c("gene_symbol", "sample_id", "cnv")] %>%
    pivot_wider(names_from = sample_id, values_from = cnv) %>%
    column_to_rownames("gene_symbol") %>%
    as.matrix()
  
  cols_with_semicolon <- grep(";", colnames(cnv_matrix))
  if (length(cols_with_semicolon) > 0) {
    cnv_matrix <- cnv_matrix[, -cols_with_semicolon, drop = FALSE]
  }
  if (is.null(metadata)) {
    tumor_cols <- grep("-01", colnames(cnv_matrix))
    if(length(tumor_cols) > 0) {
      cnv_matrix <- cnv_matrix[, tumor_cols, drop = FALSE]
    }
    # Standardize colnames to 12-char patient ID
    colnames(cnv_matrix) <- substr(colnames(cnv_matrix), 1, 12)
  } else {
    # Map colnames to patient IDs from metadata
    colnames(cnv_matrix) <- get_patient_id(colnames(cnv_matrix), metadata)
  }
  cnv_msg(sprintf("Final tumor-only matrix: %d genes x %d samples", nrow(cnv_matrix), ncol(cnv_matrix)))
  
  rm(matched_data, matched_dt, matched_data_tumor, gene_sample_cnv, tumor_sample_ids, tumor_samples_per_patient, patients_with_multiple)
  gc()
  
  # ==============================================================================
  # 8. MISSING VALUE IMPUTATION
  # ==============================================================================
  cnv_step(8, "Missing Value Imputation")
  
  n_missing <- sum(is.na(cnv_matrix))
  pct_missing <- 100 * n_missing / length(cnv_matrix)
  cnv_msg(sprintf("Missing values: %s (%.2f%% of matrix)", format(n_missing, big.mark = ","), pct_missing))
  
  cnv_matrix[is.na(cnv_matrix)] <- 0
  cnv_msg("Imputed with 0 (neutral copy number)")
  
  # ==============================================================================
  # 9. VARIANCE FEATURE SELECTION & POST-PROCESSING
  # ==============================================================================
  cnv_step(9, "Feature Selection & Post-Processing")
  
  gene_var <- rowVars(cnv_matrix)
  names(gene_var) <- rownames(cnv_matrix)
  
  top_n_param <- opt$ntop
  top_n <- min(top_n_param, length(gene_var))
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:top_n]
  cnv_final <- cnv_matrix[top_genes, , drop = FALSE]
  
  cnv_msg(sprintf("Selected top %d most variable genes", top_n))
  qc_metrics <- add_cnv_qc(qc_metrics, "parameters", "n_top_genes", top_n_param)
  qc_metrics <- add_cnv_qc(qc_metrics, "genes", "final_genes_before_immune_filter", nrow(cnv_final))
  
  # Remove immune genes
  immune_genes <- grep("^IG[HKL][VDJC]|^TR[ABDG][VDJC]", rownames(cnv_final), value = TRUE)
  if (length(immune_genes) > 0) {
    cnv_msg(sprintf("Removing %d immune receptor genes...", length(immune_genes)), level = "DETAILS")
    cnv_final <- cnv_final[!rownames(cnv_final) %in% immune_genes, , drop = FALSE]
  }
  
  # Clip extremes
  n_below <- sum(cnv_final < -5)
  n_above <- sum(cnv_final > 5)
  if (n_below > 0 || n_above > 0) {
    cnv_final[cnv_final < -5] <- -5
    cnv_final[cnv_final > 5]  <-  5
    cnv_msg(sprintf("Clipped extremes: %d below -5, %d above 5", n_below, n_above), level = "DETAILS")
  }
  
  qc_metrics <- add_cnv_qc(qc_metrics, "genes", "final_genes", nrow(cnv_final))
  cnv_msg(sprintf("Final processed CNV matrix: %d genes x %d samples", nrow(cnv_final), ncol(cnv_final)))
  
  # ==============================================================================
  # 10. GENERATE QC PLOTS
  # ==============================================================================
  generate_cnv_qc_plots(cnv_final, opt$outdir)
  
  # ==============================================================================
  # 11. EXPORT RESULTS
  # ==============================================================================
  gene_var_df <- data.frame(
    gene     = names(gene_var),
    variance = gene_var,
    selected = names(gene_var) %in% top_genes,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(metadata)) {
    sample_info <- metadata[match(colnames(cnv_final), metadata$patient_id), ]
  } else {
    sample_info <- data.frame(
      sample_id        = colnames(cnv_final),
      patient_id       = colnames(cnv_final),
      sample_class     = "primary_tumor",
      batch            = "batch1",
      stringsAsFactors = FALSE
    )
  }
  export_cnv_results(cnv_final, gene_var_df, opt$outdir, sample_info)
  
  # ==============================================================================
  # 12. SAVE QC METRICS
  # ==============================================================================
  export_cnv_qc(qc_metrics, opt$outdir)
  
  cnv_banner("CNV MODULE COMPLETE \u2714")
  
}, error = function(e) {
  cnv_msg(sprintf("FATAL ERROR: %s", e$message), level = "ERROR")
  quit(status = 1)
})
