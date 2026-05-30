#!/usr/bin/env Rscript
# ==============================================================================
# preprocess_meth.R — DNA Methylation Preprocessing Pipeline
# OmicsFlow | Phase 1: Methylation Module
# Scientific source: data/full_scripts.R (lines 630-1185)
# ==============================================================================
# USAGE:
#   Rscript preprocess_meth.R \
#     --input       path/to/methylation_beta_raw.rds \
#     --outdir      results/methylation/ \
#     --clinical    path/to/clinical_data.tsv \
#     --det-pval    path/to/methylation_detection_pvals.rds \
#     --cross-react path/to/Chen_2013_cross_reactive_probes.csv
# ==============================================================================

suppressPackageStartupMessages(library(optparse))

# ------------------------------------------------------------------------------
# CLI ARGUMENT PARSING
# ------------------------------------------------------------------------------
option_list <- list(
  make_option("--input",        type="character", default=NULL,
              help="Path to input Methylation Beta-value RDS"),
  make_option("--outdir",       type="character", default="results/methylation/",
              help="Output directory [default: results/methylation/]"),
  make_option("--metadata",     type="character", default=NULL,
              help="Path to sample metadata CSV file (optional)"),
  make_option("--clinical",     type="character", default=NULL,
              help="Path to clinical data TSV (optional, for ComBat covariate protection)"),
  make_option("--det-pval",     type="character", default=NULL,
              help="Path to detection p-values RDS (optional)"),
  make_option("--cross-react",  type="character", default=NULL,
              help="Path to cross-reactive probes CSV (optional, will download if missing)"),
  make_option("--na-threshold", type="double",    default=0.20,
              help="Max NA fraction per probe before removal [default: 0.20]"),
  make_option("--n-top",        type="integer",   default=5000L,
              help="Number of top variable probes to select for MOFA [default: 5000]"),
  make_option("--knn-k",        type="integer",   default=10L,
              help="k value for KNN imputation [default: 10]"),
  make_option("--clinical_map", type="character", default=NULL,
              help="Clinical column mapping (JSON file/string) for custom cohorts [optional]")
)

parser <- OptionParser(option_list = option_list,
                       description = "OmicsFlow DNA Methylation Preprocessing Pipeline")
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
  error = function(e) "modules/methylation"
)
source(file.path(script_dir, "utils_meth.R"))
source(file.path(script_dir, "qc_meth.R"))
source(file.path(script_dir, "export_meth.R"))
source(file.path(dirname(script_dir), "utils_metadata.R"))
source(file.path(dirname(script_dir), "utils_clinical.R"))

# ------------------------------------------------------------------------------
# INITIALISE
# ------------------------------------------------------------------------------
meth_load_packages()
set.seed(42)
meth_ensure_dirs(args$outdir)
qc <- meth_init_qc()
probe_log <- meth_init_probe_log()

metadata <- load_metadata(args$metadata)

epsilon <- 1e-7

meth_banner("OmicsFlow v1.0.0 | DNA Methylation Preprocessing")
meth_msg(sprintf("Start time : %s", Sys.time()))
meth_msg(sprintf("Input      : %s", args$input))
meth_msg(sprintf("Output dir : %s", args$outdir))
if (!is.null(metadata)) {
  meth_msg(sprintf("Metadata   : %s (%d samples)", args$metadata, nrow(metadata)))
}

# ==============================================================================
# STEP 1 — DATA LOADING & DEDUPLICATION
# ==============================================================================
meth_step(1, "Loading & Quality-Aware Deduplication")
meth_check_file(args$input, "Methylation input")

ext <- tolower(tools::file_ext(args$input))
if (ext == "rds") {
  met_beta <- readRDS(args$input)
  if (inherits(met_beta, "SummarizedExperiment")) met_beta <- assay(met_beta)
} else if (ext == "csv") {
  met_beta <- as.matrix(read.csv(args$input, row.names = 1, check.names = FALSE))
} else if (ext %in% c("tsv", "txt")) {
  met_beta <- as.matrix(read.delim(args$input, row.names = 1, check.names = FALSE))
} else {
  stop("Unsupported input format.\n\nSupported formats:\n- .rds\n- .csv\n- .tsv\n- .txt\n\nRecommended format:\n- .rds", call. = FALSE)
}
storage.mode(met_beta) <- "double"

meth_msg(sprintf("Raw matrix: %s probes x %s samples",
                 format(nrow(met_beta), big.mark=","),
                 format(ncol(met_beta), big.mark=",")))

qc <- meth_update_qc(qc, "samples", "raw", ncol(met_beta))
qc <- meth_update_qc(qc, "probes", "raw", nrow(met_beta))

# Keep primary tumor or match to metadata
if (!is.null(metadata)) {
  common_samples <- intersect(colnames(met_beta), metadata$sample_id)
  if (length(common_samples) == 0) {
    stop("None of the raw methylation matrix column names match the sample_id column in metadata.")
  }
  met_beta <- met_beta[, common_samples, drop = FALSE]
  qc <- meth_update_qc(qc, "samples", "primary_tumor", ncol(met_beta))
  meth_msg(sprintf("Samples matched to metadata: %d", ncol(met_beta)))
  
  # Deduplicate: keep sample with fewest NAs per patient
  na_counts   <- colSums(is.na(met_beta))
  met_beta    <- met_beta[, order(na_counts), drop = FALSE]
  patient_ids <- get_patient_id(colnames(met_beta), metadata)
  met_beta    <- met_beta[, !duplicated(patient_ids), drop = FALSE]
  
  # Rename to patient IDs
  colnames(met_beta) <- get_patient_id(colnames(met_beta), metadata)
  meth_msg(sprintf("After deduplication: %d unique patients", ncol(met_beta)))
  qc <- meth_update_qc(qc, "samples", "after_deduplication", ncol(met_beta))
} else {
  sample_type <- meth_extract_sample_type(colnames(met_beta))
  met_beta    <- met_beta[, !is.na(sample_type) & sample_type == "01", drop = FALSE]
  meth_msg(sprintf("Primary tumor samples: %d", ncol(met_beta)))
  qc <- meth_update_qc(qc, "samples", "primary_tumor", ncol(met_beta))
  
  # Deduplicate: keep sample with fewest NAs per patient
  na_counts   <- colSums(is.na(met_beta))
  met_beta    <- met_beta[, order(na_counts), drop = FALSE]
  patient_ids <- meth_extract_patient_id(colnames(met_beta))
  met_beta    <- met_beta[, !duplicated(patient_ids), drop = FALSE]
  
  # Rename to 12-char patient IDs
  colnames(met_beta) <- meth_extract_patient_id(colnames(met_beta))
  meth_msg(sprintf("After deduplication: %d unique patients", ncol(met_beta)))
  qc <- meth_update_qc(qc, "samples", "after_deduplication", ncol(met_beta))
}

probe_log <- meth_log_probes(probe_log, "After Dedup", met_beta)
qc <- meth_update_qc(qc, "probes", "after_dedup", nrow(met_beta))

rm(sample_type, na_counts, patient_ids)
gc()

# ==============================================================================
# STEP 2 — DETECTION P-VALUE FILTERING
# ==============================================================================
meth_step(2, "Detection P-value Filtering")

if (!is.null(args[["det-pval"]]) && file.exists(args[["det-pval"]])) {
  detP <- readRDS(args[["det-pval"]])
  if (ncol(detP) != ncol(met_beta)) {
    colnames(detP) <- meth_extract_patient_id(colnames(detP))
  }
  common_s <- intersect(colnames(met_beta), colnames(detP))
  detP     <- detP[rownames(met_beta), common_s, drop = FALSE]
  
  # Remove probes with DetP > 0.01 in > 20% of samples
  high     <- rowSums(detP > 0.01, na.rm = TRUE) / ncol(detP) > 0.20
  met_beta <- met_beta[!high, , drop = FALSE]
  
  meth_msg(sprintf("Removed %d probes with poor detection P-values", sum(high)))
  qc <- meth_update_qc(qc, "filters", "detection_pval_removed", sum(high))
  
  probe_log <- meth_log_probes(probe_log, "After DetP filter", met_beta)
  rm(detP, common_s, high)
} else {
  meth_msg("Detection P-value file not provided or not found — skipping.")
  qc <- meth_update_qc(qc, "filters", "detection_pval_removed", 0)
}
qc <- meth_update_qc(qc, "probes", "after_detection_pval", nrow(met_beta))
gc()

# ==============================================================================
# STEP 3 — PROBE FILTERING
# ==============================================================================
meth_step(3, "Probe Filtering")

anno   <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
shared <- intersect(rownames(met_beta), rownames(anno))
met_beta <- met_beta[shared, , drop = FALSE]
anno     <- anno[shared, , drop = FALSE]

# 3.1 Cross-reactive probes
cross_reac <- character(0)
cross_react_path <- args[["cross-react"]]

# Try to resolve bundled path if not provided
if (is.null(cross_react_path) || !nzchar(cross_react_path) || !file.exists(cross_react_path)) {
  bundled_path <- system.file("configs", "cross_reactive_probes.csv", package = "OmicsFlow")
  if (nzchar(bundled_path) && file.exists(bundled_path)) {
    cross_react_path <- bundled_path
  }
}

if (!is.null(cross_react_path) && nzchar(cross_react_path) && file.exists(cross_react_path)) {
  cr <- read.csv(cross_react_path, stringsAsFactors = FALSE)
  col <- intersect(c("TargetID", "Probe_ID", "probe"), names(cr))
  cross_reac <- if (length(col) > 0) as.character(cr[[col[1]]]) else as.character(cr[[1]])
  meth_msg(sprintf("Cross-reactive probe list loaded: %s probes (from %s)",
                   format(length(cross_reac), big.mark = ","), cross_react_path))
} else {
  tryCatch({
    cross_url <- "https://raw.githubusercontent.com/sirselim/illumina450k_filtering/master/Chen_2013_cross_reactive_probes.csv"
    cr <- read.csv(cross_url, stringsAsFactors = FALSE)
    col <- intersect(c("TargetID", "Probe_ID", "probe"), names(cr))
    cross_reac <- if (length(col) > 0) as.character(cr[[col[1]]]) else as.character(cr[[1]])
    meth_msg(sprintf("Cross-reactive probe list downloaded: %s probes",
                     format(length(cross_reac), big.mark = ",")))
    if (!is.null(args[["cross-react"]])) {
      write.csv(cr, args[["cross-react"]], row.names = FALSE)
    }
  }, error = function(e) meth_msg("Cross-reactive download failed — skipping.", level = "WARN"))
}
# Record QC metric AFTER cross_reac is definitively populated so the count is accurate.
# If cross_reac is still empty here, the list could not be loaded — log a clear warning.
n_cross_reactive <- sum(rownames(met_beta) %in% cross_reac)
if (length(cross_reac) == 0L) {
  meth_msg(paste0("cross_reactive_removed recorded as 0: Chen 2013 list unavailable. ",
                  "Re-run with --cross-react to obtain the correct count."),
           level = "WARN")
} else {
  meth_msg(sprintf("Cross-reactive probes to remove: %s (in %s-probe matrix)",
                   format(n_cross_reactive, big.mark = ","),
                   format(nrow(met_beta), big.mark = ",")))
}
qc <- meth_update_qc(qc, "filters", "cross_reactive_removed", n_cross_reactive)

# 3.2 SNP-associated probes
probe_rs  <- meth_safe_get_col(anno, c("Probe_rs",  "probe_rs"))
probe_maf <- meth_safe_get_col(anno, c("Probe_maf", "probe_maf"))
cpg_rs    <- meth_safe_get_col(anno, c("CpG_rs",    "cpg_rs"))
cpg_maf   <- meth_safe_get_col(anno, c("CpG_maf",   "cpg_maf"))
snp_mask  <- (!is.na(probe_rs)  & probe_rs  != "") |
             (!is.na(probe_maf) & probe_maf > 0.01) |
             (!is.na(cpg_rs)    & cpg_rs    != "") |
             (!is.na(cpg_maf)   & cpg_maf   > 0.01)
snp_probes <- rownames(anno)[snp_mask]
qc <- meth_update_qc(qc, "filters", "snp_probes_removed", sum(rownames(met_beta) %in% snp_probes))

# 3.3 Sex chromosome probes
chr_col    <- meth_safe_get_col(anno, c("chr", "Chromosome", "CHR"))
sex_probes <- rownames(anno)[chr_col %in% c("chrX", "chrY", "X", "Y")]
qc <- meth_update_qc(qc, "filters", "sex_chromosome_removed", sum(rownames(met_beta) %in% sex_probes))

# 3.4 Non-CpG probes
non_cpg <- grep("^ch\\.", rownames(met_beta), value = TRUE)
qc <- meth_update_qc(qc, "filters", "non_cpg_removed", sum(rownames(met_beta) %in% non_cpg))

# Remove all bad probes
bad      <- unique(c(cross_reac, snp_probes, sex_probes, non_cpg))
met_beta <- met_beta[!rownames(met_beta) %in% bad, , drop = FALSE]
meth_msg(sprintf("Total removed: %d probes", length(bad)))
qc <- meth_update_qc(qc, "filters", "total_bad_probes", length(bad))

probe_log <- meth_log_probes(probe_log, "After Probe filter", met_beta)
qc <- meth_update_qc(qc, "probes", "after_probe_filter", nrow(met_beta))

rm(anno, shared, cross_reac, probe_rs, probe_maf, cpg_rs, cpg_maf,
   snp_mask, chr_col, sex_probes, non_cpg, bad, snp_probes)
gc()

# ==============================================================================
# STEP 4 — NA FILTERING, M-VALUE TRANSFORMATION & IMPUTATION
# ==============================================================================
meth_step(4, "NA Filtering, Transformation & Imputation")

na_rate  <- rowMeans(is.na(met_beta))
na_removed <- sum(na_rate >= args[["na-threshold"]])
met_beta <- met_beta[na_rate < args[["na-threshold"]], , drop = FALSE]
meth_msg(sprintf("Removed %d probes with >%.0f%% NA", na_removed, args[["na-threshold"]] * 100))
qc <- meth_update_qc(qc, "filters", "na_fraction_removed", na_removed)

probe_log <- meth_log_probes(probe_log, "After NA filter", met_beta)
qc <- meth_update_qc(qc, "probes", "after_na_filter", nrow(met_beta))

met_beta <- pmax(pmin(met_beta, 1 - epsilon), epsilon)
meth_m   <- log2(met_beta / (1 - met_beta))

meth_msg(sprintf("M-value range: [%.3f, %.3f]", min(meth_m, na.rm = TRUE), max(meth_m, na.rm = TRUE)))
qc <- meth_update_qc(qc, "m_value_stats", "global_min", min(meth_m, na.rm = TRUE))
qc <- meth_update_qc(qc, "m_value_stats", "global_max", max(meth_m, na.rm = TRUE))
qc <- meth_update_qc(qc, "m_value_stats", "global_mean", mean(meth_m, na.rm = TRUE))

if (anyNA(meth_m)) {
  miss_frac <- sum(is.na(meth_m)) / length(meth_m)
  qc <- meth_update_qc(qc, "imputation", "missing_fraction", miss_frac)
  qc <- meth_update_qc(qc, "imputation", "probes_imputed", sum(rowSums(is.na(meth_m)) > 0))
  
  if (miss_frac < 0.01) {
    meth_msg("Median imputation (missing < 1%)...")
    qc <- meth_update_qc(qc, "imputation", "method_used", "Median")
    row_med <- rowMedians(meth_m, na.rm = TRUE)
    na_idx  <- which(is.na(meth_m), arr.ind = TRUE)
    meth_m[na_idx] <- row_med[na_idx[, 1]]
  } else {
    meth_msg(sprintf("KNN imputation (k=%d)...", args[["knn-k"]]))
    qc <- meth_update_qc(qc, "imputation", "method_used", "KNN")
    meth_m <- impute.knn(as.matrix(meth_m), k = args[["knn-k"]])$data
  }
  met_beta <- 2^meth_m / (2^meth_m + 1)
  meth_msg("Imputation complete.")
} else {
  meth_msg("No missing values.")
  qc <- meth_update_qc(qc, "imputation", "method_used", "None")
  qc <- meth_update_qc(qc, "imputation", "missing_fraction", 0)
  qc <- meth_update_qc(qc, "imputation", "probes_imputed", 0)
}
rm(na_rate)
if (exists("na_idx")) rm(na_idx, row_med)
if (exists("miss_frac")) rm(miss_frac)
gc()

# ==============================================================================
# STEP 5 — BATCH CORRECTION WITH COMBAT
# ==============================================================================
meth_step(5, "Batch Correction (ComBat)")

# Recover full barcodes for batch extraction
# Need to redo deduplication to find original barcodes
if (ext == "rds") {
  met_temp <- readRDS(args$input)
  if (inherits(met_temp, "SummarizedExperiment")) met_temp <- assay(met_temp)
} else if (ext == "csv") {
  met_temp <- as.matrix(read.csv(args$input, row.names = 1, check.names = FALSE))
} else if (ext %in% c("tsv", "txt")) {
  met_temp <- as.matrix(read.delim(args$input, row.names = 1, check.names = FALSE))
}

if (!is.null(metadata)) {
  common_temp <- intersect(colnames(met_temp), metadata$sample_id)
  met_temp <- met_temp[, common_temp, drop = FALSE]
  
  na_counts_temp   <- colSums(is.na(met_temp))
  patient_ids_temp <- get_patient_id(colnames(met_temp), metadata)
  keep_idx <- tapply(seq_len(ncol(met_temp)), patient_ids_temp,
                     function(idx) idx[which.min(na_counts_temp[idx])])
  keep_idx <- sort(unlist(keep_idx))
  met_temp <- met_temp[, keep_idx, drop = FALSE]
  
  orig_barcodes <- colnames(met_temp)
  names(orig_barcodes) <- get_patient_id(orig_barcodes, metadata)
} else {
  st <- meth_extract_sample_type(colnames(met_temp))
  met_temp <- met_temp[, !is.na(st) & st == "01", drop = FALSE]
  
  na_counts_temp   <- colSums(is.na(met_temp))
  patient_ids_temp <- meth_extract_patient_id(colnames(met_temp))
  keep_idx <- tapply(seq_len(ncol(met_temp)), patient_ids_temp,
                     function(idx) idx[which.min(na_counts_temp[idx])])
  keep_idx <- sort(unlist(keep_idx))
  met_temp <- met_temp[, keep_idx, drop = FALSE]
  
  orig_barcodes <- colnames(met_temp)
  names(orig_barcodes) <- meth_extract_patient_id(orig_barcodes)
}

rm(met_temp)
if (exists("st")) rm(st)
rm(na_counts_temp, patient_ids_temp, keep_idx)
meth_msg(sprintf("Original barcodes recovered: %d", length(orig_barcodes)))

# Match to current samples
patient_current <- colnames(meth_m)
orig_barcodes_matched <- orig_barcodes[patient_current]

# Extract batch
batch_vec        <- get_batch(orig_barcodes_matched, metadata, omics_type = "meth")
names(batch_vec) <- patient_current
batch_vec        <- factor(batch_vec)
if (!is.null(metadata)) {
  meth_msg(sprintf("Metadata batches: %d", nlevels(batch_vec)))
} else {
  meth_msg(sprintf("Plate batches: %d", nlevels(batch_vec)))
}
qc <- meth_update_qc(qc, "batches", "total_detected", nlevels(batch_vec))

# Remove missing batch
valid     <- !is.na(batch_vec) & batch_vec != ""
meth_m    <- meth_m[, valid, drop = FALSE]
met_beta  <- met_beta[, valid, drop = FALSE]
batch_vec <- batch_vec[valid]

# Remove singletons
tab        <- table(batch_vec)
singletons <- names(tab[tab == 1])
qc <- meth_update_qc(qc, "batches", "singleton_batches", length(singletons))

if (length(singletons) > 0) {
  keep      <- !batch_vec %in% singletons
  meth_m    <- meth_m[, keep, drop = FALSE]
  met_beta  <- met_beta[, keep, drop = FALSE]
  batch_vec <- droplevels(batch_vec[keep])
  meth_msg(sprintf("Removed %d singleton batch(es)", length(singletons)))
}
qc <- meth_update_qc(qc, "samples", "after_batch_na_removal", ncol(meth_m))
qc <- meth_update_qc(qc, "samples", "after_singleton_removal", ncol(meth_m))

# Build biological model
tss <- if (is.null(metadata)) {
  meth_extract_tss(orig_barcodes_matched)
} else {
  get_center(orig_barcodes_matched, metadata)
}

final_covs <- character(0)
model_data <- data.frame(row.names = colnames(meth_m))

if (!is.null(tss)) {
  tss <- tss[match(colnames(meth_m), patient_current)]
  if (length(unique(na.omit(tss))) > 1) {
    model_data$TSS <- factor(tss)
    final_covs <- c(final_covs, "TSS")
  }
}

if (!is.null(args[["clinical"]]) && file.exists(args[["clinical"]])) {
  tryCatch({
    # Load clinical data via the universal abstraction layer
    clinical_std <- load_clinical_data(
      file       = args[["clinical"]],
      column_map = args[["clinical_map"]],
      metadata   = metadata
    )
    
    # Deduplicate and index by patient_id
    clinical_std <- clinical_std[!duplicated(clinical_std$patient_id), ]
    rownames(clinical_std) <- clinical_std$patient_id
    
    # Match clinical rows to methylation sample columns
    matched <- clinical_std[colnames(meth_m), , drop = FALSE]
    
    # Age covariate
    if ("age" %in% colnames(matched) && any(!is.na(matched$age))) {
      model_data$age <- matched$age
      final_covs <- c(final_covs, "age")
    }
    
    # Gender covariate
    if ("gender" %in% colnames(matched) && any(!is.na(matched$gender))) {
      levs <- unique(na.omit(matched$gender))
      if (length(levs) >= 2) {
        model_data$gender <- factor(matched$gender)
        final_covs <- c(final_covs, "gender")
      }
    }
    
    # Stage covariate — read from original clinical file if available
    clinical_raw <- read.delim(args[["clinical"]], stringsAsFactors = FALSE, check.names = FALSE)
    stage_col <- NULL
    # Check for stage column: custom mapping first, then common names
    cmap <- parse_clinical_mapping(args[["clinical_map"]])
    if (!is.null(cmap$stage)) {
      stage_col <- cmap$stage
    } else {
      stage_candidates <- c("ajcc_pathologic_stage", "stage", "tumor_stage", "clinical_stage", "Stage")
      stage_col <- intersect(colnames(clinical_raw), stage_candidates)[1]
    }
    if (!is.null(stage_col) && !is.na(stage_col) && stage_col %in% colnames(clinical_raw)) {
      # Build stage vector aligned to clinical_std patient order
      stage_raw <- clinical_raw[[stage_col]]
      names(stage_raw) <- clinical_std$patient_id[seq_along(stage_raw)]
      stage_matched <- stage_raw[match(colnames(meth_m), names(stage_raw))]
      stage_clean <- meth_map_stage(stage_matched)
      levs <- unique(na.omit(stage_clean))
      if (length(levs) >= 2) {
        model_data$stage_clean <- factor(stage_clean)
        final_covs <- c(final_covs, "stage_clean")
      }
    }
    
    if (length(final_covs) > length(intersect(final_covs, "TSS"))) {
      keep_idx <- complete.cases(model_data[, final_covs, drop = FALSE])
      if (sum(keep_idx) > 10) {
        meth_m <- meth_m[, keep_idx, drop = FALSE]
        met_beta <- met_beta[, keep_idx, drop = FALSE]
        batch_vec <- batch_vec[keep_idx]
        model_data <- model_data[keep_idx, , drop = FALSE]
      }
    }
    
    meth_msg(sprintf("Clinical covariates loaded: %s", paste(final_covs, collapse = ", ")))
  }, error = function(e) meth_msg(paste("Clinical data error:", e$message), level="WARN"))
}

# Remove single-level factors
for (cov in final_covs) {
  if (cov %in% names(model_data) && is.factor(model_data[[cov]]) && nlevels(model_data[[cov]]) < 2) {
    final_covs <- setdiff(final_covs, cov)
    model_data[[cov]] <- NULL
  }
}
qc <- meth_update_qc(qc, "batches", "covariates_protected", paste(final_covs, collapse=","))

# Build model matrix
if (length(final_covs) > 0) {
  mod <- model.matrix(as.formula(paste("~", paste(final_covs, collapse = " + "))), data = model_data)
} else {
  mod <- model.matrix(~ 1, data = data.frame(row.names = colnames(meth_m)))
}
rnk <- qr(mod)$rank
if (rnk < ncol(mod)) mod <- mod[, qr(mod)$pivot[seq_len(rnk)], drop = FALSE]

qc <- meth_update_qc(qc, "samples", "after_clinical_filter", ncol(meth_m))

# Apply ComBat
if (nlevels(batch_vec) > 1) {
  meth_msg(sprintf("Running ComBat on %d batches...", nlevels(batch_vec)))
  meth_m <- ComBat(dat = as.matrix(meth_m), batch = batch_vec, mod = mod)
  met_beta <- pmax(pmin(2^meth_m / (2^meth_m + 1), 1 - epsilon), epsilon)
  meth_msg("ComBat complete.")
  qc <- meth_update_qc(qc, "batches", "combat_applied", TRUE)
} else {
  meth_msg("Only one batch — ComBat skipped.")
  qc <- meth_update_qc(qc, "batches", "combat_applied", FALSE)
}

meth_msg(sprintf("Final sample count: %d", ncol(meth_m)))
qc <- meth_update_qc(qc, "samples", "final", ncol(meth_m))

rm(tss, valid, tab, singletons, keep, mod, rnk, model_data, final_covs)
if (exists("clinical_raw")) rm(clinical_raw)
if (exists("clinical_clean")) rm(clinical_clean)
if (exists("matched")) rm(matched)
gc()

# ==============================================================================
# STEP 6 — TOP VARIABLE PROBE SELECTION (MOFA/DIABLO)
# ==============================================================================
meth_step(6, "Feature Selection (MOFA/DIABLO)")

var_vec <- rowVars(meth_m)
n_top   <- min(args[["n-top"]], nrow(meth_m))
top     <- rownames(meth_m)[order(var_vec, decreasing = TRUE)][1:n_top]

meth_m_top   <- meth_m[top, , drop = FALSE]
met_beta_top <- met_beta[top, , drop = FALSE]

meth_msg(sprintf("Selected top %d probes", n_top))
meth_msg(sprintf("Variance range: [%.4f, %.4f]", min(var_vec[top]), max(var_vec[top])))

probe_log <- meth_log_probes(probe_log, "After Feature Selection", meth_m_top)
qc <- meth_update_qc(qc, "probes", "selected_top_variable", nrow(meth_m_top))

rm(var_vec, n_top, top)
gc()

# ==============================================================================
# STEP 7 — PCA & OUTLIER DETECTION
# ==============================================================================
meth_step(7, "PCA & Outlier Detection")

pca_final <- prcomp(t(meth_m_top), scale. = TRUE)
pc1 <- pca_final$x[, 1]; pc2 <- pca_final$x[, 2]
out_pc1 <- names(pc1)[abs(pc1 - mean(pc1)) > 3 * sd(pc1)]
out_pc2 <- names(pc2)[abs(pc2 - mean(pc2)) > 3 * sd(pc2)]
all_outliers <- unique(c(out_pc1, out_pc2))

qc <- meth_update_qc(qc, "outliers", "pc1_outliers", length(out_pc1))
qc <- meth_update_qc(qc, "outliers", "pc2_outliers", length(out_pc2))
qc <- meth_update_qc(qc, "outliers", "total_flagged", length(all_outliers))
qc <- meth_update_qc(qc, "outliers", "outlier_ids", as.list(all_outliers))

if (length(all_outliers) > 0) {
  meth_msg(sprintf("Flagged %d outlier(s) — see flagged_outliers.txt", length(all_outliers)), level="WARN")
  writeLines(all_outliers, file.path(args$outdir, "flagged_outliers.txt"))
} else {
  meth_msg("No outliers detected.")
}

pca_var <- summary(pca_final)$importance
meth_msg(sprintf("PC1: %.1f%% | PC2: %.1f%% | PC3: %.1f%%", pca_var[2,1]*100, pca_var[2,2]*100, pca_var[2,3]*100))

# ==============================================================================
# STEP 8 — QC VISUALIZATIONS
# ==============================================================================
meth_step(8, "QC Visualizations")

generate_meth_qc_plots(
  pca_final    = pca_final,
  pca_var      = pca_var,
  batch_vec    = batch_vec,
  met_beta_top = met_beta_top,
  outdir       = args$outdir
)

generate_meth_validation_figures(
  meth_m    = meth_m,
  meth_beta = met_beta,
  outdir    = args$outdir
)

# ==============================================================================
# STEP 9 — AUDIT TRAIL & FINAL EXPORT
# ==============================================================================
meth_step(9, "Audit Trail & Final Export")

if (!is.null(metadata)) {
  sample_info <- metadata[match(orig_barcodes[colnames(meth_m)], metadata$sample_id), ]
} else {
  sample_info <- data.frame(
    sample_id        = orig_barcodes[colnames(meth_m)],
    patient_id       = colnames(meth_m),
    sample_class     = "primary_tumor",
    batch            = as.character(batch_vec),
    stringsAsFactors = FALSE
  )
}

qc <- meth_update_qc(qc, "output_matrices", "mofa_probes", nrow(meth_m_top))
qc <- meth_update_qc(qc, "output_matrices", "mofa_samples", ncol(meth_m_top))
qc <- meth_update_qc(qc, "output_matrices", "ml_full_probes", nrow(meth_m))
qc <- meth_update_qc(qc, "output_matrices", "ml_full_samples", ncol(meth_m))

export_meth_results(
  meth_m_top   = meth_m_top,
  met_beta_top = met_beta_top,
  meth_m_full  = meth_m,
  met_beta_full= met_beta,
  sample_info  = sample_info,
  qc_metrics   = qc,
  probe_log    = probe_log,
  outdir       = args$outdir,
  metadata_supplied = !is.null(metadata)
)

# ==============================================================================
# FINAL REPORT
# ==============================================================================
meth_banner("METHYLATION PREPROCESSING COMPLETE")
meth_msg(sprintf("Completion time: %s", Sys.time()))
meth_msg(sprintf("MOFA/DIABLO    : %d probes x %d samples", nrow(meth_m_top), ncol(meth_m_top)))
meth_msg(sprintf("ML (full)      : %d probes x %d samples", nrow(meth_m), ncol(meth_m)))
meth_msg(sprintf("Batches        : %d", nlevels(batch_vec)))
meth_msg(sprintf("Outliers       : %d", length(all_outliers)))
meth_msg("Ready for multi-omics integration!")
