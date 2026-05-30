#!/usr/bin/env Rscript
# ==============================================================================
# run_integration.R — Multi-Omics Integration Pipeline (MOFA+)
# OmicsFlow | Phase 5: Integration Module
# ==============================================================================
#
# PURPOSE:
#   Executes the scientific steps for integrating multiple omics layers
#   using Multi-Omics Factor Analysis (MOFA+).
#   
# METHODOLOGY:
#   1. Load standardized omics matrices (RNA, Meth, CNV, SNV).
#   2. Find common patients across all views.
#   3. Build MOFA object and configure options (Gaussian/Bernoulli).
#   4. Train MOFA model (up to maxiter).
#   5. Identify active factors (variance explained > 1%).
#   6. Generate validation and discovery plots.
#
# ==============================================================================
# USAGE:
#   Rscript run_integration.R \
#     --rna     results/rna/rna_processed_matrix.rds \
#     --meth    results/methylation/methylation_processed_matrix.rds \
#     --cnv     results/cnv/cnv_processed_matrix.rds \
#     --snv     results/snv/snv_processed_matrix.rds \
#     --outdir  results/integration/
# ==============================================================================

suppressPackageStartupMessages(library(optparse))

# ------------------------------------------------------------------------------
# COMMAND-LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--rna"), type = "character", default = NULL,
              help = "Path to RNA matrix (.rds)"),
  make_option(c("--meth"), type = "character", default = NULL,
              help = "Path to Methylation matrix (.rds)"),
  make_option(c("--cnv"), type = "character", default = NULL,
              help = "Path to CNV matrix (.rds)"),
  make_option(c("--snv"), type = "character", default = NULL,
              help = "Path to SNV matrix (.rds)"),
  make_option(c("-o", "--outdir"), type = "character", default = "results/integration/",
              help = "Output directory [default= %default]"),
  make_option(c("-f", "--factors"), type = "integer", default = 15,
              help = "Number of latent factors to learn [default= %default]"),
  make_option(c("-i", "--iter"), type = "integer", default = 1000,
              help = "Max iterations [default= %default]"),
  make_option(c("--metadata"), type = "character", default = NULL,
              help = "Path to sample metadata CSV file (optional)"),
  make_option(c("-s", "--seed"), type = "integer", default = 42,
              help = "Random seed for reproducibility [default= %default]")
)

opt_parser <- OptionParser(
  usage = "Usage: %prog [options]",
  option_list = option_list,
  description = "OmicsFlow Phase 5: MOFA+ Multi-Omics Integration"
)
opt <- parse_args(opt_parser)

# ------------------------------------------------------------------------------
# SOURCE MODULE FILES
# ------------------------------------------------------------------------------
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) "modules/integration"
)
source(file.path(script_dir, "utils_integration.R"))
source(file.path(script_dir, "qc_integration.R"))
source(file.path(script_dir, "export_integration.R"))

# ------------------------------------------------------------------------------
# INITIALISE
# ------------------------------------------------------------------------------
tryCatch({
  load_int_packages()
  set.seed(opt$seed)
  
  if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
  }
  
  qc_metrics <- init_int_qc()
  
  int_banner("OmicsFlow v1.0.0 | MOFA+ Integration")
  int_msg(sprintf("Start time : %s", Sys.time()))
  int_msg(sprintf("Output dir : %s", opt$outdir))
  
  # Load metadata if supplied
  metadata <- NULL
  if (!is.null(opt$metadata) && opt$metadata != "") {
    meta_utils_path <- file.path(dirname(script_dir), "utils_metadata.R")
    if (file.exists(meta_utils_path)) {
      source(meta_utils_path)
    } else if (file.exists("modules/utils_metadata.R")) {
      source("modules/utils_metadata.R")
    }
    
    if (exists("load_metadata")) {
      int_msg(sprintf("Loading metadata: %s", opt$metadata))
      metadata <- load_metadata(opt$metadata)
    } else {
      int_msg("utils_metadata.R not found. Running without metadata layer.", level = "WARN")
    }
  }
  
  # ==============================================================================
  # 1. LOAD MATRICES
  # ==============================================================================
  int_step(1, "Loading Preprocessed Matrices")
  
  mofa_data <- list()
  
  if (!is.null(opt$rna) && opt$rna != "" && file.exists(opt$rna)) {
    rna  <- load_omics_matrix(opt$rna, "RNA", metadata)
    int_msg(sprintf("RNA: %d genes x %d samples", nrow(rna), ncol(rna)))
    qc_metrics <- add_int_qc(qc_metrics, "features", "RNA", nrow(rna))
    mofa_data$RNA <- rna
  }
  if (!is.null(opt$meth) && opt$meth != "" && file.exists(opt$meth)) {
    meth <- load_omics_matrix(opt$meth, "Methylation", metadata)
    int_msg(sprintf("Methylation: %d probes x %d samples", nrow(meth), ncol(meth)))
    qc_metrics <- add_int_qc(qc_metrics, "features", "Methylation", nrow(meth))
    mofa_data$Methylation <- meth
  }
  if (!is.null(opt$cnv) && opt$cnv != "" && file.exists(opt$cnv)) {
    cnv  <- load_omics_matrix(opt$cnv, "CNV", metadata)
    int_msg(sprintf("CNV: %d genes x %d samples", nrow(cnv), ncol(cnv)))
    qc_metrics <- add_int_qc(qc_metrics, "features", "CNV", nrow(cnv))
    mofa_data$CNV <- cnv
  }
  if (!is.null(opt$snv) && opt$snv != "" && file.exists(opt$snv)) {
    snv  <- load_omics_matrix(opt$snv, "SNV", metadata)
    int_msg(sprintf("SNV: %d genes x %d samples", nrow(snv), ncol(snv)))
    qc_metrics <- add_int_qc(qc_metrics, "features", "SNV", nrow(snv))
    mofa_data$SNV <- snv
  }
  
  if (length(mofa_data) < 2) {
    stop("Multi-omics integration requires at least 2 views.")
  }
  
  # ==============================================================================
  # 2. FIND COMMON PATIENTS
  # ==============================================================================
  int_step(2, "Identifying Common Patients")
  
  common_patients <- Reduce(intersect, lapply(mofa_data, colnames))
  
  int_msg(sprintf("Common patients across provided views: %d", length(common_patients)))
  qc_metrics <- add_int_qc(qc_metrics, "samples", "common_patients", length(common_patients))
  
  if (length(common_patients) < 50) {
    int_msg("Less than 50 common patients — MOFA may underfit!", level = "WARN")
  }
  
  mofa_data <- lapply(mofa_data, function(x) x[, common_patients, drop = FALSE])
  
  # Expose individual variables if needed by other components
  if ("RNA" %in% names(mofa_data)) rna <- mofa_data$RNA
  if ("Methylation" %in% names(mofa_data)) meth <- mofa_data$Methylation
  if ("CNV" %in% names(mofa_data)) cnv <- mofa_data$CNV
  if ("SNV" %in% names(mofa_data)) snv <- mofa_data$SNV
  
  int_msg("All matrices subsetted to common patients", level = "DETAILS")
  
  # ==============================================================================
  # 3. BUILD MOFA OBJECT
  # ==============================================================================
  int_step(3, "Constructing MOFA Object")
  
  mofa <- create_mofa(mofa_data)
  
  int_msg("MOFA object created")
  int_msg(sprintf("Views: %d (%s)", length(mofa_data), paste(names(mofa_data), collapse = ", ")), level = "DETAILS")
  int_msg(sprintf("Samples: %d", length(common_patients)), level = "DETAILS")
  int_msg(sprintf("Total features: %s", format(sum(sapply(mofa_data, nrow)), big.mark = ",")), level = "DETAILS")
  
  # Export views dimensions to QC metrics
  for (vname in names(mofa_data)) {
    qc_metrics <- add_int_qc(qc_metrics, "views", vname, nrow(mofa_data[[vname]]))
  }
  
  # ==============================================================================
  # 4. CONFIGURE OPTIONS
  # ==============================================================================
  int_step(4, "Configuring MOFA Options")
  
  data_opts <- get_default_data_options(mofa)
  data_opts$scale_views <- TRUE
  int_msg("scale_views = TRUE (prevents methylation dominance)")
  
  model_opts <- get_default_model_options(mofa)
  model_opts$num_factors <- opt$factors
  
  all_likelihoods <- c(
    RNA         = "gaussian",
    Methylation = "gaussian",
    CNV         = "gaussian",
    SNV         = "bernoulli"
  )
  model_opts$likelihoods <- all_likelihoods[names(mofa_data)]
  
  int_msg(sprintf("num_factors = %d", opt$factors))
  int_msg(sprintf("likelihoods: %s", paste(sprintf("%s = %s", names(model_opts$likelihoods), model_opts$likelihoods), collapse = " | ")))
  
  train_opts <- get_default_training_options(mofa)
  train_opts$maxiter <- opt$iter
  train_opts$seed    <- 42
  train_opts$convergence_mode <- "slow"   # stricter ELBO tolerance (1e-6 vs 1e-4)
  int_msg("seed = 42 (reproducible)")
  int_msg(sprintf("maxiter = %d", opt$iter))
  int_msg("convergence_mode = 'slow' (publication-quality tolerance)")
  
  qc_metrics <- add_int_qc(qc_metrics, "model", "num_factors", opt$factors)
  qc_metrics <- add_int_qc(qc_metrics, "model", "maxiter", opt$iter)
  
  # ==============================================================================
  # 5. PREPARE & TRAIN MOFA
  # ==============================================================================
  int_step(5, "Preparing & Training MOFA Model")
  
  mofa <- prepare_mofa(
    mofa,
    data_options     = data_opts,
    model_options    = model_opts,
    training_options = train_opts
  )
  int_msg("MOFA prepared \u2014 data validated, options locked in")
  
  int_msg("Training MOFA Model... (typically takes 10-30 minutes)")
  temp_hdf5 <- file.path(opt$outdir, "mofa_temp.hdf5")
  mofa <- run_mofa(mofa, outfile = temp_hdf5, use_basilisk = TRUE)
  int_msg("MOFA training complete!")
  
  # ==============================================================================
  # 6. IDENTIFY ACTIVE FACTORS
  # ==============================================================================
  int_step(6, "Identifying Active Factors")
  
  var_exp        <- get_variance_explained(mofa)$r2_per_factor[[1]]
  active_factors <- which(rowMeans(var_exp) > 0.01)
  
  int_msg(sprintf("Active factors (mean R\u00b2 > 1%%): %d / %d", length(active_factors), nrow(var_exp)))
  int_msg(sprintf("Factors retained: %s", paste(active_factors, collapse = ", ")))
  
  qc_metrics <- add_int_qc(qc_metrics, "model", "active_factors_count", length(active_factors))
  
  # Summary of Variance Explained per layer
  var_per_view <- get_variance_explained(mofa)$r2_total[[1]]
  for (v in names(var_per_view)) {
    int_msg(sprintf("  %s R\u00b2 = %.3f", v, var_per_view[[v]]), level = "DETAILS")
    qc_metrics <- add_int_qc(qc_metrics, "model", sprintf("r2_total_%s", v), var_per_view[[v]])
  }
  
  # ==============================================================================
  # 7. GENERATE QC PLOTS & EXPORT
  # ==============================================================================
  generate_int_qc_plots(mofa, opt$outdir)
  export_int_results(mofa, active_factors, opt$outdir)
  export_int_qc(qc_metrics, opt$outdir)
  
  int_banner("INTEGRATION MODULE COMPLETE \u2714")
  
}, error = function(e) {
  int_msg(sprintf("FATAL ERROR: %s", e$message), level = "ERROR")
  quit(status = 1)
})
