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
  make_option(c("--rna"), type = "character", default = "results/rna/rna_processed_matrix.rds",
              help = "Path to RNA matrix (.rds)"),
  make_option(c("--meth"), type = "character", default = "results/methylation/methylation_processed_matrix.rds",
              help = "Path to Methylation matrix (.rds)"),
  make_option(c("--cnv"), type = "character", default = "results/cnv/cnv_processed_matrix.rds",
              help = "Path to CNV matrix (.rds)"),
  make_option(c("--snv"), type = "character", default = "results/snv/snv_processed_matrix.rds",
              help = "Path to SNV matrix (.rds)"),
  make_option(c("-o", "--outdir"), type = "character", default = "results/integration/",
              help = "Output directory [default= %default]"),
  make_option(c("-f", "--factors"), type = "integer", default = 15,
              help = "Number of latent factors to learn [default= %default]"),
  make_option(c("-i", "--iter"), type = "integer", default = 1000,
              help = "Max iterations [default= %default]")
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
  set.seed(42)
  
  if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
  }
  
  qc_metrics <- init_int_qc()
  
  int_banner("OmicsFlow v1.0.0 | MOFA+ Integration")
  int_msg(sprintf("Start time : %s", Sys.time()))
  int_msg(sprintf("Output dir : %s", opt$outdir))
  
  # ==============================================================================
  # 1. LOAD MATRICES
  # ==============================================================================
  int_step(1, "Loading Preprocessed Matrices")
  
  rna  <- load_omics_matrix(opt$rna, "RNA")
  int_msg(sprintf("RNA: %d genes x %d samples", nrow(rna), ncol(rna)))
  qc_metrics <- add_int_qc(qc_metrics, "features", "RNA", nrow(rna))
  
  meth <- load_omics_matrix(opt$meth, "Methylation")
  int_msg(sprintf("Methylation: %d probes x %d samples", nrow(meth), ncol(meth)))
  qc_metrics <- add_int_qc(qc_metrics, "features", "Methylation", nrow(meth))
  
  cnv  <- load_omics_matrix(opt$cnv, "CNV")
  int_msg(sprintf("CNV: %d genes x %d samples", nrow(cnv), ncol(cnv)))
  qc_metrics <- add_int_qc(qc_metrics, "features", "CNV", nrow(cnv))
  
  snv  <- load_omics_matrix(opt$snv, "SNV")
  int_msg(sprintf("SNV: %d genes x %d samples", nrow(snv), ncol(snv)))
  qc_metrics <- add_int_qc(qc_metrics, "features", "SNV", nrow(snv))
  
  # ==============================================================================
  # 2. FIND COMMON PATIENTS
  # ==============================================================================
  int_step(2, "Identifying Common Patients")
  
  common_patients <- Reduce(intersect, list(
    colnames(meth),
    colnames(cnv),
    colnames(rna),
    colnames(snv)
  ))
  
  int_msg(sprintf("Common patients across ALL 4 omics: %d", length(common_patients)))
  qc_metrics <- add_int_qc(qc_metrics, "samples", "common_patients", length(common_patients))
  
  if (length(common_patients) < 50) {
    int_msg("Less than 50 common patients \u2014 MOFA may underfit!", level = "WARN")
  }
  
  meth <- meth[, common_patients, drop = FALSE]
  cnv  <- cnv[,  common_patients, drop = FALSE]
  rna  <- rna[,  common_patients, drop = FALSE]
  snv  <- snv[,  common_patients, drop = FALSE]
  
  int_msg("All matrices subsetted to common patients", level = "DETAILS")
  
  # ==============================================================================
  # 3. BUILD MOFA OBJECT
  # ==============================================================================
  int_step(3, "Constructing MOFA Object")
  
  mofa_data <- list(
    RNA         = rna,
    Methylation = meth,
    CNV         = cnv,
    SNV         = snv
  )
  
  mofa <- create_mofa(mofa_data)
  
  int_msg("MOFA object created")
  int_msg(sprintf("Views: %d", length(mofa_data)), level = "DETAILS")
  int_msg(sprintf("Samples: %d", length(common_patients)), level = "DETAILS")
  int_msg(sprintf("Total features: %s", format(sum(sapply(mofa_data, nrow)), big.mark = ",")), level = "DETAILS")
  
  # ==============================================================================
  # 4. CONFIGURE OPTIONS
  # ==============================================================================
  int_step(4, "Configuring MOFA Options")
  
  data_opts <- get_default_data_options(mofa)
  data_opts$scale_views <- TRUE
  int_msg("scale_views = TRUE (prevents methylation dominance)")
  
  model_opts <- get_default_model_options(mofa)
  model_opts$num_factors <- opt$factors
  model_opts$likelihoods <- c(
    RNA         = "gaussian",
    Methylation = "gaussian",
    CNV         = "gaussian",
    SNV         = "bernoulli"
  )
  int_msg(sprintf("num_factors = %d", opt$factors))
  int_msg("likelihoods: RNA/Meth/CNV = gaussian | SNV = bernoulli")
  
  train_opts <- get_default_training_options(mofa)
  train_opts$maxiter <- opt$iter
  train_opts$seed    <- 42
  int_msg("seed = 42 (reproducible)")
  int_msg(sprintf("maxiter = %d", opt$iter))
  
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
