# ==============================================================================
# R/omicsflow.R — RStudio Orchestration Wrapper
# ==============================================================================

#' Run the full OmicsFlow pipeline from R
#'
#' Orchestrates the OmicsFlow analysis pipeline, dynamically executing preprocessing,
#' MOFA+ integration, machine learning survival analysis, and pathway enrichment
#' based on the provided subset of omics modalities. At least two modalities are
#' required for successful execution.
#'
#' @param rna Path to RNA-seq matrix file (.rds or .csv), or NULL
#' @param meth Path to DNA methylation matrix file (.rds or .csv), or NULL
#' @param cnv Path to CNV data file (.rds, .csv, or .tsv), or NULL
#' @param snv Path to SNV data file (.rds, .csv, or .tsv), or NULL
#' @param metadata Path to sample metadata CSV file
#' @param clinical Path to clinical TSV/CSV file, or NULL
#' @param clinical_map Path to clinical mapping JSON file, or NULL
#' @param outdir Base output directory for results
#' @param n_top_rna Number of top variable genes/features to select for RNA
#' @param n_top_meth Number of top variable features to select for Methylation
#' @param n_top_cnv Number of top variable features to select for CNV
#' @param mofa_factors Number of factors to train in MOFA
#' @param mofa_iter Number of iterations for MOFA training
#' @param mofa_prefilter Number of top features based on MOFA weights to retain for survival prefiltering
#' @param final_features Number of final features to retain in clinical model
#' @param meth_cross_react Path to cross-reactive probes list. Defaults to bundled Ensembl reference.
#' @param cache Path to CNV gene coordinates cache file. If NULL and configurations exist, defaults to configurations.
#' @param render_report Logical, whether to render the Quarto HTML report at the end
#' @param verbose Logical, whether to print detailed logs to console
#'
#' @return An S3 object of class \code{omicsflow_result} containing:
#'   \describe{
#'     \item{status}{Character, either "success" or "partial/failed"}
#'     \item{report_path}{Absolute path to the final Quarto report (if generated), else NULL}
#'     \item{output_dirs}{Named list of absolute paths to module output directories}
#'     \item{executed_modalities}{Character vector of executed components}
#'     \item{skipped_modalities}{Character vector of skipped components}
#'     \item{warnings}{Character vector of execution warnings}
#'     \item{errors}{Character vector of execution errors}
#'     \item{runtime}{Numeric execution time in seconds}
#'   }
#'
#' @examples
#' \dontrun{
#' res <- omicsflow(
#'   rna = "data/rna.rds",
#'   meth = "data/meth.rds",
#'   metadata = "data/sample_metadata.csv",
#'   outdir = "results/subset_run"
#' )
#' }
#'
#' @export
omicsflow <- function(rna = NULL, meth = NULL, cnv = NULL, snv = NULL,
                      metadata = NULL, clinical = NULL, clinical_map = NULL,
                      outdir = "results/myrun",
                      n_top_rna = 200, n_top_meth = 500, n_top_cnv = 500,
                      mofa_factors = 5, mofa_iter = 100,
                      mofa_prefilter = 50, final_features = 20,
                      meth_cross_react = NULL,
                      cache = NULL,
                      render_report = TRUE,
                      verbose = TRUE) {
  
  start_time <- Sys.time()
  
  if (verbose) msg_info("Starting OmicsFlow Pipeline Orchestration...")
  
  warnings_list <- character(0)
  errors_list <- character(0)
  
  # Override msg_warn to track warnings internally
  orig_msg_warn <- msg_warn
  msg_warn <- function(text) {
    warnings_list <<- c(warnings_list, text)
    orig_msg_warn(text)
  }
  
  # Override msg_fail to track errors internally
  orig_msg_fail <- msg_fail
  msg_fail <- function(text) {
    errors_list <<- c(errors_list, text)
    orig_msg_fail(text)
  }
  
  # Normalize input paths to absolute paths so they survive a setwd()
  safe_norm <- function(p, mustWork = TRUE) if (!is.null(p)) normalizePath(p, mustWork = mustWork) else NULL
  
  rna <- safe_norm(rna)
  meth <- safe_norm(meth)
  cnv <- safe_norm(cnv)
  snv <- safe_norm(snv)
  metadata <- safe_norm(metadata)
  clinical <- safe_norm(clinical)
  clinical_map <- safe_norm(clinical_map)
  meth_cross_react <- safe_norm(meth_cross_react, mustWork = FALSE)
  outdir <- normalizePath(outdir, mustWork = FALSE)
  
  # Default cache logic for CNV:
  if (is.null(cache) && !is.null(cnv)) {
    try({
      candidate_cache <- omicsflow_path("configs", "gene_coordinates.rds")
      candidate_cache2 <- omicsflow_path("realistic_cache.rds")
      if (file.exists(candidate_cache)) {
        cache <- candidate_cache
      } else if (file.exists(candidate_cache2)) {
        cache <- candidate_cache2
      }
    }, silent = TRUE)
    if (!is.null(cache) && verbose) {
      msg_info(sprintf("Using auto-resolved CNV cache: %s", cache))
    }
  }

  # Ensure meth_cross_react is set correctly
  if (is.null(meth_cross_react) || !nzchar(meth_cross_react)) {
    meth_cross_react <- omicsflow_path("configs", "cross_reactive_probes.csv")
  }

  # Step 1: Pre-flight Validation
  val_res <- validate_inputs(
    rna = rna, meth = meth, cnv = cnv, snv = snv,
    metadata = metadata, clinical = clinical, clinical_map = clinical_map,
    cnv_cache = cache
  )
  
  if (!val_res$valid) {
    stop("Pre-flight input validation failed. Please check the logs above.")
  }

  # Lazy Dependency Checks for Suggests
  if (!requireNamespace("MOFA2", quietly = TRUE)) {
    stop("MOFA2 is required. Run install_omicsflow_dependencies().")
  }
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler is required. Run install_omicsflow_dependencies().")
  }
  

  # Setup directory structure matching pipeline outputs
  dirs <- list(
    rna     = file.path(outdir, "output_rna"),
    meth    = file.path(outdir, "output_meth"),
    cnv     = file.path(outdir, "output_cnv"),
    snv     = file.path(outdir, "output_snv"),
    integ   = file.path(outdir, "output_integration"),
    ml      = file.path(outdir, "output_ml"),
    enrich  = file.path(outdir, "output_enrichment")
  )
  
  for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  
  status <- list(
    RNA = FALSE,
    METH = FALSE,
    CNV = FALSE,
    SNV = FALSE,
    INTEGRATION = FALSE,
    ML = FALSE,
    ENRICHMENT = FALSE
  )
  
  # Helper to run commands in the project root
  run_cmd <- function(name, cmd) {
    if (verbose) {
      cat(sprintf("\nRunning %s...\nCommand: %s\n", name, cmd))
    }
    
    # Execute command from the project root so scientific modules can resolve relative source() paths
    old_wd <- getwd()
    setwd(omicsflow_project_root())
    on.exit(setwd(old_wd), add = TRUE)
    
    exit_code <- system(cmd)
    
    setwd(old_wd)
    
    if (exit_code == 0) {
      if (verbose) msg_ok(sprintf("%s completed successfully.", name))
      return(TRUE)
    } else {
      if (verbose) msg_fail(sprintf("%s failed with exit code %d.", name, exit_code))
      return(FALSE)
    }
  }
  
  # --- 1. RNA Preprocessing ---
  if (!is.null(rna)) {
    cmd_rna <- paste(
      "Rscript", shQuote(omicsflow_path("modules", "rna", "preprocess_rna.R")),
      sprintf("--input %s", shQuote(rna)),
      sprintf("--metadata %s", shQuote(metadata)),
      sprintf("--outdir %s/", shQuote(dirs$rna)),
      sprintf("--n-top %d", n_top_rna),
      "--cor-low -1.0",
      "--cor-high 1.0"
    )
    status$RNA <- run_cmd("RNA Preprocessing", cmd_rna)
  } else {
    status$RNA <- TRUE # Treated as skipped/not-failed if not provided
  }
  
  # --- 2. Methylation Preprocessing ---
  if (!is.null(meth)) {
    cmd_meth <- paste(
      "Rscript", shQuote(omicsflow_path("modules", "methylation", "preprocess_meth.R")),
      sprintf("--input %s", shQuote(meth)),
      sprintf("--metadata %s", shQuote(metadata)),
      sprintf("--outdir %s/", shQuote(dirs$meth)),
      sprintf("--n-top %d", n_top_meth),
      sprintf("--cross-react %s", shQuote(meth_cross_react)),
      "--knn-k 2"
    )
    if (!is.null(clinical)) {
      cmd_meth <- paste(cmd_meth, sprintf("--clinical %s", shQuote(clinical)))
    }
    if (!is.null(clinical_map)) {
      cmd_meth <- paste(cmd_meth, sprintf("--clinical_map %s", shQuote(clinical_map)))
    }
    status$METH <- run_cmd("Methylation Preprocessing", cmd_meth)
  } else {
    status$METH <- TRUE
  }
  
  # --- 3. CNV Preprocessing ---
  if (!is.null(cnv)) {
    # Early validation already ensured cache exists
    cmd_cnv <- paste(
      "Rscript", shQuote(omicsflow_path("modules", "cnv", "preprocess_cnv.R")),
      sprintf("--input %s", shQuote(cnv)),
      sprintf("--metadata %s", shQuote(metadata)),
      sprintf("--outdir %s/", shQuote(dirs$cnv)),
      sprintf("--ntop %d", n_top_cnv),
      sprintf("--cache %s", shQuote(cache))
    )
    status$CNV <- run_cmd("CNV Preprocessing", cmd_cnv)
  } else {
    status$CNV <- TRUE
  }
  
  # --- 4. SNV Preprocessing ---
  if (!is.null(snv)) {
    cmd_snv <- paste(
      "Rscript", shQuote(omicsflow_path("modules", "snv", "preprocess_snv.R")),
      sprintf("--input %s", shQuote(snv)),
      sprintf("--metadata %s", shQuote(metadata)),
      sprintf("--outdir %s/", shQuote(dirs$snv))
    )
    status$SNV <- run_cmd("SNV Preprocessing", cmd_snv)
  } else {
    status$SNV <- TRUE
  }
  
  # --- 5. Integration ---
  # Only run integration if all preprocessed inputs exist and were successful
  provided_count <- sum(c(!is.null(rna), !is.null(meth), !is.null(cnv), !is.null(snv)))
  actual_success_count <- sum(c(
    !is.null(rna) && isTRUE(status$RNA),
    !is.null(meth) && isTRUE(status$METH),
    !is.null(cnv) && isTRUE(status$CNV),
    !is.null(snv) && isTRUE(status$SNV)
  ))
  
  preproc_success <- (actual_success_count == provided_count) && (actual_success_count >= 2)
  if (preproc_success) {
    cmd_integ <- paste(
      "Rscript", shQuote(omicsflow_path("modules", "integration", "run_integration.R")),
      sprintf("--metadata %s", shQuote(metadata)),
      sprintf("--outdir %s/", shQuote(dirs$integ)),
      sprintf("--factors %d", mofa_factors),
      sprintf("--iter %d", mofa_iter)
    )
    # Conditionally add inputs if they were provided
    if (!is.null(rna)) cmd_integ <- paste(cmd_integ, sprintf("--rna %s/rna_processed_matrix.rds", shQuote(dirs$rna)))
    if (!is.null(meth)) cmd_integ <- paste(cmd_integ, sprintf("--meth %s/methylation_processed_matrix.rds", shQuote(dirs$meth)))
    if (!is.null(cnv)) cmd_integ <- paste(cmd_integ, sprintf("--cnv %s/cnv_processed_matrix.rds", shQuote(dirs$cnv)))
    if (!is.null(snv)) cmd_integ <- paste(cmd_integ, sprintf("--snv %s/snv_processed_matrix.rds", shQuote(dirs$snv)))
    
    status$INTEGRATION <- run_cmd("MOFA+ Integration", cmd_integ)
  } else {
    if (verbose) {
      if (actual_success_count < 2) {
        msg_warn(sprintf("Skipping MOFA+ Integration: Requires at least 2 successful preprocessing views (got %d).", actual_success_count))
      } else {
        msg_warn("Skipping MOFA+ Integration due to preprocessing failures in one or more views.")
      }
    }
  }
  
  # --- 6. ML Survival ---
  if (isTRUE(status$INTEGRATION)) {
    if (!is.null(rna)) {
      cmd_ml <- paste(
        "Rscript", shQuote(omicsflow_path("modules", "ml", "run_ml.R")),
        sprintf("--mofa %s/mofa_model.rds", shQuote(dirs$integ)),
        sprintf("--metadata %s", shQuote(metadata)),
        sprintf("--outdir %s/", shQuote(dirs$ml)),
        sprintf("--mofa_prefilter %d", mofa_prefilter),
        sprintf("--final_features %d", final_features)
      )
      if (!is.null(rna)) {
        cmd_ml <- paste(cmd_ml, sprintf("--rna %s/rna_ml.rds", shQuote(dirs$rna)))
      }
      if (!is.null(clinical)) {
        cmd_ml <- paste(cmd_ml, sprintf("--clinical %s", shQuote(clinical)))
      }
      if (!is.null(clinical_map)) {
        cmd_ml <- paste(cmd_ml, sprintf("--clinical_map %s", shQuote(clinical_map)))
      }
      status$ML <- run_cmd("ML Survival Analysis", cmd_ml)
    } else {
      if (verbose) msg_warn("Skipping ML Survival Analysis: RNA modality is required to run ML survival models.")
      status$ML <- FALSE
    }
  } else {
    if (verbose) msg_warn("Skipping ML Survival Analysis since Integration was skipped or failed.")
  }
  
  # --- 7. Enrichment ---
  if (isTRUE(status$ML)) {
    cmd_enrich <- paste(
      "Rscript", shQuote(omicsflow_path("modules", "enrichment", "run_enrichment.R")),
      sprintf("--mofa %s/mofa_top_genes.rds", shQuote(dirs$ml)),
      sprintf("--rf %s/rf_top_genes.rds", shQuote(dirs$ml)),
      sprintf("--lasso %s/lasso_selected_genes.rds", shQuote(dirs$ml)),
      sprintf("--rna %s/rna_for_pathway.rds", shQuote(dirs$ml)),
      sprintf("--outdir %s/", shQuote(dirs$enrich))
    )
    status$ENRICHMENT <- run_cmd("Pathway Enrichment", cmd_enrich)
  } else {
    if (verbose) msg_warn("Skipping Pathway Enrichment since ML Survival was skipped or failed.")
  }
  
  # --- 8. Render Report ---
  if (render_report && isTRUE(status$INTEGRATION)) {
    if (verbose) msg_info("Rendering Quarto HTML Report...")
    # Map the relative results directory from the reports/ directory
    # If outdir is relative, we should adjust or pass absolute path
    abs_outdir <- normalizePath(outdir, mustWork = FALSE)
    report_qmd <- omicsflow_path("reports", "OmicsFlow_Report.qmd")
    # We render reports/OmicsFlow_Report.qmd using the generated outdir as results_dir
    cmd_report <- sprintf(
      "quarto render %s --to html -P results_dir:\"%s\" -P version:\"RStudio Runner\" -P subtitle:\"OmicsFlow Analytical Run\"",
      shQuote(report_qmd),
      gsub("\\\\", "/", abs_outdir)
    )
    exit_report <- system(cmd_report)
    if (exit_report == 0) {
      if (verbose) msg_ok("HTML Report successfully generated at reports/OmicsFlow_Report.html")
    } else {
      if (verbose) msg_fail("HTML Report generation failed.")
    }
  }
  
  
  # Prepare output object
  res <- list(
    status = if (isTRUE(status$INTEGRATION)) "success" else "partial/failed",
    report_path = if (render_report && isTRUE(status$INTEGRATION) && exit_report == 0) normalizePath(file.path(abs_outdir, "reports", "OmicsFlow_Report.html"), mustWork = FALSE) else NULL,
    output_dirs = dirs,
    executed_modalities = names(status)[status == TRUE],
    skipped_modalities = names(status)[status == FALSE],
    warnings = warnings_list,
    errors = errors_list,
    runtime = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  )
  class(res) <- "omicsflow_result"
  
  return(invisible(res))
}
