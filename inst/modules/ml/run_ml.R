#!/usr/bin/env Rscript
# ==============================================================================
# run_ml.R — Machine Learning Pipeline
# OmicsFlow | Phase 6: ML Module
# ==============================================================================
#
# PURPOSE:
#   Extracts top features from the MOFA+ Integration layer and builds
#   machine learning survival models (RF, XGBoost, LASSO).
#   
# METHODOLOGY:
#   1. Load MOFA model, RNA matrix, and Clinical Data.
#   2. Train/Test split (80/20).
#   3. Hybrid Feature Selection:
#      - MOFA weights -> biological pre-filter (Top N, unsupervised)
#      - Variance -> Final selection on training set only (zero leakage)
#   4. Train Models (Random Survival Forest, Cox XGBoost, LASSO Cox)
#   5. Evaluate on test set (Concordance Index).
#   6. Generate Kaplan-Meier curves and validation plots.
#
# ==============================================================================
# USAGE:
#   Rscript run_ml.R \
#     --mofa      results/integration/mofa_model.rds \
#     --rna       results/rna/rna_ml.rds \
#     --clinical  data/clinical_data.tsv \
#     --outdir    results/ml/
# ==============================================================================

suppressPackageStartupMessages(library(optparse))

# ------------------------------------------------------------------------------
# COMMAND-LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--mofa"), type = "character", default = "results/integration/mofa_model.rds",
              help = "Path to MOFA model (.rds)"),
  make_option(c("--rna"), type = "character", default = "results/rna/rna_ml.rds",
              help = "Path to RNA ML matrix (.rds)"),
  make_option(c("--clinical"), type = "character", default = "data/clinical_data.tsv",
              help = "Path to Clinical TSV data"),
  make_option(c("-o", "--outdir"), type = "character", default = "results/ml/",
              help = "Output directory [default= %default]"),
  make_option(c("--mofa_prefilter"), type = "integer", default = 300,
              help = "Top genes by MOFA weight [default= %default]"),
  make_option(c("--final_features"), type = "integer", default = 100,
              help = "Top genes by variance on training set [default= %default]"),
  make_option(c("--train_ratio"), type = "numeric", default = 0.80,
              help = "Train/Test split ratio [default= %default]"),
  make_option(c("--clinical_map"), type = "character", default = NULL,
              help = "Path to clinical column mapping JSON file (optional)"),
  make_option(c("--metadata"), type = "character", default = NULL,
              help = "Path to sample metadata CSV file (optional)"),
  make_option(c("-s", "--seed"), type = "integer", default = 42,
              help = "Random seed for reproducibility [default= %default]")
)

opt_parser <- OptionParser(
  usage = "Usage: %prog [options]",
  option_list = option_list,
  description = "OmicsFlow Phase 6: ML Prediction Pipeline"
)
opt <- parse_args(opt_parser)

# ------------------------------------------------------------------------------
# SOURCE MODULE FILES
# ------------------------------------------------------------------------------
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) "modules/ml"
)
source(file.path(script_dir, "utils_ml.R"))
source(file.path(script_dir, "qc_ml.R"))
source(file.path(script_dir, "export_ml.R"))

# ------------------------------------------------------------------------------
# INITIALISE
# ------------------------------------------------------------------------------
tryCatch({
  load_ml_packages()
  
  if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
  }
  
  set.seed(opt$seed)
  qc_metrics <- init_ml_qc()
  
  ml_banner("OmicsFlow v2.0.1 | ML Survival Analysis")
  ml_msg(sprintf("Start time : %s", Sys.time()))
  ml_msg(sprintf("Output dir : %s", opt$outdir))
  
  # ==============================================================================
  # 1. LOAD DATA
  # ==============================================================================
  ml_step(1, "Loading Data")
  
  loaded <- readRDS(opt$mofa)
  if (!"model" %in% names(loaded)) stop("MOFA object must contain '$model'")
  mofa <- loaded$model
  
  rna <- readRDS(opt$rna)
  
  # Source and load metadata if supplied
  metadata <- NULL
  if (!is.null(opt$metadata) && opt$metadata != "") {
    meta_utils_path <- file.path(dirname(script_dir), "utils_metadata.R")
    if (file.exists(meta_utils_path)) {
      source(meta_utils_path)
    } else if (file.exists("modules/utils_metadata.R")) {
      source("modules/utils_metadata.R")
    }
    if (exists("load_metadata")) {
      ml_msg(sprintf("Loading metadata: %s", opt$metadata))
      metadata <- load_metadata(opt$metadata)
    }
  }
  
  # Source clinical abstraction layer
  clinical_utils_path <- file.path(dirname(script_dir), "utils_clinical.R")
  if (file.exists(clinical_utils_path)) {
    source(clinical_utils_path)
  } else if (file.exists("modules/utils_clinical.R")) {
    source("modules/utils_clinical.R")
  }
  
  clinical_standardized <- load_clinical_data(opt$clinical, opt$clinical_map, metadata)
  
  rownames(rna) <- str_replace(rownames(rna), "_RNA$", "")
  
  ml_msg("All data loaded")
  ml_msg(sprintf("RNA matrix: %d genes x %d samples", nrow(rna), ncol(rna)), level="DETAILS")
  ml_msg(sprintf("Clinical rows: %d", nrow(clinical_standardized)), level="DETAILS")
  
  # ==============================================================================
  # 2. PREPARE SURVIVAL DATA
  # ==============================================================================
  ml_step(2, "Preparing Survival Data")
  
  clinical <- clinical_standardized
  if ("gender" %in% colnames(clinical)) {
    clinical$gender <- as.factor(clinical$gender)
  }
  
  clinical <- clinical %>% filter(!is.na(os_time) & os_time > 0)
  
  ml_msg(sprintf("Valid patients: %d", nrow(clinical)))
  ml_msg(sprintf("Events (deceased): %d (%.1f%%)", sum(clinical$os_event), 100 * mean(clinical$os_event)), level="DETAILS")
  ml_msg(sprintf("Median survival: %.0f days", median(clinical$os_time)), level="DETAILS")
  
  qc_metrics <- add_ml_qc(qc_metrics, "samples", "total_valid_patients", nrow(clinical))
  
  # ==============================================================================
  # 3. TRAIN / TEST SPLIT
  # ==============================================================================
  ml_step(3, sprintf("Train/Test Split (%.0f/%.0f)", opt$train_ratio*100, (1-opt$train_ratio)*100))
  
  train_idx <- createDataPartition(clinical$os_event, p = opt$train_ratio, list = FALSE)
  patients_train <- clinical$patient_id[train_idx]
  patients_test  <- clinical$patient_id[-train_idx]
  
  ml_msg(sprintf("Training patients: %d", length(patients_train)))
  ml_msg(sprintf("Test patients:     %d", length(patients_test)))
  
  qc_metrics <- add_ml_qc(qc_metrics, "samples", "train_patients", length(patients_train))
  qc_metrics <- add_ml_qc(qc_metrics, "samples", "test_patients", length(patients_test))
  
  # ==============================================================================
  # 4. HYBRID FEATURE SELECTION
  # ==============================================================================
  ml_step(4, "Hybrid Feature Selection")
  
  ml_msg(sprintf("Step A: MOFA biological pre-filter (top %d genes)...", opt$mofa_prefilter))
  rna_weights <- get_weights(mofa, views = "RNA", as.data.frame = TRUE)
  
  mofa_top_df <- rna_weights %>%
    mutate(feature = str_replace(feature, "_RNA$", "")) %>%
    group_by(feature) %>%
    summarise(max_abs_weight = max(abs(value)), .groups = "drop") %>%
    arrange(desc(max_abs_weight)) %>%
    head(opt$mofa_prefilter)
  
  mofa_candidate_genes <- intersect(mofa_top_df$feature, rownames(rna))
  ml_msg(sprintf("MOFA candidates: %d genes", length(mofa_candidate_genes)))
  
  ml_msg(sprintf("Step B: Variance-based final selection on training set (top %d)...", opt$final_features))
  common_train <- intersect(patients_train, colnames(rna))
  rna_train    <- rna[mofa_candidate_genes, common_train, drop = FALSE]
  
  gene_var       <- apply(rna_train, 1, var)
  selected_genes <- names(sort(gene_var, decreasing = TRUE))[1:opt$final_features]
  
  ml_msg(sprintf("Final feature set: %d genes", length(selected_genes)))
  
  qc_metrics <- add_ml_qc(qc_metrics, "features", "final_features_count", length(selected_genes))
  
  # ==============================================================================
  # 5. BUILD TRAINING MATRIX
  # ==============================================================================
  ml_step(5, "Building Training Matrix")
  
  X_train_df <- as.data.frame(t(rna_train[selected_genes, ]))
  X_train_df$patient_id <- rownames(X_train_df)
  
  train_data <- X_train_df %>%
    inner_join(clinical %>% select(patient_id, os_time, os_event, age, gender), by = "patient_id")
  
  X_train_mat <- as.matrix(train_data[, selected_genes])
  surv_train  <- Surv(time = train_data$os_time, event = train_data$os_event)
  
  ml_msg(sprintf("Training matrix: %d patients x %d features", nrow(X_train_mat), ncol(X_train_mat)))
  
  # ==============================================================================
  # 6. TRAIN SURVIVAL MODELS
  # ==============================================================================
  ml_step(6, "Training Survival Models")
  
  # -- Random Forest --
  ml_msg("Training Random Survival Forest (randomForestSRC)...", level = "DETAILS")
  rf_train_df <- train_data %>% select(all_of(selected_genes), os_time, os_event) %>% as.data.frame()
  rf_model <- rfsrc(Surv(os_time, os_event) ~ ., data = rf_train_df, ntree = 500, nsplit = 10, importance = TRUE, seed = 42)
  rf_oob_c <- 1 - tail(rf_model$err.rate, 1)
  ml_msg(sprintf("RF trained \u2014 OOB C-index: %.3f", rf_oob_c))
  
  # -- XGBoost --
  ml_msg("Training Cox XGBoost...", level = "DETAILS")
  label_cox <- ifelse(train_data$os_event == 1, train_data$os_time, -train_data$os_time)
  dtrain <- xgb.DMatrix(data = X_train_mat, label = label_cox)
  xgb_params <- list(objective = "survival:cox", eval_metric = "cox-nloglik", max_depth = 4, eta = 0.03, subsample = 0.8, colsample_bytree = 0.8)
  
  xgb_cv <- xgb.cv(params = xgb_params, data = dtrain, nrounds = 500, nfold = 5, early_stopping_rounds = 20, verbose = 0)
  best_rounds <- xgb_cv$evaluation_log$iter[which.min(xgb_cv$evaluation_log$test_cox_nloglik_mean)]
  xgb_model <- xgb.train(params = xgb_params, data = dtrain, nrounds = best_rounds, evals = list(train = dtrain), verbose = 0)
  ml_msg(sprintf("XGBoost trained \u2014 Best rounds: %d", best_rounds))
  
  # -- LASSO --
  ml_msg("Training LASSO Cox Regression...", level = "DETAILS")
  cv_lasso <- cv.glmnet(x = X_train_mat, y = surv_train, family = "cox", alpha = 1, nfolds = 5)
  lasso_coefs <- coef(cv_lasso, s = "lambda.min")
  lasso_genes <- rownames(lasso_coefs)[as.numeric(lasso_coefs) != 0]
  ml_msg(sprintf("LASSO trained \u2014 %d genes selected", length(lasso_genes)))
  
  # ==============================================================================
  # 7. EVALUATE ON TEST SET
  # ==============================================================================
  ml_step(7, "Evaluating on Test Set")
  
  common_test <- intersect(patients_test, colnames(rna))
  rna_test    <- rna[selected_genes, common_test, drop = FALSE]
  
  X_test_df             <- as.data.frame(t(rna_test))
  X_test_df$patient_id  <- rownames(X_test_df)
  
  test_data <- X_test_df %>%
    inner_join(clinical %>% select(patient_id, os_time, os_event), by = "patient_id")
  
  X_test_mat <- as.matrix(test_data[, selected_genes])
  surv_test  <- Surv(time = test_data$os_time, event = test_data$os_event)
  
  rf_pred   <- randomForestSRC::predict.rfsrc(rf_model, newdata = as.data.frame(test_data[, selected_genes]))
  rf_risk   <- -rf_pred$predicted
  rf_c_test <- concordance(surv_test ~ rf_risk)$concordance
  
  xgb_risk   <- -stats::predict(xgb_model, newdata = xgb.DMatrix(data = X_test_mat))
  xgb_c_test <- concordance(surv_test ~ xgb_risk)$concordance
  
  lasso_risk   <- -as.numeric(stats::predict(cv_lasso, newx = X_test_mat, s = "lambda.min", type = "link"))
  lasso_c_test <- concordance(surv_test ~ lasso_risk)$concordance
  
  ml_msg("Model Performance Summary:")
  ml_msg(sprintf("RF      C-index: %.3f", rf_c_test), level = "DETAILS")
  ml_msg(sprintf("XGBoost C-index: %.3f", xgb_c_test), level = "DETAILS")
  ml_msg(sprintf("LASSO   C-index: %.3f", lasso_c_test), level = "DETAILS")
  
  qc_metrics <- add_ml_qc(qc_metrics, "performance", "rf_c_index", rf_c_test)
  qc_metrics <- add_ml_qc(qc_metrics, "performance", "xgb_c_index", xgb_c_test)
  qc_metrics <- add_ml_qc(qc_metrics, "performance", "lasso_c_index", lasso_c_test)
  
  # ==============================================================================
  # 8. QC & EXPORT
  # ==============================================================================
  
  rf_imp_df <- data.frame(gene = names(rf_model$importance), importance = rf_model$importance, row.names = NULL) %>% arrange(desc(importance))
  xgb_imp_df <- xgb.importance(model = xgb_model) %>% rename(gene = Feature)
  
  results_summary <- data.frame(
    model       = c("Random Forest", "XGBoost Cox", "LASSO Cox"),
    c_index     = round(c(rf_c_test, xgb_c_test, lasso_c_test), 3),
    top_5_genes = c(
      paste(head(rf_imp_df$gene,   5), collapse = ", "),
      paste(head(xgb_imp_df$gene,  5), collapse = ", "),
      paste(head(lasso_genes,      5), collapse = ", ")
    )
  )
  
  top_features_export <- mofa_top_df %>%
    filter(feature %in% selected_genes) %>%
    mutate(train_variance = gene_var[feature]) %>%
    arrange(desc(train_variance))
  
  generate_ml_qc_plots(results_summary, top_features_export, test_data, xgb_risk, "os_time", "os_event", opt$outdir)
  export_ml_results(results_summary, top_features_export, rf_model, xgb_model, cv_lasso, mofa_top_df, selected_genes, lasso_genes, rna, opt$outdir)
  export_ml_qc(qc_metrics, opt$outdir)
  
  ml_banner("ML MODULE COMPLETE \u2714")
  
}, error = function(e) {
  ml_msg(sprintf("FATAL ERROR: %s", e$message), level = "ERROR")
  quit(status = 1)
})
