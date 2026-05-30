#!/usr/bin/env Rscript
# ==============================================================================
# run_realistic_validation.R — Full OmicsFlow Validation on Realistic Cohort
# ==============================================================================
#
# PURPOSE:
#   Run the entire OmicsFlow pipeline on the realistic oncology cohort:
#     1. RNA preprocessing
#     2. Methylation preprocessing
#     3. CNV preprocessing
#     4. SNV preprocessing
#     5. MOFA+ Integration
#     6. ML Survival Analysis
#     7. Pathway Enrichment
#
#   Then produce a comprehensive PASS/FAIL validation report with:
#     - Cohort summary
#     - Survival statistics
#     - Enriched pathways
#     - Latent factor summary
#     - Final verdict
#
# ==============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   OmicsFlow | Realistic Cohort End-to-End Validation        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

start_time <- Sys.time()

# ==============================================================================
# CONFIGURATION
# ==============================================================================
data_dir    <- "data/realistic_cohort"
results_dir <- "results/realistic_validation"
cache_file  <- "realistic_cache.rds"

# Verify data directory exists
if (!dir.exists(data_dir)) {
  stop("Data directory not found: ", data_dir,
       "\nRun generate_realistic_cohort.R first.")
}

# ==============================================================================
# MODULE RUNNER
# ==============================================================================
run_module <- function(name, cmd) {
  cat(sprintf("\n┌─────────────────────────────────────────────────────────────┐\n"))
  cat(sprintf("│  Running: %-49s│\n", name))
  cat(sprintf("└─────────────────────────────────────────────────────────────┘\n"))
  cat("CMD:", cmd, "\n\n")
  t0 <- Sys.time()
  res <- system(cmd, intern = FALSE)
  elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
  if (res != 0) {
    cat(sprintf("✗ %s FAILED (exit code %d, %.1fs)\n", name, res, elapsed))
    return(FALSE)
  }
  cat(sprintf("✓ %s PASSED (%.1fs)\n", name, elapsed))
  return(TRUE)
}

# ==============================================================================
# DEFINE PIPELINE COMMANDS
# ==============================================================================

# Output directories
dirs <- list(
  rna     = file.path(results_dir, "output_rna"),
  meth    = file.path(results_dir, "output_meth"),
  cnv     = file.path(results_dir, "output_cnv"),
  snv     = file.path(results_dir, "output_snv"),
  integ   = file.path(results_dir, "output_integration"),
  ml      = file.path(results_dir, "output_ml"),
  enrich  = file.path(results_dir, "output_enrichment")
)

for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Build command strings
cmds <- list()

cmds$RNA <- paste(
  "Rscript modules/rna/preprocess_rna.R",
  sprintf("--input %s/rna.rds", data_dir),
  sprintf("--metadata %s/sample_metadata.csv", data_dir),
  sprintf("--clinical %s/custom_clinical.tsv", data_dir),
  sprintf("--clinical_map %s/clinical_map.json", data_dir),
  sprintf("--outdir %s/", dirs$rna),
  "--n-top 200",
  "--cor-low -1.0",
  "--cor-high 1.0"
)

cmds$METH <- paste(
  "Rscript modules/methylation/preprocess_meth.R",
  sprintf("--input %s/meth.rds", data_dir),
  sprintf("--metadata %s/sample_metadata.csv", data_dir),
  sprintf("--clinical %s/custom_clinical.tsv", data_dir),
  sprintf("--clinical_map %s/clinical_map.json", data_dir),
  sprintf("--outdir %s/", dirs$meth),
  "--n-top 500",
  "--knn-k 2"
)

cmds$CNV <- paste(
  "Rscript modules/cnv/preprocess_cnv.R",
  sprintf("--input %s/cnv.rds", data_dir),
  sprintf("--metadata %s/sample_metadata.csv", data_dir),
  sprintf("--outdir %s/", dirs$cnv),
  "--ntop 500",
  sprintf("--cache %s", cache_file)
)

cmds$SNV <- paste(
  "Rscript modules/snv/preprocess_snv.R",
  sprintf("--input %s/snv.rds", data_dir),
  sprintf("--metadata %s/sample_metadata.csv", data_dir),
  sprintf("--outdir %s/", dirs$snv)
)

cmds$INTEGRATION <- paste(
  "Rscript modules/integration/run_integration.R",
  sprintf("--rna %s/rna_processed_matrix.rds", dirs$rna),
  sprintf("--meth %s/methylation_processed_matrix.rds", dirs$meth),
  sprintf("--cnv %s/cnv_processed_matrix.rds", dirs$cnv),
  sprintf("--snv %s/snv_processed_matrix.rds", dirs$snv),
  sprintf("--metadata %s/sample_metadata.csv", data_dir),
  sprintf("--outdir %s/", dirs$integ),
  "--factors 15",
  "--iter 1000"
)

cmds$ML <- paste(
  "Rscript modules/ml/run_ml.R",
  sprintf("--mofa %s/mofa_model.rds", dirs$integ),
  sprintf("--rna %s/rna_ml.rds", dirs$rna),
  sprintf("--clinical %s/custom_clinical.tsv", data_dir),
  sprintf("--metadata %s/sample_metadata.csv", data_dir),
  sprintf("--clinical_map %s/clinical_map.json", data_dir),
  sprintf("--outdir %s/", dirs$ml),
  "--mofa_prefilter 50",
  "--final_features 20"
)

cmds$ENRICHMENT <- paste(
  "Rscript modules/enrichment/run_enrichment.R",
  sprintf("--mofa %s/mofa_top_genes.rds", dirs$ml),
  sprintf("--rf %s/rf_top_genes.rds", dirs$ml),
  sprintf("--lasso %s/lasso_selected_genes.rds", dirs$ml),
  sprintf("--rna %s/rna_for_pathway.rds", dirs$ml),
  sprintf("--rna_full %s/rna_all_expressed_genes.rds", dirs$rna),
  sprintf("--outdir %s/", dirs$enrich)
)

# ==============================================================================
# EXECUTE PIPELINE
# ==============================================================================
status <- list()

# Stage 1-4: Preprocessing (independent)
status$RNA  <- run_module("RNA Preprocessing", cmds$RNA)
status$METH <- run_module("Methylation Preprocessing", cmds$METH)
status$CNV  <- run_module("CNV Preprocessing", cmds$CNV)
status$SNV  <- run_module("SNV Preprocessing", cmds$SNV)

# Stage 5: Integration (requires all 4 preprocessing modules)
preproc_ok <- all(unlist(status))
if (preproc_ok) {
  status$INTEGRATION <- run_module("MOFA+ Integration", cmds$INTEGRATION)
} else {
  cat("\n⚠ Skipping INTEGRATION — preprocessing failures detected\n")
  status$INTEGRATION <- FALSE
}

# Stage 6: ML (requires integration)
if (isTRUE(status$INTEGRATION)) {
  status$ML <- run_module("ML Survival Analysis", cmds$ML)
} else {
  cat("\n⚠ Skipping ML — integration not available\n")
  status$ML <- FALSE
}

# Stage 7: Enrichment (requires ML)
if (isTRUE(status$ML)) {
  status$ENRICHMENT <- run_module("Pathway Enrichment", cmds$ENRICHMENT)
} else {
  cat("\n⚠ Skipping ENRICHMENT — ML not available\n")
  status$ENRICHMENT <- FALSE
}

# ==============================================================================
# COLLECT RESULTS & BUILD REPORT
# ==============================================================================
cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║            VALIDATION REPORT                                ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

report_lines <- character()
add_line <- function(line) {
  report_lines <<- c(report_lines, line)
  cat(line, "\n")
}

add_line("═══════════════════════════════════════════════════════════════")
add_line("  OmicsFlow Realistic Cohort Validation Report")
add_line(sprintf("  Generated: %s", Sys.time()))
add_line("═══════════════════════════════════════════════════════════════")
add_line("")

# --- Module Status ---
add_line("┌─────────────────────────────────────────────────────────────┐")
add_line("│  MODULE STATUS                                             │")
add_line("└─────────────────────────────────────────────────────────────┘")
for (mod in names(status)) {
  icon <- if (isTRUE(status[[mod]])) "✓ PASS" else "✗ FAIL"
  add_line(sprintf("  %-20s %s", mod, icon))
}
add_line("")

# --- Cohort Summary ---
add_line("┌─────────────────────────────────────────────────────────────┐")
add_line("│  COHORT SUMMARY                                            │")
add_line("└─────────────────────────────────────────────────────────────┘")
tryCatch({
  meta <- read.csv(file.path(data_dir, "sample_metadata.csv"), stringsAsFactors = FALSE)
  clin <- read.delim(file.path(data_dir, "custom_clinical.tsv"), stringsAsFactors = FALSE)
  
  add_line(sprintf("  Total patients:   %d", nrow(meta)))
  add_line(sprintf("  Total batches:    %d (%s)", length(unique(meta$batch)),
                    paste(unique(meta$batch), collapse = ", ")))
  add_line(sprintf("  Total centers:    %d (%s)", length(unique(meta$center)),
                    paste(unique(meta$center), collapse = ", ")))
  
  if ("subtype" %in% colnames(clin)) {
    st_tab <- table(clin$subtype)
    for (st in names(st_tab)) {
      add_line(sprintf("  Subtype %-20s: n=%d", st, st_tab[st]))
    }
  }
  add_line("")
}, error = function(e) {
  add_line(sprintf("  [Error reading cohort data: %s]", e$message))
  add_line("")
})

# --- Survival Statistics ---
add_line("┌─────────────────────────────────────────────────────────────┐")
add_line("│  SURVIVAL STATISTICS                                       │")
add_line("└─────────────────────────────────────────────────────────────┘")
tryCatch({
  clin <- read.delim(file.path(data_dir, "custom_clinical.tsv"), stringsAsFactors = FALSE)
  
  add_line(sprintf("  Overall event rate:   %.1f%%", 100 * mean(clin$OS_event)))
  add_line(sprintf("  Overall median OS:    %.0f days", median(clin$OS_time)))
  add_line(sprintf("  OS range:             [%.0f, %.0f] days", min(clin$OS_time), max(clin$OS_time)))
  
  if ("subtype" %in% colnames(clin)) {
    add_line("")
    add_line("  Per-subtype survival:")
    for (st in unique(clin$subtype)) {
      idx <- clin$subtype == st
      add_line(sprintf("    %-20s: median=%4.0f days, events=%d/%d (%.0f%%)",
                        st, median(clin$OS_time[idx]),
                        sum(clin$OS_event[idx]), sum(idx),
                        100 * mean(clin$OS_event[idx])))
    }
  }
  add_line("")
}, error = function(e) {
  add_line(sprintf("  [Error computing survival stats: %s]", e$message))
  add_line("")
})

# --- Latent Factor Summary ---
add_line("┌─────────────────────────────────────────────────────────────┐")
add_line("│  LATENT FACTOR SUMMARY (MOFA+)                             │")
add_line("└─────────────────────────────────────────────────────────────┘")
n_active_factors <- 0
tryCatch({
  qc_int_file <- file.path(dirs$integ, "qc_metrics.json")
  if (file.exists(qc_int_file)) {
    qc_int <- jsonlite::fromJSON(qc_int_file)
    n_active <- qc_int$model$active_factors_count
    n_active_factors <- if (!is.null(n_active)) n_active else 0
    add_line(sprintf("  Active factors (R² > 1%%): %d", n_active_factors))
    
    # Report R² per view
    for (nm in names(qc_int$model)) {
      if (grepl("^r2_total_", nm)) {
        view_name <- sub("^r2_total_", "", nm)
        add_line(sprintf("    %-15s R² = %.4f", view_name, qc_int$model[[nm]]))
      }
    }
  } else {
    add_line("  [Integration QC not available]")
  }
  add_line("")
}, error = function(e) {
  add_line(sprintf("  [Error reading MOFA results: %s]", e$message))
  add_line("")
})

# --- ML Performance ---
add_line("┌─────────────────────────────────────────────────────────────┐")
add_line("│  ML PERFORMANCE                                            │")
add_line("└─────────────────────────────────────────────────────────────┘")
best_c_index <- 0
tryCatch({
  qc_ml_file <- file.path(dirs$ml, "qc_metrics.json")
  if (file.exists(qc_ml_file)) {
    qc_ml <- jsonlite::fromJSON(qc_ml_file)
    perf <- qc_ml$performance
    if (!is.null(perf)) {
      for (nm in names(perf)) {
        add_line(sprintf("  %-20s: %.3f", nm, perf[[nm]]))
        if (perf[[nm]] > best_c_index) best_c_index <- perf[[nm]]
      }
    }
  } else {
    add_line("  [ML QC not available]")
  }
  add_line("")
}, error = function(e) {
  add_line(sprintf("  [Error reading ML results: %s]", e$message))
  add_line("")
})

# --- Enriched Pathways ---
add_line("┌─────────────────────────────────────────────────────────────┐")
add_line("│  ENRICHED PATHWAYS                                         │")
add_line("└─────────────────────────────────────────────────────────────┘")
total_pathways <- 0
tryCatch({
  qc_path_file <- file.path(dirs$enrich, "qc_metrics.json")
  if (file.exists(qc_path_file)) {
    qc_path <- jsonlite::fromJSON(qc_path_file)
    enr <- qc_path$enrichment
    if (!is.null(enr)) {
      for (nm in names(enr)) {
        cnt <- enr[[nm]]
        add_line(sprintf("  %-20s: %d", nm, cnt))
        total_pathways <- total_pathways + cnt
      }
    }
  } else {
    add_line("  [Enrichment QC not available]")
  }
  
  # Try to show top pathways if available
  go_bp_file <- file.path(dirs$enrich, "go_bp_results.rds")
  if (file.exists(go_bp_file)) {
    go_bp <- readRDS(go_bp_file)
    if (is.data.frame(go_bp) && nrow(go_bp) > 0) {
      add_line("")
      add_line("  Top GO Biological Processes:")
      n_show <- min(5, nrow(go_bp))
      for (i in 1:n_show) {
        desc <- if ("Description" %in% colnames(go_bp)) go_bp$Description[i] else go_bp[i, 2]
        padj <- if ("p.adjust" %in% colnames(go_bp)) go_bp$p.adjust[i] else NA
        add_line(sprintf("    %d. %s (p.adj=%.2e)", i, substr(desc, 1, 50), padj))
      }
    }
  }
  add_line("")
}, error = function(e) {
  add_line(sprintf("  [Error reading enrichment results: %s]", e$message))
  add_line("")
})

# --- Common Patients Check ---
add_line("┌─────────────────────────────────────────────────────────────┐")
add_line("│  INTEGRATION SAMPLE OVERLAP                                │")
add_line("└─────────────────────────────────────────────────────────────┘")
n_common <- 0
tryCatch({
  qc_int_file <- file.path(dirs$integ, "qc_metrics.json")
  if (file.exists(qc_int_file)) {
    qc_int <- jsonlite::fromJSON(qc_int_file)
    n_common <- qc_int$samples$common_patients
    if (!is.null(n_common)) {
      add_line(sprintf("  Common patients across all 4 omics: %d", n_common))
    }
  } else {
    add_line("  [Integration QC not available]")
  }
  add_line("")
}, error = function(e) {
  add_line(sprintf("  [Error: %s]", e$message))
  add_line("")
})

# ==============================================================================
# FINAL VERDICT
# ==============================================================================
add_line("╔══════════════════════════════════════════════════════════════╗")
add_line("║  FINAL VERDICT                                             ║")
add_line("╚══════════════════════════════════════════════════════════════╝")

checks <- list()
checks[["All 7 modules passed"]]        <- all(unlist(status))
checks[["≥100 common patients"]]        <- !is.null(n_common) && n_common >= 100
checks[["≥1 active MOFA factor"]]       <- n_active_factors >= 1
checks[["≥1 ML C-index > 0.5"]]         <- best_c_index > 0.5
checks[["≥1 enriched pathway"]]         <- total_pathways >= 1

all_pass <- all(unlist(checks))

for (ck_name in names(checks)) {
  icon <- if (isTRUE(checks[[ck_name]])) "✓" else "✗"
  add_line(sprintf("  %s %s", icon, ck_name))
}

add_line("")
if (all_pass) {
  add_line("  ╔════════════════════════════════════════╗")
  add_line("  ║     ✓ OVERALL VERDICT: PASS            ║")
  add_line("  ╚════════════════════════════════════════╝")
} else {
  add_line("  ╔════════════════════════════════════════╗")
  add_line("  ║     ✗ OVERALL VERDICT: FAIL            ║")
  add_line("  ╚════════════════════════════════════════╝")
}

elapsed_total <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
add_line(sprintf("\n  Total runtime: %s minutes", elapsed_total))
add_line(sprintf("  Report time:   %s", Sys.time()))

# Save report to file
report_file <- file.path(results_dir, "validation_report.txt")
writeLines(report_lines, report_file)
cat(sprintf("\n\nReport saved to: %s\n", report_file))

# ==============================================================================
# HTML REPORT GENERATION
# ==============================================================================
cat("\nGenerating Quarto HTML Report...\n")
cmd_report <- "quarto render reports/OmicsFlow_Report.qmd --to html -P results_dir:../results/realistic_validation -P version:\"Realistic Validation\" -P subtitle:\"Multi-Omics Validation Cohort\""
exit_code <- system(cmd_report)
if (exit_code == 0) {
  cat("✓ HTML Report successfully generated in reports/OmicsFlow_Report.html\n")
} else {
  cat("✗ HTML Report generation failed.\n")
}
