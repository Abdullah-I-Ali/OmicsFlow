cat("--- Simulating End-to-End Execution ---\n")

run_module <- function(name, cmd) {
  cat("\nRunning", name, "...\n")
  cat("CMD:", cmd, "\n")
  res <- system(cmd, intern=FALSE)
  if (res != 0) {
    cat(name, "FAILED with exit code", res, "\n")
    return(FALSE)
  }
  cat(name, "PASSED\n")
  return(TRUE)
}

dir.create("results/synthetic_validation/output_rna", recursive=TRUE, showWarnings=FALSE)
cmd_rna <- "Rscript modules/rna/preprocess_rna.R --input data/synthetic_cohort/rna.rds --metadata data/synthetic_cohort/sample_metadata.csv --outdir results/synthetic_validation/output_rna/ --n-top 200 --cor-low -1.0 --cor-high 1.0"

dir.create("results/synthetic_validation/output_meth", recursive=TRUE, showWarnings=FALSE)
cmd_meth <- "Rscript modules/methylation/preprocess_meth.R --input data/synthetic_cohort/meth.rds --metadata data/synthetic_cohort/sample_metadata.csv --clinical data/synthetic_cohort/custom_clinical.tsv --clinical_map data/synthetic_cohort/clinical_map.json --outdir results/synthetic_validation/output_meth/ --n-top 500 --knn-k 2"

dir.create("results/synthetic_validation/output_cnv", recursive=TRUE, showWarnings=FALSE)
cmd_cnv <- "Rscript modules/cnv/preprocess_cnv.R --input data/synthetic_cohort/cnv.rds --metadata data/synthetic_cohort/sample_metadata.csv --outdir results/synthetic_validation/output_cnv/ --ntop 500 --cache fake_cache.rds"

dir.create("results/synthetic_validation/output_snv", recursive=TRUE, showWarnings=FALSE)
cmd_snv <- "Rscript modules/snv/preprocess_snv.R --input data/synthetic_cohort/snv.rds --metadata data/synthetic_cohort/sample_metadata.csv --outdir results/synthetic_validation/output_snv/"

dir.create("results/synthetic_validation/output_integration", recursive=TRUE, showWarnings=FALSE)
cmd_int <- "Rscript modules/integration/run_integration.R --rna results/synthetic_validation/output_rna/rna_processed_matrix.rds --meth results/synthetic_validation/output_meth/methylation_processed_matrix.rds --cnv results/synthetic_validation/output_cnv/cnv_processed_matrix.rds --snv results/synthetic_validation/output_snv/snv_processed_matrix.rds --metadata data/synthetic_cohort/sample_metadata.csv --outdir results/synthetic_validation/output_integration/ --factors 5 --iter 100"

dir.create("results/synthetic_validation/output_ml", recursive=TRUE, showWarnings=FALSE)
cmd_ml <- "Rscript modules/ml/run_ml.R --mofa results/synthetic_validation/output_integration/mofa_model.rds --rna results/synthetic_validation/output_rna/rna_ml.rds --clinical data/synthetic_cohort/custom_clinical.tsv --metadata data/synthetic_cohort/sample_metadata.csv --clinical_map data/synthetic_cohort/clinical_map.json --outdir results/synthetic_validation/output_ml/ --mofa_prefilter 50 --final_features 20"

dir.create("results/synthetic_validation/output_enrichment", recursive=TRUE, showWarnings=FALSE)
cmd_enrich <- "Rscript modules/enrichment/run_enrichment.R --mofa results/synthetic_validation/output_ml/mofa_top_genes.rds --rf results/synthetic_validation/output_ml/rf_top_genes.rds --lasso results/synthetic_validation/output_ml/lasso_selected_genes.rds --rna results/synthetic_validation/output_ml/rna_for_pathway.rds --outdir results/synthetic_validation/output_enrichment/"

status <- list()
status$RNA <- run_module("RNA", cmd_rna)
status$METH <- run_module("METH", cmd_meth)
status$CNV <- run_module("CNV", cmd_cnv)
status$SNV <- run_module("SNV", cmd_snv)
status$INTEGRATION <- if(all(unlist(status))) run_module("INTEGRATION", cmd_int) else FALSE
status$ML <- if(status$INTEGRATION) run_module("ML", cmd_ml) else FALSE
status$ENRICHMENT <- if(status$ML) run_module("ENRICHMENT", cmd_enrich) else FALSE

# Report generation
cat("\nRunning REPORTING ...\n")
# We just call quarto render. If quarto is not installed, it will fail.
cmd_report <- "quarto render reports/OmicsFlow_Report.qmd --to html -P results_dir:results/synthetic_validation"
cat("CMD:", cmd_report, "\n")
res <- system(cmd_report, intern=FALSE)
if (res != 0) {
    cat("REPORTING FAILED with exit code", res, "\n")
    status$REPORTING <- FALSE
} else {
    cat("REPORTING PASSED\n")
    status$REPORTING <- TRUE
}

cat("\n--- FINAL STATUS ---\n")
print(status)
