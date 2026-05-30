#!/usr/bin/env Rscript
# ==============================================================================
# diagnostic_deep_investigation.R
# Deep quantitative investigation of the three priority scientific issues
# ==============================================================================

suppressPackageStartupMessages({
  library(MOFA2)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(dplyr)
  library(jsonlite)
})

set.seed(42)

results_dir <- "results/realistic_validation"
data_dir    <- "data/realistic_cohort"

cat("\n")
cat("================================================================\n")
cat("  OmicsFlow — Deep Scientific Investigation\n")
cat("================================================================\n\n")

# ==============================================================================
# INVESTIGATION 1: ENRICHMENT BACKGROUND UNIVERSE
# ==============================================================================
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  INVESTIGATION 1: Enrichment Background Universe           ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# 1A. Show the EXACT universe currently used
rna_ml <- readRDS(file.path(results_dir, "output_ml/rna_for_pathway.rds"))
cat(sprintf("── Current Universe ──\n"))
cat(sprintf("  Matrix used for universe: output_ml/rna_for_pathway.rds\n"))
cat(sprintf("  Dimensions: %d genes × %d samples\n", nrow(rna_ml), ncol(rna_ml)))
cat(sprintf("  Gene names (first 20): %s\n", paste(head(rownames(rna_ml), 20), collapse=", ")))

# Map current universe to ENTREZ
current_universe_mapped <- bitr(rownames(rna_ml), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
cat(sprintf("  ENTREZ-mapped universe genes: %d / %d\n", nrow(current_universe_mapped), nrow(rna_ml)))

# 1B. Show what the correct universe SHOULD be
# The correct universe = all genes that passed RNA preprocessing BEFORE top-N selection
# Load the full processed RNA matrix (pre-feature-selection)
rna_processed <- readRDS(file.path(results_dir, "output_rna/rna_processed_matrix.rds"))
cat(sprintf("\n── Recommended Universe ──\n"))
cat(sprintf("  Full processed RNA matrix (Z-scored): %d genes × %d samples\n", nrow(rna_processed), ncol(rna_processed)))

# The pre-selection gene count: check QC metrics
rna_qc <- fromJSON(file.path(results_dir, "output_rna/qc_metrics.json"))
cat(sprintf("  Genes after symbol dedup (pre-selection): %d\n", rna_qc$genes$after_symbol_dedup))
cat(sprintf("  Genes after variance filter (pre-selection): %d\n", rna_qc$genes$after_variance_filter))
cat(sprintf("  Genes selected (top variable): %d\n", rna_qc$genes$selected_top_variable))

# For this synthetic run, the pre-selection pool = 1933 genes
# For real TCGA-LIHC, it would be ~15,000-18,000 genes
# The universe SHOULD be all 1933 genes, not just the 200 selected

# Ideally we'd have the pre-selection matrix saved. Let's reconstruct from raw.
cat(sprintf("\n── Universe Size Comparison ──\n"))
cat(sprintf("  Current universe:      %d genes (post top-N selection)\n", nrow(rna_ml)))
cat(sprintf("  Correct universe:      %d genes (pre top-N selection, all expressed)\n", 
            rna_qc$genes$after_variance_filter))
cat(sprintf("  Inflation factor:      %.1fx smaller than correct\n", 
            rna_qc$genes$after_variance_filter / nrow(rna_ml)))

# 1C. Estimate impact: re-run enrichment with a larger simulated universe
cat(sprintf("\n── Impact Estimation ──\n"))

# Load candidate genes
mofa_top_df    <- readRDS(file.path(results_dir, "output_ml/mofa_top_genes.rds"))
selected_genes <- readRDS(file.path(results_dir, "output_ml/rf_top_genes.rds"))
lasso_genes    <- readRDS(file.path(results_dir, "output_ml/lasso_selected_genes.rds"))

mofa_genes <- gsub("_RNA$", "", mofa_top_df$feature[1:min(50, nrow(mofa_top_df))])
rf_genes   <- selected_genes[1:min(50, length(selected_genes))]
all_candidate_genes <- unique(c(mofa_genes, rf_genes, lasso_genes))

gene_entrez <- bitr(all_candidate_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
cat(sprintf("  Candidate genes: %d unique, %d mapped to ENTREZ\n", 
            length(all_candidate_genes), nrow(gene_entrez)))

# Run enrichment with CURRENT (small) universe
cat("\n  Running GO-BP with CURRENT universe (200 genes)...\n")
go_bp_small <- enrichGO(gene = gene_entrez$ENTREZID, 
                         universe = current_universe_mapped$ENTREZID,
                         OrgDb = org.Hs.eg.db, ont = "BP", 
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05, qvalueCutoff = 0.1, 
                         readable = TRUE)
n_small <- if(!is.null(go_bp_small)) nrow(go_bp_small) else 0
cat(sprintf("  GO-BP terms (small universe): %d\n", n_small))

# Reconstruct a realistic universe: all genes that could have been expressed
# We'll map the 1933 pre-selection symbols. Since we don't have that matrix saved,
# approximate by loading the raw data and running the ID mapping pipeline.
cat("\n  Reconstructing full pre-selection gene list...\n")
rna_raw <- readRDS(file.path(data_dir, "rna.rds"))
ensembl_ids <- gsub("\\..*", "", rownames(rna_raw))
gene_symbols_full <- mapIds(org.Hs.eg.db, keys = ensembl_ids, 
                            column = "SYMBOL", keytype = "ENSEMBL", 
                            multiVals = "asNA")
valid_symbols <- unique(na.omit(gene_symbols_full))
cat(sprintf("  Full expressed gene pool: %d unique symbols\n", length(valid_symbols)))

full_universe_mapped <- bitr(valid_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
cat(sprintf("  ENTREZ-mapped full universe: %d genes\n", nrow(full_universe_mapped)))

# Run enrichment with FULL (correct) universe
cat("  Running GO-BP with CORRECT universe (~1900 genes)...\n")
go_bp_full <- enrichGO(gene = gene_entrez$ENTREZID, 
                        universe = full_universe_mapped$ENTREZID,
                        OrgDb = org.Hs.eg.db, ont = "BP", 
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05, qvalueCutoff = 0.1, 
                        readable = TRUE)
n_full <- if(!is.null(go_bp_full)) nrow(go_bp_full) else 0
cat(sprintf("  GO-BP terms (full universe): %d\n", n_full))

# Also try genome-wide universe (NULL = all annotated genes)
cat("  Running GO-BP with GENOME-WIDE universe (default ~20,000 genes)...\n")
go_bp_genome <- enrichGO(gene = gene_entrez$ENTREZID, 
                          universe = NULL,  # genome-wide default
                          OrgDb = org.Hs.eg.db, ont = "BP", 
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05, qvalueCutoff = 0.1, 
                          readable = TRUE)
n_genome <- if(!is.null(go_bp_genome)) nrow(go_bp_genome) else 0
cat(sprintf("  GO-BP terms (genome-wide): %d\n", n_genome))

# Compare top terms across universes
cat("\n── Side-by-Side Comparison ──\n")
cat(sprintf("  %-40s  %-12s  %-12s  %-12s\n", "GO Term", "Small(200)", "Full(~1900)", "Genome-wide"))
cat(sprintf("  %s\n", strrep("-", 80)))

# Get all unique terms
all_terms <- character()
if (n_small > 0) all_terms <- c(all_terms, as.data.frame(go_bp_small)$ID)
if (n_full > 0)  all_terms <- c(all_terms, as.data.frame(go_bp_full)$ID)
if (n_genome > 0) all_terms <- c(all_terms, as.data.frame(go_bp_genome)$ID)
all_terms <- unique(all_terms)

# Build comparison table
for (tid in head(all_terms, 15)) {
  desc <- ""
  p_small <- p_full <- p_genome <- NA
  
  if (n_small > 0) {
    df <- as.data.frame(go_bp_small)
    idx <- match(tid, df$ID)
    if (!is.na(idx)) { desc <- substr(df$Description[idx], 1, 38); p_small <- df$p.adjust[idx] }
  }
  if (n_full > 0) {
    df <- as.data.frame(go_bp_full)
    idx <- match(tid, df$ID)
    if (!is.na(idx)) { if (desc == "") desc <- substr(df$Description[idx], 1, 38); p_full <- df$p.adjust[idx] }
  }
  if (n_genome > 0) {
    df <- as.data.frame(go_bp_genome)
    idx <- match(tid, df$ID)
    if (!is.na(idx)) { if (desc == "") desc <- substr(df$Description[idx], 1, 38); p_genome <- df$p.adjust[idx] }
  }
  
  fmt_s <- if(is.na(p_small)) "  ---" else sprintf("  %.2e", p_small)
  fmt_f <- if(is.na(p_full))  "  ---" else sprintf("  %.2e", p_full)
  fmt_g <- if(is.na(p_genome)) "  ---" else sprintf("  %.2e", p_genome)
  
  cat(sprintf("  %-40s  %-12s  %-12s  %-12s\n", desc, fmt_s, fmt_f, fmt_g))
}

cat(sprintf("\n── Summary ──\n"))
cat(sprintf("  Universe=200 genes:    %d significant GO-BP terms\n", n_small))
cat(sprintf("  Universe=~1900 genes:  %d significant GO-BP terms\n", n_full))
cat(sprintf("  Universe=genome-wide:  %d significant GO-BP terms\n", n_genome))
if (n_small > n_full) {
  cat(sprintf("  CONFIRMED: Small universe INFLATES results by %d terms (%.0f%% inflation)\n",
              n_small - n_full, 100 * (n_small - n_full) / max(n_full, 1)))
}

# KEGG comparison
cat("\n  Running KEGG with CORRECT universe...\n")
kegg_full <- tryCatch({
  enrichKEGG(gene = gene_entrez$ENTREZID, universe = full_universe_mapped$ENTREZID,
             organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1)
}, error = function(e) NULL)
n_kegg_full <- if(!is.null(kegg_full)) nrow(kegg_full) else 0
cat(sprintf("  KEGG pathways (full universe): %d\n", n_kegg_full))


# ==============================================================================
# INVESTIGATION 2: BATCH-SUBTYPE CONFOUNDING
# ==============================================================================
cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  INVESTIGATION 2: Batch–Subtype Confounding                ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Load metadata and clinical
meta <- read.csv(file.path(data_dir, "sample_metadata.csv"), stringsAsFactors = FALSE)
clin <- read.delim(file.path(data_dir, "custom_clinical.tsv"), stringsAsFactors = FALSE)

# Merge
merged <- merge(meta, clin, by.x = "patient_id", by.y = "patient")

# 2A. Cross-tabulation: Batch × Subtype
cat("── Batch × Subtype Cross-Tabulation ──\n")
ct <- table(Batch = merged$batch, Subtype = merged$subtype)
print(ct)

# Chi-squared test
chi_test <- chisq.test(ct)
cat(sprintf("\n  Chi-squared statistic: %.2f\n", chi_test$statistic))
cat(sprintf("  Degrees of freedom:   %d\n", chi_test$parameter))
cat(sprintf("  p-value:              %.2e\n", chi_test$p.value))

# Cramér's V
n <- sum(ct)
k <- min(nrow(ct), ncol(ct))
cramers_v <- sqrt(chi_test$statistic / (n * (k - 1)))
cat(sprintf("  Cramér's V:           %.3f\n", cramers_v))
cat(sprintf("  Interpretation:       %s\n", 
            if(cramers_v < 0.1) "Negligible" 
            else if(cramers_v < 0.3) "Weak"
            else if(cramers_v < 0.5) "Moderate"
            else "Strong"))

# 2B. Center × Subtype
cat("\n── Center × Subtype Cross-Tabulation ──\n")
ct2 <- table(Center = merged$center, Subtype = merged$subtype)
print(ct2)
chi2 <- chisq.test(ct2)
cramers_v2 <- sqrt(chi2$statistic / (sum(ct2) * (min(nrow(ct2), ncol(ct2)) - 1)))
cat(sprintf("  Cramér's V: %.3f (p = %.2e)\n", cramers_v2, chi2$p.value))

# 2C. Quantify biological signal removal
cat("\n── Biological Signal Before/After Batch Correction ──\n")

# Load the raw RNA, run minimal preprocessing, then compare marker gene variance
# before and after batch correction
rna_raw <- readRDS(file.path(data_dir, "rna.rds"))
# The first 150 genes are subtype markers (50 prolif + 50 mesench + 50 immune)
# per generate_realistic_cohort.R lines 214-216

# Load the processed (post-batch-correction) matrix
rna_processed <- readRDS(file.path(results_dir, "output_rna/rna_processed_matrix.rds"))

# We need to check: do marker genes retain subtype separation after batch correction?
# Get patient-subtype mapping
subtype_map <- setNames(merged$subtype, merged$patient_id)

# For the processed matrix, compute per-gene ANOVA F-statistic (subtype effect)
common_patients <- intersect(colnames(rna_processed), names(subtype_map))
rna_sub <- rna_processed[, common_patients, drop = FALSE]
subtypes_ordered <- subtype_map[common_patients]

# F-test for each gene: expression ~ subtype
cat(sprintf("  Testing %d genes × %d patients for subtype signal...\n", 
            nrow(rna_sub), length(common_patients)))

f_stats <- apply(rna_sub, 1, function(x) {
  tryCatch({
    res <- summary(aov(x ~ subtypes_ordered))
    res[[1]]$`F value`[1]
  }, error = function(e) NA)
})

cat(sprintf("  Genes with significant subtype effect (p < 0.05): %d / %d\n",
            sum(f_stats > qf(0.95, 2, length(common_patients) - 3), na.rm = TRUE),
            length(f_stats)))
cat(sprintf("  Median F-statistic: %.2f\n", median(f_stats, na.rm = TRUE)))
cat(sprintf("  Max F-statistic: %.2f\n", max(f_stats, na.rm = TRUE)))

# Top 10 genes by subtype effect
top_f <- sort(f_stats, decreasing = TRUE)[1:10]
cat("\n  Top 10 genes by subtype effect (F-statistic):\n")
for (i in seq_along(top_f)) {
  cat(sprintf("    %2d. %-15s F = %.2f\n", i, names(top_f)[i], top_f[i]))
}

# Compare: what fraction of the original 200 genes still show subtype signal?
n_sig <- sum(f_stats > qf(0.95, 2, length(common_patients) - 3), na.rm = TRUE)
cat(sprintf("\n  Verdict: %d / %d genes (%.0f%%) retain significant subtype signal\n",
            n_sig, nrow(rna_sub), 100 * n_sig / nrow(rna_sub)))
cat(sprintf("  (after Z-scoring + batch correction)\n"))

# 2D. Batch effect on subtype marker genes specifically
# We can check if batch correction differentially reduces variance in marker vs. non-marker genes
# Since this is post-processed, let's look at inter-subtype effect sizes
cat("\n── Subtype Separation Assessment (Cohen's d for marker genes) ──\n")

# Top F-stat genes = likely marker genes. Bottom = likely noise.
top_markers <- names(sort(f_stats, decreasing = TRUE))[1:50]
bottom_genes <- names(sort(f_stats, decreasing = FALSE))[1:50]

# Mean absolute Cohen's d between subtypes for top vs. bottom genes
cohens_d <- function(x, groups) {
  levs <- unique(groups)
  ds <- numeric()
  for (i in 1:(length(levs)-1)) {
    for (j in (i+1):length(levs)) {
      g1 <- x[groups == levs[i]]
      g2 <- x[groups == levs[j]]
      pooled_sd <- sqrt(((length(g1)-1)*var(g1) + (length(g2)-1)*var(g2)) / (length(g1)+length(g2)-2))
      ds <- c(ds, abs(mean(g1) - mean(g2)) / pooled_sd)
    }
  }
  mean(ds)
}

d_top <- mean(sapply(top_markers, function(g) {
  cohens_d(rna_sub[g, ], subtypes_ordered)
}))
d_bottom <- mean(sapply(bottom_genes, function(g) {
  cohens_d(rna_sub[g, ], subtypes_ordered)
}))

cat(sprintf("  Mean |Cohen's d| (top 50 marker genes):   %.3f\n", d_top))
cat(sprintf("  Mean |Cohen's d| (bottom 50 noise genes): %.3f\n", d_bottom))
cat(sprintf("  Ratio (signal/noise):                     %.1fx\n", d_top / d_bottom))
cat(sprintf("  Interpretation: %s\n",
            if(d_top / d_bottom > 3) "Strong biological signal retained ✓"
            else if(d_top / d_bottom > 1.5) "Moderate signal retained ⚠"
            else "Signal may be attenuated ✗"))

# 2E. Recommendation
cat("\n── Recommendation for Dataset-Agnostic Batch Correction ──\n")
cat("  Current:  removeBatchEffect(rna_mat, batch = batch_info)\n")
cat("  Problem:  No design matrix → biological covariates unprotected\n")
cat("\n  Generic solution (dataset-agnostic):\n")
cat("    1. If metadata has 'subtype' or 'sample_class' → use as covariate\n")
cat("    2. If no known covariates → use design = ~1 (intercept only)\n")
cat("    3. Never pass subtype-confounded variables as batch\n")
cat("  Code:\n")
cat("    mod <- model.matrix(~1, data = data.frame(row.names = colnames(rna_mat)))\n")
cat("    rna_mat <- removeBatchEffect(rna_mat, batch = batch_info, design = mod)\n")
cat("  Note: ~1 is equivalent to current behavior. The real fix is:\n")
cat("    mod <- model.matrix(~subtype, data = sample_info)  # if subtype known\n")
cat("  For unknown subtypes, use PCA-based surrogate:\n")
cat("    pca <- prcomp(t(rna_mat))\n")
cat("    mod <- model.matrix(~pca$x[,1:2])  # protect top 2 PCs\n")


# ==============================================================================
# INVESTIGATION 3: MOFA CONVERGENCE + CNV DOMINANCE
# ==============================================================================
cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  INVESTIGATION 3: MOFA Convergence & CNV Dominance         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# 3A. Load MOFA model
mofa_obj <- readRDS(file.path(results_dir, "output_integration/mofa_model.rds"))
mofa <- mofa_obj$model

# 3B. Convergence diagnostics
cat("── MOFA Convergence Diagnostics ──\n")

# Extract ELBO trace
elbo <- NULL
tryCatch({
  # MOFA2 stores training stats
  training_stats <- mofa@training_stats
  if (!is.null(training_stats$elbo)) {
    elbo <- training_stats$elbo
    cat(sprintf("  ELBO values recovered: %d iterations\n", length(elbo)))
    
    # Remove NA/NULL entries
    elbo <- elbo[!is.na(elbo)]
    cat(sprintf("  Valid ELBO entries: %d\n", length(elbo)))
    
    if (length(elbo) >= 3) {
      cat(sprintf("  Initial ELBO:  %.2f\n", elbo[1]))
      cat(sprintf("  Final ELBO:    %.2f\n", tail(elbo, 1)))
      cat(sprintf("  Total change:  %.2f\n", tail(elbo, 1) - elbo[1]))
      
      # Relative change in last 10% of iterations
      n_elbo <- length(elbo)
      last_10pct <- max(1, floor(n_elbo * 0.9)):n_elbo
      if (length(last_10pct) >= 2) {
        delta_last <- diff(elbo[last_10pct])
        abs_deltas <- abs(delta_last)
        cat(sprintf("\n  Last 10%% of iterations (iter %d-%d):\n", 
                    last_10pct[1], tail(last_10pct, 1)))
        cat(sprintf("    Mean |ΔELBO|:  %.6f\n", mean(abs_deltas)))
        cat(sprintf("    Max  |ΔELBO|:  %.6f\n", max(abs_deltas)))
        cat(sprintf("    Min  |ΔELBO|:  %.6f\n", min(abs_deltas)))
        
        # Relative change
        rel_change <- abs_deltas / abs(elbo[last_10pct[-length(last_10pct)]])
        cat(sprintf("    Mean relative change: %.2e\n", mean(rel_change)))
        
        # Standard convergence criterion: |ΔELBO/ELBO| < 1e-5
        converged <- all(rel_change < 1e-5)
        cat(sprintf("    Convergence criterion (|ΔELBO/ELBO| < 1e-5): %s\n",
                    if(converged) "MET ✓" else "NOT MET ✗"))
        
        # Alternative: check if ELBO is still changing meaningfully
        almost_converged <- all(rel_change < 1e-4)
        cat(sprintf("    Relaxed criterion (|ΔELBO/ELBO| < 1e-4): %s\n",
                    if(almost_converged) "MET ✓" else "NOT MET ✗"))
      }
      
      # Plot ELBO trajectory summary as text
      cat("\n  ELBO Trajectory (sampled):\n")
      sample_idx <- unique(c(1, round(seq(1, n_elbo, length.out = 10)), n_elbo))
      for (idx in sample_idx) {
        pct <- 100 * idx / n_elbo
        cat(sprintf("    Iter %3d (%5.1f%%): ELBO = %.2f\n", idx, pct, elbo[idx]))
      }
    }
  } else {
    cat("  ELBO trace not found in training_stats\n")
  }
}, error = function(e) {
  cat(sprintf("  Error accessing ELBO: %s\n", e$message))
})

# Try alternative ELBO access
if (is.null(elbo)) {
  tryCatch({
    cat("  Trying alternative ELBO access...\n")
    # Check slots
    cat(sprintf("  MOFA slots: %s\n", paste(slotNames(mofa), collapse=", ")))
    if ("training_stats" %in% slotNames(mofa)) {
      ts <- slot(mofa, "training_stats")
      cat(sprintf("  training_stats names: %s\n", paste(names(ts), collapse=", ")))
    }
  }, error = function(e) {
    cat(sprintf("  Error: %s\n", e$message))
  })
}

# 3C. Variance explained — detailed per-factor per-view
cat("\n── Variance Explained: Detailed Breakdown ──\n")
var_exp <- get_variance_explained(mofa)

# Per-factor R² for each view
r2_per_factor <- var_exp$r2_per_factor[[1]]
cat("\n  Per-Factor R² (%):\n")
cat(sprintf("  %-10s", ""))
for (v in colnames(r2_per_factor)) cat(sprintf("  %-12s", v))
cat("\n")
for (i in 1:nrow(r2_per_factor)) {
  cat(sprintf("  Factor %-2d", i))
  for (j in 1:ncol(r2_per_factor)) {
    cat(sprintf("  %10.4f%%", r2_per_factor[i, j]))
  }
  cat(sprintf("  | mean = %.4f%%", mean(r2_per_factor[i, ])))
  cat("\n")
}

# Total R² per view
r2_total <- var_exp$r2_total[[1]]
cat("\n  Total R² per view:\n")
for (v in names(r2_total)) {
  cat(sprintf("    %-15s: %.4f%%\n", v, r2_total[[v]]))
}

# 3D. CNV DOMINANCE ROOT CAUSE ANALYSIS
cat("\n\n── CNV Dominance Root Cause Analysis ──\n")

# Load all four processed matrices
rna_mat  <- readRDS(file.path(results_dir, "output_rna/rna_processed_matrix.rds"))
meth_mat <- readRDS(file.path(results_dir, "output_meth/methylation_processed_matrix.rds"))
cnv_mat  <- readRDS(file.path(results_dir, "output_cnv/cnv_processed_matrix.rds"))
snv_mat  <- readRDS(file.path(results_dir, "output_snv/snv_processed_matrix.rds"))

cat("\n  Input matrix properties:\n")
cat(sprintf("    %-15s: %5d features × %3d samples | class=%s\n", 
            "RNA", nrow(rna_mat), ncol(rna_mat), class(rna_mat)[1]))
cat(sprintf("    %-15s: %5d features × %3d samples | class=%s\n", 
            "Methylation", nrow(meth_mat), ncol(meth_mat), class(meth_mat)[1]))
cat(sprintf("    %-15s: %5d features × %3d samples | class=%s\n", 
            "CNV", nrow(cnv_mat), ncol(cnv_mat), class(cnv_mat)[1]))
cat(sprintf("    %-15s: %5d features × %3d samples | class=%s\n", 
            "SNV", nrow(snv_mat), ncol(snv_mat), class(snv_mat)[1]))

# Hypothesis 1: Feature variance distribution differences
cat("\n  [Hypothesis 1] Feature variance distributions:\n")
var_rna  <- apply(rna_mat, 1, var)
var_meth <- apply(meth_mat, 1, var)
var_cnv  <- apply(cnv_mat, 1, var)
var_snv  <- apply(snv_mat, 1, var)

cat(sprintf("    RNA:   median=%.4f, mean=%.4f, max=%.4f\n", 
            median(var_rna), mean(var_rna), max(var_rna)))
cat(sprintf("    Meth:  median=%.4f, mean=%.4f, max=%.4f\n", 
            median(var_meth), mean(var_meth), max(var_meth)))
cat(sprintf("    CNV:   median=%.4f, mean=%.4f, max=%.4f\n", 
            median(var_cnv), mean(var_cnv), max(var_cnv)))
cat(sprintf("    SNV:   median=%.4f, mean=%.4f, max=%.4f\n", 
            median(var_snv), mean(var_snv), max(var_snv)))

# Hypothesis 2: RNA is Z-scored (var~1), CNV is not
cat("\n  [Hypothesis 2] Pre-scaling status:\n")
cat(sprintf("    RNA  global var: %.4f (expected ~1.0 if Z-scored)\n", var(as.numeric(rna_mat))))
cat(sprintf("    Meth global var: %.4f\n", var(as.numeric(meth_mat))))
cat(sprintf("    CNV  global var: %.4f\n", var(as.numeric(cnv_mat))))
cat(sprintf("    SNV  global var: %.4f (expected ~0.05 if binary)\n", var(as.numeric(snv_mat))))

# Value range analysis
cat("\n  Value ranges:\n")
cat(sprintf("    RNA:   [%.4f, %.4f]\n", min(rna_mat), max(rna_mat)))
cat(sprintf("    Meth:  [%.4f, %.4f]\n", min(meth_mat), max(meth_mat)))
cat(sprintf("    CNV:   [%.4f, %.4f]\n", min(cnv_mat), max(cnv_mat)))
cat(sprintf("    SNV:   [%.0f, %.0f] (binary)\n", min(snv_mat), max(snv_mat)))

# Hypothesis 3: Total sum-of-squares per view (what MOFA sees)
# MOFA with scale_views=TRUE divides each view by its total variance
# R² = fraction of that view-specific variance captured
cat("\n  [Hypothesis 3] Total sum-of-squares per view:\n")

# Find common patients
common_p <- Reduce(intersect, list(colnames(rna_mat), colnames(meth_mat), 
                                    colnames(cnv_mat), colnames(snv_mat)))
cat(sprintf("    Common patients: %d\n", length(common_p)))

ss_rna  <- sum(rna_mat[, common_p]^2)
ss_meth <- sum(meth_mat[, common_p]^2)
ss_cnv  <- sum(cnv_mat[, common_p]^2)
ss_snv  <- sum(snv_mat[, common_p]^2)

total_ss <- ss_rna + ss_meth + ss_cnv + ss_snv
cat(sprintf("    RNA  SS = %.2e (%.1f%% of total)\n", ss_rna, 100*ss_rna/total_ss))
cat(sprintf("    Meth SS = %.2e (%.1f%% of total)\n", ss_meth, 100*ss_meth/total_ss))
cat(sprintf("    CNV  SS = %.2e (%.1f%% of total)\n", ss_cnv, 100*ss_cnv/total_ss))
cat(sprintf("    SNV  SS = %.2e (%.1f%% of total)\n", ss_snv, 100*ss_snv/total_ss))

# Hypothesis 4: Effective dimensionality (PCA)
cat("\n  [Hypothesis 4] Effective dimensionality (PCA on each view):\n")
for (view_name in c("RNA", "Meth", "CNV", "SNV")) {
  mat <- switch(view_name,
                "RNA" = rna_mat[, common_p],
                "Meth" = meth_mat[, common_p],
                "CNV" = cnv_mat[, common_p],
                "SNV" = snv_mat[, common_p])
  
  pca_res <- tryCatch({
    prcomp(t(mat), scale. = FALSE, rank. = min(10, ncol(mat)-1, nrow(mat)-1))
  }, error = function(e) NULL)
  
  if (!is.null(pca_res)) {
    var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
    pc1_var <- var_explained[1] * 100
    pc2_var <- if(length(var_explained) >= 2) var_explained[2] * 100 else 0
    cum_5 <- if(length(var_explained) >= 5) sum(var_explained[1:5]) * 100 else sum(var_explained) * 100
    cat(sprintf("    %-5s: PC1=%.1f%%, PC2=%.1f%%, Top5=%.1f%%\n", 
                view_name, pc1_var, pc2_var, cum_5))
  }
}

# Hypothesis 5: Subtype signal strength per view
cat("\n  [Hypothesis 5] Subtype signal strength per view (MANOVA-like):\n")
for (view_name in c("RNA", "Meth", "CNV", "SNV")) {
  mat <- switch(view_name,
                "RNA" = rna_mat[, common_p],
                "Meth" = meth_mat[, common_p],
                "CNV" = cnv_mat[, common_p],
                "SNV" = snv_mat[, common_p])
  
  subs <- subtype_map[common_p]
  
  # Per-feature F-statistics
  f_vals <- apply(mat, 1, function(x) {
    tryCatch({
      summary(aov(x ~ subs))[[1]]$`F value`[1]
    }, error = function(e) NA)
  })
  
  n_sig <- sum(f_vals > qf(0.95, 2, length(common_p) - 3), na.rm = TRUE)
  cat(sprintf("    %-5s: %d / %d features with significant subtype effect (%.0f%%)\n",
              view_name, n_sig, nrow(mat), 100 * n_sig / nrow(mat)))
  cat(sprintf("           Median F = %.2f, Mean F = %.2f\n", 
              median(f_vals, na.rm=TRUE), mean(f_vals, na.rm=TRUE)))
}

# 3E. Root cause verdict
cat("\n── CNV Dominance: Root Cause Verdict ──\n")
cat("  The R² reported by MOFA is VIEW-SPECIFIC (% of each view's variance\n")
cat("  explained by the factors). It does NOT mean CNV drives 70% of the\n")
cat("  integrated model. It means 70% of the CNV view's internal variance\n")
cat("  is captured by the 5 factors.\n\n")
cat("  This is expected when:\n")
cat("  1. CNV has strong structured variation (subtype-specific amplifications)\n")
cat("  2. CNV has fewer features with higher per-feature variance\n")
cat("  3. CNV features are more correlated (fewer independent dimensions)\n")
cat("  4. scale_views=TRUE normalizes total view variance but does not\n")
cat("     change internal structure\n")


# ==============================================================================
# INVESTIGATION 4: MOFA RECOMMENDED SETTINGS
# ==============================================================================
cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  Recommended MOFA Settings for Publication                 ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("  Current settings:\n")
cat("    num_factors = 5\n")
cat("    maxiter     = 100\n")
cat("    seed        = 42\n")
cat("    scale_views = TRUE\n\n")

cat("  Recommended settings:\n")
cat("    num_factors = 15    (let ARD prune; starting too low may miss structure)\n")
cat("    maxiter     = 1000  (MOFA2 default; ensures convergence)\n")
cat("    seed        = 42    (keep for reproducibility)\n")
cat("    scale_views = TRUE  (keep; critical for multi-scale data)\n")
cat("    convergence_mode = 'slow'  (stricter convergence)\n\n")

cat("  Justification:\n")
cat("  - num_factors=15: MOFA2 uses Automatic Relevance Determination (ARD)\n")
cat("    to prune unused factors. Starting with 15 allows discovery of\n")
cat("    fine-grained structure while inactive factors are automatically zeroed.\n")
cat("  - maxiter=1000: The MOFA2 default. 100 iterations is only appropriate\n")
cat("    for quick exploratory runs, not publication.\n")
cat("  - convergence_mode='slow': Uses tighter ELBO tolerance (1e-6 vs 1e-4).\n")

cat("\n\n════════════════════════════════════════════════════════════════\n")
cat("  Investigation Complete\n")
cat("════════════════════════════════════════════════════════════════\n\n")
