#!/usr/bin/env Rscript
# ==============================================================================
# run_enrichment.R — Pathway Enrichment Analysis
# OmicsFlow | Phase 7: Enrichment Module
# ==============================================================================
#
# PURPOSE:
#   Executes the scientific steps for Gene Ontology (GO) and KEGG pathway 
#   enrichment on the multi-omics prognostic gene candidates.
#   
# METHODOLOGY:
#   1. Load predictive features from ML & MOFA (RF, LASSO, MOFA Top Weights).
#   2. Load RNA matrix to serve as the background universe (liver-specific).
#   3. Map SYMBOL -> ENTREZID via bitr.
#   4. Run GO Enrichment (BP, MF, CC) using clusterProfiler.
#   5. Run KEGG Enrichment.
#   6. Verify Extracellular Matrix (ECM) signature.
#   7. Generate Dotplots, Gene-Concept Networks, and Enrichment Maps.
#
# ==============================================================================
# USAGE:
#   Rscript run_enrichment.R \
#     --mofa      results/ml/mofa_top_genes.rds \
#     --rf        results/ml/rf_top_genes.rds \
#     --lasso     results/ml/lasso_selected_genes.rds \
#     --rna       results/ml/rna_for_pathway.rds \
#     --outdir    results/enrichment/
# ==============================================================================

suppressPackageStartupMessages(library(optparse))

# ------------------------------------------------------------------------------
# COMMAND-LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--mofa"), type = "character", default = "results/ml/mofa_top_genes.rds",
              help = "Path to MOFA top genes (.rds)"),
  make_option(c("--rf"), type = "character", default = "results/ml/rf_top_genes.rds",
              help = "Path to RF top genes (.rds)"),
  make_option(c("--lasso"), type = "character", default = "results/ml/lasso_selected_genes.rds",
              help = "Path to LASSO genes (.rds)"),
  make_option(c("--rna"), type = "character", default = "results/ml/rna_for_pathway.rds",
              help = "Path to RNA matrix for background universe (.rds)"),
  make_option(c("-o", "--outdir"), type = "character", default = "results/enrichment/",
              help = "Output directory [default= %default]")
)

opt_parser <- OptionParser(
  usage = "Usage: %prog [options]",
  option_list = option_list,
  description = "OmicsFlow Phase 7: Enrichment Pipeline"
)
opt <- parse_args(opt_parser)

# ------------------------------------------------------------------------------
# SOURCE MODULE FILES
# ------------------------------------------------------------------------------
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) "modules/enrichment"
)
source(file.path(script_dir, "utils_enrichment.R"))
source(file.path(script_dir, "qc_enrichment.R"))
source(file.path(script_dir, "export_enrichment.R"))

# ------------------------------------------------------------------------------
# INITIALISE
# ------------------------------------------------------------------------------
tryCatch({
  load_pathway_packages()
  set.seed(42)
  
  if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
  }
  
  qc_metrics <- init_pathway_qc()
  
  path_banner("OmicsFlow v1.0.0 | Pathway Enrichment")
  path_msg(sprintf("Start time : %s", Sys.time()))
  path_msg(sprintf("Output dir : %s", opt$outdir))
  
  # ==============================================================================
  # 1. LOAD CANDIDATE GENES
  # ==============================================================================
  path_step(1, "Collecting Candidate Genes")
  
  mofa_top_df    <- readRDS(opt$mofa)
  selected_genes <- readRDS(opt$rf)
  lasso_genes    <- readRDS(opt$lasso)
  rna            <- readRDS(opt$rna)
  
  path_msg("Variables loaded successfully")
  
  mofa_genes       <- str_replace(mofa_top_df$feature[1:min(50, nrow(mofa_top_df))], "_RNA$", "")
  rf_genes         <- selected_genes[1:min(50, length(selected_genes))]
  lasso_genes_list <- lasso_genes
  
  all_candidate_genes <- unique(c(mofa_genes, rf_genes, lasso_genes_list))
  
  path_msg(sprintf("MOFA top genes:   %d", length(mofa_genes)), level = "DETAILS")
  path_msg(sprintf("RF top genes:     %d", length(rf_genes)), level = "DETAILS")
  path_msg(sprintf("LASSO genes:      %d", length(lasso_genes_list)), level = "DETAILS")
  path_msg(sprintf("Total unique:     %d", length(all_candidate_genes)))
  
  qc_metrics <- add_pathway_qc(qc_metrics, "genes", "total_candidates", length(all_candidate_genes))
  
  # ==============================================================================
  # 2. GENE ID MAPPING
  # ==============================================================================
  path_step(2, "Gene ID Mapping (Symbol -> Entrez)")
  
  gene_entrez <- bitr(all_candidate_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  path_msg(sprintf("Candidate genes mapped: %d / %d", nrow(gene_entrez), length(all_candidate_genes)))
  
  universe_entrez <- bitr(rownames(rna), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  path_msg(sprintf("Universe genes mapped:  %d / %d", nrow(universe_entrez), nrow(rna)))
  
  qc_metrics <- add_pathway_qc(qc_metrics, "genes", "mapped_candidates", nrow(gene_entrez))
  qc_metrics <- add_pathway_qc(qc_metrics, "genes", "mapped_universe", nrow(universe_entrez))
  
  # ==============================================================================
  # 3. GO ENRICHMENT ANALYSIS
  # ==============================================================================
  path_step(3, "GO Enrichment Analysis")
  
  go_bp <- enrichGO(gene = gene_entrez$ENTREZID, universe = universe_entrez$ENTREZID, 
                    OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", 
                    pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE)
  path_msg(sprintf("GO-BP terms found: %d", if(!is.null(go_bp)) nrow(go_bp) else 0))
  
  go_mf <- enrichGO(gene = gene_entrez$ENTREZID, universe = universe_entrez$ENTREZID, 
                    OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", 
                    pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE)
  path_msg(sprintf("GO-MF terms found: %d", if(!is.null(go_mf)) nrow(go_mf) else 0))
  
  go_cc <- enrichGO(gene = gene_entrez$ENTREZID, universe = universe_entrez$ENTREZID, 
                    OrgDb = org.Hs.eg.db, ont = "CC", pAdjustMethod = "BH", 
                    pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE)
  path_msg(sprintf("GO-CC terms found: %d", if(!is.null(go_cc)) nrow(go_cc) else 0))
  
  qc_metrics <- add_pathway_qc(qc_metrics, "enrichment", "go_bp_terms", if(!is.null(go_bp)) nrow(go_bp) else 0)
  qc_metrics <- add_pathway_qc(qc_metrics, "enrichment", "go_mf_terms", if(!is.null(go_mf)) nrow(go_mf) else 0)
  qc_metrics <- add_pathway_qc(qc_metrics, "enrichment", "go_cc_terms", if(!is.null(go_cc)) nrow(go_cc) else 0)
  
  # ==============================================================================
  # 4. KEGG PATHWAY ANALYSIS
  # ==============================================================================
  path_step(4, "KEGG Pathway Analysis")
  
  kegg <- enrichKEGG(gene = gene_entrez$ENTREZID, universe = universe_entrez$ENTREZID, 
                     organism = "hsa", pAdjustMethod = "BH", 
                     pvalueCutoff = 0.05, qvalueCutoff = 0.1)
  path_msg(sprintf("KEGG pathways found: %d", if(!is.null(kegg)) nrow(kegg) else 0))
  
  qc_metrics <- add_pathway_qc(qc_metrics, "enrichment", "kegg_pathways", if(!is.null(kegg)) nrow(kegg) else 0)
  
  # ==============================================================================
  # 5. DISPLAY RESULTS
  # ==============================================================================
  path_step(5, "Top Results")
  
  if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
    path_msg("Top GO Biological Processes:", level="INFO")
    bp_df <- as.data.frame(go_bp)
    for (i in 1:min(10, nrow(bp_df))) {
      path_msg(sprintf("%2d. %-40s (p.adj = %.2e)", i, substr(bp_df$Description[i], 1, 40), bp_df$p.adjust[i]), level="DETAILS")
    }
  } else {
    path_msg("No significant GO-BP terms", level="WARN")
  }
  
  if (!is.null(kegg) && nrow(kegg@result) > 0) {
    path_msg("Top KEGG Pathways:", level="INFO")
    kegg_df <- as.data.frame(kegg)
    for (i in 1:min(10, nrow(kegg_df))) {
      path_msg(sprintf("%2d. %-40s (p.adj = %.2e)", i, substr(kegg_df$Description[i], 1, 40), kegg_df$p.adjust[i]), level="DETAILS")
    }
  } else {
    path_msg("No significant KEGG pathways", level="WARN")
  }
  
  # ==============================================================================
  # 6. ECM VALIDATION
  # ==============================================================================
  path_step(6, "ECM / Tumor Microenvironment Validation")
  ecm_keywords <- c("extracellular", "collagen", "matrix", "fibrosis", "cirrhosis", "ECM", "integrin", "laminin", "basement")
  
  for (kw in ecm_keywords) {
    found <- FALSE
    if (!is.null(go_bp) && nrow(go_bp) > 0) {
      hits <- grep(kw, go_bp@result$Description, ignore.case = TRUE)
      if (length(hits) > 0) {
        path_msg(sprintf("'%s' - %d GO-BP terms", kw, length(hits)))
        found <- TRUE
      }
    }
    if (!is.null(kegg) && nrow(kegg) > 0) {
      hits <- grep(kw, kegg@result$Description, ignore.case = TRUE)
      if (length(hits) > 0) {
        path_msg(sprintf("'%s' - %d KEGG pathways", kw, length(hits)))
        found <- TRUE
      }
    }
  }
  
  # ==============================================================================
  # 7. QC & EXPORT
  # ==============================================================================
  generate_pathway_qc_plots(if(!is.null(go_bp)) go_bp else data.frame(),
                            if(!is.null(go_mf)) go_mf else data.frame(),
                            if(!is.null(go_cc)) go_cc else data.frame(),
                            if(!is.null(kegg)) kegg else data.frame(),
                            opt$outdir)
                            
  export_pathway_results(if(!is.null(go_bp)) go_bp else data.frame(),
                         if(!is.null(go_mf)) go_mf else data.frame(),
                         if(!is.null(go_cc)) go_cc else data.frame(),
                         if(!is.null(kegg)) kegg else data.frame(),
                         opt$outdir)
                         
  export_pathway_qc(qc_metrics, opt$outdir)
  
  path_banner("PATHWAY ENRICHMENT COMPLETE \u2714")
  
}, error = function(e) {
  path_msg(sprintf("FATAL ERROR: %s", e$message), level = "ERROR")
  quit(status = 1)
})
