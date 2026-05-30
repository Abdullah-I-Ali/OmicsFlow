# ==============================================================================
# qc_enrichment.R — Quality Control & Visualization for Pathway Enrichment
# OmicsFlow — Enrichment Module
# ==============================================================================

source("modules/enrichment/utils_enrichment.R")

generate_pathway_qc_plots <- function(go_bp, go_mf, go_cc, kegg, outdir) {
  path_step("QC", "Generating Enrichment Dotplots & Networks")
  
  plot_dir <- file.path(outdir, "plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  # 1. GO-BP Dotplot
  if (nrow(go_bp) > 0) {
    p_bp <- enrichplot::dotplot(go_bp, showCategory = 10) +
      ggtitle("GO Biological Process - Prognostic Genes") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 12))
    
    ggsave(file.path(plot_dir, "pathway_go_bp.png"), p_bp, width = 10, height = 7, dpi = 300)
    path_msg("Saved: plots/pathway_go_bp.png")
    
    # 2. Gene-Concept Network
    tryCatch({
      p_cnet <- cnetplot(go_bp, showCategory = 5, colorEdge = TRUE, node_label = "all") +
        ggtitle("Gene-Concept Network (Top 5 GO-BP Terms)") +
        theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
      ggsave(file.path(plot_dir, "pathway_cnetplot.png"), p_cnet, width = 12, height = 10, dpi = 300)
      path_msg("Saved: plots/pathway_cnetplot.png")
    }, error = function(e) path_msg("Could not generate cnetplot", level="WARN"))
    
    # 3. Enrichment Map
    tryCatch({
      go_bp_sim <- pairwise_termsim(go_bp)
      p_emap <- emapplot(go_bp_sim, showCategory = 15, color = "p.adjust", layout = "nicely") +
        ggtitle("Enrichment Map of GO-BP Pathways") +
        theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
      ggsave(file.path(plot_dir, "pathway_emapplot.png"), p_emap, width = 10, height = 8, dpi = 300)
      path_msg("Saved: plots/pathway_emapplot.png")
    }, error = function(e) path_msg("Could not generate emapplot", level="WARN"))
  }
  
  # 4. GO-MF Dotplot
  if (nrow(go_mf) > 0) {
    p_mf <- enrichplot::dotplot(go_mf, showCategory = 10) +
      ggtitle("GO Molecular Function - Prognostic Genes") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 12))
    ggsave(file.path(plot_dir, "pathway_go_mf.png"), p_mf, width = 10, height = 7, dpi = 300)
    path_msg("Saved: plots/pathway_go_mf.png")
  }
  
  # 5. GO-CC Dotplot
  if (nrow(go_cc) > 0) {
    p_cc <- enrichplot::dotplot(go_cc, showCategory = 10) +
      ggtitle("GO Cellular Component - Prognostic Genes") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 12))
    ggsave(file.path(plot_dir, "pathway_go_cc.png"), p_cc, width = 10, height = 7, dpi = 300)
    path_msg("Saved: plots/pathway_go_cc.png")
  }
  
  # 6. KEGG Dotplot
  if (nrow(kegg) > 0) {
    p_kegg <- enrichplot::dotplot(kegg, showCategory = 10) +
      ggtitle("KEGG Pathways - Prognostic Genes") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 12))
    
    ggsave(file.path(plot_dir, "pathway_kegg.png"), p_kegg, width = 10, height = 7, dpi = 300)
    path_msg("Saved: plots/pathway_kegg.png")
  }
}
