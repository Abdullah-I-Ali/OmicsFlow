# ==============================================================================
# qc_snv.R — Quality Control & Visualization for SNV
# OmicsFlow — SNV Module
# ==============================================================================

source("modules/snv/utils_snv.R")

generate_snv_qc_plots <- function(snv_matrix, maf_filtered, outdir) {
  snv_step("QC", "Generating QC Visualizations")
  
  plot_dir <- file.path(outdir, "plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  # 1. Oncoplot (Top 20 mutated genes)
  png(file.path(plot_dir, "SNV_Oncoplot.png"), width = 1000, height = 800)
  oncoplot(maf = maf_filtered, top = 20)
  dev.off()
  snv_msg("Saved: plots/SNV_Oncoplot.png")
  
  # 2. Top 20 Mutated Genes Frequency Plot
  mut_freq <- rowSums(snv_matrix) / ncol(snv_matrix) * 100
  top_genes <- sort(mut_freq, decreasing = TRUE)[1:min(20, length(mut_freq))]
  
  df_top_genes <- data.frame(
    Gene = factor(names(top_genes), levels = rev(names(top_genes))), 
    Frequency = top_genes
  )
  
  fig_mut_freq <- ggplot(df_top_genes, aes(x = Gene, y = Frequency)) +
    geom_bar(stat = "identity", fill = "#c0392b", alpha = 0.8) +
    coord_flip() +
    theme_bw() +
    labs(title = "Top Mutated Genes (Mutation Frequency %)",
         subtitle = paste0("Cohort Size: ", ncol(snv_matrix), " Patients"),
         x = "Gene Symbol", y = "Percentage of Patients Mutated (%)")
  
  # 3. Tumor Mutation Burden (TMB) Distribution
  mut_per_patient <- colSums(snv_matrix)
  df_tmb <- data.frame(Patient = names(mut_per_patient), Mutations = mut_per_patient)
  
  fig_tmb <- ggplot(df_tmb, aes(x = Mutations)) +
    geom_histogram(binwidth = 1, fill = "#2980b9", color = "white", alpha = 0.9) +
    theme_classic() +
    labs(title = "Tumor Mutation Burden (TMB) Distribution",
         subtitle = "Hypermutated patients (>99th percentile) excluded",
         x = paste0("Number of Functional Mutations (in selected ", nrow(snv_matrix), " genes)"), 
         y = "Number of Patients")
  
  # 4. Gene Co-occurrence Heatmap (Top 15 Genes)
  top_n_cooccur <- min(15, nrow(snv_matrix))
  if (top_n_cooccur > 1) {
    top15 <- names(sort(rowSums(snv_matrix), decreasing = TRUE)[1:top_n_cooccur])
    snv_top15 <- snv_matrix[top15, ]
    
    cor_mat <- cor(t(snv_top15), method = "pearson")
    cor_melt <- reshape2::melt(cor_mat)
    
    fig_co_occur <- ggplot(cor_melt, aes(Var1, Var2, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 0, limit = c(-0.5, 0.5), space = "Lab", 
                           name="Correlation\n(Pearson)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      labs(title = sprintf("Co-mutation Matrix (Top %d Genes)", top_n_cooccur),
           subtitle = "Red = Co-occurrence | Blue = Mutual Exclusivity",
           x = "", y = "")
           
    final_snv_plot <- grid.arrange(fig_mut_freq, fig_tmb, fig_co_occur, 
                                   layout_matrix = rbind(c(1, 2), c(3, 3)), 
                                   heights = c(1, 1.2))
  } else {
    final_snv_plot <- grid.arrange(fig_mut_freq, fig_tmb, ncol = 1)
  }
  
  # Save statistical validation figures
  ggsave(file.path(plot_dir, "SNV_Research_Validation_Figures.png"), 
         plot = final_snv_plot, width = 11, height = 10, dpi = 300)
  snv_msg("Saved: plots/SNV_Research_Validation_Figures.png")
}
