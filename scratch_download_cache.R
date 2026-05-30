suppressPackageStartupMessages(library(biomaRt))

cat("Downloading from Ensembl (1-3 minutes)...\n")
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror  = "useast"
)
gene_coords <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters    = "chromosome_name",
  values     = c(1:22, "X", "Y"),
  mart       = mart
)
gene_coords <- gene_coords[gene_coords$hgnc_symbol != "", ]
saveRDS(gene_coords, "gene_coords_hg38.rds")
cat(sprintf("Saved cache with %d genes. Size: %s\n", nrow(gene_coords), format(object.size(gene_coords), units="MB")))
