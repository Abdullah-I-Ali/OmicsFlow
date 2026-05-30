# ==============================================================================
# R/utils_cnv.R — CNV Utility Functions
# ==============================================================================

#' Generate a CNV Gene Coordinate Cache from Ensembl
#'
#' Downloads the latest hg38 gene coordinates (start and end positions) for all 
#' human genes using the \code{biomaRt} package. This cache is required by the CNV 
#' preprocessing module to map genomic segments to gene-level features.
#' 
#' Note: Queries to Ensembl may take 1-3 minutes and can occasionally fail due
#' to network timeouts or server unresponsiveness.
#'
#' @param destfile Character. Path where the resulting \code{.rds} cache file will be saved.
#' @param verbose Logical. If TRUE, prints progress messages.
#'
#' @return The normalized path to the created cache file (invisibly).
#' @export
generate_cnv_cache <- function(destfile = "gene_coords_hg38.rds", verbose = TRUE) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' is required to generate the CNV cache. Please install it or run install_omicsflow_dependencies().")
  }
  
  if (verbose) {
    message("Downloading gene coordinates from Ensembl (this may take 1-3 minutes)...")
  }
  
  tryCatch({
    mart <- biomaRt::useEnsembl(
      biomart = "genes",
      dataset = "hsapiens_gene_ensembl",
      mirror  = "useast" # "useast" mirror is often more reliable
    )
    
    gene_coords <- biomaRt::getBM(
      attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
      filters    = "chromosome_name",
      values     = c(1:22, "X", "Y"),
      mart       = mart
    )
    
    # Filter out empty symbols
    gene_coords <- gene_coords[gene_coords$hgnc_symbol != "", ]
    
    saveRDS(gene_coords, destfile)
    
    if (verbose) {
      message(sprintf("Successfully saved cache with %s genes to: %s", 
                      format(nrow(gene_coords), big.mark = ","), destfile))
    }
    
    return(invisible(normalizePath(destfile)))
    
  }, error = function(e) {
    stop("Failed to download gene coordinates from Ensembl: ", e$message, 
         "\n\nEnsembl servers are frequently unresponsive. Please try again later, or use the bundled default cache.")
  })
}
