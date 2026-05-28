# ==============================================================================
# R/detect_sample_ids.R — Extract Sample IDs from Omics Inputs
# ==============================================================================

#' Extract sample identifiers from heterogeneous omics input files
#'
#' Reads each supplied omics file and extracts sample identifiers using
#' format-appropriate strategies: column names for matrix files (RNA,
#' methylation), unique values from identifier columns for tabular files
#' (CNV, SNV). All parameters are optional; only supplied modalities are
#' processed.
#'
#' @param rna Path to RNA-seq matrix file (.rds or .csv), or NULL
#' @param meth Path to DNA methylation matrix file (.rds or .csv), or NULL
#' @param cnv Path to CNV data file (.rds, .csv, or .tsv), or NULL
#' @param snv Path to SNV data file (.rds, .csv, or .tsv), or NULL
#' @return Named list with elements \code{$rna}, \code{$meth}, \code{$cnv},
#'   \code{$snv} (character vectors or NULL)
#' @keywords internal
detect_sample_ids <- function(rna = NULL, meth = NULL, cnv = NULL, snv = NULL) {
  res <- list(rna = NULL, meth = NULL, cnv = NULL, snv = NULL)
  
  # 1. RNA (Matrix colnames)
  if (!is.null(rna) && file.exists(rna)) {
    mat_rna <- read_matrix_safe(rna)
    if (!is.null(mat_rna)) {
      res$rna <- colnames(mat_rna)
    }
  }
  
  # 2. Methylation (Matrix colnames)
  if (!is.null(meth) && file.exists(meth)) {
    mat_meth <- read_matrix_safe(meth)
    if (!is.null(mat_meth)) {
      res$meth <- colnames(mat_meth)
    }
  }
  
  # 3. CNV (Typically unique values from Sample column)
  if (!is.null(cnv) && file.exists(cnv)) {
    df_cnv <- read_table_safe(cnv)
    if (!is.null(df_cnv)) {
      # Try "Sample" or first column
      if ("Sample" %in% colnames(df_cnv)) {
        res$cnv <- unique(as.character(df_cnv$Sample))
      } else if (ncol(df_cnv) > 0) {
        res$cnv <- unique(as.character(df_cnv[[1]]))
      }
    }
  }
  
  # 4. SNV (Unique values from Tumor_Sample_Barcode column)
  if (!is.null(snv) && file.exists(snv)) {
    df_snv <- read_table_safe(snv)
    if (!is.null(df_snv)) {
      # Try "Tumor_Sample_Barcode" or first column
      if ("Tumor_Sample_Barcode" %in% colnames(df_snv)) {
        res$snv <- unique(as.character(df_snv$Tumor_Sample_Barcode))
      } else if (ncol(df_snv) > 0) {
        res$snv <- unique(as.character(df_snv[[1]]))
      }
    }
  }
  
  return(res)
}
