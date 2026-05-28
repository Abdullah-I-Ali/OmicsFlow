# ==============================================================================
# R/install_dependencies.R — Dedicated Installer for Scientific Dependencies
# ==============================================================================

#' Install required scientific dependencies for OmicsFlow
#'
#' Installs the heavy scientific frameworks (CRAN and Bioconductor packages)
#' required by OmicsFlow's analytical modules. This keeps the core package
#' installation lightweight while providing a convenient helper to configure
#' the full computational environment.
#'
#' @param force Logical, whether to force re-installation of already installed
#'   packages. Defaults to \code{FALSE}.
#'
#' @return Invisible \code{TRUE} if all installations complete successfully.
#'
#' @examples
#' \dontrun{
#' install_omicsflow_dependencies()
#' }
#'
#' @export
install_omicsflow_dependencies <- function(force = FALSE) {
  
  msg_info("Checking CRAN dependencies...")
  
  cran_pkgs <- c("survminer", "xgboost", "randomForestSRC", "optparse")
  cran_to_install <- if (force) cran_pkgs else cran_pkgs[!cran_pkgs %in% installed.packages()[, "Package"]]
  
  if (length(cran_to_install) > 0) {
    msg_info(sprintf("Installing %d CRAN package(s): %s", 
                     length(cran_to_install), paste(cran_to_install, collapse = ", ")))
    install.packages(cran_to_install, repos = "https://cloud.r-project.org", quiet = TRUE)
  } else {
    msg_ok("All CRAN dependencies are already installed.")
  }
  
  msg_info("Checking Bioconductor dependencies...")
  
  bioc_pkgs <- c("MOFA2", "clusterProfiler", "edgeR", "limma", "sva", "org.Hs.eg.db")
  bioc_to_install <- if (force) bioc_pkgs else bioc_pkgs[!bioc_pkgs %in% installed.packages()[, "Package"]]
  
  if (length(bioc_to_install) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      msg_warn("BiocManager not found. Installing BiocManager...")
      install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
    }
    
    msg_info(sprintf("Installing %d Bioconductor package(s): %s", 
                     length(bioc_to_install), paste(bioc_to_install, collapse = ", ")))
    
    # Suppress update prompts during automated install
    BiocManager::install(bioc_to_install, update = FALSE, ask = FALSE, quiet = TRUE)
  } else {
    msg_ok("All Bioconductor dependencies are already installed.")
  }
  
  # Final verification check
  all_pkgs <- c(cran_pkgs, bioc_pkgs)
  missing_pkgs <- all_pkgs[!all_pkgs %in% installed.packages()[, "Package"]]
  
  if (length(missing_pkgs) > 0) {
    msg_fail(sprintf("Failed to install the following packages: %s", paste(missing_pkgs, collapse = ", ")))
    stop("Dependency installation incomplete. Please check error logs.")
  } else {
    msg_ok("All OmicsFlow scientific dependencies successfully installed and verified!")
    return(invisible(TRUE))
  }
}
