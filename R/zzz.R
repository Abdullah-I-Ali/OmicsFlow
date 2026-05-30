# ==============================================================================
# R/zzz.R — Package Load Hooks for OmicsFlow
# ==============================================================================

.onAttach <- function(libname, pkgname) {
  ver <- utils::packageVersion("OmicsFlow")
  packageStartupMessage(sprintf(
    "OmicsFlow v%s — Multi-Omics Integration Framework\nType ?omicsflow for help getting started.",
    ver
  ))
}

.onLoad <- function(libname, pkgname) {
  # Initialize the root cache as NULL (lazy detection via system.file())
  .omicsflow_env$root <- NULL
}
