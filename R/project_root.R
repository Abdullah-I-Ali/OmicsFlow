# ==============================================================================
# R/project_root.R — OmicsFlow Project Root Detection and Resolution
# ==============================================================================

# Internal environment to cache the resolved project root
.omicsflow_env <- new.env(parent = emptyenv())

#' Set the OmicsFlow project root directory
#'
#' Explicitly configure the path to the OmicsFlow repository root where
#' scientific modules (\code{modules/}), reports (\code{reports/}), and
#' configuration files reside. This is useful when calling OmicsFlow functions
#' from a working directory other than the repository root.
#'
#' @param path Absolute or relative path to the OmicsFlow repository root.
#'   Must contain \code{modules/} and \code{reports/} subdirectories.
#'
#' @return The normalized path (invisibly).
#'
#' @examples
#' \dontrun{
#' set_omicsflow_root("~/projects/OmicsFlow")
#' }
#'
#' @export
set_omicsflow_root <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  if (!dir.exists(path)) {
    stop("OmicsFlow root directory does not exist: ", path)
  }
  # Validate that this looks like an OmicsFlow project
  modules_dir <- file.path(path, "modules")
  reports_dir <- file.path(path, "reports")
  if (!dir.exists(modules_dir)) {
    stop("Invalid OmicsFlow root: 'modules/' directory not found in: ", path)
  }
  if (!dir.exists(reports_dir)) {
    stop("Invalid OmicsFlow root: 'reports/' directory not found in: ", path)
  }
  .omicsflow_env$root <- path
  msg_ok(sprintf("OmicsFlow project root set to: %s", path))
  invisible(path)
}

#' Resolve the OmicsFlow project root directory
#'
#' Internal function that returns the project root, using (in order):
#' \enumerate{
#'   \item An explicitly set root via \code{\link{set_omicsflow_root}}
#'   \item Auto-detection by walking upward from the current working directory
#'   \item Fallback to the current working directory (with a warning)
#' }
#'
#' @return Absolute path to the OmicsFlow project root.
#'
#' @keywords internal
omicsflow_project_root <- function() {
  # 1. Return cached root if explicitly set

  if (!is.null(.omicsflow_env$root)) {
    return(.omicsflow_env$root)
  }

  # 2. Auto-detect by walking upward from getwd()
  candidate <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:10) {
    if (.is_omicsflow_root(candidate)) {
      .omicsflow_env$root <- candidate
      return(candidate)
    }
    parent <- dirname(candidate)
    if (parent == candidate) break # filesystem root reached
    candidate <- parent
  }

  # 3. Fallback to getwd() with warning
  msg_warn(paste(
    "Could not auto-detect OmicsFlow project root.",
    "Using current working directory:", getwd(),
    "\nHint: call set_omicsflow_root('/path/to/OmicsFlow') to set explicitly."
  ))
  return(normalizePath(getwd(), mustWork = FALSE))
}

#' Check whether a directory looks like the OmicsFlow project root
#'
#' @param path Directory path to check
#' @return Logical
#' @keywords internal
.is_omicsflow_root <- function(path) {
  dir.exists(file.path(path, "modules")) &&
    dir.exists(file.path(path, "reports")) &&
    (file.exists(file.path(path, "DESCRIPTION")) ||
     file.exists(file.path(path, "OmicsFlow.Rproj")))
}

#' Resolve a path relative to the OmicsFlow project root
#'
#' @param ... Path components to join after the project root
#' @return Absolute path
#' @keywords internal
omicsflow_path <- function(...) {
  file.path(omicsflow_project_root(), ...)
}
