# ==============================================================================
# R/project_root.R — OmicsFlow Resource Resolution
#
# Provides path resolution for OmicsFlow resources (modules, reports, configs).
# Resolution order:
#   1. Explicit override via set_omicsflow_root()
#   2. Installed package resources via system.file()
#   3. Development: walk up from getwd() to find repository root
# ==============================================================================

# Internal environment to cache the resolved project root
.omicsflow_env <- new.env(parent = emptyenv())

#' Set the OmicsFlow project root directory
#'
#' Explicitly configure the path to the OmicsFlow resource root where
#' scientific modules (\code{modules/}), reports (\code{reports/}), and
#' configuration files reside. This is primarily useful for developers
#' working directly in the OmicsFlow repository or users who need to
#' override the default resolution.
#'
#' For most users who install OmicsFlow via \code{install_github()}, this
#' function is \strong{not needed} — resource paths are resolved automatically
#' from the installed package.
#'
#' @param path Absolute or relative path to a directory containing a
#'   \code{modules/} subdirectory.
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
  # Validate that this looks like an OmicsFlow resource root
  modules_dir <- file.path(path, "modules")
  if (!dir.exists(modules_dir)) {
    stop("Invalid OmicsFlow root: 'modules/' directory not found in: ", path)
  }
  .omicsflow_env$root <- path
  msg_ok(sprintf("OmicsFlow project root set to: %s", path))
  invisible(path)
}

#' Resolve the OmicsFlow resource root directory
#'
#' Internal function that returns the root directory containing OmicsFlow
#' resources (modules, reports, configs). Uses the following resolution
#' order:
#' \enumerate{
#'   \item An explicitly set root via \code{\link{set_omicsflow_root}}
#'   \item The installed package directory (via \code{system.file()})
#'   \item Auto-detection by walking upward from the current working directory
#' }
#'
#' After a standard \code{install_github()} installation, resolution (2) is
#' used and no user configuration is needed.
#'
#' @return Absolute path to the OmicsFlow resource root.
#'
#' @keywords internal
omicsflow_root <- function() {
  # 1. Return cached/explicit root if set
  if (!is.null(.omicsflow_env$root)) {
    return(.omicsflow_env$root)
  }

  # 2. Installed package path via system.file()
  pkg_modules <- system.file("modules", package = "OmicsFlow")
  if (nzchar(pkg_modules) && dir.exists(pkg_modules)) {
    pkg_root <- system.file(package = "OmicsFlow")
    .omicsflow_env$root <- pkg_root
    return(pkg_root)
  }

  # 3. Development: walk upward from getwd()
  candidate <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:10) {
    if (.is_omicsflow_root(candidate)) {
      .omicsflow_env$root <- candidate
      return(candidate)
    }
    parent <- dirname(candidate)
    if (parent == candidate) break  # filesystem root reached
    candidate <- parent
  }

  stop(
    "Cannot locate OmicsFlow resources (modules/, reports/).\n",
    "If you installed via install_github(), please reinstall the package.\n",
    "If working from a source checkout, call set_omicsflow_root('/path/to/OmicsFlow')."
  )
}

#' Resolve the OmicsFlow project root directory
#'
#' @return Absolute path to the OmicsFlow project root.
#' @keywords internal
omicsflow_project_root <- function() {
  omicsflow_root()
}

#' Check whether a directory looks like the OmicsFlow project root
#'
#' @param path Directory path to check
#' @return Logical
#' @keywords internal
.is_omicsflow_root <- function(path) {
  dir.exists(file.path(path, "modules")) &&
    (dir.exists(file.path(path, "reports")) ||
     file.exists(file.path(path, "DESCRIPTION")) ||
     file.exists(file.path(path, "OmicsFlow.Rproj")))
}

#' Resolve a path relative to the OmicsFlow resource root
#'
#' @param ... Path components to join after the root
#' @return Absolute path
#' @keywords internal
omicsflow_path <- function(...) {
  file.path(omicsflow_root(), ...)
}
