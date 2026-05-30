# ==============================================================================
# R/utils_io.R — Shared I/O Utilities for OmicsFlow Usability Layer
# ==============================================================================

#' Safely read a matrix from .rds or .csv
#'
#' @param path Path to an .rds or .csv file
#' @return Matrix or data.frame, or NULL on failure
#' @keywords internal
read_matrix_safe <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    msg_warn(sprintf("File not found: %s", path))
    return(NULL)
  }
  ext <- tolower(tools::file_ext(path))
  tryCatch({
    if (ext == "rds") {
      readRDS(path)
    } else if (ext == "csv") {
      read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    } else if (ext %in% c("tsv", "txt")) {
      read.delim(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    } else {
      msg_warn(sprintf("Unsupported input format: .%s\nSupported formats:\n- .rds\n- .csv\n- .tsv\n- .txt\n\nRecommended format:\n- .rds", ext))
      NULL
    }
  }, error = function(e) {
    msg_warn(sprintf("Failed to read %s: %s", basename(path), e$message))
    NULL
  })
}

#' Safely read a tabular file (.csv, .tsv, or .rds)
#'
#' @param path Path to a .csv, .tsv, or .rds file
#' @return Data frame, or NULL on failure
#' @keywords internal
read_table_safe <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    msg_warn(sprintf("File not found: %s", path))
    return(NULL)
  }
  ext <- tolower(tools::file_ext(path))
  tryCatch({
    if (ext == "rds") {
      readRDS(path)
    } else if (ext == "csv") {
      read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (ext %in% c("tsv", "txt")) {
      read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      msg_warn(sprintf("Unsupported input format: .%s\nSupported formats:\n- .rds\n- .csv\n- .tsv\n- .txt\n\nRecommended format:\n- .rds", ext))
      NULL
    }
  }, error = function(e) {
    msg_warn(sprintf("Failed to read %s: %s", basename(path), e$message))
    NULL
  })
}

# --- Formatted console output helpers ---

#' Print a success message to the console
#' @param text Message text
#' @keywords internal
msg_ok <- function(text) {
  cat(sprintf("[OK]      %s\n", text))
}

#' Print a warning message to the console
#' @param text Message text
#' @keywords internal
msg_warn <- function(text) {
  cat(sprintf("[WARNING] %s\n", text))
}

#' Print a failure message to the console
#' @param text Message text
#' @keywords internal
msg_fail <- function(text) {
  cat(sprintf("[FAIL]    %s\n", text))
}

#' Print an info message to the console
#' @param text Message text
#' @keywords internal
msg_info <- function(text) {
  cat(sprintf("[INFO]    %s\n", text))
}

#' Write a vector of formatted status lines to a text file
#'
#' @param lines Character vector of report lines
#' @param path Output file path
#' @keywords internal
write_validation_report <- function(lines, path) {
  header <- c(
    "==============================================================================",
    sprintf("  OmicsFlow — Input Validation Report"),
    sprintf("  Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "==============================================================================",
    ""
  )
  writeLines(c(header, lines), con = path)
  msg_ok(sprintf("Validation report written to: %s", path))
}
