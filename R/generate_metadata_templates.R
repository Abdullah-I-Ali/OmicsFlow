# ==============================================================================
# R/generate_metadata_templates.R — Orchestrate Metadata Template Generation
# ==============================================================================

#' Generate OmicsFlow metadata and clinical mapping templates
#'
#' Scans the supplied omics data files, detects sample identifiers, infers
#' patient IDs, and generates standardized metadata templates that can be
#' reviewed and edited before pipeline execution. Accepts any subset of
#' omics modalities — no modality is mandatory.
#'
#' @param rna Path to RNA-seq matrix file (.rds, .csv, .tsv, .txt), or NULL
#' @param meth Path to DNA methylation matrix file (.rds, .csv, .tsv, .txt), or NULL
#' @param cnv Path to CNV data file (.rds, .csv, .tsv, .txt), or NULL
#' @param snv Path to SNV data file (.rds, .csv, .tsv, .txt), or NULL
#' @param clinical Path to clinical TSV/CSV file (optional). If NULL, an
#'   empty clinical template and default column mapping are generated.
#' @param output_dir Path to output templates directory (created if absent)
#'
#' @return An S3 object of class \code{omicsflow_templates} containing:
#'   \describe{
#'     \item{status}{Character, "success" or "failed"}
#'     \item{output_dir}{Absolute path to the output directory}
#'     \item{templates_generated}{List of paths to generated files (metadata, clinical, clinical_map, validation)}
#'     \item{warnings}{Character vector of warnings}
#'     \item{errors}{Character vector of errors}
#'     \item{runtime}{Numeric execution time in seconds}
#'   }
#'
#' @examples
#' \dontrun{
#' # Generate templates from RNA-seq data only
#' generate_metadata_templates(rna = "data/rna.rds")
#'
#' # Generate templates from RNA + methylation with clinical data
#' generate_metadata_templates(
#'   rna = "data/rna.rds",
#'   meth = "data/meth.rds",
#'   clinical = "data/clinical.tsv"
#' )
#' }
#'
#' @export
generate_metadata_templates <- function(rna = NULL, meth = NULL, cnv = NULL, snv = NULL,
                                        clinical = NULL, output_dir = "templates") {
  
  start_time <- Sys.time()
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  report_lines <- c()
  report_lines <- c(report_lines, "### 1. File Summary")
  
  warnings_list <- character(0)
  
  # Override msg_warn to also track warnings
  orig_msg_warn <- msg_warn
  msg_warn <- function(text) {
    warnings_list <<- c(warnings_list, text)
    orig_msg_warn(text)
  }
  
  # Log files received
  check_file_log <- function(name, path) {
    if (!is.null(path)) {
      if (file.exists(path)) {
        sprintf("  - %s: %s (Found)", name, path)
      } else {
        sprintf("  - %s: %s (Not found)", name, path)
      }
    } else {
      sprintf("  - %s: Not provided", name)
    }
  }
  
  report_lines <- c(report_lines, check_file_log("RNA-seq", rna))
  report_lines <- c(report_lines, check_file_log("Methylation", meth))
  report_lines <- c(report_lines, check_file_log("CNV", cnv))
  report_lines <- c(report_lines, check_file_log("SNV", snv))
  report_lines <- c(report_lines, check_file_log("Clinical", clinical))
  report_lines <- c(report_lines, "")
  
  # Step 1: Detect sample IDs
  sample_ids_list <- detect_sample_ids(rna = rna, meth = meth, cnv = cnv, snv = snv)
  
  report_lines <- c(report_lines, "### 2. Sample Identification")
  for (m in names(sample_ids_list)) {
    count <- if (is.null(sample_ids_list[[m]])) 0 else length(sample_ids_list[[m]])
    report_lines <- c(report_lines, sprintf("  - %s: %d sample(s) detected", toupper(m), count))
  }
  report_lines <- c(report_lines, "")
  
  # Union sample IDs
  all_samples <- unique(unlist(sample_ids_list))
  
  # Step 2: Infer patient IDs
  patient_ids_df <- infer_patient_ids(all_samples)
  
  report_lines <- c(report_lines, "### 3. Patient ID Inference Summary")
  if (nrow(patient_ids_df) > 0) {
    methods_freq <- table(patient_ids_df$method)
    for (method_name in names(methods_freq)) {
      report_lines <- c(report_lines, sprintf("  - Method '%s': %d sample(s) mapped", method_name, methods_freq[method_name]))
    }
  } else {
    report_lines <- c(report_lines, "  - No sample IDs to map.")
  }
  report_lines <- c(report_lines, "")
  
  # Step 3: Build sample metadata table
  metadata_df <- build_metadata_table(sample_ids_list, patient_ids_df)
  
  # Write sample metadata CSV
  metadata_path <- file.path(output_dir, "sample_metadata.csv")
  write.csv(metadata_df, file = metadata_path, row.names = FALSE, quote = TRUE)
  msg_ok(sprintf("Sample metadata template written to: %s", metadata_path))
  
  # Step 4: Clinical column mapping if clinical is provided
  clinical_map_path <- NULL
  clinical_template_path <- NULL
  
  report_lines <- c(report_lines, "### 4. Clinical Column Mapping")
  if (!is.null(clinical) && file.exists(clinical)) {
    det_clinical <- detect_clinical_columns(clinical)
    
    # Write JSON mapping file
    clinical_map_path <- file.path(output_dir, "clinical_map.json")
    jsonlite::write_json(det_clinical$mapping, path = clinical_map_path, auto_unbox = TRUE, pretty = TRUE)
    msg_ok(sprintf("Clinical column mapping written to: %s", clinical_map_path))
    
    # Write summary
    for (field in names(det_clinical$mapping)) {
      report_lines <- c(report_lines, sprintf("  - Field '%s' -> Column '%s' (Confidence: %s)",
                                              field, det_clinical$mapping[[field]], det_clinical$confidence[[field]]))
    }
    
    if (length(det_clinical$warnings) > 0) {
      report_lines <- c(report_lines, "")
      report_lines <- c(report_lines, "Warnings:")
      for (w in det_clinical$warnings) {
        report_lines <- c(report_lines, sprintf("  [WARNING] %s", w))
        msg_warn(w)
      }
    }
  } else {
    report_lines <- c(report_lines, "  - [WARNING] No clinical file provided. Generating empty template.")
    msg_warn("No clinical file provided. Generating empty clinical template.")
    
    clinical_template_path <- file.path(output_dir, "custom_clinical_template.tsv")
    df_empty_clin <- data.frame(
      patient_id = character(0),
      os_time = numeric(0),
      os_event = integer(0),
      age = numeric(0),
      gender = character(0),
      stringsAsFactors = FALSE
    )
    write.table(df_empty_clin, file = clinical_template_path, sep = "\t", row.names = FALSE, quote = FALSE)
    msg_ok(sprintf("Empty clinical template written to: %s", clinical_template_path))
    
    clinical_map_path <- file.path(output_dir, "clinical_map.json")
    default_mapping <- list(
      patient_id = "patient_id",
      os_time = "os_time",
      os_event = "os_event",
      age = "age",
      gender = "gender"
    )
    jsonlite::write_json(default_mapping, path = clinical_map_path, auto_unbox = TRUE, pretty = TRUE)
    msg_ok(sprintf("Default clinical mapping written to: %s", clinical_map_path))
  }
  report_lines <- c(report_lines, "")
  
  # Write validation report txt
  validation_path <- file.path(output_dir, "validation_report.txt")
  write_validation_report(report_lines, validation_path)
  
  res <- list(
    status = "success",
    output_dir = normalizePath(output_dir, mustWork = FALSE),
    templates_generated = list(
      metadata = metadata_path,
      clinical = clinical_template_path,
      clinical_map = clinical_map_path,
      validation = validation_path
    ),
    warnings = warnings_list,
    errors = character(0),
    runtime = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  )
  class(res) <- "omicsflow_templates"
  
  return(invisible(res))
}
