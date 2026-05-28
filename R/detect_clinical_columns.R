# ==============================================================================
# R/detect_clinical_columns.R — Dictionary-based Clinical Mapping Auto-detection
# ==============================================================================

#' Dictionary-based auto-detection of clinical column semantics
#'
#' Inspects column names in a clinical data file and attempts to map them
#' to standardized OmicsFlow fields (\code{patient_id}, \code{os_time},
#' \code{os_event}, \code{age}, \code{gender}) using exact, case-insensitive,
#' fuzzy, and partial substring matching strategies.
#'
#' @param clinical_path Path to clinical TSV/CSV file
#' @return A list with elements:
#'   \describe{
#'     \item{mapping}{Named list of detected column mappings}
#'     \item{confidence}{Named list of match confidence levels
#'       (\code{"exact"}, \code{"fuzzy"}, or \code{"unresolved"})}
#'     \item{warnings}{Character vector of warnings for unresolved required fields}
#'   }
#' @keywords internal
detect_clinical_columns <- function(clinical_path) {
  warnings_vec <- character(0)
  
  if (is.null(clinical_path) || !file.exists(clinical_path)) {
    return(list(
      mapping = list(),
      confidence = list(),
      warnings = c(sprintf("Clinical file not found: %s", clinical_path))
    ))
  }
  
  df_clinical <- read_table_safe(clinical_path)
  if (is.null(df_clinical) || ncol(df_clinical) == 0) {
    return(list(
      mapping = list(),
      confidence = list(),
      warnings = c(sprintf("Failed to read clinical file or file is empty: %s", clinical_path))
    ))
  }
  
  col_names <- colnames(df_clinical)
  
  # Define candidate dictionaries (exact/high priority first)
  dict <- list(
    patient_id = list(
      exact = c("patient_id", "PatientID", "patient", "ID", "bcr_patient_barcode"),
      fuzzy = c("case_id", "SampleID", "sample_id", "subject", "subject_id", "barcode", "patientbarcode")
    ),
    os_time = list(
      exact = c("os_time", "survival_time", "time", "overall_survival_time", "days_to_death"),
      fuzzy = c("days", "months", "years", "duration", "survival", "overall_survival", "os_months", "os_days")
    ),
    os_event = list(
      exact = c("os_event", "event", "vital_status", "status", "censored"),
      fuzzy = c("death", "deceased", "os_status", "survival_status", "dead_or_alive")
    ),
    age = list(
      exact = c("age", "age_at_diagnosis", "Age"),
      fuzzy = c("years_old", "age_years", "birth_year")
    ),
    gender = list(
      exact = c("gender", "sex", "Gender", "Sex"),
      fuzzy = c("biological_sex", "male_female")
    )
  )
  
  mapping <- list()
  confidence <- list()
  
  for (field in names(dict)) {
    # 1. Look for exact matches (case-sensitive first)
    exact_candidates <- dict[[field]]$exact
    found_exact <- col_names[col_names %in% exact_candidates]
    
    if (length(found_exact) > 0) {
      mapping[[field]] <- found_exact[1]
      confidence[[field]] <- "exact"
      next
    }
    
    # 2. Look for case-insensitive matches in exact candidates
    col_names_lower <- tolower(col_names)
    found_exact_ci <- col_names[col_names_lower %in% tolower(exact_candidates)]
    if (length(found_exact_ci) > 0) {
      mapping[[field]] <- found_exact_ci[1]
      confidence[[field]] <- "fuzzy"
      next
    }
    
    # 3. Look for fuzzy/ci candidates
    fuzzy_candidates <- dict[[field]]$fuzzy
    found_fuzzy <- col_names[col_names_lower %in% tolower(fuzzy_candidates)]
    if (length(found_fuzzy) > 0) {
      mapping[[field]] <- found_fuzzy[1]
      confidence[[field]] <- "fuzzy"
      next
    }
    
    # 4. Partial substring matches (e.g. if field is "age" and column is "age_at_dx")
    partial_matches <- col_names[grepl(field, col_names, ignore.case = TRUE)]
    if (length(partial_matches) > 0) {
      mapping[[field]] <- partial_matches[1]
      confidence[[field]] <- "fuzzy"
      next
    }
    
    # 5. Unresolved
    confidence[[field]] <- "unresolved"
    # For required columns in survival pipeline, warn if missing
    if (field %in% c("patient_id", "os_time", "os_event")) {
      warnings_vec <- c(warnings_vec, sprintf("Could not auto-detect column for required clinical field: %s", field))
    }
  }
  
  return(list(
    mapping = mapping,
    confidence = confidence,
    warnings = warnings_vec
  ))
}
