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
  col_names_lower <- tolower(col_names)
  
  mapping <- list()
  confidence <- list()
  
  # ---------------------------------------------------------------------------
  # 1. Dual-column Survival Time Detection (Pattern-based & Dataset-Agnostic)
  # ---------------------------------------------------------------------------
  # Cohorts like TCGA use two columns for survival time:
  # - A column indicating event/death time (e.g. days_to_death)
  # - A column indicating censored/follow-up time (e.g. days_to_last_follow_up)
  event_time_col <- NULL
  censor_time_col <- NULL
  
  # Event/death time pattern: (death, dead, deceased, event) AND (day, month, year, time, survival)
  for (i in seq_along(col_names)) {
    col <- col_names[i]
    cl <- col_names_lower[i]
    if ((grepl("death", cl) || grepl("dead", cl) || grepl("deceased", cl) || grepl("event", cl)) &&
        (grepl("day", cl) || grepl("month", cl) || grepl("year", cl) || grepl("time", cl) || grepl("survival", cl))) {
      event_time_col <- col
      break
    }
  }
  
  # Censored/follow-up time pattern: (follow, censor, alive) AND (day, month, year, time, survival)
  for (i in seq_along(col_names)) {
    col <- col_names[i]
    cl <- col_names_lower[i]
    if ((grepl("follow", cl) || grepl("censor", cl) || grepl("alive", cl)) &&
        (grepl("day", cl) || grepl("month", cl) || grepl("year", cl) || grepl("time", cl) || grepl("survival", cl))) {
      censor_time_col <- col
      break
    }
  }
  
  dual_detected <- FALSE
  if (!is.null(event_time_col) && !is.null(censor_time_col) && (event_time_col != censor_time_col)) {
    mapping[["os_time"]] <- paste0(event_time_col, ",", censor_time_col)
    confidence[["os_time"]] <- "exact"
    dual_detected <- TRUE
  }
  
  # ---------------------------------------------------------------------------
  # 2. Extensible Dictionary-based Matching
  # ---------------------------------------------------------------------------
  dict <- list(
    patient_id = list(
      exact = c("patient_id", "patientid", "patient", "id", "bcr_patient_barcode", "patient_barcode", "sample", "sample_id", "sampleid", "subject", "subject_id", "case_id", "barcode", "patientbarcode"),
      fuzzy = c("individual", "individual_id", "unique_id", "uid", "uuid", "sample_barcode", "usubjid", "submitter_id")
    ),
    os_time = list(
      exact = c("os_time", "survival_time", "overall_survival_time", "time", "followup_time", "survival_months", "survival_days", "duration", "followup_months", "survival"),
      fuzzy = c("days", "months", "years", "overall_survival", "os_months", "os_days", "os", "followup", "follow_up", "t_os", "time_to_event")
    ),
    os_event = list(
      exact = c("os_event", "event", "vital_status", "status", "censored", "deceased", "death_status", "survival_status", "os_status", "dead_or_alive"),
      fuzzy = c("vital", "censor", "censor_status", "is_dead", "death", "dead", "alive", "d_os")
    ),
    age = list(
      exact = c("age", "age_at_diagnosis", "age_years", "age_at_dx", "age_at_index"),
      fuzzy = c("years_old", "birth_year", "dob", "birth", "age_at_start")
    ),
    gender = list(
      exact = c("gender", "sex"),
      fuzzy = c("biological_sex", "male_female", "gender_code")
    )
  )
  
  for (field in names(dict)) {
    if (field == "os_time" && dual_detected) {
      next
    }
    
    # 1. Exact case-insensitive matches in exact candidates
    exact_candidates <- tolower(dict[[field]]$exact)
    found_exact <- col_names[col_names_lower %in% exact_candidates]
    if (length(found_exact) > 0) {
      mapping[[field]] <- found_exact[1]
      confidence[[field]] <- "exact"
      next
    }
    
    # 2. Fuzzy/ci candidates
    fuzzy_candidates <- tolower(dict[[field]]$fuzzy)
    found_fuzzy <- col_names[col_names_lower %in% fuzzy_candidates]
    if (length(found_fuzzy) > 0) {
      mapping[[field]] <- found_fuzzy[1]
      confidence[[field]] <- "fuzzy"
      next
    }
    
    # 3. Partial substring/pattern matches
    pattern_field <- gsub("_", "", field)
    partial_matches <- col_names[grepl(pattern_field, col_names_lower) | grepl(field, col_names_lower)]
    if (length(partial_matches) > 0) {
      mapping[[field]] <- partial_matches[1]
      confidence[[field]] <- "fuzzy"
      next
    }
    
    # 4. Check if candidates are partial substrings of column names
    candidate_matched <- FALSE
    for (candidate in c(dict[[field]]$exact, dict[[field]]$fuzzy)) {
      cand_lower <- tolower(candidate)
      matches <- col_names[grepl(cand_lower, col_names_lower)]
      if (length(matches) > 0) {
        mapping[[field]] <- matches[1]
        confidence[[field]] <- "fuzzy"
        candidate_matched <- TRUE
        break
      }
    }
    if (candidate_matched) next
    
    # 5. Unresolved
    confidence[[field]] <- "unresolved"
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
