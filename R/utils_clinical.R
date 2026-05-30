# ==============================================================================
# R/utils_clinical.R — Clinical Data Utilities (Package-Internal)
#
# These functions were originally in modules/utils_clinical.R and loaded via
# source() at runtime.  They are now proper package-internal functions so that
# the R API layer can call them without requiring the repository source tree.
# ==============================================================================

#' Load and parse clinical mapping configuration
#'
#' @param mapping_input Can be a JSON file path, a JSON string, or a named list
#' @return A named list of column mappings
#' @keywords internal
parse_clinical_mapping <- function(mapping_input) {
  if (is.null(mapping_input) || mapping_input == "") {
    return(NULL)
  }

  if (is.list(mapping_input)) {
    return(mapping_input)
  }

  # If it is a string
  if (is.character(mapping_input)) {
    # Check if it points to an existing file
    if (file.exists(mapping_input)) {
      tryCatch({
        return(jsonlite::fromJSON(mapping_input))
      }, error = function(e) {
        # If JSON failed, try parsing it as key-value pairs (e.g. patient_id:ID,os_time:time)
        lines <- readLines(mapping_input, warn = FALSE)
        parsed <- list()
        for (line in lines) {
          line <- trimws(line)
          if (line == "" || grepl("^#", line)) next
          parts <- strsplit(line, "[:=]")[[1]]
          if (length(parts) >= 2) {
            key <- trimws(parts[1])
            val <- trimws(paste(parts[-1], collapse = ":"))
            parsed[[key]] <- val
          }
        }
        if (length(parsed) > 0) return(parsed)
        stop(sprintf("Failed to parse clinical mapping file %s: %s", mapping_input, e$message))
      })
    } else {
      # Try parsing direct JSON string
      tryCatch({
        return(jsonlite::fromJSON(mapping_input))
      }, error = function(e) {
        # Try parsing comma-separated key=value or key:value
        parsed <- list()
        pairs <- strsplit(mapping_input, ",")[[1]]
        for (pair in pairs) {
          parts <- strsplit(pair, "[:=]")[[1]]
          if (length(parts) >= 2) {
            key <- trimws(parts[1])
            val <- trimws(paste(parts[-1], collapse = ":"))
            parsed[[key]] <- val
          }
        }
        if (length(parsed) > 0) return(parsed)
        stop(sprintf("Invalid clinical mapping format. Must be JSON file path, JSON string, or key=value list: %s", mapping_input))
      })
    }
  }

  stop("Unknown clinical mapping format")
}

#' Load, standardize and validate clinical data
#'
#' @param file Path to clinical data (TSV/CSV)
#' @param column_map Optional clinical mapping (file path, JSON string, or list)
#' @param metadata Optional metadata dataframe (Phase 1 sample metadata) to align patients
#' @return Standardized clinical data frame with columns: patient_id, os_time, os_event, age, gender
#' @keywords internal
load_clinical_data <- function(file, column_map = NULL, metadata = NULL) {
  suppressPackageStartupMessages(library(dplyr))
  if (is.null(file) || file == "") {
    stop("Clinical data file path must be provided")
  }
  if (!file.exists(file)) {
    stop("Clinical data file not found: ", file)
  }

  # Read clinical data
  # Auto-detect separator
  ext <- tolower(tools::file_ext(file))
  if (ext == "csv") {
    clinical_raw <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    clinical_raw <- read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
  }

  # Parse mapping configuration
  mapping <- parse_clinical_mapping(column_map)

  # Determine if we fall back to TCGA auto-detection
  is_tcga <- FALSE
  if (is.null(mapping)) {
    if ("bcr_patient_barcode" %in% colnames(clinical_raw)) {
      is_tcga <- TRUE
    } else {
      # Custom non-TCGA default mapping: look for standard columns
      mapping <- list(
        patient_id = intersect(colnames(clinical_raw), c("patient_id", "PatientID", "patient", "ID", "SampleID", "sample_id"))[1],
        os_time    = intersect(colnames(clinical_raw), c("os_time", "survival_time", "time", "days_to_death", "days"))[1],
        os_event   = intersect(colnames(clinical_raw), c("os_event", "event", "vital_status", "status"))[1],
        age        = intersect(colnames(clinical_raw), c("age", "age_at_diagnosis", "Age"))[1],
        gender     = intersect(colnames(clinical_raw), c("gender", "sex", "Gender", "Sex"))[1]
      )
      # Clean up NULL mappings
      mapping <- mapping[!sapply(mapping, is.null)]
      mapping <- mapping[!is.na(mapping)]
    }
  }

  # Standardize dataframe
  standardized <- data.frame(row.names = NULL)

  if (is_tcga) {
    # Extract using hardcoded TCGA logic
    required_tcga <- c("bcr_patient_barcode", "vital_status")
    missing_tcga <- setdiff(required_tcga, colnames(clinical_raw))
    if (length(missing_tcga) > 0) {
      stop("TCGA format auto-detected but missing expected columns: ", paste(missing_tcga, collapse = ", "))
    }

    # TCGA Event parsing
    vital <- clinical_raw$vital_status
    os_event <- ifelse(vital == "Dead", 1, ifelse(vital == "Alive" | vital == "Censored", 0, NA))
    # If they are NA because they don't match "Dead"/"Alive", try case-insensitive or default to "Dead" = 1, others = 0
    os_event[is.na(os_event)] <- ifelse(tolower(vital[is.na(os_event)]) == "dead", 1, 0)

    # TCGA Time parsing (days_to_death for dead, days_to_last_follow_up for alive)
    days_death <- if ("days_to_death" %in% colnames(clinical_raw)) as.numeric(clinical_raw$days_to_death) else rep(NA_real_, nrow(clinical_raw))
    days_follow <- if ("days_to_last_follow_up" %in% colnames(clinical_raw)) as.numeric(clinical_raw$days_to_last_follow_up) else rep(NA_real_, nrow(clinical_raw))

    os_time <- ifelse(vital == "Dead", days_death, days_follow)

    # Fallback if both might have values or one is missing
    os_time[is.na(os_time) & vital == "Dead"] <- days_death[is.na(os_time) & vital == "Dead"]
    os_time[is.na(os_time)] <- days_follow[is.na(os_time)]

    # Age (TCGA is in days)
    age_days <- if ("age_at_diagnosis" %in% colnames(clinical_raw)) as.numeric(clinical_raw$age_at_diagnosis) else NA_real_
    age_years <- age_days / 365.25

    # Gender
    gender_val <- if ("gender" %in% colnames(clinical_raw)) clinical_raw$gender else rep(NA_character_, nrow(clinical_raw))

    standardized <- data.frame(
      patient_id = clinical_raw$bcr_patient_barcode,
      os_time    = os_time,
      os_event   = os_event,
      age        = age_years,
      gender     = gender_val,
      stringsAsFactors = FALSE
    )
  } else {
    # Extract using the parsed custom mapping
    patient_col <- mapping$patient_id
    time_col    <- mapping$os_time
    event_col   <- mapping$os_event
    age_col     <- mapping$age
    gender_col  <- mapping$gender

    if (is.null(patient_col) || !patient_col %in% colnames(clinical_raw)) {
      stop("patient_id column mapping not found in clinical columns: ", paste(colnames(clinical_raw), collapse = ", "))
    }

    # patient_id
    patient_ids <- as.character(clinical_raw[[patient_col]])

    # os_time
    os_time <- rep(NA_real_, nrow(clinical_raw))
    if (!is.null(time_col)) {
      if (length(time_col) == 1 && grepl(",", time_col)) {
        time_cols <- trimws(strsplit(time_col, ",")[[1]])
      } else {
        time_cols <- time_col
      }

      if (length(time_cols) == 2) {
        # Dual-column format specified manually (e.g. days_to_death, days_to_last_follow_up)
        # We need an event column to decide
        t1 <- as.numeric(clinical_raw[[time_cols[1]]])
        t2 <- as.numeric(clinical_raw[[time_cols[2]]])

        # Check event_status
        raw_events <- if (!is.null(event_col) && event_col %in% colnames(clinical_raw)) clinical_raw[[event_col]] else rep(NA, nrow(clinical_raw))
        # If event is dead/1, use t1, else t2
        is_dead <- rep(FALSE, length(raw_events))
        if (!is.null(event_col)) {
          is_dead <- tolower(as.character(raw_events)) %in% c("dead", "1", "true", "deceased")
        }
        os_time <- ifelse(is_dead, t1, t2)
        os_time[is.na(os_time)] <- ifelse(!is.na(t1[is.na(os_time)]), t1[is.na(os_time)], t2[is.na(os_time)])
      } else if (length(time_cols) == 1 && time_cols %in% colnames(clinical_raw)) {
        os_time <- as.numeric(clinical_raw[[time_cols]])
      }
    }

    # os_event
    os_event <- rep(NA_integer_, nrow(clinical_raw))
    if (!is.null(event_col) && event_col %in% colnames(clinical_raw)) {
      raw_event <- clinical_raw[[event_col]]
      # standard conversions
      char_event <- tolower(trimws(as.character(raw_event)))
      os_event <- ifelse(char_event %in% c("1", "true", "dead", "deceased", "yes"), 1,
                         ifelse(char_event %in% c("0", "false", "alive", "censored", "no"), 0, NA_integer_))
      # For anything remaining NA, if it is numeric and > 0, treat as 1
      is_num <- !is.na(suppressWarnings(as.numeric(raw_event)))
      if (any(is_num)) {
        num_vals <- as.numeric(raw_event)
        os_event[is_num & is.na(os_event)] <- ifelse(num_vals[is_num & is.na(os_event)] > 0, 1, 0)
      }
      # Default any remaining NAs to 0
      os_event[is.na(os_event)] <- 0
    }

    # age
    age_years <- rep(NA_real_, nrow(clinical_raw))
    if (!is.null(age_col) && age_col %in% colnames(clinical_raw)) {
      raw_age <- as.numeric(clinical_raw[[age_col]])
      # Heuristic for age unit (days vs years) based on median
      valid_ages <- raw_age[!is.na(raw_age)]
      if (length(valid_ages) > 0 && median(valid_ages, na.rm = TRUE) > 365) {
        age_years <- raw_age / 365.25
      } else {
        age_years <- raw_age
      }
    }

    # gender
    gender_val <- rep(NA_character_, nrow(clinical_raw))
    if (!is.null(gender_col) && gender_col %in% colnames(clinical_raw)) {
      gender_val <- as.character(clinical_raw[[gender_col]])
    }

    standardized <- data.frame(
      patient_id = patient_ids,
      os_time    = os_time,
      os_event   = os_event,
      age        = age_years,
      gender     = gender_val,
      stringsAsFactors = FALSE
    )
  }

  # Standardize patient_id formatting if metadata is available
  if (!is.null(metadata)) {
    clinical_pids <- standardized$patient_id
    meta_pids <- unique(metadata$patient_id)
    meta_sids <- unique(metadata$sample_id)

    match_direct <- clinical_pids %in% meta_pids
    match_sample <- clinical_pids %in% meta_sids

    if (sum(match_direct) < sum(match_sample)) {
      mapped_idx <- match(clinical_pids, metadata$sample_id)
      standardized$patient_id[!is.na(mapped_idx)] <- metadata$patient_id[mapped_idx[!is.na(mapped_idx)]]
    }
  }

  # Remove duplicate rows for same patient_id, keeping the one with valid survival time if possible
  if (nrow(standardized) > 0) {
    standardized <- standardized %>%
      arrange(patient_id, desc(!is.na(os_time)), desc(os_time)) %>%
      distinct(patient_id, .keep_all = TRUE)
  }

  return(standardized)
}
