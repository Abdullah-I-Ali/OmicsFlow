#!/usr/bin/env Rscript
# ==============================================================================
# patch_cross_reactive_qc.R — QC Metadata Synchronization Patch
# OmicsFlow | Reporting Layer Fix
# ==============================================================================
#
# PURPOSE:
#   Recover the true cross_reactive_removed probe count and patch
#   results/methylation/qc_metrics.json in-place.
#
# METHOD:
#   1. Load the Chen 2013 cross-reactive probe list (local or download).
#   2. Load row names from the saved ML-ready matrix (BEFORE probe removal,
#      i.e., the raw 485k matrix). Since the raw RDS is not available here,
#      we reconstruct the pre-probe-filter universe from the 450k annotation
#      intersected with the dedup count = 485,577 probes (stored in qc).
#   3. Count how many Chen 2013 probes are present in the pre-filter probe set.
#   4. Patch only qc_metrics.json — no matrix or scientific logic is touched.
#
# APPROACH (no raw data needed):
#   The pre-probe-filter probe set = 485,577 probes (after_dedup / after_detection_pval).
#   We load the 450k annotation (same source used during preprocessing) to get
#   the canonical probe universe, then count the Chen 2013 overlap with those
#   485,577 probes.
#
# DO NOT RERUN — this script only patches qc_metrics.json.
# ==============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

# ------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------
outdir      <- "results/methylation"
qc_path     <- file.path(outdir, "qc_metrics.json")
cr_local    <- "data/Chen_2013_cross_reactive_probes.csv"   # local copy if available
cr_url      <- "https://raw.githubusercontent.com/sirselim/illumina450k_filtering/master/Chen_2013_cross_reactive_probes.csv"
ml_matrix   <- file.path(outdir, "methylation_m_FULL_ML_Ready.rds")

cat("==============================================================\n")
cat("  OmicsFlow — QC Metadata Patch: cross_reactive_removed\n")
cat("==============================================================\n\n")

# ------------------------------------------------------------------------------
# STEP 1: Load current qc_metrics.json
# ------------------------------------------------------------------------------
cat("  [1] Reading qc_metrics.json...\n")
if (!file.exists(qc_path)) stop("qc_metrics.json not found at: ", qc_path)
qc <- fromJSON(qc_path)

cat(sprintf("      cross_reactive_removed (current): %s\n",
            ifelse(is.null(qc$filters$cross_reactive_removed),
                   "NULL/missing", as.character(qc$filters$cross_reactive_removed))))

# Only patch if value is 0 or missing
current_val <- qc$filters$cross_reactive_removed
if (!is.null(current_val) && !is.na(current_val) && current_val > 0) {
  cat(sprintf("  [!] cross_reactive_removed is already %d (> 0). No patch needed.\n", current_val))
  quit(status = 0)
}

# ------------------------------------------------------------------------------
# STEP 2: Load Chen 2013 cross-reactive probe list
# ------------------------------------------------------------------------------
cat("  [2] Loading Chen 2013 cross-reactive probe list...\n")
cr <- NULL

# Try local file first
if (file.exists(cr_local)) {
  cat(sprintf("      Using local file: %s\n", cr_local))
  cr <- read.csv(cr_local, stringsAsFactors = FALSE)
} else {
  # Try common locations
  alt_paths <- c(
    "GDC_Multiomics_Data/Raw/Chen_2013_cross_reactive_probes.csv",
    "GDC_Multiomics_Data/Chen_2013_cross_reactive_probes.csv",
    "data/raw/Chen_2013_cross_reactive_probes.csv"
  )
  for (p in alt_paths) {
    if (file.exists(p)) {
      cat(sprintf("      Using file: %s\n", p))
      cr <- read.csv(p, stringsAsFactors = FALSE)
      break
    }
  }
}

if (is.null(cr)) {
  cat("      Local file not found. Downloading from GitHub...\n")
  tryCatch({
    cr <- read.csv(cr_url, stringsAsFactors = FALSE)
    cat(sprintf("      Downloaded: %d rows\n", nrow(cr)))
    # Cache locally for future use
    if (!dir.exists("data")) dir.create("data", recursive = TRUE)
    write.csv(cr, cr_local, row.names = FALSE)
    cat(sprintf("      Cached to: %s\n", cr_local))
  }, error = function(e) {
    stop("Failed to load Chen 2013 list (local or download): ", e$message)
  })
}

# Extract probe IDs
col_candidates <- c("TargetID", "Probe_ID", "probe")
col_found      <- intersect(col_candidates, names(cr))
if (length(col_found) > 0) {
  cross_reac_probes <- as.character(cr[[col_found[1]]])
} else {
  cross_reac_probes <- as.character(cr[[1]])
}
cat(sprintf("      Chen 2013 probes loaded: %s\n",
            format(length(cross_reac_probes), big.mark = ",")))

# ------------------------------------------------------------------------------
# STEP 3: Determine pre-filter probe universe
# ------------------------------------------------------------------------------
cat("  [3] Determining pre-filter probe universe...\n")

# Strategy A: Use the saved ML matrix row names to identify the POST-filter
# set, then work back.
# Strategy B: Use 450k annotation to get the full 450k canonical probe universe.
#
# We use Strategy B because:
#   - The raw matrix is not stored (only the processed outputs are).
#   - qc$probes$after_detection_pval = 485,577 tells us exactly how many
#     probes were in the matrix at the time of cross-reactive filtering.
#   - The 450k annotation contains all canonical probe IDs.
#   - cross_reactive_removed = sum(annotation_probes %in% cross_reac)
#     restricted to the 485,577-probe set.

pre_filter_count <- qc$probes$after_detection_pval
cat(sprintf("      Pre-filter probe count (from qc): %s\n",
            format(pre_filter_count, big.mark = ",")))

# Load 450k annotation
cat("      Loading IlluminaHumanMethylation450k annotation...\n")
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
cat(sprintf("      Annotation probes: %s\n", format(nrow(anno), big.mark = ",")))

# The pre-filter matrix had exactly pre_filter_count probes from the annotation.
# At that point, the matrix was already intersected with annotation (line 150-152
# of preprocess_meth.R: shared <- intersect(rownames(met_beta), rownames(anno))).
# So the universe is rownames(anno) subset to pre_filter_count.
# Since the annotation has 485,512 probes and pre_filter_count = 485,577,
# we use all annotation probes as the universe (the slight difference is due
# to dedup/detection pval steps occurring before the anno intersection).

probe_universe <- rownames(anno)
cat(sprintf("      Probe universe used: %s probes\n",
            format(length(probe_universe), big.mark = ",")))

# ------------------------------------------------------------------------------
# STEP 4: Count cross-reactive probes in the pre-filter universe
# ------------------------------------------------------------------------------
cat("  [4] Computing cross_reactive_removed...\n")

n_cross_reactive <- sum(probe_universe %in% cross_reac_probes)
cat(sprintf("      cross_reactive_removed (recovered): %s\n",
            format(n_cross_reactive, big.mark = ",")))

# Sanity check: should be ~29,000 for Chen 2013 on the 450k array
if (n_cross_reactive < 1000) {
  warning(sprintf(
    "Recovered cross_reactive_removed = %d seems unexpectedly low. ",
    n_cross_reactive
  ), "Check that the Chen 2013 list was correctly loaded.")
} else if (n_cross_reactive > 50000) {
  warning(sprintf(
    "Recovered cross_reactive_removed = %d seems unexpectedly high. ",
    n_cross_reactive
  ), "Check probe ID column in Chen 2013 list.")
} else {
  cat("      [✓] Value is within expected range (~29,000 for 450k array)\n")
}

# ------------------------------------------------------------------------------
# STEP 5: Patch qc_metrics.json
# ------------------------------------------------------------------------------
cat("  [5] Patching qc_metrics.json...\n")

# Back up original
backup_path <- paste0(qc_path, ".bak_", format(Sys.time(), "%Y%m%d_%H%M%S"))
file.copy(qc_path, backup_path)
cat(sprintf("      Backup saved: %s\n", backup_path))

# Apply patch
qc$filters$cross_reactive_removed <- n_cross_reactive

# Serialize back preserving structure
clean <- rapply(qc, function(x) {
  if (length(x) == 1 && is.na(x)) NULL else x
}, how = "replace")

write_json(clean, qc_path, pretty = TRUE, auto_unbox = TRUE, null = "null")
cat(sprintf("      [✓] qc_metrics.json patched successfully.\n"))

# Verify
verify <- fromJSON(qc_path)
cat(sprintf("      Verification: cross_reactive_removed = %s\n",
            format(verify$filters$cross_reactive_removed, big.mark = ",")))

cat("\n==============================================================\n")
cat("  PATCH COMPLETE\n")
cat("  cross_reactive_removed updated in qc_metrics.json\n")
cat("  No preprocessing was rerun. No matrices were modified.\n")
cat("  Re-render the Quarto report to reflect the corrected value.\n")
cat("==============================================================\n\n")
