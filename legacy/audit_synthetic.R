suppressPackageStartupMessages(library(jsonlite))

cat("--- SYNTHETIC COHORT INTEGRITY AUDIT ---\n")
dir <- "data/synthetic_cohort"
errors <- character()

# Load files
rna <- readRDS(file.path(dir, "rna.rds"))
meth <- readRDS(file.path(dir, "meth.rds"))
cnv <- readRDS(file.path(dir, "cnv.rds"))
snv <- readRDS(file.path(dir, "snv.rds"))
meta <- read.csv(file.path(dir, "sample_metadata.csv"), stringsAsFactors=FALSE)
clin <- read.table(file.path(dir, "custom_clinical.tsv"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
cmap <- read_json(file.path(dir, "clinical_map.json"))

# 1. Modality sample counts
cat("1. Sample Counts & Dimensions\n")
cat("RNA:", nrow(rna), "x", ncol(rna), "\n")
cat("METH:", nrow(meth), "x", ncol(meth), "\n")
cat("CNV:", nrow(cnv), "rows (", length(unique(cnv$Sample)), "unique samples )\n")
cat("SNV:", nrow(snv), "rows (", length(unique(snv$Tumor_Sample_Barcode)), "unique samples )\n")
cat("Metadata:", nrow(meta), "x", ncol(meta), "\n")
cat("Clinical:", nrow(clin), "x", ncol(clin), "\n")

if (ncol(rna) != 100) errors <- c(errors, paste("RNA sample count", ncol(rna), "!= 100"))
if (ncol(meth) != 100) errors <- c(errors, paste("METH sample count", ncol(meth), "!= 100"))
if (length(unique(cnv$Sample)) != 100) errors <- c(errors, "CNV unique samples != 100")
if (length(unique(snv$Tumor_Sample_Barcode)) != 100) errors <- c(errors, "SNV unique samples != 100")

# Check dimensions are sufficient (e.g. > 100 features for ML)
if (nrow(rna) < 100) errors <- c(errors, "RNA genes < 100")
if (nrow(meth) < 100) errors <- c(errors, "METH probes < 100")

# 2. Patient IDs match metadata
cat("\n2. ID matching\n")
rna_samples <- colnames(rna)
meth_samples <- colnames(meth)
cnv_samples <- unique(cnv$Sample)
snv_samples <- unique(snv$Tumor_Sample_Barcode)

if (!all(rna_samples %in% meta$sample_id)) errors <- c(errors, "RNA samples not in metadata")
if (!all(meth_samples %in% meta$sample_id)) errors <- c(errors, "METH samples not in metadata")
if (!all(cnv_samples %in% meta$sample_id)) errors <- c(errors, "CNV samples not in metadata")
if (!all(snv_samples %in% meta$sample_id)) errors <- c(errors, "SNV samples not in metadata")

# 3. Metadata IDs match Clinical IDs
clin_id_col <- cmap$patient_id
if (!all(meta$patient_id %in% clin[[clin_id_col]])) errors <- c(errors, "Metadata patient IDs not in clinical data")
if (!all(clin[[clin_id_col]] %in% meta$patient_id)) errors <- c(errors, "Clinical patient IDs not in metadata")

# 4. No TCGA-style identifiers
tcga_regex <- "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}"
if (any(grepl(tcga_regex, meta$sample_id))) errors <- c(errors, "TCGA-style identifiers found in sample_id")
if (any(grepl(tcga_regex, meta$patient_id))) errors <- c(errors, "TCGA-style identifiers found in patient_id")

# 5. Required OmicsFlow columns
req_meta_cols <- c("sample_id", "patient_id", "sample_class", "batch", "center")
if (!all(req_meta_cols %in% colnames(meta))) errors <- c(errors, "Missing required metadata columns")

req_clin_mapped_cols <- c(cmap$patient_id, cmap$survival_time, cmap$survival_event, cmap$age, cmap$gender)
req_clin_mapped_cols <- req_clin_mapped_cols[!sapply(req_clin_mapped_cols, is.null)]
if (!all(req_clin_mapped_cols %in% colnames(clin))) errors <- c(errors, "Clinical map references missing columns in custom_clinical.tsv")

cat("\nErrors detected:\n")
if (length(errors) == 0) {
  cat("NONE\n")
} else {
  for (e in errors) cat("- ", e, "\n")
}
