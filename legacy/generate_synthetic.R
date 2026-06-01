suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(jsonlite))

out_dir <- "data/synthetic_cohort"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 1. IDENTIFIERS
n_samples <- 100
patients <- sprintf("PAT_%03d", 1:n_samples)
samples <- sprintf("SAMP_%03d", 1:n_samples)

# 2. METADATA
metadata <- data.frame(
  sample_id = samples,
  patient_id = patients,
  sample_class = rep("Tumor", n_samples),
  batch = sample(c("B1", "B2", "B3"), n_samples, replace = TRUE),
  center = sample(c("Center_X", "Center_Y"), n_samples, replace = TRUE),
  stringsAsFactors = FALSE
)
write.csv(metadata, file.path(out_dir, "sample_metadata.csv"), row.names = FALSE)

# 3. CLINICAL DATA
clinical <- data.frame(
  patient = patients,
  OS_time = round(runif(n_samples, 10, 3000)),
  OS_event = sample(c(0, 1), n_samples, replace = TRUE, prob = c(0.4, 0.6)),
  age = round(rnorm(n_samples, 60, 10)),
  sex = sample(c("Male", "Female"), n_samples, replace = TRUE),
  stringsAsFactors = FALSE
)
write.table(clinical, file.path(out_dir, "custom_clinical.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# 4. CLINICAL MAP JSON
clin_map <- list(
  patient_id = "patient",
  os_time = "OS_time",
  os_event = "OS_event",
  age = "age",
  gender = "sex"
)
write_json(clin_map, file.path(out_dir, "clinical_map.json"), auto_unbox = TRUE, pretty = TRUE)

# 5. RNA MATRIX
cat("Generating RNA...\n")
# Get some valid Ensembl IDs
all_genes <- keys(org.Hs.eg.db, keytype="ENSEMBL")
# some might not have symbols, let's pick 2000 random ones
set.seed(123)
selected_ensembl <- sample(all_genes, 2000)

rna_mat <- matrix(rnbinom(n_samples * 2000, size=1, mu=500), nrow=2000, ncol=n_samples)
rownames(rna_mat) <- selected_ensembl
colnames(rna_mat) <- samples
saveRDS(rna_mat, file.path(out_dir, "rna.rds"))

# 6. METHYLATION MATRIX
cat("Generating METH...\n")
probes <- sprintf("cg%08d", 1:5000)
meth_mat <- matrix(runif(n_samples * 5000, 0, 1), nrow=5000, ncol=n_samples)
rownames(meth_mat) <- probes
colnames(meth_mat) <- samples
saveRDS(meth_mat, file.path(out_dir, "meth.rds"))

# 7. CNV MATRIX
cat("Generating CNV...\n")
cnv_list <- lapply(1:n_samples, function(i) {
  n_seg <- sample(10:50, 1)
  data.frame(
    Sample = samples[i],
    Chromosome = sample(paste0("chr", 1:22), n_seg, replace=TRUE),
    Start = sample(10000:50000000, n_seg),
    End = sample(50000001:100000000, n_seg),
    Segment_Mean = rnorm(n_seg, mean=0, sd=1.5),
    Num_Probes = sample(10:500, n_seg, replace=TRUE),
    stringsAsFactors = FALSE
  )
})
cnv_df <- do.call(rbind, cnv_list)
saveRDS(cnv_df, file.path(out_dir, "cnv.rds"))

# 8. SNV MATRIX
cat("Generating SNV...\n")
# Get valid symbols for SNV
all_symbols <- keys(org.Hs.eg.db, keytype="SYMBOL")
selected_symbols <- sample(all_symbols, 500)

snv_list <- lapply(1:n_samples, function(i) {
  n_mut <- sample(5:30, 1)
  data.frame(
    Tumor_Sample_Barcode = samples[i],
    Hugo_Symbol = sample(selected_symbols, n_mut, replace=TRUE),
    Chromosome = sample(1:22, n_mut, replace=TRUE),
    Start_Position = sample(10000:50000000, n_mut),
    End_Position = sample(10000:50000000, n_mut), # dummy
    Variant_Classification = sample(c("Missense_Mutation", "Nonsense_Mutation", "Silent"), n_mut, replace=TRUE),
    Variant_Type = "SNP",
    Reference_Allele = "A",
    Tumor_Seq_Allele2 = "T",
    stringsAsFactors = FALSE
  )
})
snv_df <- do.call(rbind, snv_list)
saveRDS(snv_df, file.path(out_dir, "snv.rds"))

cat("All synthetic data generated successfully.\n")
