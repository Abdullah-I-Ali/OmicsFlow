#!/usr/bin/env Rscript
# ==============================================================================
# generate_realistic_cohort.R — Realistic Oncology Cohort Generator
# OmicsFlow Validation | Biologically Structured Multi-Omics Dataset
# ==============================================================================
#
# PURPOSE:
#   Generate a realistic multi-omics oncology cohort with:
#     - 180 patients across 3 biologically distinct subtypes
#     - Subtype-specific molecular markers
#     - Survival-associated molecular features
#     - Batch effects and center effects
#     - Clinically realistic survival distributions
#
# SUBTYPES:
#   A — "Proliferative"   (60 pts): TP53-mutant, MYC-amplified, worst prognosis
#   B — "Mesenchymal"     (60 pts): ECM signature, intermediate prognosis
#   C — "Immune-enriched" (60 pts): Immune markers, best prognosis
#
# OUTPUT: data/realistic_cohort/
# ==============================================================================

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(jsonlite)
})

set.seed(2024)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  OmicsFlow | Realistic Oncology Cohort Generator\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

out_dir <- "data/realistic_cohort"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==============================================================================
# 1. COHORT DESIGN
# ==============================================================================
cat("[1/8] Defining cohort structure...\n")

n_per_subtype <- 60
n_samples     <- n_per_subtype * 3
subtypes      <- rep(c("Proliferative", "Mesenchymal", "Immune_enriched"),
                     each = n_per_subtype)

patients <- sprintf("PAT_%03d", 1:n_samples)
samples  <- sprintf("SAMP_%03d", 1:n_samples)

# Assign batches with partial subtype confounding (realistic)
# B1, B2 → Center_Alpha; B3 → Center_Beta
# Proliferative enriched in B1/B2, Immune enriched in B3, Mesenchymal spread
batch_assignments <- character(n_samples)
center_assignments <- character(n_samples)
for (i in 1:n_samples) {
  if (subtypes[i] == "Proliferative") {
    # 80% in B1/B2, 20% in B3
    batch_assignments[i] <- sample(c("B1", "B2", "B3"), 1, prob = c(0.40, 0.40, 0.20))
  } else if (subtypes[i] == "Mesenchymal") {
    # Evenly distributed
    batch_assignments[i] <- sample(c("B1", "B2", "B3"), 1, prob = c(0.33, 0.33, 0.34))
  } else {
    # 20% in B1/B2, 80% in B3
    batch_assignments[i] <- sample(c("B1", "B2", "B3"), 1, prob = c(0.10, 0.10, 0.80))
  }
  # Distribute batches across centers so they are not perfectly confounded
  center_assignments[i] <- sample(c("Center_Alpha", "Center_Beta"), 1, prob = c(0.6, 0.4))
}

cat(sprintf("  Patients: %d (3 subtypes × %d)\n", n_samples, n_per_subtype))
cat(sprintf("  Batches: B1=%d, B2=%d, B3=%d\n",
            sum(batch_assignments == "B1"),
            sum(batch_assignments == "B2"),
            sum(batch_assignments == "B3")))
cat(sprintf("  Centers: Alpha=%d, Beta=%d\n",
            sum(center_assignments == "Center_Alpha"),
            sum(center_assignments == "Center_Beta")))

# ==============================================================================
# 2. METADATA
# ==============================================================================
cat("[2/8] Generating sample metadata...\n")

metadata <- data.frame(
  sample_id    = samples,
  patient_id   = patients,
  sample_class = rep("Tumor", n_samples),
  batch        = batch_assignments,
  center       = center_assignments,
  stringsAsFactors = FALSE
)
write.csv(metadata, file.path(out_dir, "sample_metadata.csv"), row.names = FALSE)
cat(sprintf("  Saved: sample_metadata.csv (%d rows)\n", nrow(metadata)))

# ==============================================================================
# 3. SURVIVAL DATA (Weibull per subtype + molecular modulation)
# ==============================================================================
cat("[3/8] Generating clinically realistic survival data...\n")

# Age distribution — slightly different per subtype
ages <- numeric(n_samples)
ages[subtypes == "Proliferative"]   <- round(rnorm(sum(subtypes == "Proliferative"), 58, 11))
ages[subtypes == "Mesenchymal"]     <- round(rnorm(sum(subtypes == "Mesenchymal"), 63, 10))
ages[subtypes == "Immune_enriched"] <- round(rnorm(sum(subtypes == "Immune_enriched"), 55, 12))
ages <- pmax(25, pmin(85, ages))

sex <- sample(c("Male", "Female"), n_samples, replace = TRUE, prob = c(0.55, 0.45))

# Stage distribution — subtype-dependent
stages <- character(n_samples)
for (i in 1:n_samples) {
  if (subtypes[i] == "Proliferative") {
    stages[i] <- sample(c("I", "II", "III", "IV"), 1, prob = c(0.05, 0.15, 0.40, 0.40))
  } else if (subtypes[i] == "Mesenchymal") {
    stages[i] <- sample(c("I", "II", "III", "IV"), 1, prob = c(0.10, 0.30, 0.35, 0.25))
  } else {
    stages[i] <- sample(c("I", "II", "III", "IV"), 1, prob = c(0.25, 0.35, 0.25, 0.15))
  }
}

# Weibull survival with subtype-specific parameters
# shape < 1 → decreasing hazard, shape > 1 → increasing hazard
weibull_shape <- ifelse(subtypes == "Proliferative", 0.9,
                 ifelse(subtypes == "Mesenchymal", 1.0, 1.1))
weibull_scale <- ifelse(subtypes == "Proliferative", 500,
                 ifelse(subtypes == "Mesenchymal", 1100, 2000))

# Individual hazard modulation based on age & stage
hazard_mod <- rep(1.0, n_samples)
hazard_mod <- hazard_mod * (1 + (ages - 60) * 0.005)  # older → slightly worse
hazard_mod[stages == "III"] <- hazard_mod[stages == "III"] * 1.2
hazard_mod[stages == "IV"]  <- hazard_mod[stages == "IV"]  * 1.5

# Generate true survival times (Weibull)
true_time <- rweibull(n_samples,
                      shape = weibull_shape,
                      scale = weibull_scale / hazard_mod)
true_time <- pmax(10, round(true_time))

# Administrative censoring at 2500 days + random censoring
admin_censor <- 2500
random_censor <- runif(n_samples, 200, 3000)
censor_time <- pmin(admin_censor, random_censor)

os_event <- as.integer(true_time <= censor_time)
os_time  <- ifelse(os_event == 1, true_time, pmin(true_time, censor_time))
os_time  <- pmax(10, round(os_time))

cat(sprintf("  Overall event rate: %.1f%%\n", 100 * mean(os_event)))
for (st in c("Proliferative", "Mesenchymal", "Immune_enriched")) {
  idx <- subtypes == st
  cat(sprintf("  %s: median OS=%.0f days, event rate=%.1f%%\n",
              st, median(os_time[idx]), 100 * mean(os_event[idx])))
}

clinical <- data.frame(
  patient  = patients,
  OS_time  = os_time,
  OS_event = os_event,
  age      = ages,
  sex      = sex,
  subtype  = subtypes,
  stage    = stages,
  stringsAsFactors = FALSE
)
write.table(clinical, file.path(out_dir, "custom_clinical.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Clinical map JSON
clin_map <- list(
  patient_id = "patient",
  os_time    = "OS_time",
  os_event   = "OS_event",
  age        = "age",
  gender     = "sex"
)
write_json(clin_map, file.path(out_dir, "clinical_map.json"),
           auto_unbox = TRUE, pretty = TRUE)
cat("  Saved: custom_clinical.tsv + clinical_map.json\n")

# ==============================================================================
# 4. RNA EXPRESSION (2000 genes × 180 samples)
# ==============================================================================
cat("[4/8] Generating RNA expression with biological signal...\n")

# Use real biological pathways so enrichment analysis finds true signal
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))

# GO:0007049 (cell cycle), GO:0030198 (ECM), GO:0006955 (immune response)
go_prolif <- na.omit(unique(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys="GO:0007049", columns="ENSEMBL", keytype="GOALL")$ENSEMBL)))
go_mesench <- na.omit(unique(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys="GO:0030198", columns="ENSEMBL", keytype="GOALL")$ENSEMBL)))
go_immune <- na.omit(unique(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys="GO:0006955", columns="ENSEMBL", keytype="GOALL")$ENSEMBL)))

marker_block_size <- 50
prolif_markers_genes <- sample(go_prolif, marker_block_size)
mesench_markers_genes <- sample(go_mesench, marker_block_size)
immune_markers_genes <- sample(go_immune, marker_block_size)

all_genes <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
pool <- setdiff(all_genes, c(prolif_markers_genes, mesench_markers_genes, immune_markers_genes))
random_genes <- sample(pool, 2000 - 3 * marker_block_size)

selected_ensembl <- c(prolif_markers_genes, mesench_markers_genes, immune_markers_genes, random_genes)

# Base expression: negative binomial (realistic count data)
rna_mat <- matrix(rnbinom(n_samples * 2000, size = 2, mu = 300),
                  nrow = 2000, ncol = n_samples)
rownames(rna_mat) <- selected_ensembl
colnames(rna_mat) <- samples

# --- Subtype-specific marker genes ---
# These now map to real Ensembl IDs from their respective GO pathways
prolif_markers <- 1:marker_block_size
mesench_markers <- (marker_block_size + 1):(2 * marker_block_size)
immune_markers  <- (2 * marker_block_size + 1):(3 * marker_block_size)

for (i in 1:n_samples) {
  if (subtypes[i] == "Proliferative") {
    rna_mat[prolif_markers, i] <- rna_mat[prolif_markers, i] * 
      sample(seq(3, 6, by = 0.5), marker_block_size, replace = TRUE)
  } else if (subtypes[i] == "Mesenchymal") {
    rna_mat[mesench_markers, i] <- rna_mat[mesench_markers, i] * 
      sample(seq(3, 6, by = 0.5), marker_block_size, replace = TRUE)
  } else {
    rna_mat[immune_markers, i] <- rna_mat[immune_markers, i] * 
      sample(seq(3, 6, by = 0.5), marker_block_size, replace = TRUE)
  }
}

# --- Survival-associated genes ---
# Top 20 genes (indices 151-170) have expression correlated with survival
surv_gene_idx <- 151:170
# Create a survival risk score based on clinical data
risk_score <- scale(-os_time + rnorm(n_samples, 0, 200))  # higher = worse
for (g in surv_gene_idx) {
  # Add expression component correlated with risk
  signal <- as.numeric(risk_score) * runif(1, 100, 300)
  rna_mat[g, ] <- pmax(0, rna_mat[g, ] + signal)
}

# --- Batch effects (additive, moderate) ---
batch_shift <- list(
  B1 = rnorm(2000, mean = 0,   sd = 30),
  B2 = rnorm(2000, mean = 50,  sd = 30),
  B3 = rnorm(2000, mean = -40, sd = 30)
)
for (i in 1:n_samples) {
  rna_mat[, i] <- pmax(0, rna_mat[, i] + batch_shift[[batch_assignments[i]]])
}

# Ensure integer counts
storage.mode(rna_mat) <- "integer"
rna_mat[rna_mat < 0] <- 0L

saveRDS(rna_mat, file.path(out_dir, "rna.rds"))
cat(sprintf("  RNA matrix: %d genes × %d samples\n", nrow(rna_mat), ncol(rna_mat)))
cat(sprintf("  Expression range: [%d, %d]\n", min(rna_mat), max(rna_mat)))

# ==============================================================================
# 5. DNA METHYLATION (5000 probes × 180 samples)
# ==============================================================================
cat("[5/8] Generating DNA methylation with subtype patterns...\n")

# Use REAL Illumina 450k probe names from annotation
# This is essential so preprocess_meth.R can match against the 450k annotation
suppressPackageStartupMessages({
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
real_probes <- rownames(anno)
# Select 5000 random real probes (excluding sex chromosomes and non-CpG)
chr_col <- anno$chr
valid_idx <- which(!chr_col %in% c("chrX", "chrY") & grepl("^cg", real_probes))
probes <- sample(real_probes[valid_idx], 5000)

# Base beta values: bimodal distribution (realistic)
meth_mat <- matrix(NA_real_, nrow = 5000, ncol = n_samples)
for (j in 1:5000) {
  # Bimodal: ~40% probes hypomethylated, ~60% hypermethylated
  if (runif(1) < 0.4) {
    meth_mat[j, ] <- rbeta(n_samples, shape1 = 2, shape2 = 8)  # low beta
  } else {
    meth_mat[j, ] <- rbeta(n_samples, shape1 = 8, shape2 = 2)  # high beta
  }
}
rownames(meth_mat) <- probes
colnames(meth_mat) <- samples

# --- Subtype-specific methylation patterns ---
# Block 1 (probes 1-200): hypomethylated in Proliferative
# Block 2 (probes 201-400): intermediate in Mesenchymal
# Block 3 (probes 401-600): hypermethylated in Immune_enriched

for (i in 1:n_samples) {
  if (subtypes[i] == "Proliferative") {
    meth_mat[1:200, i] <- rbeta(200, shape1 = 1.5, shape2 = 10)    # hypo
    meth_mat[201:400, i] <- rbeta(200, shape1 = 4, shape2 = 4)     # mid
    meth_mat[401:600, i] <- rbeta(200, shape1 = 6, shape2 = 3)     # mid-high
  } else if (subtypes[i] == "Mesenchymal") {
    meth_mat[1:200, i] <- rbeta(200, shape1 = 4, shape2 = 4)       # mid
    meth_mat[201:400, i] <- rbeta(200, shape1 = 2, shape2 = 8)     # hypo (ECM)
    meth_mat[401:600, i] <- rbeta(200, shape1 = 5, shape2 = 3)     # mid-high
  } else {
    meth_mat[1:200, i] <- rbeta(200, shape1 = 5, shape2 = 3)       # mid-high
    meth_mat[201:400, i] <- rbeta(200, shape1 = 5, shape2 = 3)     # mid-high
    meth_mat[401:600, i] <- rbeta(200, shape1 = 10, shape2 = 1.5)  # hyper
  }
}

# --- Survival-linked probes (probes 601-630) ---
surv_probe_idx <- 601:630
for (p in surv_probe_idx) {
  # Probes correlated with survival: high beta → worse prognosis
  meth_mat[p, ] <- pmin(1, pmax(0,
    0.5 + as.numeric(risk_score) * runif(1, 0.05, 0.15) + rnorm(n_samples, 0, 0.08)
  ))
}

# --- Batch effects (small beta offset) ---
batch_meth_shift <- list(
  B1 = 0.00,
  B2 = 0.03,
  B3 = -0.02
)
for (i in 1:n_samples) {
  meth_mat[, i] <- meth_mat[, i] + batch_meth_shift[[batch_assignments[i]]]
}

# Clamp to [0, 1] while preserving matrix dimensions
meth_mat[meth_mat < 0.001] <- 0.001
meth_mat[meth_mat > 0.999] <- 0.999

saveRDS(meth_mat, file.path(out_dir, "meth.rds"))
cat(sprintf("  Methylation matrix: %d probes × %d samples\n", nrow(meth_mat), ncol(meth_mat)))
cat(sprintf("  Beta range: [%.3f, %.3f]\n", min(meth_mat), max(meth_mat)))

# ==============================================================================
# 6. COPY NUMBER VARIATION (segment data)
# ==============================================================================
cat("[6/8] Generating CNV segments with subtype-specific amplifications...\n")

# Gene symbols for CNV mapping
all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
# Pick a stable set for the cache
cnv_gene_symbols <- sample(all_symbols, 5000)

cnv_list <- lapply(1:n_samples, function(i) {
  # Base: 20-60 random segments
  n_seg <- sample(20:60, 1)

  chrs <- sample(paste0("chr", 1:22), n_seg, replace = TRUE)
  starts <- sample(10000:150000000, n_seg)
  ends <- starts + sample(50000:5000000, n_seg, replace = TRUE)
  seg_mean <- rnorm(n_seg, mean = 0, sd = 0.5)  # background noise
  n_probes <- sample(20:500, n_seg, replace = TRUE)

  # --- Subtype-specific focal events ---
  if (subtypes[i] == "Proliferative") {
    # MYC amplification (chr8q), CDKN2A deletion (chr9p)
    extra <- data.frame(
      Sample    = samples[i],
      Chromosome = c("chr8", "chr9"),
      Start     = c(127000000, 21900000),
      End       = c(128500000, 22100000),
      Segment_Mean = c(runif(1, 1.5, 3.0), runif(1, -3.0, -1.5)),
      Num_Probes   = c(sample(100:300, 1), sample(50:150, 1)),
      stringsAsFactors = FALSE
    )
  } else if (subtypes[i] == "Mesenchymal") {
    # chr1q gain, chr16q loss
    extra <- data.frame(
      Sample    = samples[i],
      Chromosome = c("chr1", "chr16"),
      Start     = c(145000000, 46000000),
      End       = c(155000000, 90000000),
      Segment_Mean = c(runif(1, 0.5, 1.5), runif(1, -1.5, -0.5)),
      Num_Probes   = c(sample(100:400, 1), sample(100:300, 1)),
      stringsAsFactors = FALSE
    )
  } else {
    # Immune: fewer CNV events, smaller magnitude
    extra <- data.frame(
      Sample    = samples[i],
      Chromosome = c("chr6"),
      Start     = c(29000000),
      End       = c(33000000),
      Segment_Mean = c(runif(1, -0.5, 0.5)),
      Num_Probes   = c(sample(50:200, 1)),
      stringsAsFactors = FALSE
    )
  }

  base_df <- data.frame(
    Sample    = samples[i],
    Chromosome = chrs,
    Start     = starts,
    End       = ends,
    Segment_Mean = seg_mean,
    Num_Probes   = n_probes,
    stringsAsFactors = FALSE
  )

  rbind(base_df, extra)
})

cnv_df <- do.call(rbind, cnv_list)

# Survival-linked: patients with high risk → larger CNV amplitude
for (i in 1:n_samples) {
  idx <- cnv_df$Sample == samples[i]
  amp_factor <- 1 + as.numeric(risk_score[i]) * 0.15
  cnv_df$Segment_Mean[idx] <- cnv_df$Segment_Mean[idx] * max(0.3, amp_factor)
}

saveRDS(cnv_df, file.path(out_dir, "cnv.rds"))
cat(sprintf("  CNV segments: %d total across %d patients\n",
            nrow(cnv_df), length(unique(cnv_df$Sample))))

# ==============================================================================
# 7. SOMATIC MUTATIONS (SNV, MAF-like)
# ==============================================================================
cat("[7/8] Generating SNV mutations with driver gene enrichment...\n")

# Define driver genes per subtype
driver_genes <- list(
  Proliferative   = c("TP53", "MYC", "CCND1", "CDK4", "RB1",
                       "CDKN2A", "E2F1", "AURKA", "PLK1", "BUB1"),
  Mesenchymal     = c("COL1A1", "FN1", "TGFBR2", "SMAD4", "VIM",
                       "SNAI1", "TWIST1", "ZEB1", "MMP2", "MMP9"),
  Immune_enriched = c("CD8A", "PDCD1", "CTLA4", "LAG3", "HAVCR2",
                       "IDO1", "CD274", "STAT1", "JAK2", "IRF1")
)

# Passenger genes — random background
passenger_genes <- sample(setdiff(all_symbols,
                                  unlist(driver_genes)), 400)

snv_list <- lapply(1:n_samples, function(i) {
  # TMB varies by subtype: Proliferative = high, Immune = moderate, Mesenchymal = low
  if (subtypes[i] == "Proliferative") {
    n_driver <- sample(3:7, 1)
    n_passenger <- sample(10:25, 1)
    drivers <- driver_genes$Proliferative
  } else if (subtypes[i] == "Mesenchymal") {
    n_driver <- sample(2:5, 1)
    n_passenger <- sample(5:15, 1)
    drivers <- driver_genes$Mesenchymal
  } else {
    n_driver <- sample(2:6, 1)
    n_passenger <- sample(8:20, 1)
    drivers <- driver_genes$Immune_enriched
  }

  # Cross-subtype TP53 mutations (seen across subtypes but enriched in Proliferative)
  if (subtypes[i] != "Proliferative" && runif(1) < 0.15) {
    drivers <- c(drivers, "TP53")
  }

  selected_drivers <- sample(drivers, min(n_driver, length(drivers)))
  selected_passengers <- sample(passenger_genes, n_passenger)
  all_muts <- c(selected_drivers, selected_passengers)

  n_mut <- length(all_muts)

  data.frame(
    Tumor_Sample_Barcode = samples[i],
    Hugo_Symbol = all_muts,
    Chromosome = sample(1:22, n_mut, replace = TRUE),
    Start_Position = sample(10000:150000000, n_mut),
    End_Position = sample(10000:150000000, n_mut),
    Variant_Classification = sample(
      c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
        "Frame_Shift_Ins", "Splice_Site"),
      n_mut, replace = TRUE,
      prob = c(0.65, 0.15, 0.08, 0.07, 0.05)
    ),
    Variant_Type = "SNP",
    Reference_Allele = sample(c("A", "C", "G", "T"), n_mut, replace = TRUE),
    Tumor_Seq_Allele2 = sample(c("A", "C", "G", "T"), n_mut, replace = TRUE),
    stringsAsFactors = FALSE
  )
})

snv_df <- do.call(rbind, snv_list)

# TP53 mutations associated with worse survival — ensure presence in worst-prognosis patients
tp53_patients <- snv_df$Tumor_Sample_Barcode[snv_df$Hugo_Symbol == "TP53"]
cat(sprintf("  TP53 mutated in %d patients\n", length(unique(tp53_patients))))

saveRDS(snv_df, file.path(out_dir, "snv.rds"))
cat(sprintf("  SNV mutations: %d total across %d patients\n",
            nrow(snv_df), length(unique(snv_df$Tumor_Sample_Barcode))))

# ==============================================================================
# 8. SUMMARY
# ==============================================================================
cat("\n[8/8] Generating cohort summary...\n")

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║          REALISTIC COHORT GENERATION COMPLETE               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat(sprintf("\nOutput directory: %s/\n", out_dir))
cat(sprintf("  sample_metadata.csv  : %d samples\n", nrow(metadata)))
cat(sprintf("  custom_clinical.tsv  : %d patients\n", nrow(clinical)))
cat(sprintf("  clinical_map.json    : %d mappings\n", length(clin_map)))
cat(sprintf("  rna.rds              : %d × %d\n", nrow(rna_mat), ncol(rna_mat)))
cat(sprintf("  meth.rds             : %d × %d\n", nrow(meth_mat), ncol(meth_mat)))
cat(sprintf("  cnv.rds              : %d segments\n", nrow(cnv_df)))
cat(sprintf("  snv.rds              : %d mutations\n", nrow(snv_df)))

cat("\nSubtype distribution:\n")
print(table(subtypes))

cat("\nSurvival summary by subtype:\n")
for (st in c("Proliferative", "Mesenchymal", "Immune_enriched")) {
  idx <- subtypes == st
  cat(sprintf("  %-20s: n=%d, median_OS=%.0f, events=%d (%.0f%%)\n",
              st, sum(idx), median(os_time[idx]),
              sum(os_event[idx]), 100 * mean(os_event[idx])))
}

cat("\nBatch × Subtype cross-table:\n")
print(table(Batch = batch_assignments, Subtype = subtypes))

cat("\nAll files generated successfully.\n")
