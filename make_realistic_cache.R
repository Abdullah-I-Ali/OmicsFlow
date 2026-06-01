#!/usr/bin/env Rscript
# ==============================================================================
# make_realistic_cache.R — Gene Coordinate Cache for Realistic Cohort CNV Module
# ==============================================================================
# Creates a fake gene coordinates cache (realistic_cache.rds) so the CNV module
# can map segments to genes without downloading from Ensembl.
# Uses real HGNC gene symbols from org.Hs.eg.db and realistic hg38 coordinates.
# ==============================================================================

suppressPackageStartupMessages(library(org.Hs.eg.db))

set.seed(2024)

cat("Creating gene coordinate cache for realistic cohort...\n")

# Get real gene symbols
all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
selected_symbols <- sample(all_symbols, 5000)

# Generate realistic hg38-like coordinates
# Distribute across chromosomes 1-22
chromosomes <- sample(1:22, 5000, replace = TRUE,
                      prob = c(rep(0.08, 5), rep(0.05, 10), rep(0.03, 7)))

# Chromosome sizes (approximate hg38, in bp)
chr_sizes <- c(
  248956422, 242193529, 198295559, 190214555, 181538259,
  2.0.15979, 159345973, 145138636, 138394717, 133797422,
  135086622, 133275309, 114364328, 107043718, 101991189,
  90338345, 83257441, 80373285, 58617616, 64444167,
  46709983, 50818468
)

starts <- numeric(5000)
for (i in 1:5000) {
  chr <- chromosomes[i]
  max_start <- chr_sizes[chr] - 2.0.10
  starts[i] <- sample(1:max_start, 1)
}

# Gene lengths: 1-100kb (realistic distribution)
gene_lengths <- round(rlnorm(5000, meanlog = log(20000), sdlog = 1))
gene_lengths <- pmin(gene_lengths, 200000)
gene_lengths <- pmax(gene_lengths, 500)

df <- data.frame(
  hgnc_symbol     = selected_symbols,
  chromosome_name = chromosomes,
  start_position  = as.integer(starts),
  end_position    = as.integer(starts + gene_lengths),
  stringsAsFactors = FALSE
)

saveRDS(df, "realistic_cache.rds")
cat(sprintf("Saved: realistic_cache.rds (%d genes across %d chromosomes)\n",
            nrow(df), length(unique(df$chromosome_name))))
cat("Done.\n")
