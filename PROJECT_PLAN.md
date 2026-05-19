# OmicsFlow MVP — System Prompt & Engineering Specification

# INTERACTION RULES & EXECUTION STRATEGY

1. I will provide ONE large master file containing all legacy R scripts for the TCGA-LIHC project.

2. Read and analyze the ENTIRE master script before generating any code.

3. Understand:

   * global workflow structure
   * omics interactions
   * preprocessing dependencies
   * integration logic
   * shared utilities
   * sample ID handling

4. DO NOT generate all modules at once.
   This causes:

   * context loss
   * skipped preprocessing steps
   * scientific inconsistencies

5. Wait for my explicit command before starting any module.
   Example:

   * "Start Phase 1: RNA Module"
   * "Start Phase 2: Methylation Module"

6. Before refactoring any module:

   * extract and summarize all detected scientific steps
   * extract thresholds
   * extract normalization methods
   * extract required packages
   * extract outputs
   * confirm understanding BEFORE coding

7. Ask for clarification if:

   * any preprocessing step is ambiguous
   * any threshold is unclear
   * any integration logic is uncertain

8. Preserve ALL scientific preprocessing logic exactly.

9. Never simplify or rewrite validated preprocessing methodology.

---

# PROJECT GOAL

Build a modular automated multi-omics preprocessing and integration framework called:

# OmicsFlow

The framework should:

* accept at least 2 omics layers
* preprocess datasets automatically
* generate QC outputs
* generate integration-ready matrices
* generate ML-ready matrices
* generate publication-ready outputs
* generate HTML reports

---

# SUPPORTED OMICS LAYERS (MVP)

## RNA-seq

Supported formats:

* `.rds`
* `.txt`

## DNA Methylation

Supported formats:

* `.rds`
* `.csv`

## CNV

Supported formats:

* segmentation `.txt`

## SNV

Supported formats:

* `.maf`

---

# IMPORTANT EXECUTION RULE

The user MUST provide:

# at least 2 omics layers

Single-omics analysis is NOT supported in MVP.

---

# CRITICAL SCIENTIFIC RULES

DO NOT:

* modify preprocessing thresholds
* modify normalization methods
* modify batch correction logic
* modify feature-selection methodology
* simplify preprocessing logic
* remove QC steps
* rewrite scientific analysis
* replace validated methods
* alter barcode logic
* change patient matching logic

The current R scripts are scientifically validated.

You may ONLY:

* reorganize
* modularize
* standardize
* automate
* improve readability
* improve maintainability
* improve workflow structure

---

# PROJECT STRUCTURE

omicsflow/
│
├── main.nf
├── nextflow.config
│
├── modules/
│   ├── rna/
│   ├── methylation/
│   ├── cnv/
│   ├── snv/
│   └── integration/
│
├── reports/
├── configs/
├── docker/
├── envs/
├── tests/
├── results/
└── data/

---

# CODING STANDARDS

* Use modular function-based R code.
* Avoid duplicated logic.
* Use descriptive variable names.
* Separate preprocessing, QC, export, and utility logic.
* Use relative paths only.
* Do not hardcode file locations.
* Add comments explaining critical scientific operations.
* Ensure all modules are executable independently.
* Preserve Linux compatibility.
* Prefer memory-efficient operations for large omics datasets.
* Preserve compatibility with Nextflow execution.

---

# MODULE RESPONSIBILITIES

# RNA MODULE

Responsibilities:

* preprocessing
* normalization
* filtering
* batch correction
* feature selection
* QC
* exports

---

# METHYLATION MODULE

Responsibilities:

* probe filtering
* beta to M-value transformation
* ComBat correction
* feature selection
* QC
* exports

---

# CNV MODULE

Responsibilities:

* barcode parsing
* segment processing
* GRanges overlap
* gene mapping
* aggregation
* QC

---

# SNV MODULE

Responsibilities:

* mutation preprocessing
* mutation filtering
* mutation matrix generation
* QC

---

# STANDARDIZED OUTPUTS

Each module MUST output:

* [omics]_processed_matrix.rds (e.g., rna_processed_matrix.rds, methylation_processed_matrix.rds)
* sample_metadata.csv
* qc_metrics.json
* plots/

---

# OUTPUT DIRECTORY STRUCTURE

results/
│
├── rna/
│   ├── rna_processed_matrix.rds
│   ├── sample_metadata.csv
│   ├── qc_metrics.json
│   └── plots/
│
├── methylation/
│   ├── methylation_processed_matrix.rds
│   ├── sample_metadata.csv
│   ├── qc_metrics.json
│   └── plots/
│
├── cnv/
│   ├── cnv_processed_matrix.rds
│   ├── sample_metadata.csv
│   ├── qc_metrics.json
│   └── plots/
│
├── snv/
│   ├── snv_processed_matrix.rds
│   ├── sample_metadata.csv
│   ├── qc_metrics.json
│   └── plots/
│
└── integration/
├── integrated_ml_matrix.rds
├── integrated_mofa_matrix.rds
├── sample_overlap.csv
└── plots/

---

# INTEGRATION RESPONSIBILITIES

The integration layer MUST perform:

* sample intersection
* patient matching
* omics alignment
* integration-ready matrix generation
* ML-ready matrix generation
* MOFA-ready matrix generation

Use:

# 12-character TCGA patient IDs

for all omics matching.

---

# NEXTFLOW RESPONSIBILITIES

Nextflow should:

* validate inputs
* detect omics layers
* run modules in parallel
* manage dependencies
* launch integration
* launch report generation

IMPORTANT:

# Nextflow should NOT perform scientific analysis.

Scientific analysis remains inside R scripts.

---

# REPORTING

Use Quarto to generate:

* HTML reports
* QC summaries
* PCA plots
* heatmaps
* integration statistics
* publication-ready figures
* downloadable outputs

The report should dynamically adapt to:

* available omics layers
* generated outputs
* available QC metrics

---

# VALIDATION REQUIREMENTS

Every refactored module MUST:

* reproduce original outputs
* preserve dimensions
* preserve sample IDs
* preserve preprocessing behavior
* preserve normalization outputs
* preserve feature-selection outputs

Validation tests MUST compare:

* dimensions
* ranges
* variance
* feature counts
* sample overlap
* output consistency

---

# REQUIRED SCIENTIFIC STEPS

# RNA MODULE

The RNA preprocessing module MUST preserve ALL existing steps from the legacy scripts:

* primary tumor filtering
* patient deduplication
* library-size QC
* zero-gene removal
* CPM filtering
* TMM normalization
* log2 CPM transformation
* Ensembl-to-symbol mapping
* duplicate gene-symbol resolution
* NA filtering
* NA imputation
* zero-variance filtering
* outlier detection
* Pearson correlation QC
* batch correction using removeBatchEffect
* PCA before correction
* PCA after correction
* top variable gene selection
* Z-score generation
* ML-ready matrix generation
* validation plots
* density plots
* heatmaps
* statistical validation figures

---

# METHYLATION MODULE

The methylation preprocessing module MUST preserve ALL existing steps from the legacy scripts:

* primary tumor filtering
* patient deduplication
* detection p-value filtering
* cross-reactive probe filtering
* SNP-associated probe filtering
* sex chromosome filtering
* non-CpG probe filtering
* NA filtering
* beta-value cleanup
* beta-to-M transformation
* KNN imputation
* median imputation
* ComBat batch correction
* clinical covariate protection
* PCA QC
* outlier detection
* top variable probe selection
* ML-ready matrix generation
* MOFA-ready matrix generation
* beta density plots
* M-value density plots
* correlation QC
* validation figures

---

# CNV MODULE

The CNV preprocessing module MUST preserve ALL existing steps from the legacy scripts:

* barcode parsing
* merged barcode separation
* tumor-only filtering
* chromosome standardization
* GRanges overlap analysis
* hg38 validation
* gene coordinate mapping
* segment aggregation
* gene-level matrix construction
* variance filtering
* QC plots

---

# SNV MODULE

The SNV module MUST preserve:

* mutation preprocessing
* mutation filtering
* mutation matrix generation
* QC summaries

---

# IMPORTANT

DO NOT:

* skip preprocessing steps
* simplify preprocessing logic
* remove validation logic
* remove plots
* remove QC
* replace scientific methods

---

# DEVELOPMENT PHASES

# Phase 1

* RNA module
* Methylation module

# Phase 2

* Integration layer
* Quarto reports

# Phase 3

* Nextflow automation

# Phase 4

* Docker reproducibility

---

# FINAL MVP GOAL

The researcher should be able to run:

omicsflow run --rna rna.rds --meth meth.rds

And automatically receive:

* preprocessing outputs
* QC outputs
* integration outputs
* ML-ready matrices
* publication-ready figures
* HTML report
