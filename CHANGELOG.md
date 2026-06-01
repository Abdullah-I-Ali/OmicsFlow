# Changelog

All notable changes to this project will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.1] - 2026-06-01

### Added
- Zenodo DOI integration and citation framework.

### Fixed
- Fixed typo in `validate_inputs` logging (`log_warning` -> `log_warn`).
- GitHub Actions CI now correctly runs realistic cohort generation prior to testing.

## [1.0.0] - 2026-05-19

### Added
- **Orchestration Layer:** Implemented end-to-end multi-omics execution workflow using Nextflow v23+.
- **RNA-Seq Module:** TMM library size normalization, outlier rejection, log2 CPM conversion, and batch effect correction using `limma::removeBatchEffect`.
- **DNA Methylation Module:** Missingness-based probe filtering, invariant CpG deletion, cross-reactive probe elimination, and high-variance feature selection.
- **CNV Module:** Genomic mapping translating segment-level files to gene-level data matrices.
- **SNV Module:** Occurrence matrix parser filtering mutations by custom frequency thresholds.
- **MOFA+ Integration Module:** Automated training, latent factor variance explanation maps, and factor weight extraction.
- **Survival Machine Learning Module:** Performance comparison of LASSO Cox regression, Random Survival Forest, and XGBoost using cross-validation to isolate prognostic signatures.
- **Pathway Enrichment Module:** Functional pathway analysis across Gene Ontology (BP, MF, CC) and KEGG databases via `clusterProfiler`.
- **Automated Quarto Reporting:** Professional CSS layout generating standalone HTML reports with interactive Lightbox overlays for high-resolution graphics.
- **Reproducibility Layer:** Version-pinned Conda environment configurations and multi-platform SLURM/Docker execution profiles.

### Fixed
- **Reporting Data-Source Bug:** Synchronized report enrichment counters with actual pipeline run outputs by correcting pathways configuration (`qc_metrics.json` path mapping).
- **Metric Formatting Errors:** Corrected grammatical and trailing punctuation formatting in dynamic multi-view text generation within report chunks.
- **Summary Status Rows:** Updated final pipeline table to list all active modules and correctly report integrated metrics.
