# Implementation Plan - Phase 1: Universal Metadata Layer

Implement Phase 1: Universal Metadata Layer for OmicsFlow to allow flexible, metadata-driven sample tracking and batch correction while maintaining absolute backward compatibility with existing TCGA workflows.

## User Review Required

> [!IMPORTANT]
> - **Zero-Impact Fallback**: When `--metadata` is not supplied, the pipeline will fallback 100% to the current TCGA parsing rules. No existing clinical, TCGA barcode, or plate ID behavior is altered.
> - **Unified Metadata Schema**: Standardized output `sample_metadata.csv` will now conform to the universal schema (`sample_id`, `patient_id`, `sample_class`, `batch`, `center` [optional]) for both TCGA mode and metadata mode. For TCGA mode, fields are auto-filled using TCGA parsing, which keeps existing columns `patient_id` and `batch` intact to ensure backward compatibility with current validation tests.

## Open Questions

There are no open questions; all requirements are fully specified and the design aligns precisely with the TCGA audit report.

## Proposed Changes

### Shared Metadata Component

Standardize the metadata schema and loading utilities.

#### [NEW] [utils_metadata.R](file:///c:/Users/ABDULLAH/Desktop/omicsflow/modules/utils_metadata.R)
- Implement schema parsing and validation for the required columns: `sample_id`, `patient_id`, `sample_class`, `batch`, and `center` (optional).
- Add helper functions `get_patient_id`, `get_sample_class`, `get_batch`, and `get_center` supporting both metadata lookup and TCGA fallback parsing.

---

### RNA Module

#### [MODIFY] [preprocess_rna.R](file:///c:/Users/ABDULLAH/Desktop/omicsflow/modules/rna/preprocess_rna.R)
- Add `--metadata` argument.
- Replace manual substring operations (`substr`) with shared helper functions.
- Update sample selection, deduplication, and ComBat batch correction to utilize metadata values.
- Build and export the standardized metadata schema matching the universal specification.

#### [MODIFY] [utils_rna.R](file:///c:/Users/ABDULLAH/Desktop/omicsflow/modules/rna/utils_rna.R)
- Update `rna_validate_barcodes` to bypass TCGA-specific length/prefix checks when metadata is supplied.

---

### Methylation Module

#### [MODIFY] [preprocess_meth.R](file:///c:/Users/ABDULLAH/Desktop/omicsflow/modules/methylation/preprocess_meth.R)
- Add `--metadata` argument.
- Refactor primary tumor filtering, deduplication, and ComBat plate ID mapping to leverage metadata values.
- Extract biological covariate `tss` from `center` when metadata is supplied, falling back to TCGA TSS parsing.
- Build and export the standardized metadata schema.

#### [MODIFY] [utils_meth.R](file:///c:/Users/ABDULLAH/Desktop/omicsflow/modules/methylation/utils_meth.R)
- Update `meth_validate_barcodes` to bypass TCGA-specific validation when metadata is supplied.

---

### CNV Module

#### [MODIFY] [preprocess_cnv.R](file:///c:/Users/ABDULLAH/Desktop/omicsflow/modules/cnv/preprocess_cnv.R)
- Add `--metadata` argument.
- Refactor CNV sample identification, tumor selection, patient mapping, and deduplication to use metadata lookup.

#### [MODIFY] [export_cnv.R](file:///c:/Users/ABDULLAH/Desktop/omicsflow/modules/cnv/export_cnv.R)
- Refactor `export_cnv_results` to accept an optional `sample_info` data frame conforming to the standardized metadata schema.

---

### SNV Module

#### [MODIFY] [preprocess_snv.R](file:///c:/Users/ABDULLAH/Desktop/omicsflow/modules/snv/preprocess_snv.R)
- Add `--metadata` argument.
- Refactor barcode standardization in step 4 to use the universal metadata helper.

#### [MODIFY] [export_snv.R](file:///c:/Users/ABDULLAH/Desktop/omicsflow/modules/snv/export_snv.R)
- Refactor `export_snv_results` to accept an optional `sample_info` data frame conforming to the standardized metadata schema.

---

### Nextflow Orchestration Layer

#### [MODIFY] [nextflow.config](file:///c:/Users/ABDULLAH/Desktop/omicsflow/nextflow.config)
- Add `params.metadata = null` default value.

#### [MODIFY] [main.nf](file:///c:/Users/ABDULLAH/Desktop/omicsflow/main.nf)
- Update input validation in `validateInputs()` to verify metadata file existence.
- Print `Metadata` in the startup banner.
- Update module triggers to create a `ch_metadata` channel and pass it to `PREPROCESS_RNA`, `PREPROCESS_METH`, `PREPROCESS_CNV`, and `PREPROCESS_SNV`.

#### [MODIFY] [preprocess_rna.nf](file:///c:/Users/ABDULLAH/Desktop/omicsflow/modules_nf/preprocess_rna.nf)
- Add `metadata_file` to process input channel and pass `--metadata` to the underlying R execution.

#### [MODIFY] [preprocess_meth.nf](file:///c:/Users/ABDULLAH/Desktop/omicsflow/modules_nf/preprocess_meth.nf)
- Add `metadata_file` to process input channel and pass `--metadata` to the R execution.

#### [MODIFY] [preprocess_cnv.nf](file:///c:/Users/ABDULLAH/Desktop/omicsflow/modules_nf/preprocess_cnv.nf)
- Add `metadata_file` to process input channel and pass `--metadata` to the R execution.

#### [MODIFY] [preprocess_snv.nf](file:///c:/Users/ABDULLAH/Desktop/omicsflow/modules_nf/preprocess_snv.nf)
- Add `metadata_file` to process input channel and pass `--metadata` to the R execution.

---

### Verification and Test Suite

#### [NEW] [test_metadata_mode.R](file:///c:/Users/ABDULLAH/Desktop/omicsflow/tests/test_metadata_mode.R)
- Create a script that programmatically builds a valid `sample_metadata.csv` from existing TCGA input RDS barcodes.
- Run preprocess pipelines for all four modules (RNA, Methylation, CNV, SNV) in both TCGA mode and Metadata mode.
- Verify that matrix dimensions, ranges, and processed outputs are structurally identical or scientifically equivalent.

#### [MODIFY] [test_rna.R](file:///c:/Users/ABDULLAH/Desktop/omicsflow/tests/test_rna.R)
- Add `--metadata` argument.
- Adjust validation constraints (length & TCGA prefix checks) if metadata is provided.

#### [MODIFY] [test_meth.R](file:///c:/Users/ABDULLAH/Desktop/omicsflow/tests/test_meth.R)
- Add `--metadata` argument.
- Adjust validation constraints if metadata is provided.

---

## Verification Plan

### Automated Tests
1. Run all standard validation tests in legacy TCGA fallback mode to ensure absolute compatibility:
   ```bash
   Rscript tests/test_rna.R --outdir results/rna/
   Rscript tests/test_meth.R --outdir results/methylation/
   Rscript tests/test_cnv.R --outdir results/cnv/
   Rscript tests/test_snv.R --outdir results/snv/
   ```
2. Run the new universal metadata validation script to verify equivalence of outputs under metadata-driven mode:
   ```bash
   Rscript tests/test_metadata_mode.R
   ```
3. Run Nextflow pipeline verification with and without `--metadata` parameter:
   ```bash
   nextflow run main.nf --rna data/rna_expression_raw.rds --meth data/methylation_beta_raw.rds ...
   ```
