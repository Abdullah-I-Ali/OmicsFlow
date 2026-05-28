# CNV Module Failure Investigation Report

## 1. Trace & Exact Code Path
The failure occurs during the very first step of the CNV preprocessing execution (`STEP 1: Data Loading & Barcode Parsing`). 

The pipeline initializes by calling the data loading function:
- **Trigger**: `preprocess_cnv.R` calls `load_cnv_data(opt$input, opt$cache)` (Line 102).
- **Execution**: The `load_cnv_data` function loads the raw `.rds` object and immediately validates its column structure.
- **Failure**: The script checks the data frame against a hardcoded list of `required_cols`. Because the synthetic dataset does not contain `GDC_Aliquot`, the script halts execution via `stop()`.

## 2. Location Identification
- **File**: `modules/cnv/utils_cnv.R`
- **Function**: `load_cnv_data()`
- **Line Numbers**: 76–80

```r
# utils_cnv.R (Lines 76-80)
required_cols <- c("GDC_Aliquot", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean", "Sample")
missing_cols <- setdiff(required_cols, colnames(cnv_data))
if (length(missing_cols) > 0) {
  stop("Missing required columns in CNV data: ", paste(missing_cols, collapse = ", "))
}
```

## 3. Dependency Classification
**Verdict: True TCGA Dependency Leak**
This is a leaked assumption within the source code, not a flaw in the synthetic data generation. `GDC_Aliquot` is an abstraction specific to the NCI Genomic Data Commons (GDC) format used for TCGA cohorts. 
A codebase search confirms that `GDC_Aliquot` is never actually referenced or used downstream in `preprocess_cnv.R`; all downstream aggregation relies entirely on the `Sample` column. Therefore, strictly requiring it is a residual artifact from an older TCGA-only version of the pipeline.

## 4. Why Metadata Mode Did Not Bypass This
The column validation is executed unconditionally inside the `load_cnv_data` function, which acts purely as an I/O loader. The `metadata` object is not passed to this function, nor is there any conditional logic inside `load_cnv_data` to check if the pipeline is operating in legacy or custom metadata mode. The hardcoded TCGA requirement acts as a structural gatekeeper, crashing the pipeline before the dynamic metadata mappings are ever evaluated.

## 5. Proposed Minimal Code Change
The most robust and minimal fix is to safely remove `"GDC_Aliquot"` from the `required_cols` list inside `utils_cnv.R`.

**From:**
```r
required_cols <- c("GDC_Aliquot", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean", "Sample")
```

**To:**
```r
required_cols <- c("Chromosome", "Start", "End", "Num_Probes", "Segment_Mean", "Sample")
```

This single line modification will permit non-TCGA, metadata-driven datasets to successfully pass the ingestion step while preserving all structural requirements needed for downstream genomic mapping.
