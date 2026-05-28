# Realistic Oncology Validation Cohort Generation & OmicsFlow Pipeline Verification

## Objective
Generate a realistic multi-omics oncology validation cohort and run it through the full OmicsFlow pipeline to evaluate its correctness, robustness, and ability to detect biological signals under realistic conditions such as missing data, confounders, and batch effects.

## What Was Accompliished

### 1. Realistic Cohort Generation (`generate_realistic_cohort.R`)
Generated a synthetic dataset of 180 patients that simulates a heterogeneous oncology cohort. Instead of purely random noise, we injected realistic biological signals linked to three distinct subtypes: **Proliferative**, **Mesenchymal**, and **Immune_enriched**.

**Clinical & Survival Data:**
- 74.4% overall event rate with distinct survival distributions.
- **Immune_enriched** patients had the most favorable prognosis (median OS 958 days).
- **Proliferative** patients had the worst prognosis (median OS 278 days).

**Data Layers Generated:**
- **RNA Expression**: 2000 genes mapped to genuine Ensembl IDs. Subtype-specific markers were exclusively drawn from biologically relevant GO pathways (Cell Cycle `GO:0007049`, ECM `GO:0030198`, and Immune Response `GO:0006955`).
- **DNA Methylation**: 5000 valid `cg` probes matching the true Illumina 450k annotation matrix, displaying bimodal (hypo/hyper-methylated) density distributions and realistic subtype-specific patterns.
- **Copy Number Variation (CNV)**: 7685 genomic segments mapped precisely to human genome build hg38. 
- **Single Nucleotide Variants (SNV)**: 3271 mutations targeting realistic drivers (e.g., *TP53*, *MYC*, *PIK3CA*).

**Confounding Factors Included:**
- Two simulated study centers (Alpha, Beta)
- Three experimental batches with baseline additive batch effects. 

### 2. End-to-End Pipeline Orchestration (`run_realistic_validation.R`)
We orchestrated a full, automated execution of the 7-module OmicsFlow pipeline.

- **RNA Preprocessing**: Successfully mapped symbols, filtered low expression, standardized via Z-score, and mitigated batch effects.
- **Methylation Preprocessing**: Matched probes effectively to the true 450k annotation package, eliminated cross-reactive probes, imputed missing data, and successfully corrected for batch effects using `ComBat`.
- **CNV Preprocessing**: Standardized segments, built `GRanges` objects, retrieved genomic coordinates from our generated cache, and correctly merged them into a clean sample × gene matrix.
- **SNV Preprocessing**: Parsed MAF files, aggregated at the patient level, and removed hypermutated subjects.
- **MOFA+ Integration**: Merged all 4 data modalities across 176 common patients, discovering 5 major latent active factors.
- **Machine Learning**: Extracted top MOFA+ features and trained Random Forest (C-index: 0.820), XGBoost (C-index: 0.830), and LASSO Cox survival models.
- **Pathway Enrichment**: Conducted clusterProfiler `enrichGO` against the selected features.

## Validation Results

**Final Verdict**: `PASS` :white_check_mark:

The system successfully recovered the injected biological signals and successfully completed all steps of the pipeline.

```text
  ✓ All 7 modules passed 
  ✓ ≥100 common patients (176)
  ✓ ≥1 active MOFA factor (5)
  ✓ ≥1 ML C-index > 0.5 (Max: 0.830)
  ✓ ≥1 enriched pathway (24 total GO terms)
```

**Significant Enriched Pathways Found**:
The Random Forest/MOFA analysis recovered 23 GO Biological Processes and 1 GO Molecular Function linked to the exact marker genes we injected (e.g. immune response, extracellular matrix organization, cell cycle processes).

## Visual Artifacts
The pipeline automatically produced a highly detailed `validation_report.txt` and various multi-omics visualizations in `results/realistic_validation/`, including:
- Sample correlation heatmaps
- Density charts and PCA plots per modality
- MOFA+ variance explained plots
- Final Pathway Enrichment Network (`cnetplot`, `emapplot`)
