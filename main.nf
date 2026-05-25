#!/usr/bin/env nextflow

// =============================================================================
// main.nf — OmicsFlow: Multi-Omics Integration Pipeline
// =============================================================================
//
// USAGE:
//   nextflow run main.nf \
//     --rna   data/rna_expression_raw.rds \
//     --meth  data/methylation_beta_raw.rds \
//     --cnv   data/cnv_segment_raw.rds \
//     --snv   data/snv_mutation_raw.rds \
//     --clinical data/clinical_data.tsv
//
// RESUME:
//   nextflow run main.nf -resume
//
// =============================================================================

nextflow.enable.dsl = 2

// =============================================================================
// IMPORT PROCESS MODULES
// =============================================================================

include { PREPROCESS_RNA  } from './modules_nf/preprocess_rna'
include { PREPROCESS_METH } from './modules_nf/preprocess_meth'
include { PREPROCESS_CNV  } from './modules_nf/preprocess_cnv'
include { PREPROCESS_SNV  } from './modules_nf/preprocess_snv'
include { RUN_INTEGRATION } from './modules_nf/run_integration'
include { RUN_ML          } from './modules_nf/run_ml'
include { RUN_ENRICHMENT  } from './modules_nf/run_enrichment'
include { GENERATE_REPORT } from './modules_nf/generate_report'

// =============================================================================
// PARAMETER VALIDATION
// =============================================================================

def validateInputs() {

    // --- Count supplied omics layers ---
    def omics_count = 0
    def omics_names = []

    if (params.rna)  { omics_count++; omics_names << 'RNA'         }
    if (params.meth) { omics_count++; omics_names << 'Methylation' }
    if (params.cnv)  { omics_count++; omics_names << 'CNV'         }
    if (params.snv)  { omics_count++; omics_names << 'SNV'         }

    // --- Minimum omics validation ---
    if (omics_count < 2) {
        error """
        ╔══════════════════════════════════════════════════════════════╗
        ║  INPUT VALIDATION FAILED                                    ║
        ╠══════════════════════════════════════════════════════════════╣
        ║  At least 2 omics layers are required.                      ║
        ║  Provided: ${omics_count} (${omics_names.join(', ') ?: 'none'})
        ║                                                              ║
        ║  Supply inputs via:                                          ║
        ║    --rna   <path>   RNA-seq expression RDS                   ║
        ║    --meth  <path>   Methylation beta-value RDS               ║
        ║    --cnv   <path>   CNV segment data RDS                     ║
        ║    --snv   <path>   SNV mutation data RDS                    ║
        ╚══════════════════════════════════════════════════════════════╝
        """.stripIndent()
    }

    // --- File existence validation ---
    if (params.rna  && !file(params.rna).exists())  { error "RNA input not found: ${params.rna}"   }
    if (params.meth && !file(params.meth).exists()) { error "Meth input not found: ${params.meth}" }
    if (params.cnv  && !file(params.cnv).exists())  { error "CNV input not found: ${params.cnv}"   }
    if (params.snv  && !file(params.snv).exists())  { error "SNV input not found: ${params.snv}"   }

    if (params.clinical && !file(params.clinical).exists()) {
        error "Clinical data not found: ${params.clinical}"
    }

    if (params.metadata && !file(params.metadata).exists()) {
        error "Metadata file not found: ${params.metadata}"
    }

    if (params.clinical_map && !file(params.clinical_map).exists()) {
        error "Clinical map file not found: ${params.clinical_map}"
    }

    if (params.validation_keywords && (params.validation_keywords.endsWith(".json") || params.validation_keywords.endsWith(".txt")) && !file(params.validation_keywords).exists()) {
        error "Validation keywords file not found: ${params.validation_keywords}"
    }

    def csv_pattern = ~/.*\.csv$/
    if (params.metadata && !(params.metadata =~ csv_pattern)) {
        error "Metadata input must be .csv format: ${params.metadata}"
    }

    // --- Format validation ---
    def rds_pattern = ~/.*\.rds$/
    if (params.rna  && !(params.rna  =~ rds_pattern)) { error "RNA input must be .rds format: ${params.rna}"   }
    if (params.meth && !(params.meth =~ rds_pattern)) { error "Meth input must be .rds format: ${params.meth}" }
    if (params.cnv  && !(params.cnv  =~ rds_pattern)) { error "CNV input must be .rds format: ${params.cnv}"   }
    if (params.snv  && !(params.snv  =~ rds_pattern)) { error "SNV input must be .rds format: ${params.snv}"   }

    // --- Full pipeline requirement ---
    def full_pipeline = (omics_count == 4)
    if (!full_pipeline) {
        log.warn """
        ┌──────────────────────────────────────────────────────────────┐
        │  PARTIAL MODE: ${omics_count}/4 omics provided                          │
        │  Only preprocessing will run for supplied layers.            │
        │  Integration / ML / Enrichment require all 4 omics.          │
        └──────────────────────────────────────────────────────────────┘
        """.stripIndent()
    }

    return [omics_count: omics_count, omics_names: omics_names, full_pipeline: full_pipeline]
}

// =============================================================================
// STARTUP BANNER
// =============================================================================

def printBanner() {
    log.info """
    ╔══════════════════════════════════════════════════════════════╗
    ║                                                              ║
    ║    ██████  ███    ███ ██  ██████ ███████                     ║
    ║   ██    ██ ████  ████ ██ ██      ██                          ║
    ║   ██    ██ ██ ████ ██ ██ ██      ███████                     ║
    ║   ██    ██ ██  ██  ██ ██ ██           ██                     ║
    ║    ██████  ██      ██ ██  ██████ ███████                     ║
    ║                                                              ║
    ║   ███████ ██       ██████  ██     ██                         ║
    ║   ██      ██      ██    ██ ██     ██                         ║
    ║   █████   ██      ██    ██ ██  █  ██                         ║
    ║   ██      ██      ██    ██ ██ ███ ██                         ║
    ║   ██      ███████  ██████   ███ ███                          ║
    ║                                                              ║
    ║   OmicsFlow v1.0.0 — Stable Core Release                    ║
    ║   Author: Abdullah Ibrahim Ali                               ║
    ║                                                              ║
    ╠══════════════════════════════════════════════════════════════╣
    ║  RNA          : ${params.rna  ?: 'not provided'}
    ║  Methylation  : ${params.meth ?: 'not provided'}
    ║  CNV          : ${params.cnv  ?: 'not provided'}
    ║  SNV          : ${params.snv  ?: 'not provided'}
    ║  Clinical     : ${params.clinical ?: 'not provided'}
    ║  Metadata     : ${params.metadata ?: 'not provided'}
    ║  Output       : ${params.outdir}
    ╚══════════════════════════════════════════════════════════════╝
    """.stripIndent()
}

// =============================================================================
// MAIN WORKFLOW
// =============================================================================

workflow {

    // --- Banner & Validation ---
    printBanner()
    def validation = validateInputs()

    // =========================================================================
    // PHASE 1: PREPROCESSING (parallel execution for all supplied omics)
    // =========================================================================

    // --- Metadata Channel ---
    ch_metadata = params.metadata
        ? Channel.fromPath(params.metadata, checkIfExists: true)
        : Channel.fromPath("${projectDir}/assets/NO_FILE_METADATA", checkIfExists: false).ifEmpty(file("NO_FILE_METADATA"))

    // --- RNA ---
    if (params.rna) {
        ch_rna_input = Channel.fromPath(params.rna, checkIfExists: true)
        PREPROCESS_RNA(ch_rna_input, ch_metadata)
    }

    // --- Methylation ---
    if (params.meth) {
        ch_meth_input = Channel.fromPath(params.meth, checkIfExists: true)

        // Optional auxiliary files (use placeholder if not provided)
        ch_clinical = params.clinical
            ? Channel.fromPath(params.clinical, checkIfExists: true)
            : Channel.fromPath("${projectDir}/assets/NO_FILE_CLINICAL", checkIfExists: false).ifEmpty(file("NO_FILE_CLINICAL"))

        ch_cross_react = params.cross_react
            ? Channel.fromPath(params.cross_react, checkIfExists: true)
            : Channel.fromPath("${projectDir}/assets/NO_FILE_CROSSREACT", checkIfExists: false).ifEmpty(file("NO_FILE_CROSSREACT"))

        PREPROCESS_METH(ch_meth_input, ch_clinical, ch_cross_react, ch_metadata)
    }

    // --- CNV ---
    if (params.cnv) {
        ch_cnv_input = Channel.fromPath(params.cnv, checkIfExists: true)

        ch_gene_coords = params.gene_coords
            ? Channel.fromPath(params.gene_coords, checkIfExists: true)
            : Channel.fromPath("${projectDir}/assets/NO_FILE_CACHE", checkIfExists: false).ifEmpty(file("NO_FILE_CACHE"))

        PREPROCESS_CNV(ch_cnv_input, ch_gene_coords, ch_metadata)
    }

    // --- SNV ---
    if (params.snv) {
        ch_snv_input = Channel.fromPath(params.snv, checkIfExists: true)
        PREPROCESS_SNV(ch_snv_input, ch_metadata)
    }

    // =========================================================================
    // PHASE 2: INTEGRATION → ML → ENRICHMENT (sequential, requires all 4)
    // =========================================================================

    if (validation.full_pipeline) {

        // --- MOFA+ Integration ---
        // Waits for all 4 preprocessing outputs automatically via channels
        RUN_INTEGRATION(
            PREPROCESS_RNA.out.processed_matrix,
            PREPROCESS_METH.out.processed_matrix,
            PREPROCESS_CNV.out.processed_matrix,
            PREPROCESS_SNV.out.processed_matrix
        )

        // --- ML Survival Analysis ---
        // Requires: MOFA model + RNA ML matrix + Clinical data
        ch_clinical_ml = Channel.fromPath(params.clinical, checkIfExists: true)

        RUN_ML(
            RUN_INTEGRATION.out.mofa_model,
            PREPROCESS_RNA.out.ml_matrix,
            ch_clinical_ml
        )

        // --- Pathway Enrichment ---
        // Requires: ML downstream outputs
        RUN_ENRICHMENT(
            RUN_ML.out.mofa_top_genes,
            RUN_ML.out.rf_top_genes,
            RUN_ML.out.lasso_genes,
            RUN_ML.out.rna_for_pathway
        )

        // --- Automated Report Generation ---
        // Collects all result directories and generates the final HTML report
        ch_results = RUN_ENRICHMENT.out.all_outputs
            .map { it -> file(params.outdir) }

        GENERATE_REPORT(ch_results)
    }
}

// =============================================================================
// COMPLETION HANDLER
// =============================================================================

workflow.onComplete {
    def duration = workflow.duration
    def status   = workflow.success ? 'COMPLETED SUCCESSFULLY' : 'FAILED'

    log.info """
    ╔══════════════════════════════════════════════════════════════╗
    ║  OmicsFlow Pipeline ${status}
    ╠══════════════════════════════════════════════════════════════╣
    ║  Duration    : ${duration}
    ║  Completed   : ${workflow.complete}
    ║  Exit status : ${workflow.exitStatus}
    ║  Work dir    : ${workflow.workDir}
    ║  Results     : ${params.outdir}
    ╚══════════════════════════════════════════════════════════════╝
    """.stripIndent()
}

workflow.onError {
    log.error """
    ╔══════════════════════════════════════════════════════════════╗
    ║  OmicsFlow Pipeline ERROR                                    ║
    ╠══════════════════════════════════════════════════════════════╣
    ║  Error message: ${workflow.errorMessage}
    ║  Exit status  : ${workflow.exitStatus}
    ╚══════════════════════════════════════════════════════════════╝
    """.stripIndent()
}
