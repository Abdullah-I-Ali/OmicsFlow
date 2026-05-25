// =============================================================================
// preprocess_rna.nf — Nextflow Process: RNA-seq Preprocessing
// OmicsFlow Orchestration Layer
// =============================================================================

process PREPROCESS_RNA {
    tag "RNA"
    label 'process_high'

    publishDir "${params.outdir}/rna", mode: 'copy', overwrite: true

    input:
    path rna_input
    path metadata_file

    output:
    path "output_rna/rna_processed_matrix.rds", emit: processed_matrix
    path "output_rna/rna_ml.rds",               emit: ml_matrix
    path "output_rna/*",                         emit: all_outputs

    script:
    def metadata_arg = metadata_file.name != 'NO_FILE_METADATA' ? "--metadata ${metadata_file}" : ''
    """
    Rscript ${projectDir}/modules/rna/preprocess_rna.R \\
        --input   ${rna_input} \\
        --outdir  output_rna/ \\
        --n-top   ${params.rna_ntop} \\
        --cor-low ${params.rna_cor_low} \\
        --cor-high ${params.rna_cor_high} \\
        ${metadata_arg}
    """
}
