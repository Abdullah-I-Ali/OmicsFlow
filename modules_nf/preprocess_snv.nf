// =============================================================================
// preprocess_snv.nf — Nextflow Process: SNV Mutation Preprocessing
// OmicsFlow Orchestration Layer
// =============================================================================

process PREPROCESS_SNV {
    tag "SNV"
    label 'process_medium'

    publishDir "${params.outdir}/snv", mode: 'copy', overwrite: true

    input:
    path snv_input
    path metadata_file

    output:
    path "output_snv/snv_processed_matrix.rds", emit: processed_matrix
    path "output_snv/*",                         emit: all_outputs

    script:
    def metadata_arg = metadata_file.name != 'NO_FILE_METADATA' ? "--metadata ${metadata_file}" : ''
    """
    Rscript ${projectDir}/modules/snv/preprocess_snv.R \\
        --input   ${snv_input} \\
        --outdir  output_snv/ \\
        --freq    ${params.snv_freq} \\
        --topn    ${params.snv_topn} \\
        ${metadata_arg}
    """
}
