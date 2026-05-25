// =============================================================================
// preprocess_meth.nf — Nextflow Process: DNA Methylation Preprocessing
// OmicsFlow Orchestration Layer
// =============================================================================

process PREPROCESS_METH {
    tag "Methylation"
    label 'process_high_memory'

    publishDir "${params.outdir}/methylation", mode: 'copy', overwrite: true

    input:
    path meth_input
    path clinical_file
    path cross_react_file
    path metadata_file

    output:
    path "output_meth/methylation_processed_matrix.rds", emit: processed_matrix
    path "output_meth/*",                                 emit: all_outputs

    script:
    def clinical_arg   = clinical_file.name   != 'NO_FILE_CLINICAL'   ? "--clinical ${clinical_file}"       : ''
    def cross_arg      = cross_react_file.name != 'NO_FILE_CROSSREACT' ? "--cross-react ${cross_react_file}" : ''
    def metadata_arg   = metadata_file.name   != 'NO_FILE_METADATA'   ? "--metadata ${metadata_file}"       : ''
    def clinical_map_arg = params.clinical_map ? "--clinical_map ${params.clinical_map}" : ''
    """
    Rscript ${projectDir}/modules/methylation/preprocess_meth.R \\
        --input   ${meth_input} \\
        --outdir  output_meth/ \\
        --n-top   ${params.meth_ntop} \\
        ${clinical_arg} \\
        ${cross_arg} \\
        ${metadata_arg} \\
        ${clinical_map_arg}
    """
}
