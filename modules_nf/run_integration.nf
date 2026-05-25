// =============================================================================
// run_integration.nf — Nextflow Process: MOFA+ Multi-Omics Integration
// OmicsFlow Orchestration Layer
// =============================================================================

process RUN_INTEGRATION {
    tag "MOFA+"
    label 'process_high'

    publishDir "${params.outdir}/integration", mode: 'copy', overwrite: true

    input:
    path rna_matrix
    path meth_matrix
    path cnv_matrix
    path snv_matrix

    output:
    path "output_integration/mofa_model.rds", emit: mofa_model
    path "output_integration/*",               emit: all_outputs

    script:
    def metadata_arg = params.metadata ? "--metadata ${params.metadata}" : ""
    """
    Rscript ${projectDir}/modules/integration/run_integration.R \\
        --rna     ${rna_matrix} \\
        --meth    ${meth_matrix} \\
        --cnv     ${cnv_matrix} \\
        --snv     ${snv_matrix} \\
        --outdir  output_integration/ \\
        --factors ${params.mofa_factors} \\
        --iter    ${params.mofa_iter} \\
        ${metadata_arg}
    """
}
