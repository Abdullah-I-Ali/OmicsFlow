// =============================================================================
// preprocess_cnv.nf — Nextflow Process: CNV Segment Preprocessing
// OmicsFlow Orchestration Layer
// =============================================================================

process PREPROCESS_CNV {
    tag "CNV"
    label 'process_medium'

    publishDir "${params.outdir}/cnv", mode: 'copy', overwrite: true

    input:
    path cnv_input
    path gene_coords_cache
    path metadata_file

    output:
    path "output_cnv/cnv_processed_matrix.rds", emit: processed_matrix
    path "output_cnv/*",                         emit: all_outputs

    script:
    def cache_arg    = gene_coords_cache.name != 'NO_FILE_CACHE' ? "--cache ${gene_coords_cache}" : ''
    def metadata_arg = metadata_file.name != 'NO_FILE_METADATA' ? "--metadata ${metadata_file}" : ''
    """
    Rscript ${projectDir}/modules/cnv/preprocess_cnv.R \\
        --input   ${cnv_input} \\
        --outdir  output_cnv/ \\
        --ntop    ${params.cnv_ntop} \\
        ${cache_arg} \\
        ${metadata_arg}
    """
}
