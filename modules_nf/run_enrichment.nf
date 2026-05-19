// =============================================================================
// run_enrichment.nf — Nextflow Process: Pathway Enrichment Analysis
// OmicsFlow Orchestration Layer
// =============================================================================

process RUN_ENRICHMENT {
    tag "Enrichment"
    label 'process_medium'

    publishDir "${params.outdir}/enrichment", mode: 'copy', overwrite: true

    input:
    path mofa_top_genes
    path rf_top_genes
    path lasso_genes
    path rna_for_pathway

    output:
    path "output_enrichment/enrichment_results.rds", emit: enrichment_results
    path "output_enrichment/*",                       emit: all_outputs

    script:
    """
    Rscript ${projectDir}/modules/enrichment/run_enrichment.R \\
        --mofa    ${mofa_top_genes} \\
        --rf      ${rf_top_genes} \\
        --lasso   ${lasso_genes} \\
        --rna     ${rna_for_pathway} \\
        --outdir  output_enrichment/
    """
}
