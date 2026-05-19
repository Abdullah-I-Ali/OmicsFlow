// =============================================================================
// run_ml.nf — Nextflow Process: Machine Learning Survival Pipeline
// OmicsFlow Orchestration Layer
// =============================================================================

process RUN_ML {
    tag "ML"
    label 'process_high'

    publishDir "${params.outdir}/ml", mode: 'copy', overwrite: true

    input:
    path mofa_model
    path rna_ml_matrix
    path clinical_file

    output:
    path "output_ml/mofa_top_genes.rds",        emit: mofa_top_genes
    path "output_ml/rf_top_genes.rds",           emit: rf_top_genes
    path "output_ml/lasso_selected_genes.rds",   emit: lasso_genes
    path "output_ml/rna_for_pathway.rds",        emit: rna_for_pathway
    path "output_ml/*",                           emit: all_outputs

    script:
    """
    Rscript ${projectDir}/modules/ml/run_ml.R \\
        --mofa            ${mofa_model} \\
        --rna             ${rna_ml_matrix} \\
        --clinical        ${clinical_file} \\
        --outdir          output_ml/ \\
        --mofa_prefilter  ${params.ml_mofa_prefilter} \\
        --final_features  ${params.ml_final_features} \\
        --train_ratio     ${params.ml_train_ratio}
    """
}
