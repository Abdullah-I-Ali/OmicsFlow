// =============================================================================
// generate_report.nf — Nextflow Process: Quarto Report Generation
// OmicsFlow Orchestration Layer
// =============================================================================

process GENERATE_REPORT {
    tag "Report"
    label 'process_medium'

    publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true

    input:
    path results_dir

    output:
    path "OmicsFlow_Report.html", emit: report_html

    script:
    """
    # Copy report template and stylesheet into working directory
    cp ${projectDir}/reports/OmicsFlow_Report.qmd .
    cp ${projectDir}/reports/styles.css .

    # Render the Quarto report with results path parameter
    quarto render OmicsFlow_Report.qmd \
        --to html \
        -P results_dir:${results_dir}
    """
}
