
process RescueCompoundHet {
    container params.container

    publishDir params.output_dir, mode: 'copy'

    input:
        path results_with_sv_json
        path full_report_json
        path panelapp_json
        path rescue_script

    def timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm')

    output:
        path "${params.cohort}_results_final_${timestamp}.json"

    script:
    """
    echo "=== Rescuing Compound Heterozygous Variants (SNV+SV) ===" >&2
    echo "Input: ${results_with_sv_json}" >&2
    echo "Full report: ${full_report_json}" >&2
    echo "PanelApp: ${panelapp_json}" >&2

    python3 ${rescue_script} \\
        ${results_with_sv_json} \\
        ${full_report_json} \\
        ${panelapp_json} \\
        ${params.cohort}_results_final_${timestamp}.json

    echo "=== Compound Het Rescue Complete ===" >&2
    echo "Output: ${params.cohort}_results_final_${timestamp}.json" >&2
    """
}
