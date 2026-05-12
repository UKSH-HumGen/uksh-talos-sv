process TransferAnnotationsToMatrixTable {
    container params.container

    input:
        path annotations
        tuple path(vcf), path(tbi)

    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        path "${params.cohort}.mt"
        path "${params.cohort}_merged_with_samples.vcf.bgz"
        path "${params.cohort}_merged_with_samples.vcf.bgz.tbi"

    script:
        """
        set -ex

        TransferAnnotationsToMatrixTable \
            --input ${vcf} \
            --annotations ${annotations} \
            --output ${params.cohort}.mt

        # Copy the original merged VCF with sample data for use by Talos
        cp ${vcf} ${params.cohort}_merged_with_samples.vcf.bgz
        cp ${tbi} ${params.cohort}_merged_with_samples.vcf.bgz.tbi
        """
}
