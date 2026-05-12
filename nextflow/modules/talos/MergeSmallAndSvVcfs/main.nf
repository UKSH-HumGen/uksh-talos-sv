process MergeSmallAndSvVcfs {
    container params.container

    // Merge small variants and structural variants VCFs for unified MOI processing
    publishDir params.output_dir, mode: 'copy'

    input:
        tuple path(small_vcf), path(small_vcf_tbi)
        tuple path(sv_vcf), path(sv_vcf_tbi)

    output:
        tuple \
            path("${params.cohort}_all_variants_labelled.vcf.bgz"), \
            path("${params.cohort}_all_variants_labelled.vcf.bgz.tbi")

    script:
    """
    echo "=== Merging Small Variants and Structural Variants ==="

    # Check if SV VCF has any variants
    sv_count=\$(bcftools view --no-header ${sv_vcf} | wc -l)
    small_count=\$(bcftools view --no-header ${small_vcf} | wc -l)

    echo "Small variants: \$small_count"
    echo "Structural variants: \$sv_count"

    if [ \$sv_count -gt 0 ]; then
        echo "Merging both small variants and SVs..."

        # Check if sample counts match
        small_samples=\$(bcftools query -l ${small_vcf} | wc -l)
        sv_samples=\$(bcftools query -l ${sv_vcf} | wc -l)

        echo "Small VCF samples: \$small_samples"
        echo "SV VCF samples: \$sv_samples"

        if [ \$small_samples -ne \$sv_samples ]; then
            echo "Sample counts differ - subsetting SV VCF to match small variants VCF samples"

            # Get samples from small VCF
            bcftools query -l ${small_vcf} > samples_to_keep.txt

            # Subset SV VCF to match
            bcftools view -S samples_to_keep.txt -O z -o ${sv_vcf}.subset.vcf.bgz ${sv_vcf}
            bcftools index -t ${sv_vcf}.subset.vcf.bgz

            sv_vcf_to_use="${sv_vcf}.subset.vcf.bgz"
            echo "Subsetted SV VCF to \$(bcftools query -l \$sv_vcf_to_use | wc -l) samples"
        else
            sv_vcf_to_use="${sv_vcf}"
        fi

        # Fix missing FILTER definitions in small variants VCF header
        echo "Adding missing FILTER definitions to small variants VCF..."
        cat > filter_defs.txt <<EOF
##FILTER=<ID=DRAGENSnpHardQUAL,Description="DRAGEN SNP hard quality filter">
##FILTER=<ID=DRAGENIndelHardQUAL,Description="DRAGEN Indel hard quality filter">
##FILTER=<ID=LowDepth,Description="Low depth filter">
EOF

        bcftools annotate -h filter_defs.txt -O z -o ${small_vcf}.fixed.vcf.bgz ${small_vcf}
        bcftools index -t ${small_vcf}.fixed.vcf.bgz
        small_vcf_to_use="${small_vcf}.fixed.vcf.bgz"

        # Concatenate small variants and SVs, then sort
        # Use --allow-overlaps since SVs can overlap with SNVs
        # Sort to ensure proper ordering
        bcftools concat \
            --allow-overlaps \
            -O u \
            \$small_vcf_to_use \$sv_vcf_to_use | \
        bcftools sort \
            -O z \
            -o ${params.cohort}_all_variants_labelled.vcf.bgz

        bcftools index -t ${params.cohort}_all_variants_labelled.vcf.bgz

        merged_count=\$(bcftools view --no-header ${params.cohort}_all_variants_labelled.vcf.bgz | wc -l)
        echo "Merged total: \$merged_count variants"
        echo "✓ Successfully merged small variants and SVs for unified MOI processing"
    else
        echo "No SVs found, using small variants only..."

        # Just copy the small variants VCF
        cp ${small_vcf} ${params.cohort}_all_variants_labelled.vcf.bgz
        cp ${small_vcf_tbi} ${params.cohort}_all_variants_labelled.vcf.bgz.tbi

        echo "✓ Using small variants only (\$small_count variants)"
    fi

    echo "=== Merge Complete ==="
    """
}
