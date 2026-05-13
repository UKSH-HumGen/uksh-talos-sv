process RunHailFilteringSv {
    container params.container

    // runs the structural variant categorization
    publishDir params.output_dir, mode: 'copy'

    input:
        path merged_vcf
        path panelapp_data
        path pedigree
        path mane_json
        path talos_config
        path ensembl_gff

    output:
        tuple \
            path("${params.cohort}_structural_variants_labelled.vcf.bgz"), \
            path("${params.cohort}_structural_variants_labelled.vcf.bgz.tbi")

    """
    export TALOS_CONFIG=${talos_config}

    # Extract structural variants from merged VCF (DRAGEN format with LOSS/GAIN IDs)
    echo "Extracting structural variants from merged VCF..."
    bcftools view \
        --include 'ID~"DRAGEN:LOSS" || ID~"DRAGEN:GAIN" || ALT~"<DEL>" || ALT~"<DUP>" || ALT~"<INV>" || ALT~"<BND>" || ALT~"<INS>"' \
        -O z \
        -o structural_variants.vcf.bgz \
        ${merged_vcf}

    # Index the SV VCF
    bcftools index -t structural_variants.vcf.bgz

    # Check if structural variants file has any variants (more than just header)
    variant_count=\$(bcftools view --no-header structural_variants.vcf.bgz | wc -l)

    if [ \$variant_count -gt 0 ]; then
        echo "Found \$variant_count structural variants"
        echo "Categorizing structural variants based on gene/exon overlap..."

        # Set basic IDs first
        bcftools annotate \
            --set-id 'SV_%CHROM:%POS' \
            -O z \
            -o structural_variants_with_ids.vcf.bgz \
            structural_variants.vcf.bgz

        bcftools index -t structural_variants_with_ids.vcf.bgz

        # Categorize SVs using Python script
        # This adds gene_id, categorybooleanhighimpact, and other fields
        CategorizeSVs \
            --input structural_variants_with_ids.vcf.bgz \
            --output ${params.cohort}_structural_variants_labelled.vcf.bgz \
            --gff3 ${ensembl_gff} \
            --panelapp ${panelapp_data}

        # Index the categorized output
        bcftools index -t ${params.cohort}_structural_variants_labelled.vcf.bgz

        echo "Successfully categorized \$variant_count structural variants"
    else
        echo "No structural variants found, creating empty output file"
        # Create empty output
        bcftools view -O z -o ${params.cohort}_structural_variants_labelled.vcf.bgz structural_variants.vcf.bgz
        bcftools index -t ${params.cohort}_structural_variants_labelled.vcf.bgz
    fi
    """
}
