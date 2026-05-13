process MergeStructuralVariants {
    container params.container

    publishDir params.output_dir, mode: 'copy'

    input:
        path talos_result_json
        tuple path(sv_vcf), path(sv_vcf_idx)
        path pedigree
        path panelapp_data
        path gene_to_phenotype
        path talos_config

    output:
        path "${params.cohort}_results_with_sv_*.json"

    script:
    def timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm')
    """
    export TALOS_CONFIG=${talos_config}
    export GFF3_FILE=${projectDir}/../large_files/Homo_sapiens.GRCh38.113.gff3.gz

    echo "=== Starting SV Merge with Phenotype Matching ==="

    # Convert VCF to uncompressed format for Python script
    bcftools view ${sv_vcf} > structural_variants.vcf

    # Extract all sample IDs from pedigree file (no header row)
    SAMPLE_IDS=\$(cut -f2 ${pedigree})
    echo "Samples in pedigree: \$SAMPLE_IDS"

    # Start with the base Talos JSON
    cp ${talos_result_json} intermediate.json

    # Process each sample separately to ensure SVs are assigned to correct samples
    for SAMPLE_ID in \$SAMPLE_IDS; do
        echo ""
        echo "Processing sample: \$SAMPLE_ID"

        # Run enhanced Python script with config-driven phenotype matching
        # This script will:
        # 1. Parse VCF header to get sample names
        # 2. Filter SVs based on genotype (only include variants this sample has)
        # 3. Annotate SVs with overlapping genes from GFF3
        # 4. Perform HPO phenotype matching between patient and gene
        # 5. Calculate phenotype similarity scores
        # 6. Apply filters from config.toml [StructuralVariantFiltering] section

        python3 ${projectDir}/../scripts/merge_sv_with_phenotype_matching.py \\
            structural_variants.vcf \\
            intermediate.json \\
            output_temp.json \\
            \$SAMPLE_ID \\
            --pedigree ${pedigree} \\
            --gen2phen ${gene_to_phenotype} \\
            --config ${talos_config}

        # Update intermediate JSON for next sample
        mv output_temp.json intermediate.json
    done

    # Rename final output
    mv intermediate.json ${params.cohort}_results_with_sv_${timestamp}.json

    echo ""
    echo "Structural variant integration with phenotype matching completed"
    echo "Output: ${params.cohort}_results_with_sv_${timestamp}.json"

    # Clean up temporary files
    rm -f structural_variants.vcf
    """
}
