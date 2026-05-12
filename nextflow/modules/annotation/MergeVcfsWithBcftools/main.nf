
process MergeVcfsWithBcftools {
    container params.container

    input:
        path vcfs
        path tbis
        path regions
        path ref_genome

    // merge the VCFs into a single VCF
    publishDir params.cohort_output_dir

    output:
        tuple \
            path("${params.cohort}_merged.vcf.bgz"), \
            path("${params.cohort}_merged.vcf.bgz.tbi")

    script:

        def input = (vcfs.collect().size() > 1) ? vcfs.sort{ it.name } : vcfs
        """
        set -e
        # https://github.com/samtools/bcftools/issues/1189
        # -m none means don't merge multi-allelic sites, keep everything atomic
        # "-m none" means don't merge multi-allelic sites, keep everything atomic, we're splitting in the next step
        # -0 to set all missing genotypes to HomWT - gap-filling with Missing (default) reduces the AN, so callset
        # frequency filters can appear to show an inflated AC/AN ratio

        echo "=== Starting VCF Merge with SV Preservation ===" >&2
        echo "Input VCFs: $input" >&2

        # Pre-filter: Remove DRAGEN reference blocks (REF=N with no ALT) to avoid merge conflicts
        echo "Step 0: Removing DRAGEN reference blocks from input VCFs..." >&2
        for vcf in $input; do
            base=\$(basename "\$vcf" .vcf.gz)
            bcftools view -e 'REF="N" && ALT="."' -Oz "\$vcf" > "\${base}.filtered.vcf.gz"
            bcftools index -t "\${base}.filtered.vcf.gz"
        done

        # First merge all variants
        echo "Step 1: Merging all variants..." >&2
        bcftools merge \
        	--force-single \
        	-m none \
        	-R ${regions} \
        	-0 \
        	-Ou \
        	--no-version \
        	*.filtered.vcf.gz \
        	-o merged_all.bcf

        # Split into SVs and non-SVs
        # SVs: anything with SVTYPE in INFO or symbolic alleles like <DEL>, <DUP>, <INS>, <INV>, <CNV>, <BND>
        echo "Step 2: Separating SVs from non-SVs..." >&2
        bcftools view -i 'INFO/SVTYPE!="." || ALT~"<"' -Oz merged_all.bcf -o svs.vcf.gz
        bcftools index -t svs.vcf.gz
        SV_COUNT=\$(bcftools view -H svs.vcf.gz | wc -l)
        echo "  Found \$SV_COUNT structural variants" >&2

        # Non-SVs: exclude anything with SVTYPE or symbolic alleles
        bcftools view -e 'INFO/SVTYPE!="." || ALT~"<"' -Oz merged_all.bcf -o non_svs.vcf.gz
        bcftools index -t non_svs.vcf.gz
        NON_SV_COUNT=\$(bcftools view -H non_svs.vcf.gz | wc -l)
        echo "  Found \$NON_SV_COUNT non-SV variants (SNVs/indels)" >&2

        # Normalize only non-SVs (SNVs and small indels)
        # This preserves SV coordinates while properly normalizing SNVs
        echo "Step 3: Normalizing non-SVs only (preserving SV coordinates)..." >&2
        bcftools norm \
        	-m -any \
        	-c ws \
        	-Oz \
        	--no-version \
        	-f ${ref_genome} \
        	non_svs.vcf.gz \
        	-o non_svs_norm.vcf.gz
        bcftools index -t non_svs_norm.vcf.gz

        # Decompose multiallelic SVs (preserves coordinates, just splits alleles)
        echo "Step 4a: Decomposing multiallelic SVs..." >&2
        bcftools norm -m -any -Oz svs.vcf.gz -o svs_decomposed.vcf.gz
        bcftools index -t svs_decomposed.vcf.gz
        SV_DECOMP_COUNT=\$(bcftools view -H svs_decomposed.vcf.gz | wc -l)
        echo "  After decomposition: \$SV_DECOMP_COUNT SV records" >&2

        # Concatenate decomposed SVs with normalized non-SVs
        echo "Step 4b: Concatenating SVs (decomposed) with non-SVs (normalized)..." >&2
        bcftools concat \
        	--allow-overlaps \
        	-Ou \
        	svs_decomposed.vcf.gz non_svs_norm.vcf.gz | \
        bcftools sort \
            -Ou \
            - | \
        bcftools +fill-tags \
            --no-version \
            -Oz \
            -o "${params.cohort}_merged.vcf.bgz" \
            -W=tbi \
            - -- \
            -t AC,AF,AN

        # Validation
        FINAL_COUNT=\$(bcftools view -H "${params.cohort}_merged.vcf.bgz" | wc -l)
        FINAL_SV_COUNT=\$(bcftools view -H -i 'INFO/SVTYPE!="." || ALT~"<"' "${params.cohort}_merged.vcf.bgz" | wc -l)
        echo "Step 5: Validation complete" >&2
        echo "  Total variants in merged VCF: \$FINAL_COUNT" >&2
        echo "  SVs in merged VCF: \$FINAL_SV_COUNT" >&2

        # Clean up temporary files
        rm -f merged_all.bcf svs.vcf.gz* svs_decomposed.vcf.gz* non_svs.vcf.gz* non_svs_norm.vcf.gz* *.filtered.vcf.gz*

        echo "=== VCF Merge Complete - SV Coordinates Preserved ===" >&2
        """
}
