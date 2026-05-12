process AnnotateCsqWithBcftools {
    container params.container

    input:
        path(vcf)
        path(gff3)
        path(reference)

    // annotate this VCF with gnomAD data
    publishDir params.cohort_output_dir, mode: 'copy'

    output:
        path("${params.cohort}_csq.vcf.bgz")

    script:
    """
    # Remove 'chr' prefix from chromosome names to match GFF3 format
    bcftools annotate --rename-chrs <(echo 'chr1 1'; echo 'chr2 2'; echo 'chr3 3'; echo 'chr4 4'; echo 'chr5 5'; echo 'chr6 6'; echo 'chr7 7'; echo 'chr8 8'; echo 'chr9 9'; echo 'chr10 10'; echo 'chr11 11'; echo 'chr12 12'; echo 'chr13 13'; echo 'chr14 14'; echo 'chr15 15'; echo 'chr16 16'; echo 'chr17 17'; echo 'chr18 18'; echo 'chr19 19'; echo 'chr20 20'; echo 'chr21 21'; echo 'chr22 22'; echo 'chrX X'; echo 'chrY Y'; echo 'chrM MT') \
        -Ov ${vcf} | \\
    bcftools view -e 'ALT~"^<"' -Ov | \\
    bcftools view -e '(POS=248361635 && CHROM="1") || (POS=11155134 && CHROM="2") || (POS=16323553 && CHROM="6")' -Ov | \\
    bcftools csq --force -f "${reference}" \
        --local-csq \
        -g ${gff3} \
        -B 20 \
        -Ov | \\
    bcftools annotate --rename-chrs <(echo '1 chr1'; echo '2 chr2'; echo '3 chr3'; echo '4 chr4'; echo '5 chr5'; echo '6 chr6'; echo '7 chr7'; echo '8 chr8'; echo '9 chr9'; echo '10 chr10'; echo '11 chr11'; echo '12 chr12'; echo '13 chr13'; echo '14 chr14'; echo '15 chr15'; echo '16 chr16'; echo '17 chr17'; echo '18 chr18'; echo '19 chr19'; echo '20 chr20'; echo '21 chr21'; echo '22 chr22'; echo 'X chrX'; echo 'Y chrY'; echo 'MT chrM') \
        -Oz -o "${params.cohort}_csq.vcf.bgz"
    """
}
