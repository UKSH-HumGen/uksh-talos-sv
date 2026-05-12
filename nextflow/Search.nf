nextflow.enable.dsl=2

params.region         = params.region ?: ''
params.vcfs           = params.vcfs   ?: ''
params.sample_map     = params.sample_map ?: ''  // 2-column TSV: SAMPLE<TAB>PATIENT_ID
params.out            = params.out ?: 'region_hits'

if (!params.region || !params.vcfs || !params.sample_map) {
    log.error """
    Required:
      --region       chr:start-end   (e.g. chr7:55000000-56000000)
      --vcfs         path/pattern to annotated VCF/BCF files (bgzipped & indexed)
      --sample_map   TSV mapping: SAMPLE<TAB>PATIENT_ID
    """
    System.exit(1)
}

Channel
  .fromPath(params.vcfs.split(','))
  .ifEmpty { error "No VCF/BCF matched: ${params.vcfs}" }
  .set { VCF_FILES }

SAMPLE_MAP = file(params.sample_map)

process ExtractRegionVariants {
    tag { "${vcf.simpleName}" }
    publishDir "${params.out}/per_vcf", mode: 'copy'

    input:
    path vcf
    val region from params.region
    path SAMPLE_MAP

    output:
    path "${vcf.simpleName}.region_hits.tsv"

    container 'biocontainers/bcftools:v1.20-1-deb_cv1'  // or your env with bcftools >=1.10

    script:
    // Parse region into CHR, START, END for SV overlap logic
    """
    set -euo pipefail

    REGION='${region}'
    CHR=\$(echo "\$REGION" | cut -d: -f1)
    S=\$(echo "\$REGION" | cut -d: -f2 | cut -d- -f1)
    E=\$(echo "\$REGION" | cut -d- -f2)

    # Ensure index exists (no-op if already present)
    if [ ! -f "${vcf}.csi" ] && [ ! -f "${vcf}.tbi" ]; then
      bcftools index -f "${vcf}"
    fi

    # Filter records overlapping the interval.
    # Logic:
    # - Point/short variants (no INFO/END): POS in [S,E]
    # - SV/CNV with END: interval [POS, END] intersects [S,E]
    bcftools view -i '((!INFO/END) && CHROM==\$CHR && POS>=\$S && POS<=\$E) || (INFO/END && CHROM==\$CHR && INFO/END>=\$S && POS<=\$E)' \
                  -Oz -o region.subset.vcf.gz "${vcf}"
    bcftools index -f region.subset.vcf.gz

    # From the subset, expand to long format with per-sample non-reference carriers.
    # Keep only samples where GT is not 0/0 and not ./.
    # Output columns:
    # CHROM POS END REF ALT SVTYPE VARIANT_ID SAMPLE GT GQ DP
    bcftools query -f '%CHROM\\t%POS\\t[%INFO/END]\\t%REF\\t%ALT\\t%INFO/SVTYPE\\t%ID[\\t%SAMPLE=%GT,%GQ,%DP]\\n' region.subset.vcf.gz \
    | awk -F'\\t' '
        BEGIN{OFS="\\t"}
        {
          chrom=\$1; pos=\$2; end=\$3; ref=\$4; alt=\$5; svtype=\$6; vid=\$7;
          # iterate sample triplets starting at field 8
          for (i=8; i<=NF; i++){
            split(\$i, a, "="); smp=a[1]; split(a[2], b, ","); gt=b[1]; gq=b[2]; dp=b[3];
            if (gt != "0/0" && gt != "0|0" && gt != "./." && gt != ".|.") {
              print chrom, pos, (end==""?"." : end), ref, alt, (svtype==""?"." : svtype), (vid==""?chrom":"pos":"-":alt), smp, gt, gq, dp;
            }
          }
        }
      ' \
    > carriers.tmp.tsv

    # Join to patient IDs
    # SAMPLE_MAP: SAMPLE<TAB>PATIENT_ID
    # Output columns:
    # CHROM POS END REF ALT SVTYPE VARIANT_ID SAMPLE PATIENT_ID GT GQ DP
    join -t $'\\t' -1 8 -2 1 \
         <(sort -k8,8 carriers.tmp.tsv) \
         <(sort -k1,1 "${SAMPLE_MAP}") \
    | awk -F'\\t' 'BEGIN{OFS="\\t"} { print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$12,\$9,\$10,\$11 }' \
    > "${vcf.simpleName}.region_hits.tsv"
    """
}

process MergeAll {
    publishDir "${params.out}", mode: 'copy'
    input:
    path(hit) from ExtractRegionVariants.out.collect()

    output:
    path "all_region_hits.tsv"

    script:
    """
    set -euo pipefail
    header="CHROM\tPOS\tEND\tREF\tALT\tSVTYPE\tVARIANT_ID\tSAMPLE\tPATIENT_ID\tGT\tGQ\tDP"
    echo -e "\$header" > all_region_hits.tsv
    for f in ${hit.join(' ')}; do
      cat "\$f" >> all_region_hits.tsv
    done
    """
}

workflow {
    take: VCF_FILES
    main:
        ExtractRegionVariants(VCF_FILES, SAMPLE_MAP)
        MergeAll(ExtractRegionVariants.out)
    emit:
        MergeAll.out
}
