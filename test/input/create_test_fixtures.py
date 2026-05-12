#!/usr/bin/env python3
"""
Generates the synthetic VCF test fixtures that are not committed to the repo
(because cyvcf2 / bgzip are not available in the bare CI workspace).

Run once on a machine with htslib + bgzip + tabix available:

    python3 test/input/create_test_fixtures.py

Fixtures produced
-----------------
1_labelled_variant.vcf.bgz + .tbi
    Two small-variant calls on chr20 for a trio (male proband, mother_1, father_1).
    Variant 1: chr20:63406931 C>CGG  - CategoryBoolean3=1, CategorySample4=male
    Variant 2: chr20:63406991 C>CGG  - CategoryBoolean3=1, CategoryBoolean1=1, CategorySample4=male
    Both annotated against ENSG00000075043.

newphase.vcf.bgz + .tbi
    Same two loci, same trio.
    mother_1 carries both variants as phased hets (0|1) in phase-set PS=420.
    male and father_1 are WT (0/0).
    Variant 1: CategoryBoolean3=1, CategorySample4=mother_1
    Variant 2: CategoryBoolean3=1, CategorySample4=mother_1
"""

import subprocess
import sys
import os

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def write_vcf(path: str, content: str) -> None:
    """Write plain VCF text, then bgzip + tabix it."""
    plain = path.replace('.vcf.bgz', '.vcf')
    with open(plain, 'w') as fh:
        fh.write(content)
    subprocess.run(['bgzip', '-f', plain], check=True)
    subprocess.run(['tabix', '-p', 'vcf', path], check=True)
    print(f'Created: {path}')


# ---------------------------------------------------------------------------
# 1_labelled_variant.vcf.bgz
# ---------------------------------------------------------------------------
LABELLED_HEADER = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr20,length=64444167>
##INFO=<ID=Gene_id,Number=1,Type=String,Description="Ensembl gene ID">
##INFO=<ID=CategoryBoolean3,Number=1,Type=Integer,Description="High Impact category">
##INFO=<ID=CategoryBoolean1,Number=1,Type=Integer,Description="ClinVar P/LP category">
##INFO=<ID=CategorySample4,Number=.,Type=String,Description="De Novo category samples">
##INFO=<ID=gnomad_af,Number=1,Type=Float,Description="gnomAD AF">
##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmale\tmother_1\tfather_1
"""

LABELLED_ROWS = [
    # Variant 1: cat3 + cat4(de_novo) for male only
    'chr20\t63406931\t.\tC\tCGG\t100\tPASS\t'
    'Gene_id=ENSG00000075043;CategoryBoolean3=1;CategorySample4=male;gnomad_af=0.0001;AC=1;AF=0.0001\t'
    'GT:AD:DP:GQ\t0/1:20,20:40:99\t0/0:40,0:40:99\t0/0:40,0:40:99\n',
    # Variant 2: cat1 + cat3 + cat4(de_novo) for male only
    'chr20\t63406991\t.\tC\tCGG\t100\tPASS\t'
    'Gene_id=ENSG00000075043;CategoryBoolean3=1;CategoryBoolean1=1;CategorySample4=male;gnomad_af=0.0001;AC=1;AF=0.0001\t'
    'GT:AD:DP:GQ\t0/1:20,20:40:99\t0/0:40,0:40:99\t0/0:40,0:40:99\n',
]

labelled_path = os.path.join(OUT_DIR, '1_labelled_variant.vcf.bgz')
write_vcf(labelled_path, LABELLED_HEADER + ''.join(LABELLED_ROWS))


# ---------------------------------------------------------------------------
# newphase.vcf.bgz
# ---------------------------------------------------------------------------
PHASE_HEADER = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr20,length=64444167>
##INFO=<ID=Gene_id,Number=1,Type=String,Description="Ensembl gene ID">
##INFO=<ID=CategoryBoolean3,Number=1,Type=Integer,Description="High Impact category">
##INFO=<ID=CategorySample4,Number=.,Type=String,Description="De Novo category samples">
##INFO=<ID=gnomad_af,Number=1,Type=Float,Description="gnomAD AF">
##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmale\tmother_1\tfather_1
"""

PHASE_ROWS = [
    # Variant 1: mother_1 het phased (0|1) PS=420, male and father_1 WT
    'chr20\t63406931\t.\tC\tCGG\t100\tPASS\t'
    'Gene_id=ENSG00000075043;CategoryBoolean3=1;CategorySample4=mother_1;gnomad_af=0.0001;AC=1;AF=0.0001\t'
    'GT:AD:DP:GQ:PS\t0/0:40,0:40:99:.\t0|1:20,20:40:99:420\t0/0:40,0:40:99:.\n',
    # Variant 2: same phase set for mother_1
    'chr20\t63406991\t.\tC\tCGG\t100\tPASS\t'
    'Gene_id=ENSG00000075043;CategoryBoolean3=1;CategorySample4=mother_1;gnomad_af=0.0001;AC=1;AF=0.0001\t'
    'GT:AD:DP:GQ:PS\t0/0:40,0:40:99:.\t0|1:20,20:40:99:420\t0/0:40,0:40:99:.\n',
]

phase_path = os.path.join(OUT_DIR, 'newphase.vcf.bgz')
write_vcf(phase_path, PHASE_HEADER + ''.join(PHASE_ROWS))

print('\nDone. Both fixtures written and indexed.')
