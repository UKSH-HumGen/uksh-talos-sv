"""
Categorize structural variants based on gene/exon overlap and impact.

This script:
1. Reads SVs from a VCF file
2. Parses GFF3 for gene/exon coordinates
3. Determines gene overlap for each SV
4. Calculates impact based on overlap with coding regions
5. Adds category boolean fields (categorybooleanhighimpact, etc.)
6. Writes categorized SV VCF

Similar to RunHailFiltering but specialized for structural variants.
"""

import gzip
from argparse import ArgumentParser
from collections import defaultdict
from dataclasses import dataclass
from typing import IO, TYPE_CHECKING

import pysam
from loguru import logger

from talos.models import PanelApp
from talos.utils import read_json_from_path

if TYPE_CHECKING:
    from collections.abc import Callable


@dataclass
class Gene:
    """Represents a gene with its coordinates and exons."""

    gene_id: str
    gene_name: str
    chrom: str
    start: int
    end: int
    biotype: str
    exons: list[tuple[int, int]]  # List of (start, end) tuples


def parse_gff3_genes(gff3_path: str) -> dict[str, list[Gene]]:
    """
    Parse GFF3 file and extract gene/exon information.

    GFF3 hierarchy: gene → transcript → exon
    We need to map transcript_id → gene_id to associate exons with genes.

    Returns:
        dict mapping chromosome -> list of Gene objects
    """
    logger.info(f'Parsing GFF3 file: {gff3_path}')

    genes_by_chr = defaultdict(list)
    gene_dict = {}  # gene_id -> Gene object
    transcript_to_gene = {}  # transcript_id -> gene_id

    open_func: Callable[..., IO[str]] = gzip.open if gff3_path.endswith('.gz') else open  # type: ignore[assignment]

    # First pass: collect genes and transcripts
    with open_func(gff3_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, _source, feature_type, start, end, _score, _strand, _phase, attributes = fields

            # Parse attributes into dict
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

            # Process genes
            if feature_type == 'gene':
                gene_id = attr_dict.get('gene_id', attr_dict.get('ID', '').replace('gene:', ''))
                gene_name = attr_dict.get('Name', attr_dict.get('gene_name', ''))
                biotype = attr_dict.get('biotype', attr_dict.get('gene_type', 'unknown'))

                if gene_id:
                    gene = Gene(
                        gene_id=gene_id,
                        gene_name=gene_name,
                        chrom=chrom,
                        start=int(start),
                        end=int(end),
                        biotype=biotype,
                        exons=[],
                    )
                    gene_dict[gene_id] = gene
                    genes_by_chr[chrom].append(gene)

            # Process transcripts (to map transcript → gene)
            elif feature_type in ['mRNA', 'transcript', 'lnc_RNA', 'miRNA', 'ncRNA', 'pseudogenic_transcript']:
                transcript_id = attr_dict.get('ID', '').replace('transcript:', '')
                parent = attr_dict.get('Parent', '').replace('gene:', '')
                gene_id = attr_dict.get('gene_id', parent)

                if transcript_id and gene_id:
                    transcript_to_gene[transcript_id] = gene_id

    # Second pass: collect exons and associate with genes via transcripts
    with open_func(gff3_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, _source, feature_type, start, end, _score, _strand, _phase, attributes = fields

            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

            # Process exons
            if feature_type == 'exon':
                parent = attr_dict.get('Parent', '').replace('transcript:', '')
                exon_gene_id = transcript_to_gene.get(parent)

                if exon_gene_id and exon_gene_id in gene_dict:
                    gene_dict[exon_gene_id].exons.append((int(start), int(end)))

    total_genes = sum(len(genes) for genes in genes_by_chr.values())
    total_exons = sum(len(gene.exons) for gene in gene_dict.values())
    logger.info(f'Parsed {total_genes} genes with {total_exons} total exons from GFF3')

    return genes_by_chr


def calculate_overlap(start1: int, end1: int, start2: int, end2: int) -> int:
    """Calculate overlap between two intervals."""
    return max(0, min(end1, end2) - max(start1, start2))


def find_overlapping_genes(
    chrom: str, sv_start: int, sv_end: int, genes_by_chr: dict, debug: bool = False
) -> list[Gene]:
    """
    Find genes that overlap with an SV.
    Handles both chr2 and 2 chromosome naming conventions.

    Returns:
        list of Gene objects that overlap the SV
    """
    overlapping = []

    # Try both with and without 'chr' prefix
    test_chroms = [chrom, chrom.replace('chr', ''), f'chr{chrom}']

    if debug:
        logger.debug(f'  Looking for chromosome: {chrom}')
        logger.debug(f'  Will try: {test_chroms}')
        logger.debug(f'  Available chromosomes: {list(genes_by_chr.keys())[:10]}')

    for test_chrom in test_chroms:
        genes_on_chr = genes_by_chr.get(test_chrom, [])
        if debug and test_chrom in genes_by_chr:
            logger.debug(f'  Found {len(genes_on_chr)} genes for chromosome {test_chrom}')
            if genes_on_chr:
                logger.debug(
                    f'    First gene: {genes_on_chr[0].gene_name} at {genes_on_chr[0].start}-{genes_on_chr[0].end}'
                )
                logger.debug(f'    SV position: {sv_start}-{sv_end}')

        for gene in genes_on_chr:
            overlap = calculate_overlap(sv_start, sv_end, gene.start, gene.end)
            if debug and overlap > 0:
                logger.debug(f'    ✓ Overlap found with {gene.gene_name}: {overlap}bp')
            if overlap > 0:
                overlapping.append(gene)
        if overlapping:  # Stop once we find matches
            break

    return overlapping


def calculate_exon_overlap(sv_start: int, sv_end: int, gene: Gene) -> tuple[int, float]:
    """
    Calculate how much of a gene's exonic sequence overlaps with an SV.

    Returns:
        (exon_overlap_bp, exon_overlap_percentage)
    """
    if not gene.exons:
        return 0, 0.0

    total_exon_overlap = 0
    total_exon_length = 0

    for exon_start, exon_end in gene.exons:
        overlap = calculate_overlap(sv_start, sv_end, exon_start, exon_end)
        total_exon_overlap += overlap
        total_exon_length += exon_end - exon_start

    exon_overlap_pct = (total_exon_overlap / total_exon_length * 100) if total_exon_length > 0 else 0.0

    return total_exon_overlap, exon_overlap_pct


def categorize_sv(
    chrom: str,
    pos: int,
    sv_end: int,
    sv_type: str,
    sv_len: int,
    genes_by_chr: dict,
    green_genes: set,
    debug: bool = False,
) -> dict:
    """
    Categorize a structural variant based on gene/exon overlap.

    Returns:
        dict of category annotations to add to VCF INFO
    """
    categories = {
        'gene_id': '.',
        'gene_symbol': '.',
        'categorybooleanhighimpact': '0',
        'categorybooleansvdb': '1',  # All SVs are marked as SVDB
        'exon_overlap_bp': '0',
        'exon_overlap_pct': '0.0',
    }

    # Find overlapping genes
    overlapping_genes = find_overlapping_genes(chrom, pos, sv_end, genes_by_chr, debug=debug)

    if not overlapping_genes:
        return categories

    # Filter to protein-coding genes if possible
    protein_coding = [g for g in overlapping_genes if g.biotype == 'protein_coding']
    if protein_coding:
        overlapping_genes = protein_coding

    # Filter to green genes if possible
    green_overlapping = [g for g in overlapping_genes if g.gene_id in green_genes]
    if green_overlapping:
        overlapping_genes = green_overlapping

    # Use the first (or most significant) overlapping gene
    gene = overlapping_genes[0]

    # Calculate exon overlap
    exon_overlap_bp, exon_overlap_pct = calculate_exon_overlap(pos, sv_end, gene)

    # Determine if high impact
    # High impact if:
    # - Deletion/Duplication in protein-coding gene
    # - Overlaps exons (>0% exonic overlap)
    # - Or is a large SV (>1000bp) affecting any part of a green gene
    is_high_impact = False
    if sv_type in ['DEL', 'DUP', 'CNV']:
        if exon_overlap_pct > 0:  # Any exonic overlap
            is_high_impact = True
        elif abs(sv_len) > 1000 and gene.gene_id in green_genes:
            is_high_impact = True

    categories.update(
        {
            'gene_id': gene.gene_id,
            'gene_symbol': gene.gene_name or gene.gene_id,
            'categorybooleanhighimpact': '1' if is_high_impact else '0',
            'exon_overlap_bp': str(int(exon_overlap_bp)),
            'exon_overlap_pct': f'{exon_overlap_pct:.2f}',
        }
    )

    return categories


def categorize_sv_vcf(input_vcf: str, output_vcf: str, gff3_path: str, panelapp_path: str):
    """
    Main function to categorize SVs in a VCF file.
    """
    logger.info('=== Starting SV Categorization ===')

    # Parse GFF3 for gene coordinates
    genes_by_chr = parse_gff3_genes(gff3_path)

    # Load PanelApp data for green genes
    # The genes dict is keyed by Ensembl gene_id, so we can just use the keys
    logger.info(f'Loading PanelApp data: {panelapp_path}')
    panel_data = read_json_from_path(panelapp_path, return_model=PanelApp)
    green_genes = set(panel_data.genes.keys())
    logger.info(f'Loaded {len(green_genes)} green genes from PanelApp')

    # Open input/output VCFs
    vcf_in = pysam.VariantFile(input_vcf)

    # Add new INFO fields to header
    vcf_in.header.info.add('gene_id', '1', 'String', 'Ensembl gene ID')
    vcf_in.header.info.add('gene_symbol', '1', 'String', 'Gene symbol')
    vcf_in.header.info.add('categorybooleanhighimpact', '1', 'Integer', 'High impact SV category')
    vcf_in.header.info.add('categorybooleansvdb', '1', 'Integer', 'Structural variant database category')
    vcf_in.header.info.add('exon_overlap_bp', '1', 'Integer', 'Exonic basepairs overlapped')
    vcf_in.header.info.add('exon_overlap_pct', '1', 'Float', 'Percentage of exons overlapped')

    vcf_out = pysam.VariantFile(output_vcf, 'wz', header=vcf_in.header)

    # Process each variant
    total_svs = 0
    categorized_svs = 0
    high_impact_svs = 0

    # Debug: Log first few variants to diagnose chromosome matching
    for record in vcf_in:
        total_svs += 1

        # Extract SV information
        chrom = record.chrom
        pos = record.pos
        sv_type = record.info.get('SVTYPE', 'UNK')
        sv_end = record.info.get('END', pos)
        sv_len = record.info.get('SVLEN', 0)

        # Handle SVLEN as tuple (pysam returns tuples for multi-value fields)
        if isinstance(sv_len, tuple):
            sv_len = sv_len[0] if sv_len else 0

        # Calculate END from SVLEN if END is not set properly
        # For DEL/DUP: END = POS + abs(SVLEN)
        # For INS: END = POS (insertions are at a single position)
        if sv_end == pos and sv_len != 0:
            if sv_type in ['DEL', 'DUP', 'CNV']:
                sv_end = pos + abs(sv_len)

        # Enable detailed debug logging for first 3 SVs
        debug = total_svs <= 3
        if debug:
            logger.info(f'\n=== Analyzing SV #{total_svs}: {chrom}:{pos}-{sv_end} ({sv_type}, SVLEN={sv_len}) ===')

        # Categorize the SV
        categories = categorize_sv(chrom, pos, sv_end, sv_type, sv_len, genes_by_chr, green_genes, debug=debug)

        # Add categories to INFO
        if categories['gene_id'] != '.':
            categorized_svs += 1
            if categories['categorybooleanhighimpact'] == '1':
                high_impact_svs += 1

        for key, value in categories.items():
            if key in ['categorybooleanhighimpact', 'categorybooleansvdb', 'exon_overlap_bp']:
                record.info[key] = int(value)
            elif key == 'exon_overlap_pct':
                record.info[key] = float(value)
            else:
                record.info[key] = value

        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()

    logger.info(f'Processed {total_svs} structural variants')
    logger.info(f'  - {categorized_svs} overlapping genes ({categorized_svs / total_svs * 100:.1f}%)')
    logger.info(f'  - {high_impact_svs} high impact ({high_impact_svs / total_svs * 100:.1f}%)')
    logger.info(f'Output written to: {output_vcf}')
    logger.info('=== SV Categorization Complete ===')


def cli_main():
    """CLI entry point."""
    parser = ArgumentParser(description='Categorize structural variants based on gene/exon overlap')
    parser.add_argument('--input', required=True, help='Input SV VCF file')
    parser.add_argument('--output', required=True, help='Output categorized VCF file')
    parser.add_argument('--gff3', required=True, help='GFF3 gene annotation file')
    parser.add_argument('--panelapp', required=True, help='PanelApp JSON file')

    args = parser.parse_args()

    categorize_sv_vcf(
        input_vcf=args.input,
        output_vcf=args.output,
        gff3_path=args.gff3,
        panelapp_path=args.panelapp,
    )


if __name__ == '__main__':
    cli_main()
