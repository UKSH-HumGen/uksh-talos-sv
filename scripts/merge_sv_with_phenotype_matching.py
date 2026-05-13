#!/usr/bin/env python3
"""
Enhanced script to merge structural variants with comprehensive gene annotation and HPO phenotype matching.

This script:
1. Annotates SVs with overlapping genes from GFF3
2. Performs HPO phenotype matching between patient and gene
3. Calculates phenotype similarity scores
4. Filters to show only phenotype-matched SVs (optional)
"""

import gzip
import json
import os
import sys
from collections import defaultdict
from datetime import datetime

import toml

# Import the base functionality from the enhanced merge script
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def load_sv_filter_config(config_path=None):
    """Load structural variant filtering configuration from TOML file."""
    if config_path is None:
        config_path = os.environ.get('TALOS_CONFIG', 'inputs/config.toml')

    try:
        with open(config_path) as f:
            config = toml.load(f)

        sv_config = config.get('StructuralVariantFiltering', {})
        print(f'Loaded SV filtering config from {config_path}')
        print(f'  - require_gene_overlap: {sv_config.get("require_gene_overlap", False)}')
        print(f'  - require_exon_overlap: {sv_config.get("require_exon_overlap", False)}')
        print(f'  - require_phenotype_match: {sv_config.get("require_phenotype_match", False)}')
        print(f'  - min_phenotype_score: {sv_config.get("min_phenotype_score", 0.0)}')
        print(f'  - allowed_sv_types: {sv_config.get("allowed_sv_types", [])}')
        print(f'  - require_pass_only: {sv_config.get("require_pass_only", False)}')
        print(f'  - min_sv_size: {sv_config.get("min_sv_size", 0)}')

        return sv_config
    except Exception as e:
        print(f'Warning: Could not load SV filter config from {config_path}: {e}')
        print('Using default filtering parameters')
        return {}


def load_gene_to_hpo_mapping(gen2phen_file):
    """Load gene to HPO term mappings from genes_to_phenotype.txt file."""
    gene_to_hpo = defaultdict(set)

    try:
        with open(gen2phen_file) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    gene_symbol = fields[1]
                    hpo_term = fields[2]
                    gene_to_hpo[gene_symbol].add(hpo_term)

        print(f'Loaded HPO mappings for {len(gene_to_hpo)} genes')
    except Exception as e:
        print(f'Warning: Could not load gene-to-HPO mappings: {e}')

    return gene_to_hpo


def load_patient_hpo_terms(pedigree_file, sample_id):
    """Load patient HPO terms from pedigree file."""
    patient_hpo = set()

    try:
        with open(pedigree_file) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 7:
                    sample = fields[0]
                    if sample == sample_id:
                        hpo_string = fields[6]
                        # Parse comma-separated HPO terms
                        hpo_terms = hpo_string.split(',')
                        patient_hpo = {term.strip() for term in hpo_terms if term.strip().startswith('HP:')}
                        break

        print(f'Loaded {len(patient_hpo)} HPO terms for patient {sample_id}')
        print(f'Patient HPO terms: {patient_hpo}')
    except Exception as e:
        print(f'Warning: Could not load patient HPO terms: {e}')

    return patient_hpo


def calculate_phenotype_match(gene_symbols, patient_hpo, gene_to_hpo):
    """
    Calculate phenotype matching between patient HPO terms and gene-associated HPO terms.

    Returns:
        - matched_genes: list of genes with phenotype matches
        - match_details: dict with details about the matches
    """
    matched_genes = []
    match_details = {'total_overlaps': 0, 'matching_hpo_terms': set(), 'genes_with_matches': []}

    for gene_symbol in gene_symbols:
        gene_hpo = gene_to_hpo.get(gene_symbol, set())

        if not gene_hpo:
            continue

        # Find direct HPO term overlaps
        overlapping_hpo = patient_hpo & gene_hpo

        if overlapping_hpo:
            matched_genes.append(gene_symbol)
            match_details['total_overlaps'] += len(overlapping_hpo)
            match_details['matching_hpo_terms'].update(overlapping_hpo)
            match_details['genes_with_matches'].append(
                {
                    'gene': gene_symbol,
                    'matching_terms': list(overlapping_hpo),
                    'match_count': len(overlapping_hpo),
                    'total_gene_hpo': len(gene_hpo),
                    'match_percentage': (len(overlapping_hpo) / len(patient_hpo)) * 100 if patient_hpo else 0,
                }
            )

    return matched_genes, match_details


def load_gene_annotations(gff_file=None):
    """Load gene annotations from GFF3 file organized by chromosome."""
    gene_annotations = {}

    if not gff_file:
        # Try to find GFF file
        script_dir = os.path.dirname(os.path.abspath(__file__))
        possible_paths = [
            'large_files/Homo_sapiens.GRCh38.113.gff3.gz',
            '../large_files/Homo_sapiens.GRCh38.113.gff3.gz',
            '../../large_files/Homo_sapiens.GRCh38.113.gff3.gz',
            '../../../large_files/Homo_sapiens.GRCh38.113.gff3.gz',
            '../../../../large_files/Homo_sapiens.GRCh38.113.gff3.gz',
            os.path.join(script_dir, '../large_files/Homo_sapiens.GRCh38.113.gff3.gz'),
            '/Users/Daniel-K/Desktop/Reanalysis-Talos/talos-main/large_files/Homo_sapiens.GRCh38.113.gff3.gz',
        ]
        for path in possible_paths:
            if os.path.exists(path):
                gff_file = path
                break

    if not gff_file or not os.path.exists(gff_file):
        print('Warning: No GFF3 file found for gene annotation')
        return gene_annotations

    try:
        opener = gzip.open if gff_file.endswith('.gz') else open
        gene_count = 0

        with opener(gff_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                chrom, _source, feature, start, end, _score, strand, _phase, attributes = fields

                if feature != 'gene':
                    continue

                # Parse attributes
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value

                gene_id = attr_dict.get('ID', '').replace('gene:', '')
                gene_name = attr_dict.get('Name', '')
                gene_biotype = attr_dict.get('biotype', 'protein_coding')

                if gene_name and chrom:
                    chrom_clean = chrom.replace('chr', '')

                    if chrom_clean not in gene_annotations:
                        gene_annotations[chrom_clean] = []

                    gene_annotations[chrom_clean].append(
                        {
                            'gene_id': gene_id,
                            'gene_name': gene_name,
                            'biotype': gene_biotype,
                            'start': int(start),
                            'end': int(end),
                            'strand': strand,
                        }
                    )
                    gene_count += 1

        # Sort genes by start position for efficient overlap detection
        for chrom in gene_annotations:
            gene_annotations[chrom].sort(key=lambda x: x['start'])

        print(f'Loaded {gene_count} genes from GFF3 file')

        # Debug: print chromosome distribution
        chrom_counts = {chrom: len(genes) for chrom, genes in gene_annotations.items()}
        print(
            f'Gene distribution: {sorted(chrom_counts.items(), key=lambda x: int(x[0]) if x[0].isdigit() else 99)[:10]}'
        )

    except Exception as e:
        print(f'Error loading GFF3 file: {e}')
        import traceback

        traceback.print_exc()

    return gene_annotations


def load_exon_annotations():
    """Load exon annotations from GFF3 file."""
    gff3_file = os.environ.get('GFF3_FILE', '../large_files/Homo_sapiens.GRCh38.113.gff3.gz')

    exons_by_chr = defaultdict(list)

    try:
        print(f'Loading exon annotations from {gff3_file}...')
        with gzip.open(gff3_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                chrom, _source, feature_type, start, end, _score, strand, _phase, _attributes = fields

                # Only process exon features
                if feature_type != 'exon':
                    continue

                chrom_normalized = chrom.replace('chr', '')

                exons_by_chr[chrom_normalized].append({'start': int(start), 'end': int(end), 'strand': strand})

        # Sort exons by start position for each chromosome
        for chrom in exons_by_chr:
            exons_by_chr[chrom].sort(key=lambda x: x['start'])

        total_exons = sum(len(exons) for exons in exons_by_chr.values())
        print(f'Loaded {total_exons:,} exons across {len(exons_by_chr)} chromosomes')
    except Exception as e:
        print(f'Warning: Could not load exon annotations: {e}')
        return {}

    return exons_by_chr


def check_exon_overlap(chrom, start, end, exon_annotations):
    """Check if SV overlaps with any exons."""
    chrom_clean = str(chrom).replace('chr', '')

    if chrom_clean not in exon_annotations:
        return False

    # Check for overlap with any exon
    for exon in exon_annotations[chrom_clean]:
        if not (end < exon['start'] or start > exon['end']):
            return True

    return False


def find_overlapping_genes(chrom, start, end, gene_annotations):
    """Find genes that overlap with the structural variant."""
    overlapping_genes = []

    chrom_clean = str(chrom).replace('chr', '')

    if chrom_clean not in gene_annotations:
        return overlapping_genes

    for gene in gene_annotations[chrom_clean]:
        gene_start = gene['start']
        gene_end = gene['end']

        # Check for overlap
        if not (end < gene_start or start > gene_end):
            # Calculate overlap
            overlap_start = max(start, gene_start)
            overlap_end = min(end, gene_end)
            overlap_length = overlap_end - overlap_start + 1

            gene_length = gene_end - gene_start + 1
            sv_length = end - start + 1

            gene['overlap_start'] = overlap_start
            gene['overlap_end'] = overlap_end
            gene['overlap_length'] = overlap_length
            gene['gene_overlap_percentage'] = (overlap_length / gene_length) * 100 if gene_length > 0 else 0
            gene['sv_overlap_percentage'] = (overlap_length / sv_length) * 100 if sv_length > 0 else 0

            overlapping_genes.append(gene.copy())

    return overlapping_genes


def parse_vcf_header(vcf_file):
    """Parse VCF header to extract sample names."""
    sample_names = []
    with open(vcf_file) as f:
        for line in f:
            if line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                # Sample names start after FORMAT column (column 9)
                if len(fields) > 9:
                    sample_names = fields[9:]
                break
    return sample_names


def parse_vcf_line(line, sample_names, target_sample_id=None):
    """Parse a VCF line and extract structural variant information.

    Args:
        line: VCF line to parse
        sample_names: List of sample names from VCF header
        target_sample_id: If provided, only return variant if this sample has it (non-ref genotype)

    Returns:
        Variant dict if it's a valid SV and (if target_sample_id specified) the target sample has it,
        None otherwise
    """
    if line.startswith('#'):
        return None

    fields = line.strip().split('\t')
    if len(fields) < 10:
        return None

    chrom, pos, var_id, ref, alt, qual, filter_status, info, format_field = fields[:9]
    sample_data_fields = fields[9:]  # All sample columns

    # Parse INFO field
    info_dict = {}
    for item in info.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True

    # Process structural variants
    if 'SVTYPE' not in info_dict or info_dict['SVTYPE'] not in ['CNV', 'DEL', 'INS', 'DUP', 'INV', 'BND']:
        return None

    # Parse format field
    format_keys = format_field.split(':')

    # If target_sample_id specified, find which sample column it is and check genotype
    genotype = './.'
    if target_sample_id and sample_names:
        try:
            sample_idx = sample_names.index(target_sample_id)
            if sample_idx < len(sample_data_fields):
                sample_values = sample_data_fields[sample_idx].split(':')
                sample_dict = dict(zip(format_keys, sample_values))
                genotype = sample_dict.get('GT', './.')

                # Filter out variants where this sample is homozygous reference (0/0 or 0|0)
                if genotype in ['0/0', '0|0', './.', '.|.']:
                    return None  # This sample doesn't have the variant
        except (ValueError, IndexError) as e:
            print(f'Warning: Could not find sample {target_sample_id} in VCF: {e}')
            return None
    elif len(sample_data_fields) > 0:
        # No target sample specified, use first sample's genotype
        sample_values = sample_data_fields[0].split(':')
        sample_dict = dict(zip(format_keys, sample_values))
        genotype = sample_dict.get('GT', './.')

    variant = {
        'chrom': chrom.replace('chr', ''),
        'pos': int(pos),
        'id': var_id,
        'ref': ref,
        'alt': alt,
        'qual': float(qual) if qual != '.' else 0,
        'filter': filter_status,
        'svtype': info_dict.get('SVTYPE', ''),
        'svlen': int(info_dict.get('SVLEN', 0)) if info_dict.get('SVLEN', '').lstrip('-').isdigit() else 0,
        'end': int(info_dict.get('END', pos)) if info_dict.get('END', '').isdigit() else int(pos),
        'genotype': genotype,
    }

    return variant


def create_talos_variant_with_phenotype(sv_variant, sample_id, overlapping_genes, phenotype_match_details):
    """Convert structural variant to Talos format with phenotype matching information."""

    # Determine variant type
    consequence = 'structural_variant'
    if sv_variant['svtype'] == 'CNV':
        consequence = 'copy_number_loss' if sv_variant['svlen'] < 0 else 'copy_number_gain'
    elif sv_variant['svtype'] == 'DEL':
        consequence = 'copy_number_loss'
    elif sv_variant['svtype'] == 'DUP':
        consequence = 'copy_number_gain'

    # Extract gene information
    gene_id = 'unknown'
    gene_name = 'unknown'
    transcript = 'unknown'
    biotype = 'unknown'

    if overlapping_genes:
        # Prioritize protein-coding genes with highest overlap
        priority_biotypes = ['protein_coding', 'lncRNA', 'miRNA']

        def gene_priority(gene):
            biotype = gene.get('biotype', 'unknown')
            try:
                biotype_priority = priority_biotypes.index(biotype)
            except ValueError:
                biotype_priority = len(priority_biotypes)
            gene_overlap = gene.get('gene_overlap_percentage', 0)
            return (biotype_priority, -gene_overlap)

        sorted_genes = sorted(overlapping_genes, key=gene_priority)
        primary_gene = sorted_genes[0]

        gene_name = primary_gene['gene_name']
        gene_id = primary_gene['gene_id']
        biotype = primary_gene['biotype']

        if len(overlapping_genes) > 1:
            gene_name = f'{gene_name} (+{len(overlapping_genes) - 1} more)'

        overlap_pct = primary_gene.get('gene_overlap_percentage', 0)
        transcript = f'{primary_gene["gene_name"]}_transcript (affects {overlap_pct:.1f}% of gene)'

    # Build support categories as list
    support_categories = ['Structural Variant', 'SVDB', 'High Impact']

    # Add phenotype match category if there are matches
    if phenotype_match_details['total_overlaps'] > 0:
        support_categories.append('Phenotype Match')

    # Create categories dict (for Talos format) - map each category to today's date
    today = datetime.now().strftime('%Y-%m-%d')
    categories_dict = dict.fromkeys(support_categories, today)

    # Calculate phenotype match score
    # Note: phenotype_labels will be set by HPOFlagging module later in the pipeline
    phenotype_match_score = 0.0
    if phenotype_match_details['genes_with_matches']:
        # Use the highest match percentage among all overlapping genes
        phenotype_match_score = max(g['match_percentage'] for g in phenotype_match_details['genes_with_matches'])

    # Determine correct ALT allele
    # For DRAGEN CNV calls, convert to DEL/DUP based on SVLEN sign
    if sv_variant['svtype'] == 'CNV':
        alt_allele = '<DEL>' if sv_variant['svlen'] < 0 else '<DUP>'
    else:
        # Use original ALT from VCF for other SV types
        alt_allele = sv_variant['alt']

    talos_variant = {
        'sample': sample_id,
        'gene': gene_id,  # Top-level gene field for HTML display
        'var_data': {
            'coordinates': {
                'chrom': sv_variant['chrom'],
                'pos': sv_variant['pos'],
                'ref': f'{sv_variant["end"]}',  # Show end position in ref field (HTML will show pos-ref)
                'alt': alt_allele,  # Use corrected ALT allele (DEL/DUP for CNVs)
            },
            'genotype': sv_variant['genotype'],  # Add genotype field to var_data
            'info': {
                'ac': 1 if sv_variant['genotype'] in ['0/1', '1/0', '1/1'] else 0,
                'af': 0.5
                if sv_variant['genotype'] in ['0/1', '1/0']
                else (1.0 if sv_variant['genotype'] == '1/1' else 0.0),
                'an': 2,
                'gene_id': gene_id,
                'gene_symbol': gene_name,
                'var_link': f'{sv_variant["chrom"]}-{sv_variant["pos"]}-{sv_variant["ref"]}-{sv_variant["alt"]}',
                'svlen': sv_variant['svlen'],
                'svtype': sv_variant['svtype'],
                'end_pos': sv_variant['end'],
                'genotype': sv_variant['genotype'],  # Also add to info dict for consistency
                # Phenotype matching fields
                'phenotype_match': len(phenotype_match_details['matching_hpo_terms']) > 0,
                'phenotype_match_score': phenotype_match_score,
                'phenotype_matched_genes': [g['gene'] for g in phenotype_match_details['genes_with_matches']],
                'phenotype_matching_terms': list(phenotype_match_details['matching_hpo_terms']),
                'categorybooleanphenomatch': len(phenotype_match_details['matching_hpo_terms']) > 0,
                # Standard fields
                'categorybooleanhighimpact': True,
                'categorybooleansvdb': True,
                'gnomad_af': 0.0,
            },
            'support_categories': support_categories,
            'transcript_consequences': [
                {
                    'consequence': consequence,
                    'gene_id': gene_id,
                    'gene': gene_name,
                    'transcript': transcript,
                    'mane_id': '',
                    'mane': '',
                    'biotype': biotype,
                    'dna_change': (
                        f'{sv_variant["chrom"]}:{sv_variant["pos"]}-{sv_variant["end"]}'
                        f'({sv_variant["svtype"]},{abs(sv_variant["svlen"])}bp)'
                    ),
                    'amino_acid_change': '',
                    'codon': '',
                }
            ],
            'dna_change': (
                f'{sv_variant["chrom"]}:{sv_variant["pos"]}-{sv_variant["end"]}'
                f'({sv_variant["svtype"]},{abs(sv_variant["svlen"])}bp)'
            ),
            'gene_symbol': gene_name,
            'transcript': transcript,
        },
        'phenotypes': [],  # Will be filled by HPOFlagging module
        'date_of_phenotype_match': None,  # Will be set by HPOFlagging if there's a match
        'categories': categories_dict,
        'first_tagged': today,
        'labels': [],
        'reasons': (
            f'Structural variant with {len(overlapping_genes)} overlapping gene(s), '
            f'{phenotype_match_details["total_overlaps"]} HPO term matches'
        ),
        'support_vars': [],
        'found_in_current_run': True,
        'panels': {'forced': {}, 'matched': {}},
        'phenotype_labels': [],
    }

    return talos_variant


def main(
    vcf_file,
    input_json_file,
    output_json_file,
    sample_id,
    pedigree_file,
    gen2phen_file,
    phenotype_filter_only=None,
    min_phenotype_score=None,
    require_gene_overlap=None,
    sv_types=None,
    pass_only=None,
    min_size=None,
    require_exon_overlap=None,
    config_path=None,
):
    """
    Merge structural variants with gene annotation and phenotype matching.

    Args:
        vcf_file: Path to structural variants VCF
        input_json_file: Path to existing Talos JSON results
        output_json_file: Path to output JSON with merged SVs
        sample_id: Sample ID to process
        pedigree_file: Path to pedigree file with patient HPO terms
        gen2phen_file: Path to genes_to_phenotype.txt file
        phenotype_filter_only: If True, only include SVs with phenotype matches (overrides config)
        min_phenotype_score: Minimum phenotype match score (0-100) (overrides config)
        require_gene_overlap: If True, only include SVs that overlap at least one gene (overrides config)
        sv_types: List of SV types to include (overrides config)
        pass_only: If True, only include SVs with FILTER=PASS (overrides config)
        min_size: Minimum SV size in bp (overrides config)
        require_exon_overlap: If True, only include SVs that overlap at least one exon (overrides config)
        config_path: Path to config.toml file (defaults to TALOS_CONFIG env or inputs/config.toml)
    """

    print('\n=== Starting SV Merge with Phenotype Matching ===\n')

    # Load configuration from TOML file
    sv_config = load_sv_filter_config(config_path)

    # Use command-line arguments if provided, otherwise use config values
    require_gene_overlap = (
        require_gene_overlap if require_gene_overlap is not None else sv_config.get('require_gene_overlap', False)
    )
    require_exon_overlap = (
        require_exon_overlap if require_exon_overlap is not None else sv_config.get('require_exon_overlap', False)
    )
    phenotype_filter_only = (
        phenotype_filter_only if phenotype_filter_only is not None else sv_config.get('require_phenotype_match', False)
    )
    min_phenotype_score = (
        min_phenotype_score if min_phenotype_score is not None else sv_config.get('min_phenotype_score', 0.0)
    )
    sv_types = sv_types if sv_types is not None else sv_config.get('allowed_sv_types', None)
    pass_only = pass_only if pass_only is not None else sv_config.get('require_pass_only', False)
    min_size = min_size if min_size is not None else sv_config.get('min_sv_size', 0)
    max_size = sv_config.get('max_sv_size', 0)

    print('\n=== Applied Filtering Parameters ===')
    print(f'  Gene overlap required: {require_gene_overlap}')
    print(f'  Exon overlap required: {require_exon_overlap}')
    print(f'  Phenotype match required: {phenotype_filter_only}')
    print(f'  Min phenotype score: {min_phenotype_score}%')
    print(f'  Allowed SV types: {sv_types if sv_types else "All"}')
    print(f'  PASS filter only: {pass_only}')
    print(f'  Min SV size: {min_size}bp')
    if max_size > 0:
        print(f'  Max SV size: {max_size}bp')
    print()

    # Load resources
    print('Loading annotation resources...')
    gene_annotations = load_gene_annotations()

    # Load exon annotations if required
    exon_annotations = {}
    if require_exon_overlap:
        exon_annotations = load_exon_annotations()

    gene_to_hpo = load_gene_to_hpo_mapping(gen2phen_file)
    patient_hpo = load_patient_hpo_terms(pedigree_file, sample_id)

    if not patient_hpo:
        print('WARNING: No patient HPO terms found. Phenotype matching will not be performed.')

    # Load existing JSON
    with open(input_json_file) as f:
        talos_data = json.load(f)

    # Parse VCF header to get sample names
    print('\nParsing VCF header...')
    sample_names = parse_vcf_header(vcf_file)
    if sample_names:
        print(f'Found {len(sample_names)} samples in VCF: {", ".join(sample_names)}')
        if sample_id not in sample_names:
            print(f"WARNING: Target sample '{sample_id}' not found in VCF header!")
    else:
        print('WARNING: Could not parse sample names from VCF header')

    # Parse VCF and filter for target sample
    print(f'\nParsing structural variants from VCF for sample {sample_id}...')
    all_svs = []
    with open(vcf_file) as f:
        for line in f:
            sv = parse_vcf_line(line, sample_names, target_sample_id=sample_id)
            if sv:
                all_svs.append(sv)

    print(f'Found {len(all_svs)} structural variants for sample {sample_id}')

    # Process each SV
    annotated_svs = []
    phenotype_matched_count = 0
    filtered_no_gene = 0
    filtered_no_exon = 0
    filtered_sv_type = 0
    filtered_quality = 0
    filtered_size = 0

    for sv in all_svs:
        # Apply SV type filter
        if sv_types and sv['svtype'] not in sv_types:
            filtered_sv_type += 1
            continue

        # Apply quality filter
        if pass_only and sv['filter'] != 'PASS':
            filtered_quality += 1
            continue

        # Apply size filters
        sv_size = abs(sv['svlen'])
        if min_size > 0 and sv_size < min_size:
            filtered_size += 1
            continue
        if max_size > 0 and sv_size > max_size:
            filtered_size += 1
            continue

        # Apply exon overlap filter
        if require_exon_overlap:
            has_exon_overlap = check_exon_overlap(sv['chrom'], sv['pos'], sv['end'], exon_annotations)
            if not has_exon_overlap:
                filtered_no_exon += 1
                continue

        # Find overlapping genes
        overlapping_genes = find_overlapping_genes(sv['chrom'], sv['pos'], sv['end'], gene_annotations)

        # Apply gene overlap filter
        if require_gene_overlap and not overlapping_genes:
            filtered_no_gene += 1
            continue

        # Perform phenotype matching
        gene_symbols = [g['gene_name'] for g in overlapping_genes]
        matched_genes, match_details = calculate_phenotype_match(gene_symbols, patient_hpo, gene_to_hpo)

        # Calculate phenotype match score
        phenotype_score = 0.0
        if match_details['genes_with_matches']:
            phenotype_score = max(g['match_percentage'] for g in match_details['genes_with_matches'])

        # Apply phenotype filtering if requested
        if phenotype_filter_only and not matched_genes:
            continue

        if phenotype_score < min_phenotype_score:
            continue

        # Track phenotype matches
        if matched_genes:
            phenotype_matched_count += 1

        # Create Talos variant with phenotype information
        talos_variant = create_talos_variant_with_phenotype(sv, sample_id, overlapping_genes, match_details)

        annotated_svs.append(talos_variant)

    print(f'\nProcessed {len(all_svs)} structural variants:')
    print(f'  - {len(annotated_svs)} passed filtering')
    print(f'  - {phenotype_matched_count} with phenotype matches')
    if filtered_sv_type > 0:
        print(f'  - {filtered_sv_type} filtered out (SV type)')
    if filtered_quality > 0:
        print(f'  - {filtered_quality} filtered out (quality/FILTER)')
    if filtered_size > 0:
        print(f'  - {filtered_size} filtered out (size < {min_size}bp)')
    if require_exon_overlap:
        print(f'  - {filtered_no_exon} filtered out (no exon overlap)')
    if require_gene_overlap:
        print(f'  - {filtered_no_gene} filtered out (no gene overlap)')

    # Add to JSON
    if sample_id in talos_data['results']:
        existing_variants = talos_data['results'][sample_id]['variants']
        existing_variants.extend(annotated_svs)

        print(f'\nAdded {len(annotated_svs)} structural variants to sample {sample_id}')
        print(f'Total variants for {sample_id}: {len(existing_variants)}')

        # Update metadata
        if 'metadata' in talos_data['results'][sample_id]:
            talos_data['results'][sample_id]['metadata']['present_in_sv'] = True
    else:
        print(f'ERROR: Sample {sample_id} not found in JSON results')
        return 1

    # Write output
    with open(output_json_file, 'w') as f:
        json.dump(talos_data, f, indent=2)

    print(f'\n=== Merged results written to {output_json_file} ===\n')

    # Print summary of phenotype-matched SVs
    if phenotype_matched_count > 0:
        print(f'\nPhenotype-matched structural variants ({phenotype_matched_count}):')
        for var in annotated_svs:
            if var['var_data']['info'].get('phenotype_match'):
                coords = var['var_data']['coordinates']
                info = var['var_data']['info']
                print(
                    f'  {coords["chrom"]}:{coords["pos"]}-{info["end_pos"]} ({info["svtype"]}, {abs(info["svlen"])}bp)'
                )
                print(f'    Gene: {info.get("gene_symbol", "unknown")}')
                print(f'    Matched genes: {", ".join(info.get("phenotype_matched_genes", []))}')
                print(f'    HPO matches: {", ".join(info.get("phenotype_matching_terms", []))}')
                print(f'    Match score: {info.get("phenotype_match_score", 0):.1f}%')

    return 0


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Merge structural variants with gene annotation and phenotype matching',
        epilog='Note: All filter options can be configured in config.toml [StructuralVariantFiltering] section. '
        'Command-line arguments override config file settings.',
    )
    parser.add_argument('vcf_file', help='Path to structural variants VCF file')
    parser.add_argument('input_json', help='Path to input Talos JSON results')
    parser.add_argument('output_json', help='Path to output JSON file')
    parser.add_argument('sample_id', help='Sample ID to process')
    parser.add_argument('--pedigree', required=True, help='Path to pedigree file with HPO terms')
    parser.add_argument('--gen2phen', required=True, help='Path to genes_to_phenotype.txt file')
    parser.add_argument(
        '--config', dest='config_path', help='Path to config.toml file (default: $TALOS_CONFIG or inputs/config.toml)'
    )

    # Optional overrides (if not provided, uses config.toml settings)
    parser.add_argument(
        '--phenotype-only', action='store_true', default=None, help='Override: Only include SVs with phenotype matches'
    )
    parser.add_argument('--min-score', type=float, default=None, help='Override: Minimum phenotype match score (0-100)')
    parser.add_argument(
        '--require-gene-overlap',
        action='store_true',
        default=None,
        help='Override: Only include SVs that overlap at least one gene',
    )
    parser.add_argument(
        '--require-exon-overlap',
        action='store_true',
        default=None,
        help='Override: Only include SVs that overlap at least one exon',
    )
    parser.add_argument('--sv-types', nargs='+', default=None, help='Override: SV types to include (e.g., DEL CNV)')
    parser.add_argument(
        '--pass-only', action='store_true', default=None, help='Override: Only include SVs with FILTER=PASS'
    )
    parser.add_argument('--min-size', type=int, default=None, help='Override: Minimum SV size in bp')

    args = parser.parse_args()

    sys.exit(
        main(
            args.vcf_file,
            args.input_json,
            args.output_json,
            args.sample_id,
            args.pedigree,
            args.gen2phen,
            args.phenotype_only,
            args.min_score,
            args.require_gene_overlap,
            args.sv_types,
            args.pass_only,
            args.min_size,
            args.require_exon_overlap,
            args.config_path,
        )
    )
