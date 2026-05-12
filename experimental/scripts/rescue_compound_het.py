#!/usr/bin/env python3
"""
Rescue compound heterozygous variants (SNV + SV) for AR genes.

After SVs are merged, check if any genes have:
1. A heterozygous SV already in results
2. A heterozygous SNV that was filtered out by MOI validation

If both exist in the same AR gene, add the SNV back as a compound het.
"""

import json
import sys
import argparse
from pathlib import Path
from cyvcf2 import VCF
from cyvcf2 import VCF


def load_json(filepath):
    """Load JSON file."""
    with open(filepath, 'r') as f:
        return json.load(f)


def save_json(data, filepath):
    """Save JSON file."""
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)


def extract_genes_with_het_svs(sample_variants):
    """
    Extract genes that have heterozygous SVs in the current results.

    Returns:
        dict: {gene_name: [variant_ids]} for genes with het SVs
    """
    genes_with_het_svs = {}

    for variant in sample_variants:
        var_data = variant.get('var_data', {})
        info = var_data.get('info', {})

        # Check if this is a structural variant
        if info.get('svtype'):
            gene_symbol = var_data.get('gene_symbol') or info.get('gene_symbol')

            if gene_symbol:
                # Get genotype - could be in different places
                genotype = var_data.get('genotype') or info.get('genotype')

                # Check if het (0/1, 1/0, ./1, 1/.)
                if genotype and ('0/1' in genotype or '1/0' in genotype or
                                './1' in genotype or '1/.' in genotype):
                    if gene_symbol not in genes_with_het_svs:
                        genes_with_het_svs[gene_symbol] = []

                    var_id = info.get('var_link', f"{var_data['coordinates']['chrom']}-{var_data['coordinates']['pos']}")
                    genes_with_het_svs[gene_symbol].append(var_id)

    return genes_with_het_svs


def check_ar_gene(gene_symbol, panelapp_data):
    """
    Check if a gene has AR (autosomal recessive/biallelic) mode of inheritance in PanelApp.

    Returns:
        bool: True if gene has AR/Biallelic MOI
    """
    # Look up gene by symbol in the genes dict
    genes = panelapp_data.get('genes', {})

    for gene_id, gene_info in genes.items():
        if gene_info.get('symbol') == gene_symbol:
            moi = str(gene_info.get('moi', '')).upper()

            # Check for AR mode of inheritance (Biallelic, AR, Recessive)
            if 'BIALLELIC' in moi or 'AR' in moi or 'RECESSIVE' in moi:
                return True

    return False


def find_compound_het_snvs(sample_id, genes_with_het_svs, full_report, panelapp_data):
    """
    Find heterozygous SNVs that were filtered out but could form compound hets with SVs.

    Args:
        sample_id: Sample identifier
        genes_with_het_svs: Dict of genes with het SVs
        full_report: Full Talos report with all variants (including filtered)
        panelapp_data: PanelApp data for MOI checking

    Returns:
        list: SNVs to rescue as compound hets
    """
    rescued_snvs = []

    # Get all variants for this sample from full report (before final filtering)
    sample_data = full_report.get('results', {}).get(sample_id, {})
    all_variants = sample_data.get('variants', [])

    for variant in all_variants:
        var_data = variant.get('var_data', {})
        info = var_data.get('info', {})

        # Skip if this is an SV (we only want SNVs/indels here)
        if info.get('svtype'):
            continue

        gene_symbol = var_data.get('gene_symbol') or info.get('gene_symbol')

        # Check if this SNV is in a gene that has a het SV
        if gene_symbol in genes_with_het_svs:
            # Check if this is an AR gene
            if not check_ar_gene(gene_symbol, panelapp_data):
                continue

            # Get genotype
            genotype = var_data.get('genotype') or info.get('genotype')

            # Check if het SNV
            if genotype and ('0/1' in genotype or '1/0' in genotype):
                # This is a potential compound het!
                # Add metadata about compound het status
                variant['compound_het_partner'] = genes_with_het_svs[gene_symbol]
                variant['compound_het_type'] = 'SNV+SV'
                variant['rescued_reason'] = f'Compound het with SV in AR gene {gene_symbol}'

                # Update support categories
                if 'support_categories' not in variant:
                    variant['support_categories'] = []
                if 'Compound Het (SNV+SV)' not in variant['support_categories']:
                    variant['support_categories'].append('Compound Het (SNV+SV)')

                # Update categories dict
                if 'categories' not in variant:
                    variant['categories'] = {}
                variant['categories']['Compound Het (SNV+SV)'] = full_report.get('date', '2025-10-08')

                rescued_snvs.append(variant)

    return rescued_snvs


def rescue_compound_hets(results_json, full_report_json, panelapp_json, output_json):
    """
    Main function to rescue compound heterozygous SNVs.

    Args:
        results_json: Path to results JSON (after SV merge)
        full_report_json: Path to full report JSON (before final filtering)
        panelapp_json: Path to PanelApp JSON
        output_json: Path to output JSON
    """
    # Load data
    print("=== Compound Heterozygous Rescue ===")
    print(f"Loading results from: {results_json}")
    results = load_json(results_json)

    print(f"Loading full report from: {full_report_json}")
    full_report = load_json(full_report_json)

    print(f"Loading PanelApp data from: {panelapp_json}")
    panelapp_data = load_json(panelapp_json)

    # Process each sample
    total_rescued = 0

    for sample_id, sample_data in results.get('results', {}).items():
        print(f"\nProcessing sample: {sample_id}")

        current_variants = sample_data.get('variants', [])

        # Find genes with het SVs
        genes_with_het_svs = extract_genes_with_het_svs(current_variants)

        if genes_with_het_svs:
            print(f"  Found {len(genes_with_het_svs)} genes with heterozygous SVs:")
            for gene, sv_ids in genes_with_het_svs.items():
                print(f"    - {gene}: {len(sv_ids)} SV(s)")
        else:
            print("  No heterozygous SVs found")
            continue

        # Find compound het SNVs to rescue
        rescued_snvs = find_compound_het_snvs(
            sample_id,
            genes_with_het_svs,
            full_report,
            panelapp_data
        )

        if rescued_snvs:
            print(f"  Rescued {len(rescued_snvs)} compound het SNVs:")
            for snv in rescued_snvs:
                gene = snv['var_data'].get('gene_symbol', 'Unknown')
                pos = snv['var_data']['coordinates']['pos']
                chrom = snv['var_data']['coordinates']['chrom']
                print(f"    - {gene} {chrom}:{pos} (compound het with SV)")

            # Add rescued SNVs to results
            results['results'][sample_id]['variants'].extend(rescued_snvs)
            total_rescued += len(rescued_snvs)
        else:
            print("  No compound het SNVs to rescue")

    # Save results
    print(f"\n=== Summary ===")
    print(f"Total compound het SNVs rescued: {total_rescued}")
    print(f"Writing results to: {output_json}")
    save_json(results, output_json)
    print("=== Compound Het Rescue Complete ===")


def main():
    parser = argparse.ArgumentParser(
        description='Rescue compound heterozygous SNVs that pair with SVs in AR genes'
    )
    parser.add_argument(
        'results_json',
        help='Results JSON with SVs merged (from MergeStructuralVariants)'
    )
    parser.add_argument(
        'full_report_json',
        help='Full report JSON with all variants before MOI filtering (from HPOFlagging)'
    )
    parser.add_argument(
        'panelapp_json',
        help='PanelApp JSON with gene panels and MOI data'
    )
    parser.add_argument(
        'output_json',
        help='Output JSON with rescued compound hets'
    )

    args = parser.parse_args()

    rescue_compound_hets(
        args.results_json,
        args.full_report_json,
        args.panelapp_json,
        args.output_json
    )


if __name__ == '__main__':
    main()
