#!/usr/bin/env python3
"""
Generate a pedigree file from Cases.txt
"""

import re
import sys

def clean_hpo_terms(hpo_str):
    """Clean HPO terms by removing descriptions and extra formatting"""
    # Remove quotes
    hpo_str = hpo_str.strip().strip('"')

    # Split by comma or semicolon
    terms = re.split(r'[,;]', hpo_str)

    clean_terms = []
    for term in terms:
        term = term.strip()
        # Extract only HP:XXXXXXX pattern
        match = re.search(r'(HP:\d+)', term)
        if match:
            clean_terms.append(match.group(1))

    return ','.join(clean_terms) if clean_terms else ''

def parse_cases_file(filepath):
    """Parse Cases.txt and generate pedigree entries"""
    pedigree_lines = []

    with open(filepath, 'r') as f:
        # Read all content
        content = f.read()

    # Split into lines and skip header
    lines = content.strip().split('\n')[1:]

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if not line:
            i += 1
            continue

        # Handle multiline entries - continue reading until we have 4 tab-separated parts
        while '\t' in line and line.count('\t') < 3 and i + 1 < len(lines):
            i += 1
            next_line = lines[i].strip()
            if next_line and not next_line.startswith('HP:'):
                break
            line = line.rstrip() + ' ' + next_line

        parts = line.split('\t')
        if len(parts) < 3:
            i += 1
            continue

        vcf_id = parts[0].strip()
        gender = parts[1].strip().lower() if len(parts) > 1 else ''
        case_type = parts[2].strip().upper() if len(parts) > 2 else 'SINGLETON'
        hpo_raw = parts[3].strip() if len(parts) > 3 else ''

        # Clean HPO terms
        hpo_terms = clean_hpo_terms(hpo_raw)

        # Map gender: f->2 (female), m->1 (male)
        sex = '2' if gender == 'f' else '1' if gender == 'm' else '0'

        # Determine family structure based on case type
        # Use the VCF ID as both family and sample ID
        # Mark probands as affected (2), parents as unaffected (1)

        family_id = vcf_id  # Use VCF ID as family identifier

        if case_type == 'TRIO':
            # Proband (affected)
            father_id = f"{vcf_id}-FATHER"
            mother_id = f"{vcf_id}-MOTHER"
            proband_line = f"{family_id}\t{vcf_id}\t{father_id}\t{mother_id}\t{sex}\t2"
            if hpo_terms:
                proband_line += f"\t{hpo_terms}"
            pedigree_lines.append(proband_line)

            # Father (unaffected)
            pedigree_lines.append(f"{family_id}\t{father_id}\t0\t0\t1\t1")

            # Mother (unaffected)
            pedigree_lines.append(f"{family_id}\t{mother_id}\t0\t0\t2\t1")

        elif case_type == 'DUO':
            # Proband (affected)
            mother_id = f"{vcf_id}-MOTHER"
            proband_line = f"{family_id}\t{vcf_id}\t0\t{mother_id}\t{sex}\t2"
            if hpo_terms:
                proband_line += f"\t{hpo_terms}"
            pedigree_lines.append(proband_line)

            # Mother (unaffected)
            pedigree_lines.append(f"{family_id}\t{mother_id}\t0\t0\t2\t1")

        elif case_type in ['QUADA', 'QUADB']:
            # Proband (affected) with both parents and a sibling
            father_id = f"{vcf_id}-FATHER"
            mother_id = f"{vcf_id}-MOTHER"
            sibling_id = f"{vcf_id}-SIBLING"

            proband_line = f"{family_id}\t{vcf_id}\t{father_id}\t{mother_id}\t{sex}\t2"
            if hpo_terms:
                proband_line += f"\t{hpo_terms}"
            pedigree_lines.append(proband_line)

            # Father (unaffected)
            pedigree_lines.append(f"{family_id}\t{father_id}\t0\t0\t1\t1")

            # Mother (unaffected)
            pedigree_lines.append(f"{family_id}\t{mother_id}\t0\t0\t2\t1")

            # Sibling (unaffected, unknown sex)
            pedigree_lines.append(f"{family_id}\t{sibling_id}\t{father_id}\t{mother_id}\t0\t1")

        else:
            # SINGLETON or unknown - just the proband
            proband_line = f"{family_id}\t{vcf_id}\t0\t0\t{sex}\t2"
            if hpo_terms:
                proband_line += f"\t{hpo_terms}"
            pedigree_lines.append(proband_line)

        i += 1

    return pedigree_lines

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Generate a Talos-compatible pedigree file from a Cases.txt sheet.'
    )
    parser.add_argument(
        '--input', '-i', required=True,
        help='Path to Cases.txt (tab-separated: VCF_ID, Gender, CASE_TYPE, HPO).'
    )
    parser.add_argument(
        '--output', '-o', required=True,
        help='Path where the generated .ped file will be written.'
    )
    args = parser.parse_args()

    pedigree_entries = parse_cases_file(args.input)

    with open(args.output, 'w') as f:
        for entry in pedigree_entries:
            f.write(entry + '\n')

    print(f"Generated pedigree file with {len(pedigree_entries)} entries")
    print(f"Output: {args.output}")
