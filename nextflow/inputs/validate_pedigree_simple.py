#!/usr/bin/env python3
"""
Simple stand-alone validation of a generated pedigree file
(no dependency on the talos package).

Usage:
    python validate_pedigree_simple.py <path/to/pedigree.ped>
"""

import re
import sys

if len(sys.argv) != 2:
    print(__doc__)
    sys.exit(2)

pedigree_file = sys.argv[1]

errors = []
warnings = []
families = {}
all_samples = set()
probands = []

with open(pedigree_file) as f:
    for line_num, line in enumerate(f, 1):
        line = line.strip()
        if not line:
            continue

        parts = line.split('\t')

        # Check minimum required columns
        if len(parts) < 6:
            errors.append(f'Line {line_num}: Expected at least 6 columns, got {len(parts)}')
            continue

        family_id = parts[0]
        sample_id = parts[1]
        father_id = parts[2]
        mother_id = parts[3]
        sex = parts[4]
        affected = parts[5]
        hpo_terms = parts[6] if len(parts) > 6 else None

        # Validate family ID
        if family_id == '0':
            errors.append(f"Line {line_num}: Family ID cannot be '0'")

        # Validate sample ID
        if sample_id == '0':
            errors.append(f"Line {line_num}: Sample ID cannot be '0'")

        # Validate sex
        if sex not in ['0', '1', '2']:
            errors.append(f"Line {line_num}: Sex must be 0, 1, or 2, got '{sex}'")

        # Validate affected status
        if affected not in ['0', '1', '2']:
            errors.append(f"Line {line_num}: Affected status must be 0, 1, or 2, got '{affected}'")

        # Validate HPO terms if present
        if hpo_terms:
            hpo_list = hpo_terms.split(',')
            for hpo in hpo_list:
                if not re.match(r'^HP:\d+$', hpo.strip()):
                    errors.append(f"Line {line_num}: Invalid HPO term format '{hpo}'")

        # Track samples and families
        all_samples.add(sample_id)
        if family_id not in families:
            families[family_id] = []
        families[family_id].append(
            {
                'sample_id': sample_id,
                'father_id': father_id,
                'mother_id': mother_id,
                'sex': sex,
                'affected': affected,
                'hpo_terms': hpo_terms,
            }
        )

        # Track probands (affected individuals)
        if affected == '2':
            probands.append(sample_id)

# Check parent IDs reference existing samples or are '0'
for family_id, members in families.items():
    family_samples = {m['sample_id'] for m in members}
    for member in members:
        if member['father_id'] != '0' and member['father_id'] not in family_samples:
            warnings.append(f"Sample {member['sample_id']}: Father ID '{member['father_id']}' not found in family")
        if member['mother_id'] != '0' and member['mother_id'] not in family_samples:
            warnings.append(f"Sample {member['sample_id']}: Mother ID '{member['mother_id']}' not found in family")

# Print results
print('=' * 70)
print('PEDIGREE FILE VALIDATION RESULTS')
print('=' * 70)

if errors:
    print(f'\n✗ ERRORS FOUND ({len(errors)}):')
    for error in errors[:10]:  # Show first 10 errors
        print(f'  - {error}')
    if len(errors) > 10:
        print(f'  ... and {len(errors) - 10} more errors')
else:
    print('\n✓ No errors found!')

if warnings:
    print(f'\n⚠ WARNINGS ({len(warnings)}):')
    for warning in warnings[:10]:  # Show first 10 warnings
        print(f'  - {warning}')
    if len(warnings) > 10:
        print(f'  ... and {len(warnings) - 10} more warnings')
else:
    print('\n✓ No warnings!')

print('\n' + '=' * 70)
print('STATISTICS')
print('=' * 70)
print(f'Total samples:       {len(all_samples)}')
print(f'Total families:      {len(families)}')
print(f'Affected (probands): {len(probands)}')
print(f'Unaffected:          {len(all_samples) - len(probands)}')

print('\nFamily size distribution:')
family_sizes = {}
for family_id, members in families.items():
    size = len(members)
    family_sizes[size] = family_sizes.get(size, 0) + 1

for size in sorted(family_sizes.keys()):
    print(f'  {size} members: {family_sizes[size]} families')

print('\nSample families (first 5):')
for i, (family_id, members) in enumerate(list(families.items())[:5]):
    print(f'\n  Family: {family_id} ({len(members)} members)')
    for member in members:
        hpo_info = f' [{len(member["hpo_terms"].split(","))} HPO terms]' if member['hpo_terms'] else ''
        affected_str = 'Affected' if member['affected'] == '2' else 'Unaffected'
        sex_str = 'Male' if member['sex'] == '1' else 'Female' if member['sex'] == '2' else 'Unknown'
        print(f'    - {member["sample_id"]}: {sex_str}, {affected_str}{hpo_info}')

if not errors:
    print('\n' + '=' * 70)
    print('✓ PEDIGREE FILE IS VALID AND READY TO USE!')
    print('=' * 70)
else:
    print('\n' + '=' * 70)
    print('✗ PEDIGREE FILE HAS ERRORS - PLEASE FIX BEFORE USING')
    print('=' * 70)
