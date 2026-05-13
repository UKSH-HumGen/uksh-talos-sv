#!/usr/bin/env python3
"""
Validate a generated pedigree file using Talos's pedigree parser.

Usage:
    python validate_pedigree.py <path/to/pedigree.ped>

This script assumes the `talos` package is importable in the current
environment (e.g. installed via `pip install -e .` from the repo root).
"""

import sys

from talos.pedigree_parser import PedigreeParser


def main(pedigree_file: str) -> int:
    try:
        pedigree = PedigreeParser(pedigree_file)

        print('Pedigree file validated successfully.')
        print(f'  Total participants: {len(pedigree.participants)}')
        print(f'  Total families:     {len(pedigree.by_family)}')
        print(f'  Affected members:   {len(pedigree.get_affected_members())}')
        print(f'  Unaffected members: {len(pedigree.participants) - len(pedigree.get_affected_members())}')

        print('\nSample family breakdown (first 5 families):')
        for family_id, members in list(pedigree.by_family.items())[:5]:
            print(f'  Family {family_id}: {len(members)} members')
            for member in members:
                print(f'    - {member.sample_id} (Sex: {member.sex}, Affected: {member.affected})')
                if member.hpo_terms:
                    print(f'      HPO terms: {len(member.hpo_terms)} terms')

        print('\nPedigree file is valid.')
        return 0

    except Exception as exc:
        print(f'Error validating pedigree file: {exc}')
        return 1


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
