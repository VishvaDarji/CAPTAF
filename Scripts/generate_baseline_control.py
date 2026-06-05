#!/usr/bin/env python3
"""
Generate baseline control JSON file for AlphaFold3.
Baseline = POI + Target with NO PTMs (unmodified proteins).

Usage: python generate_baseline_control.py --poi_fasta <file> --target_fasta <file> --output <file>
"""

import argparse
import json
from pathlib import Path
import sys


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Generate baseline control JSON file (no PTMs)'
    )
    parser.add_argument('--poi_fasta', required=True, help='POI protein FASTA file')
    parser.add_argument('--target_fasta', required=True, help='Target protein FASTA file')
    parser.add_argument('--output', required=True, help='Output JSON file')
    parser.add_argument('--verbose', action='store_true', help='Print detailed information')
    return parser.parse_args()


def read_fasta_sequence(fasta_file):
    """Read sequence from FASTA file."""
    try:
        with open(fasta_file, 'r') as f:
            lines = f.readlines()

        sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
        
        header = lines[0].strip() if lines else ""
        if '|' in header:
            protein_id = header.split('|')[1]
        else:
            protein_id = Path(fasta_file).stem

        print(f"   Read FASTA: {fasta_file}")
        print(f"   Protein ID: {protein_id}")
        print(f"   Sequence length: {len(sequence)} residues")

        return sequence, protein_id

    except Exception as e:
        print(f"❌ Error reading FASTA file {fasta_file}: {e}")
        sys.exit(1)


def generate_baseline_json(poi_fasta, target_fasta, output_file, verbose=False):
    """Generate baseline AlphaFold3 JSON with NO PTMs."""
    
    print(f"\n🔬 Generating baseline control JSON")
    print(f"   POI FASTA: {poi_fasta}")
    print(f"   Target FASTA: {target_fasta}")
    
    # Read sequences
    poi_sequence, poi_id = read_fasta_sequence(poi_fasta)
    target_sequence, target_id = read_fasta_sequence(target_fasta)
    
    # Create AlphaFold3 JSON structure WITHOUT modifications
    af3_json = {
        "name": f"Baseline Control: {poi_id} + {target_id} (No PTMs)",
        "modelSeeds": [10, 42],
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": poi_sequence
                }
            },
            {
                "protein": {
                    "id": "B",
                    "sequence": target_sequence
                }
            }
        ],
        "dialect": "alphafold3",
        "version": 1
    }
    
    # Save JSON file
    try:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(af3_json, f, indent=2)
        
        print(f"   ✅ Baseline JSON saved: {output_path}")
        print(f"   POI: {poi_id} ({len(poi_sequence)} residues)")
        print(f"   Target: {target_id} ({len(target_sequence)} residues)")
        print(f"   Modifications: NONE (baseline)")
        
        if verbose:
            print(f"\n   JSON Preview:")
            print(json.dumps(af3_json, indent=2))
        
        return True
        
    except Exception as e:
        print(f"❌ Error writing JSON file {output_file}: {e}")
        return False


def main():
    args = parse_arguments()
    
    print("=" * 70)
    print("🎯 GENERATING BASELINE CONTROL JSON")
    print("=" * 70)
    
    success = generate_baseline_json(
        poi_fasta=args.poi_fasta,
        target_fasta=args.target_fasta,
        output_file=args.output,
        verbose=args.verbose
    )
    
    if success:
        print(f"\n✅ Baseline control JSON generated successfully!")
        sys.exit(0)
    else:
        print(f"\n❌ Baseline control JSON generation failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
