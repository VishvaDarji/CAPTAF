#!/usr/bin/env python3
"""
Generate positive control JSON file for AlphaFold3.
Positive Control = POI + Target with ALL validated PTMs from dbPTM.

Usage: python generate_positive_control.py --poi_fasta <file> --target_fasta <file> --poi_ptms <file> --target_ptms <file> --output <file>
"""

import argparse
import json
import pandas as pd
from pathlib import Path
import sys


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Generate positive control JSON file (all validated PTMs)'
    )
    parser.add_argument('--poi_fasta', required=True, help='POI protein FASTA file')
    parser.add_argument('--target_fasta', required=True, help='Target protein FASTA file')
    parser.add_argument('--poi_ptms', required=True, help='POI validated PTMs (processed)')
    parser.add_argument('--target_ptms', required=True, help='Target validated PTMs (processed)')
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

        return sequence, protein_id

    except Exception as e:
        print(f"❌ Error reading FASTA file {fasta_file}: {e}")
        sys.exit(1)


def load_ccd_mapping():
    """Load PTM to CCD code mapping."""
    mapping_file = "PTM_to_CCD_mapping.csv"
    
    if not Path(mapping_file).exists():
        print(f"❌ PTM to CCD mapping file not found: {mapping_file}")
        sys.exit(1)
    
    try:
        mapping_df = pd.read_csv(mapping_file, encoding='utf-8-sig')
        mapping_df.columns = mapping_df.columns.str.strip().str.replace('\ufeff', '').str.lower()
        
        for col in mapping_df.columns:
            if mapping_df[col].dtype == 'object':
                mapping_df[col] = (mapping_df[col].astype(str)
                                 .str.strip()
                                 .str.replace('`', '', regex=False)
                                 .str.replace("'", '', regex=False)
                                 .str.replace('"', '', regex=False)
                                 .str.strip())
        
        return mapping_df
        
    except Exception as e:
        print(f"❌ Error loading CCD mapping: {e}")
        sys.exit(1)


def get_ccd_code(ptm_type, aa, mapping_df, is_target=False):
    """Get CCD code for PTM."""
    
    ptm_lower = ptm_type.lower().strip()
    
    if 'dephosphorylation' in ptm_lower:
        return None
    if 'formylation' in ptm_lower:
        return 'FOR'
    if 'amidation' in ptm_lower:
        return 'NH2'
    
    if is_target and aa:
        target_residue_format = f"({aa})"
        exact_matches = mapping_df[
            (mapping_df['ptm_type'].str.lower() == ptm_type.strip().lower()) &
            (mapping_df['residue'].str.contains(target_residue_format, case=False, na=False, regex=False))
        ]
    else:
        exact_matches = mapping_df[mapping_df['ptm_type'].str.lower() == ptm_type.strip().lower()]
    
    if len(exact_matches) > 0:
        ccd_code = str(exact_matches.iloc[0]['ccd']).strip().replace('`', '').replace("'", '').replace('"', '')
        return ccd_code
    
    return None


def read_validated_ptms(ptm_file, mapping_df, protein_type="POI", protein_name="Unknown", missing_ptms_list=None, sequence_length=None):
    """Read validated PTMs from processed dbPTM file."""
    
    if missing_ptms_list is None:
        missing_ptms_list = []
    
    try:
        if not Path(ptm_file).exists():
            print(f"   ⚠️  PTM file not found: {ptm_file}")
            return [], missing_ptms_list
        
        ptm_df = pd.read_csv(ptm_file, encoding='utf-8-sig')
        ptm_df.columns = ptm_df.columns.str.strip().str.replace('\ufeff', '').str.lower()
        
        print(f"   Read {protein_type} PTMs: {len(ptm_df)} total")
        
        if len(ptm_df) == 0:
            return [], missing_ptms_list
        
        modifications = []
        skipped = 0
        
        for _, row in ptm_df.iterrows():
            try:
                position = row.get('position')
                ptm_type = row.get('ptm_type')
                aa = row.get('aa', None)
                
                if pd.notna(position) and pd.notna(ptm_type):
                    if sequence_length and int(position) > sequence_length:
                        print(f"   ⚠️  INVALID position {position} > sequence length {sequence_length}, skipping {ptm_type}")
                        skipped += 1
                        continue
                    
                    is_target = (protein_type == "Target")
                    ccd_code = get_ccd_code(str(ptm_type), str(aa) if aa else '', mapping_df, is_target)
                    
                    if ccd_code:
                        modifications.append({
                            "ptmType": ccd_code,
                            "ptmPosition": int(position)
                        })
                    else:
                        # Track missing CCD code with protein name
                        skipped += 1
                        missing_ptms_list.append({
                            'PTM_Type': str(ptm_type),
                            'Amino_Acid': str(aa) if pd.notna(aa) else 'N/A',
                            'Protein_Type': f'{protein_type}_PositiveControl_{protein_name}',
                            'Context': f'{ptm_type}' + (f' + {aa}' if pd.notna(aa) else '')
                        })
                        print(f"   ⚠️  Missing CCD code: {ptm_type}" + (f" + {aa}" if pd.notna(aa) else ""))
                        
            except Exception as e:
                print(f"   Warning: Error processing PTM: {e}")
                continue
        
        print(f"   {protein_type}: {len(modifications)} PTMs mapped, {skipped} skipped")
        return modifications, missing_ptms_list
        
    except Exception as e:
        print(f"   ❌ Error reading PTM file: {e}")
        return [], missing_ptms_list


def save_missing_ptms(missing_ptms_list):
    """Save or append missing PTMs to CSV file."""
    if not missing_ptms_list:
        return
    
    missing_df = pd.DataFrame(missing_ptms_list)
    output_file = Path("missing_CCD_codes.csv")
    
    # If file exists, append; otherwise create new
    if output_file.exists():
        existing_df = pd.read_csv(output_file)
        combined_df = pd.concat([existing_df, missing_df], ignore_index=True)
        # Remove duplicates
        combined_df = combined_df.drop_duplicates(subset=['PTM_Type', 'Amino_Acid', 'Protein_Type'])
        combined_df.to_csv(output_file, index=False)
        print(f"\n   📝 Appended {len(missing_df)} missing PTMs to {output_file}")
        print(f"   Total unique missing PTMs: {len(combined_df)}")
    else:
        missing_df.to_csv(output_file, index=False)
        print(f"\n   📝 Created {output_file} with {len(missing_df)} missing PTMs")


def generate_positive_control_json(poi_fasta, target_fasta, poi_ptms_file, 
                                   target_ptms_file, output_file, verbose=False):
    """Generate positive control JSON with ALL validated PTMs."""
    
    print(f"\n🔬 Generating positive control JSON")
    print(f"   POI FASTA: {poi_fasta}")
    print(f"   Target FASTA: {target_fasta}")
    print(f"   POI PTMs: {poi_ptms_file}")
    print(f"   Target PTMs: {target_ptms_file}")
    
    # Load CCD mapping
    mapping_df = load_ccd_mapping()
    
    # Read sequences
    poi_sequence, poi_id = read_fasta_sequence(poi_fasta)
    target_sequence, target_id = read_fasta_sequence(target_fasta)
    
    print(f"   POI: {poi_id} ({len(poi_sequence)} residues)")
    print(f"   Target: {target_id} ({len(target_sequence)} residues)")
    
    # Track missing PTMs
    missing_ptms_list = []
    
    # Read validated PTMs with protein names
    poi_modifications, missing_ptms_list = read_validated_ptms(
        poi_ptms_file, mapping_df, "POI", poi_id, missing_ptms_list, sequence_length=len(poi_sequence)
    )
    target_modifications, missing_ptms_list = read_validated_ptms(
        target_ptms_file, mapping_df, "Target", target_id, missing_ptms_list, sequence_length=len(target_sequence)
    )
    
    # Save missing PTMs to CSV
    if missing_ptms_list:
        save_missing_ptms(missing_ptms_list)
        print(f"   ⚠️  {len(missing_ptms_list)} PTMs skipped due to missing CCD codes")
    
    # Create JSON
    af3_json = {
        "name": f"Positive Control: {poi_id} + {target_id} (All Validated PTMs)",
        "modelSeeds": [10, 42],
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": poi_sequence,
                    "modifications": poi_modifications
                }
            },
            {
                "protein": {
                    "id": "B",
                    "sequence": target_sequence,
                    "modifications": target_modifications
                }
            }
        ],
        "dialect": "alphafold3",
        "version": 1
    }
    
    # Save JSON
    try:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(af3_json, f, indent=2)
        
        print(f"   ✅ Positive control JSON saved: {output_path}")
        print(f"   POI modifications: {len(poi_modifications)}")
        print(f"   Target modifications: {len(target_modifications)}")
        
        if verbose:
            print(f"\n   JSON Preview:")
            print(json.dumps(af3_json, indent=2)[:1000] + "...")
        
        return True
        
    except Exception as e:
        print(f"❌ Error writing JSON: {e}")
        return False


def main():
    args = parse_arguments()
    
    print("=" * 70)
    print("🎯 GENERATING POSITIVE CONTROL JSON")
    print("=" * 70)
    
    success = generate_positive_control_json(
        poi_fasta=args.poi_fasta,
        target_fasta=args.target_fasta,
        poi_ptms_file=args.poi_ptms,
        target_ptms_file=args.target_ptms,
        output_file=args.output,
        verbose=args.verbose
    )
    
    if success:
        print(f"\n✅ Positive control JSON generated successfully!")
        sys.exit(0)
    else:
        print(f"\n❌ Positive control JSON generation failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
