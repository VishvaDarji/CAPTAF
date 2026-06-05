#!/usr/bin/env python3
"""
Generate AlphaFold3 JSON input files from POI and target protein data with CCD code mapping.

Usage: python generate_af3_json.py --poi_fasta <file> --target_fasta <file> --poi_ptms <file> --target_ptms <file> --output <file>
"""

import argparse
import json
import pandas as pd
from pathlib import Path
import sys
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate AlphaFold3 JSON input files')
    parser.add_argument('--poi_fasta', required=True, help='POI protein FASTA file')
    parser.add_argument('--target_fasta', required=True, help='Target protein FASTA file')
    parser.add_argument('--poi_ptms', required=True, help='POI PTMs file (filtered)')
    parser.add_argument('--target_ptms', required=True, help='Target PTMs file (processed)')
    parser.add_argument('--output', required=True, help='Output JSON file')
    parser.add_argument('--verbose', action='store_true', help='Print detailed information')

    return parser.parse_args()

def read_fasta_sequence(fasta_file):
    """Read sequence from FASTA file."""
    try:
        with open(fasta_file, 'r') as f:
            lines = f.readlines()

        # Skip header line and concatenate sequence lines
        sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))

        # Also extract protein ID from header if available
        header = lines[0].strip() if lines else ""
        protein_id = header.split('|')[1] if '|' in header else Path(fasta_file).stem

        print(f"Read FASTA: {fasta_file}")
        print(f"   Protein ID: {protein_id}")
        print(f"   Sequence length: {len(sequence)} residues")
        print(f"   First 50 residues: {sequence[:50]}...")

        return sequence, protein_id

    except Exception as e:
        print(f"Error reading FASTA file {fasta_file}: {e}")
        sys.exit(1)

def load_ccd_mapping():
    """Load PTM to CCD code mapping."""
    mapping_file = "PTM_to_CCD_mapping.csv"

    if not Path(mapping_file).exists():
        print(f"Error: PTM to CCD mapping file not found: {mapping_file}")
        sys.exit(1)

    try:
        mapping_df = pd.read_csv(mapping_file, encoding='utf-8-sig')
        # NORMALIZE MAPPING COLUMNS 
        mapping_df.columns = mapping_df.columns.str.strip().str.replace('\ufeff', '').str.lower()
        
        # CLEAN THE DATA: remove extra whitespace, backticks, and quotes from all string columns
        for col in mapping_df.columns:
            if mapping_df[col].dtype == 'object':
                # Clean step by step to handle nested quotes and backticks
                mapping_df[col] = (mapping_df[col].astype(str)
                                 .str.strip()                    # Remove leading/trailing whitespace
                                 .str.replace('`', '', regex=False)  # Remove backticks
                                 .str.replace("'", '', regex=False)  # Remove single quotes
                                 .str.replace('"', '', regex=False)  # Remove double quotes
                                 .str.strip())                   # Strip again after cleaning
        
        print(f"Loaded CCD mapping: {mapping_file}")
        print(f"   Total mappings: {len(mapping_df)}")
        print(f"   Normalized columns: {list(mapping_df.columns)}")

        # Show sample of the CLEANED mapping data
        print(f"   Sample cleaned mapping data:")
        sample_data = mapping_df[['ptm_type', 'ccd']].head(3)
        for _, row in sample_data.iterrows():
            print(f"     '{row['ptm_type']}' -> CCD: '{row['ccd']}'")

        return mapping_df

    except Exception as e:
        print(f"Error loading CCD mapping file {mapping_file}: {e}")
        sys.exit(1)

def normalize_ptm_name(ptm_name):
    """Normalize PTM names for robust matching by removing common formatting variations."""
    if pd.isna(ptm_name):
        return ""
    
    normalized = str(ptm_name).strip()
    # Convert to lowercase
    normalized = normalized.lower()
    # Replace common separators with nothing to standardize
    normalized = normalized.replace('_', '').replace('-', '').replace(' ', '')
    # Remove common prefixes/suffixes that might vary
    normalized = normalized.replace('modification', '').replace('residue', '')
    
    return normalized

def get_ccd_code_for_target(ptm_type, aa, mapping_df, missing_ptms_list, protein_name="Unknown"):
    """Get CCD code for target PTM based on PTM_Type and AA (target residue in parentheses)."""
    # Check for special cases first
    ptm_lower = ptm_type.lower().strip()
    
    # Skip dephosphorylation
    if 'dephosphorylation' in ptm_lower:
        print(f"   Target: Skipping dephosphorylation PTM '{ptm_type}'")
        return None
    
    # Handle special cases
    if 'formylation' in ptm_lower:
        print(f"   Target: Special case - PTM '{ptm_type}' -> CCD 'FOR'")
        return 'FOR'
    
    if 'amidation' in ptm_lower:
        print(f"   Target: Special case - PTM '{ptm_type}' -> CCD 'NH2'")
        return 'NH2'
    
    # Normalize the input PTM type
    normalized_input_ptm = normalize_ptm_name(ptm_type)
    target_residue_format = f"({aa})"
    
    print(f"   Looking for target PTM: '{ptm_type}' (normalized: '{normalized_input_ptm}') + AA '{aa}'")

    # Try exact match first (case insensitive)
    exact_matches = mapping_df[
        (mapping_df['ptm_type'].str.lower() == ptm_type.strip().lower()) &
        (mapping_df['residue'].str.contains(target_residue_format, case=False, na=False, regex=False))
    ]
    
    if len(exact_matches) > 0:
        # Clean the CCD code thoroughly
        ccd_code = clean_ccd_code(exact_matches.iloc[0]['ccd'])
        print(f"   Target: EXACT match - PTM '{ptm_type}' + AA '{aa}' -> CCD '{ccd_code}'")
        return ccd_code
    
    # Try normalized fuzzy matching
    print(f"   No exact match found, trying normalized matching...")
    
    if 'normalized_ptm' not in mapping_df.columns:
        mapping_df['normalized_ptm'] = mapping_df['ptm_type'].apply(normalize_ptm_name)
    
    fuzzy_matches = mapping_df[
        (mapping_df['normalized_ptm'] == normalized_input_ptm) &
        (mapping_df['residue'].str.contains(target_residue_format, case=False, na=False, regex=False))
    ]
    
    if len(fuzzy_matches) > 0:
        # Clean the CCD code thoroughly
        ccd_code = clean_ccd_code(fuzzy_matches.iloc[0]['ccd'])
        original_ptm = fuzzy_matches.iloc[0]['ptm_type']
        print(f"   Target: NORMALIZED match - '{ptm_type}' matched to '{original_ptm}' + AA '{aa}' -> CCD '{ccd_code}'")
        return ccd_code
    
    # Try partial matching for cases like "Asymmetric dimethylarginine" vs "Asymmetric dimethylarginine (ADMA)"
    print(f"   No normalized match found, trying partial matching...")
    
    # Use first 8 characters for partial matching to avoid Series issues
    ptm_partial = ptm_type.strip().lower()[:8] if len(ptm_type.strip()) >= 8 else ptm_type.strip().lower()
    partial_matches = mapping_df[
        (mapping_df['ptm_type'].str.lower().str.contains(ptm_partial, case=False, na=False, regex=False)) &
        (mapping_df['residue'].str.contains(target_residue_format, case=False, na=False, regex=False))
    ]
    
    if len(partial_matches) > 0:
        # Clean the CCD code thoroughly
        ccd_code = clean_ccd_code(partial_matches.iloc[0]['ccd'])
        original_ptm = partial_matches.iloc[0]['ptm_type']
        print(f"   Target: PARTIAL match - '{ptm_type}' matched to '{original_ptm}' + AA '{aa}' -> CCD '{ccd_code}'")
        return ccd_code
    # No match found - add to missing PTMs list
    missing_ptms_list.append({
        'PTM_Type': ptm_type,
        'Amino_Acid': aa,
        'Protein_Type': f'Target_{protein_name}',
        'Context': f'{ptm_type} + {aa}'
    })
    print(f"   Warning: No CCD mapping found for target PTM '{ptm_type}' + AA '{aa}' - added to missing list")
    return None

def get_ccd_code_for_poi(ptm_type, mapping_df, missing_ptms_list, protein_name="Unknown"):
    """Get CCD code for POI PTM based on PTM_Type only."""
    # Check for special cases first
    ptm_lower = ptm_type.lower().strip()
    
    # Skip dephosphorylation
    if 'dephosphorylation' in ptm_lower:
        print(f"   POI: Skipping dephosphorylation PTM '{ptm_type}'")
        return None
    
    # Handle special cases
    if 'formylation' in ptm_lower:
        print(f"   POI: Special case - PTM '{ptm_type}' -> CCD 'FOR'")
        return 'FOR'
    
    if 'amidation' in ptm_lower:
        print(f"   POI: Special case - PTM '{ptm_type}' -> CCD 'NH2'")
        return 'NH2'
    
    # Normalize the input PTM type
    normalized_input_ptm = normalize_ptm_name(ptm_type)
    
    print(f"   Looking for POI PTM: '{ptm_type}' (normalized: '{normalized_input_ptm}')")
    
    # Try exact match first (case insensitive)
    exact_matches = mapping_df[mapping_df['ptm_type'].str.lower() == ptm_type.strip().lower()]
    
    if len(exact_matches) > 0:
        # Clean the CCD code thoroughly
        ccd_code = clean_ccd_code(exact_matches.iloc[0]['ccd'])
        print(f"   POI: EXACT match - PTM '{ptm_type}' -> CCD '{ccd_code}'")
        return ccd_code
    
    # Try normalized fuzzy matching
    print(f"   No exact match found, trying normalized matching...")
    
    if 'normalized_ptm' not in mapping_df.columns:
        mapping_df['normalized_ptm'] = mapping_df['ptm_type'].apply(normalize_ptm_name)
    
    fuzzy_matches = mapping_df[mapping_df['normalized_ptm'] == normalized_input_ptm]
    
    if len(fuzzy_matches) > 0:
        # Clean the CCD code thoroughly
        ccd_code = clean_ccd_code(fuzzy_matches.iloc[0]['ccd'])
        original_ptm = fuzzy_matches.iloc[0]['ptm_type']
        print(f"   POI: NORMALIZED match - '{ptm_type}' matched to '{original_ptm}' -> CCD '{ccd_code}'")
        return ccd_code
    
    # Try partial matching for cases like "Asymmetric dimethylarginine" vs "Asymmetric dimethylarginine (ADMA)"
    print(f"   No normalized match found, trying partial matching...")
    
    # Use first 8 characters for partial matching to avoid Series issues
    ptm_partial = ptm_type.strip().lower()[:8] if len(ptm_type.strip()) >= 8 else ptm_type.strip().lower()
    partial_matches = mapping_df[
        mapping_df['ptm_type'].str.lower().str.contains(ptm_partial, case=False, na=False, regex=False)
    ]
    
    if len(partial_matches) > 0:
        # Clean the CCD code thoroughly
        ccd_code = clean_ccd_code(partial_matches.iloc[0]['ccd'])
        original_ptm = partial_matches.iloc[0]['ptm_type']
        print(f"   POI: PARTIAL match - '{ptm_type}' matched to '{original_ptm}' -> CCD '{ccd_code}'")
        return ccd_code

    # No match found - add to missing PTMs list
    missing_ptms_list.append({
        'PTM_Type': ptm_type,
        'Amino_Acid': 'N/A',
        'Protein_Type': f'POI_{protein_name}',
        'Context': ptm_type
    })
    print(f"   Warning: No CCD mapping found for POI PTM '{ptm_type}' - added to missing list")
    return None

def save_missing_ptms(missing_ptms_list):
    """Save or append missing PTMs to CSV file for later analysis."""
    
    output_file = Path("missing_CCD_codes.csv")
    
    # Even if no new missing PTMs, ensure file exists with proper headers
    if not missing_ptms_list:
        if not output_file.exists():
            # Create empty file with headers
            empty_df = pd.DataFrame(columns=['PTM_Type', 'Amino_Acid', 'Protein_Type', 'Context'])
            empty_df.to_csv(output_file, index=False)
            print("   📝 Created empty missing_CCD_codes.csv (no missing PTMs yet)")
        else:
            print("   ℹ️  No new missing PTMs to add")
        return str(output_file)

    missing_df = pd.DataFrame(missing_ptms_list)
    
    # If file exists, load and merge intelligently
    if output_file.exists():
        try:
            existing_df = pd.read_csv(output_file)
            
            # Handle case where existing file has no data (just headers)
            if len(existing_df) == 0:
                combined_df = missing_df.copy()
            else:
                # Combine with existing data
                combined_df = pd.concat([existing_df, missing_df], ignore_index=True)
            
            # Remove duplicates based on PTM_Type, Amino_Acid, AND Protein_Type
            combined_df = combined_df.drop_duplicates(
                subset=['PTM_Type', 'Amino_Acid', 'Protein_Type'],
                keep='first'
            )
            
            # Sort for readability
            combined_df = combined_df.sort_values(['Protein_Type', 'PTM_Type', 'Amino_Acid'])
            combined_df.to_csv(output_file, index=False)
            
            new_unique = len(combined_df) - len(existing_df)
            print(f"\nMissing PTM report:")
            print(f"   Processed {len(missing_df)} PTM entries from this run")
            print(f"   Added {max(0, new_unique)} new unique missing PTMs")
            print(f"   Total unique missing PTMs: {len(combined_df)}")
        except Exception as e:
            print(f"   Warning: Could not merge with existing file: {e}")
            print(f"   Creating new file instead")
            missing_df = missing_df.drop_duplicates(subset=['PTM_Type', 'Amino_Acid', 'Protein_Type'])
            missing_df = missing_df.sort_values(['Protein_Type', 'PTM_Type', 'Amino_Acid'])
            missing_df.to_csv(output_file, index=False)
    else:
        # Create new file
        missing_df = missing_df.drop_duplicates(subset=['PTM_Type', 'Amino_Acid', 'Protein_Type'])
        missing_df = missing_df.sort_values(['Protein_Type', 'PTM_Type', 'Amino_Acid'])
        missing_df.to_csv(output_file, index=False)
        print(f"\nMissing PTM report:")
        print(f"   Created new file: {output_file}")
        print(f"   Saved {len(missing_df)} unique missing PTMs")

    # Show summary of protein types with missing PTMs
    if output_file.exists():
        try:
            final_df = pd.read_csv(output_file)
            if len(final_df) > 0:
                protein_type_counts = final_df['Protein_Type'].value_counts()
                print(f"\n   Missing PTMs by protein type:")
                for ptype, count in protein_type_counts.head(15).items():
                    print(f"     {ptype}: {count} missing PTMs")
                    
                # Show unique PTM types
                unique_ptm_types = final_df['PTM_Type'].unique()
                print(f"\n   Unique PTM types missing CCD codes ({len(unique_ptm_types)} total):")
                for ptm in unique_ptm_types[:10]:
                    print(f"     - {ptm}")
                if len(unique_ptm_types) > 10:
                    print(f"     ... and {len(unique_ptm_types) - 10} more")
        except Exception as e:
            print(f"   Could not generate summary: {e}")

    return str(output_file)

def clean_ccd_code(ccd_code):
    """Clean CCD code by removing quotes, backticks, and extra whitespace."""
    if pd.isna(ccd_code):
        return ""
    
    cleaned = str(ccd_code).strip()
    # Remove various quote types and backticks
    cleaned = cleaned.replace('`', '').replace("'", '').replace('"', '').strip()
    return cleaned

def read_target_ptm_data(ptm_file, mapping_df, missing_ptms_list, protein_type="Target"):
    """Read target PTM data from processed CSV file."""
    try:
        if not Path(ptm_file).exists():
            print(f"PTM file not found: {ptm_file}")
            return []

        ptm_df = pd.read_csv(ptm_file, encoding='utf-8-sig')
        ptm_df.columns = ptm_df.columns.str.strip().str.replace('\ufeff', '').str.lower()

        # EXTRACT PROTEIN NAME FROM FILENAME
        protein_name = Path(ptm_file).stem.replace('_processed', '').replace('_Target', '')
        
        print(f"Read target PTM file: {ptm_file}")
        print(f"   Protein name: {protein_name}")
        print(f"   Normalized columns: {list(ptm_df.columns)}")
        print(f"   Total PTMs: {len(ptm_df)}")

        if len(ptm_df) == 0:
            return []

        required_columns = ['position', 'ptm_type', 'aa']
        missing_columns = [col for col in required_columns if col not in ptm_df.columns]
        if missing_columns:
            print(f"   Error: Missing required columns: {missing_columns}")
            return []

        modifications = []
        skipped_count = 0

        for idx, row in ptm_df.iterrows():
            try:
                position = row['position']
                ptm_type = row['ptm_type']
                aa = row['aa']

                print(f"   Row {idx}: position={position}, ptm_type={ptm_type}, aa={aa}")

                if pd.notna(position) and pd.notna(ptm_type) and pd.notna(aa):
                    # PASS PROTEIN_NAME HERE
                    ccd_code = get_ccd_code_for_target(
                        str(ptm_type), str(aa), mapping_df, missing_ptms_list, protein_name
                    )

                    if ccd_code is None:
                        print(f"   Skipped PTM: Position {position}, Type '{ptm_type}' + AA '{aa}' (no CCD mapping)")
                        skipped_count += 1
                        continue

                    modifications.append({
                        "ptmType": ccd_code,
                        "ptmPosition": int(position)
                    })
                    print(f"   Added PTM: Position {position}, Type '{ptm_type}' + AA '{aa}' -> CCD '{ccd_code}'")

            except Exception as e:
                print(f"   Error processing row {idx}: {e}")
                continue

        print(f"   Generated {len(modifications)} target modifications")
        if skipped_count > 0:
            print(f"   Skipped {skipped_count} PTMs (dephosphorylation, missing mappings, etc.)")
        return modifications

    except Exception as e:
        print(f"Error reading target PTM file {ptm_file}: {e}")
        return []

def read_poi_ptm_data(ptm_file, mapping_df, missing_ptms_list, protein_type="POI"):
    """Read POI PTM data from filtered TSV file."""
    try:
        if not Path(ptm_file).exists():
            print(f"PTM file not found: {ptm_file}")
            return []

        ptm_df = pd.read_csv(ptm_file, sep='\t', encoding='utf-8-sig')
        ptm_df.columns = ptm_df.columns.str.strip().str.replace('\ufeff', '').str.lower()

        # EXTRACT PROTEIN NAME AND VARIANT FROM FILENAME
        # e.g., "P27348_v1_buried_conserved.tsv" -> protein_name="P27348", variant="v1_buried_conserved"
        filename_parts = Path(ptm_file).stem.split('_', 1)
        protein_name = filename_parts[0]
        variant_name = filename_parts[1] if len(filename_parts) > 1 else "unknown"
        protein_name_full = f"{protein_name}_{variant_name}"
        
        print(f"Read POI PTM file: {ptm_file}")
        print(f"   Protein name: {protein_name_full}")
        print(f"   Normalized columns: {list(ptm_df.columns)}")
        print(f"   Total PTMs: {len(ptm_df)}")

        if len(ptm_df) == 0:
            return []

        required_columns = ['position', 'ptm_type']
        missing_columns = [col for col in required_columns if col not in ptm_df.columns]
        if missing_columns:
            print(f"   Error: Missing required columns: {missing_columns}")
            return []

        modifications = []
        skipped_count = 0

        for idx, row in ptm_df.iterrows():
            try:
                position = row['position']
                ptm_type = row['ptm_type']

                print(f"   Row {idx}: position={position}, ptm_type={ptm_type}")

                if pd.notna(position) and pd.notna(ptm_type):
                    # PASS PROTEIN_NAME_FULL HERE
                    ccd_code = get_ccd_code_for_poi(
                        str(ptm_type), mapping_df, missing_ptms_list, protein_name_full
                    )

                    if ccd_code is None:
                        print(f"   Skipped PTM: Position {position}, Type '{ptm_type}' (no CCD mapping)")
                        skipped_count += 1
                        continue

                    modifications.append({
                        "ptmType": ccd_code,
                        "ptmPosition": int(position)
                    })
                    print(f"   Added PTM: Position {position}, Type '{ptm_type}' -> CCD '{ccd_code}'")

            except Exception as e:
                print(f"   Error processing row {idx}: {e}")
                continue

        print(f"   Generated {len(modifications)} POI modifications")
        if skipped_count > 0:
            print(f"   Skipped {skipped_count} PTMs (dephosphorylation, missing mappings, etc.)")
        return modifications

    except Exception as e:
        print(f"Error reading POI PTM file {ptm_file}: {e}")
        return []

def generate_af3_json(poi_fasta, target_fasta, poi_ptms, target_ptms, output_file, verbose=False):
    """Generate AlphaFold3 JSON input file."""

    print(f"Generating AlphaFold3 JSON")
    print(f"   POI FASTA: {poi_fasta}")
    print(f"   Target FASTA: {target_fasta}")
    print(f"   POI PTMs: {poi_ptms}")
    print(f"   Target PTMs: {target_ptms}")

    # Initialize missing PTMs list to track unmapped PTMs
    missing_ptms_list = []

    # Load CCD mapping
    mapping_df = load_ccd_mapping()

    # Read sequences
    poi_sequence, poi_id = read_fasta_sequence(poi_fasta)
    target_sequence, target_id = read_fasta_sequence(target_fasta)

    # GENERATE NEW FORMAT FILENAME: POI_Target.json
    base_dir = Path(output_file).parent
    new_filename = f"{poi_id}_{target_id}.json"
    final_output_file = Path(output_file)
    
    print(f"   Output: {final_output_file}")

    # Read PTM data with CCD mapping and track missing PTMs
    poi_modifications = read_poi_ptm_data(poi_ptms, mapping_df, missing_ptms_list, "POI")
    target_modifications = read_target_ptm_data(target_ptms, mapping_df, missing_ptms_list, "Target")

    print(f"\nModification Summary:")
    print(f"   POI modifications: {len(poi_modifications)}")
    print(f"   Target modifications: {len(target_modifications)}")

    # Save missing PTMs report
    save_missing_ptms(missing_ptms_list)

    # Validate modifications before JSON creation
    print(f"\nValidating modifications...")
    poi_modifications = validate_modifications(poi_modifications, "POI")
    target_modifications = validate_modifications(target_modifications, "Target")

    # Create AlphaFold3 JSON structure
    af3_json = {
        "name": f"POI-Target Complex: {poi_id} + {target_id}",
        "modelSeeds": [10, 42],
        "sequences": [
            {
                "protein": {
                    "id": "A",  # POI protein
                    "sequence": poi_sequence,
                    "modifications": poi_modifications
                }
            },
            {
                "protein": {
                    "id": "B",  # Target protein
                    "sequence": target_sequence,
                    "modifications": target_modifications
                }
            }
        ],
        "dialect": "alphafold3",
        "version": 1
    }

    # Save JSON file
    try:
        final_output_file.parent.mkdir(parents=True, exist_ok=True)

        # Test JSON serialization before writing
        try:
            json_string = json.dumps(af3_json, indent=2)
            print(f"   JSON validation: PASSED")
        except Exception as json_error:
            print(f"   JSON validation: FAILED - {json_error}")
            print(f"   Problematic modifications:")
            print(f"     POI: {poi_modifications}")
            print(f"     Target: {target_modifications}")
            return False

        with open(final_output_file, 'w') as f:
            f.write(json_string)

        print(f"AlphaFold3 JSON generated successfully")
        print(f"   Output file: {final_output_file}")
        print(f"   File size: {final_output_file.stat().st_size} bytes")

        if verbose:
            print(f"\nJSON Preview:")
            print(json_string[:1000] + "...")

        # Validation summary
        print(f"\nValidation Summary:")
        print(f"   POI sequence: {len(poi_sequence)} residues")
        print(f"   Target sequence: {len(target_sequence)} residues")
        print(f"   POI modifications: {len(poi_modifications)}")
        print(f"   Target modifications: {len(target_modifications)}")
        print(f"   CCD code mapping: APPLIED")
        print(f"   Missing PTMs tracked: {len(missing_ptms_list)} instances")

        return True

    except Exception as e:
        print(f"Error writing JSON file {final_output_file}: {e}")
        return False

def validate_modifications(modifications, protein_type):
    """Validate and clean modifications for JSON compatibility."""
    valid_modifications = []
    
    for i, mod in enumerate(modifications):
        try:
            # Check required fields
            if 'ptmType' not in mod or 'ptmPosition' not in mod:
                print(f"   Warning: {protein_type} modification {i} missing required fields, skipping")
                continue
            
            # Clean and validate ptmType
            ptm_type = str(mod['ptmType']).strip()
            if not ptm_type:
                print(f"   Warning: {protein_type} modification {i} has empty ptmType, skipping")
                continue
            
            # Validate ptmPosition
            try:
                position = int(mod['ptmPosition'])
                if position <= 0:
                    print(f"   Warning: {protein_type} modification {i} has invalid position {position}, skipping")
                    continue
            except (ValueError, TypeError):
                print(f"   Warning: {protein_type} modification {i} has non-numeric position, skipping")
                continue
            
            # Add validated modification
            valid_modifications.append({
                "ptmType": ptm_type,
                "ptmPosition": position
            })
            
        except Exception as e:
            print(f"   Error validating {protein_type} modification {i}: {e}")
            continue
    
    print(f"   {protein_type}: {len(valid_modifications)}/{len(modifications)} modifications validated")
    return valid_modifications

def main():
    args = parse_arguments()

    try:
        success = generate_af3_json(
            poi_fasta=args.poi_fasta,
            target_fasta=args.target_fasta,
            poi_ptms=args.poi_ptms,
            target_ptms=args.target_ptms,
            output_file=args.output,
            verbose=args.verbose
        )

        if success:
            print("\n🎉 JSON generation completed!")
            print("⚠️  Note: PTM types need CCD code mapping before AlphaFold3 execution")
            sys.exit(0)
        else:
            print("\n❌ JSON generation failed")
            sys.exit(1)
            
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
