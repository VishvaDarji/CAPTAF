#!/usr/bin/env python3
"""
Convert MusiteDeep results file to TSV format compatible with downstream filtering scripts.
Usage: python convert_musitedeep_to_tsv.py <input_results_file> <output_tsv_file>
"""

import sys
import pandas as pd
import os

def convert_musitedeep_to_tsv(input_file, output_file):
    """Convert MusiteDeep results to standardized TSV format with proper PTM parsing."""
    
    # Read the results file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Parse and convert to DataFrame format
    data = []
    
    for line in lines:
        line = line.strip()
        
        # Skip empty lines and error messages
        if not line or line.startswith('ERROR:') or line.startswith('LIBCIFPP'):
            continue
            
        # Skip protein headers (lines starting with >)
        if line.startswith('>'):
            continue
            
        # Skip the column header line
        if line.startswith('Position') and 'Residue' in line:
            continue
        
        # Parse data lines
        parts = line.split('\t')
        if len(parts) >= 3:
            try:
                position = int(parts[0].strip())
                residue = parts[1].strip()
                ptmscores_raw = parts[2].strip()  # This contains all PTM:score pairs
                cutoff_status = parts[3].strip() if len(parts) > 3 else 'None'
                
                # Parse multiple PTMs from PTMscores column
                # Example: "Phosphothreonine:0.106;O-linked_glycosylation:0.085"
                if ':' in ptmscores_raw:
                    # Split by semicolon to get individual PTM predictions
                    ptm_predictions = ptmscores_raw.split(';')
                    
                    for ptm_prediction in ptm_predictions:
                        if ':' in ptm_prediction:
                            ptm_type, score_str = ptm_prediction.split(':', 1)
                            try:
                                score = float(score_str)
                                
                                # Add separate row for each PTM prediction
                                data.append({
                                    'Position': position,
                                    'AA': residue,  # Use AA as the main column name
                                    'PTM_Type': ptm_type.strip(),
                                    'PTM_Score': score
                                })
                                
                            except ValueError:
                                print(f"Warning: Could not parse score from '{ptm_prediction}'")
                                continue
                        
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse line '{line}': {e}")
                continue
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Create DataFrame and save
    if data:
        df = pd.DataFrame(data)
        # Sort by position, then by PTM_Type for consistency
        df = df.sort_values(['Position', 'PTM_Type'])
        df.to_csv(output_file, sep='\t', index=False)
        print(f'Converted {len(data)} PTM predictions to TSV format')
        print(f'Output saved to: {output_file}')
        
        # Show preview of columns and data
        print(f'Columns: {list(df.columns)}')
        if len(df) > 0:
            print('Preview:')
            print(df.head(5).to_string())
    else:
        # Create empty file with proper headers if no data
        headers = ['Position', 'AA', 'PTM_Type', 'PTM_Score']
        empty_df = pd.DataFrame(columns=headers)
        empty_df.to_csv(output_file, sep='\t', index=False)
        print('No valid PTM predictions found - created empty file with headers')

def main():
    if len(sys.argv) != 3:
        print("Usage: python convert_musitedeep_to_tsv.py <input_results_file> <output_tsv_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(input_file):
        print(f"Error: Input file does not exist: {input_file}")
        sys.exit(1)
    
    print(f"Converting {input_file} to {output_file}")
    convert_musitedeep_to_tsv(input_file, output_file)

if __name__ == "__main__":
    main()
