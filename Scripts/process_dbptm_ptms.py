#!/usr/bin/env python3
"""
Process dbPTM target files by adding AA column and removing unnecessary columns.

This script:
1. Reads the target FASTA sequence
2. Adds an "AA" column with the amino acid at each PTM position
3. Removes unnecessary columns: Source, Substrate, Location_Raw
4. Updates the dbPTM CSV files in place

Usage: python process_dbptm_ptms.py <target_fasta> <dbptm_csv>
"""

import argparse
import pandas as pd
import sys
from pathlib import Path
from urllib.request import urlopen
import re

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process dbPTM target files')
    parser.add_argument('target_fasta', help='Target protein FASTA file')
    parser.add_argument('dbptm_csv', help='dbPTM CSV file to process')
    parser.add_argument('--verbose', action='store_true', help='Print detailed information')
    
    return parser.parse_args()

def read_fasta_sequence(fasta_file):
    """Read protein sequence from FASTA file."""
    try:
        with open(fasta_file, 'r') as f:
            lines = f.readlines()
        
        # Skip header line and concatenate sequence lines
        sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
        
        # Extract protein ID from header
        header = lines[0].strip() if lines else ""
        protein_id = header.split('|')[1] if '|' in header else Path(fasta_file).stem
        
        print(f"✅ Read FASTA: {fasta_file}")
        print(f"   Protein ID: {protein_id}")
        print(f"   Sequence length: {len(sequence)} residues")
        print(f"   First 50 residues: {sequence[:50]}...")
        
        return sequence, protein_id
        
    except Exception as e:
        print(f"❌ Error reading FASTA file {fasta_file}: {e}")
        sys.exit(1)

def get_amino_acid_at_position(sequence, position):
    """Get amino acid at specified position (1-indexed)."""
    try:
        # Convert to 0-indexed
        pos_index = int(position) - 1
        
        if pos_index < 0 or pos_index >= len(sequence):
            print(f"⚠️  Position {position} is out of range (sequence length: {len(sequence)})")
            return 'X'  # Return X for unknown/invalid positions
        
        return sequence[pos_index]
        
    except (ValueError, TypeError):
        print(f"⚠️  Invalid position: {position}")
        return 'X'

# Add these imports at the top
from urllib.request import urlopen
import re

# Add this helper function after get_amino_acid_at_position():
def fetch_full_protein_sequence(uniprot_id):
    """Fetch full protein sequence from UniProt"""
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        data = urlopen(url, timeout=10).read().decode()
        lines = data.split('\n')[1:]
        return ''.join(lines)
    except Exception as e:
        print(f"⚠️  Could not fetch full sequence for {uniprot_id}: {e}")
        return None

def filter_ptms_if_fragment(df, target_sequence, protein_id, verbose=False):
    """
    Check if target is a fragment and filter PTMs accordingly.
    
    Returns: filtered DataFrame
    """
    if len(df) == 0:
        return df
    
    # Check if this looks like a fragment:
    # If max PTM position >> sequence length, it's a fragment
    max_ptm_pos = df['Position'].max()
    seq_len = len(target_sequence)
    
    # If max position is close to sequence length, it's the full protein
    if max_ptm_pos <= seq_len * 1.2:  # Allow 20% tolerance
        print(f"   ℹ️  Target appears to be full protein (max PTM pos: {max_ptm_pos}, seq len: {seq_len})")
        return df
    
    # This is a fragment! Need to filter
    print(f"   🔍 Fragment detected: seq len={seq_len}, max PTM pos={max_ptm_pos}")
    print(f"   Filtering PTMs to match fragment region...")
    
    # Fetch full protein sequence
    full_seq = fetch_full_protein_sequence(protein_id)
    if not full_seq:
        print(f"   ⚠️  Cannot filter without full sequence - keeping all PTMs")
        return df
    
    print(f"   Full protein: {len(full_seq)} aa")
    
    # Find where fragment occurs in full protein
    start_idx = full_seq.find(target_sequence)
    if start_idx == -1:
        print(f"   ❌ Fragment not found in full protein sequence!")
        print(f"   Fragment: {target_sequence}")
        print(f"   Keeping all PTMs (may cause position mismatches)")
        return df
    
    # Convert to 1-indexed positions
    start_pos = start_idx + 1
    end_pos = start_idx + len(target_sequence)
    
    print(f"   ✅ Fragment found at positions {start_pos}-{end_pos}")
    
    # Filter PTMs to this range
    df_filtered = df[(df['Position'] >= start_pos) & (df['Position'] <= end_pos)].copy()
    
    # Remap positions to fragment-relative (1-indexed)
    df_filtered['Position'] = df_filtered['Position'] - start_pos + 1
    
    print(f"   ✅ Filtered {len(df)} → {len(df_filtered)} PTMs")
    print(f"   Remapped positions to range 1-{len(target_sequence)}")
    
    if len(df_filtered) == 0:
        print(f"   ⚠️  WARNING: No PTMs found in fragment region!")
    
    return df_filtered

def process_dbptm_file(target_fasta, dbptm_csv, verbose=False):
    """Process dbPTM file by adding AA column and removing unnecessary columns."""
    
    print(f"🔬 Processing dbPTM file")
    print(f"   Target FASTA: {target_fasta}")
    print(f"   dbPTM CSV: {dbptm_csv}")
    
    # Check if files exist
    if not Path(target_fasta).exists():
        print(f"❌ Target FASTA file not found: {target_fasta}")
        return False
    
    if not Path(dbptm_csv).exists():
        print(f"❌ dbPTM CSV file not found: {dbptm_csv}")
        return False
    
    # Read target sequence
    sequence, protein_id = read_fasta_sequence(target_fasta)
    
    # Read dbPTM CSV
    try:
        df = pd.read_csv(dbptm_csv)
        print(f"✅ Read dbPTM CSV: {dbptm_csv}")
        print(f"   Total PTMs: {len(df)}")
        print(f"   Columns: {list(df.columns)}")
        
        
        # === Check if dbPTM is empty ===
        if len(df) == 0:
            print(f"\n{'='*60}")
            print(f"❌ ERROR: No PTMs found in dbPTM for this target!")
            print(f"{'='*60}")
            print(f"Target: {target_fasta}")
            print(f"dbPTM file: {dbptm_csv}")
            print(f"\nPossible reasons:")
            print(f"  1. Protein not present in dbPTM database")
            print(f"  2. No experimentally validated PTMs for this protein")
            print(f"  3. Network/scraping error during dbPTM fetch")
            print(f"\nFor validation peptides with KNOWN PTM sites:")
            print(f"  Manually create the target PTM file:")
            print(f"  \n  Example for phospho-Ser at position 5:")
            print(f"  cat > {dbptm_csv} << 'EOF'")
            print(f"  UniProt_ID,Position,AA,PTM_Type")
            print(f"  P04049,5,S,Phosphorylation")
            print(f"  EOF")
            print(f"\n{'='*60}")
            sys.exit(1)  # Stop pipeline with error

    except Exception as e:
        print(f"❌ Error reading dbPTM CSV {dbptm_csv}: {e}")
        return False
    

    # ===  Filter PTMs if target is a fragment ===
    df = filter_ptms_if_fragment(df, sequence, protein_id, verbose)
    
    # === Check if filtering removed all PTMs ===
    if len(df) == 0:
        print(f"\n{'='*60}")
        print(f"❌ ERROR: No PTMs found in target region after filtering!")
        print(f"{'='*60}")
        print(f"Target: {target_fasta}")
        print(f"Target sequence: {sequence}")
        print(f"\nThis peptide/fragment has no PTMs in dbPTM.")
        print(f"\nFor validation with KNOWN PTM sites:")
        print(f"  Manually create both target PTM files:")
        print(f"\n  # Raw dbPTM file:")
        print(f"  cat > {dbptm_csv.replace('_processed', '')} << 'EOF'")
        print(f"  UniProt_ID,Position,PTM_Type")
        print(f"  {protein_id},5,Phosphorylation")
        print(f"  EOF")
        print(f"\n  # Processed file:")
        print(f"  cat > {dbptm_csv} << 'EOF'")
        print(f"  UniProt_ID,Position,AA,PTM_Type")
        print(f"  {protein_id},5,S,Phosphorylation")
        print(f"  EOF")
        print(f"\n{'='*60}")
        sys.exit(1)
    
    # Show sample data before processing
    if verbose and len(df) > 0:
        print(f"\n📊 Sample data before processing:")
        print(df.head().to_string(index=False))
    
    # Add AA column
    print(f"\n🧬 Adding AA column...")
    aa_list = []
    position_issues = 0
    
    for idx, row in df.iterrows():
        position = row.get('Position', None)
        if pd.isna(position):
            aa = 'X'
            position_issues += 1
        else:
            aa = get_amino_acid_at_position(sequence, position)
            if aa == 'X':
                position_issues += 1
        
        aa_list.append(aa)
        
        if verbose:
            ptm_type = row.get('PTM_Type', 'Unknown')
            print(f"   Position {position}: {aa} ({ptm_type})")
    
    # Add the AA column to the dataframe
    df['AA'] = aa_list
    
    print(f"   ✅ Added AA column")
    if position_issues > 0:
        print(f"   ⚠️  {position_issues} positions had issues (marked as 'X')")
    
    # Remove unnecessary columns
    columns_to_remove = ['Source', 'Substrate', 'Location_Raw']
    columns_removed = []
    
    for col in columns_to_remove:
        if col in df.columns:
            df = df.drop(columns=[col])
            columns_removed.append(col)
    
    if columns_removed:
        print(f"   ✅ Removed columns: {', '.join(columns_removed)}")
    else:
        print(f"   ℹ️  No unnecessary columns found to remove")
    
    # Reorder columns to put AA after Position
    desired_order = ['UniProt_ID', 'Position', 'AA', 'PTM_Type']
    
    # Keep any additional columns that might exist
    remaining_cols = [col for col in df.columns if col not in desired_order]
    final_order = desired_order + remaining_cols
    
    # Only include columns that actually exist
    final_order = [col for col in final_order if col in df.columns]
    df = df[final_order]
    
    print(f"   ✅ Reordered columns: {final_order}")
    
    # Show sample data after processing
    if verbose and len(df) > 0:
        print(f"\n📊 Sample data after processing:")
        print(df.head().to_string(index=False))
    
    # Save the processed file
    try:
        df.to_csv(dbptm_csv, index=False)
        print(f"✅ Updated dbPTM file: {dbptm_csv}")
        print(f"   Final columns: {list(df.columns)}")
        print(f"   Total PTMs: {len(df)}")
        
        # Summary statistics
        if 'AA' in df.columns:
            aa_counts = df['AA'].value_counts()
            print(f"   AA distribution: {dict(aa_counts.head())}")
        
        if 'PTM_Type' in df.columns:
            ptm_counts = df['PTM_Type'].value_counts()
            print(f"   PTM types: {dict(ptm_counts.head())}")
        
        return True
        
    except Exception as e:
        print(f"❌ Error saving processed file {dbptm_csv}: {e}")
        return False

def main():
    args = parse_arguments()
    
    success = process_dbptm_file(
        target_fasta=args.target_fasta,
        dbptm_csv=args.dbptm_csv,
        verbose=args.verbose
    )
    
    if not success:
        print(f"\n❌ Processing failed!")
        sys.exit(1)
    
    print(f"\n🎉 Processing completed successfully!")

if __name__ == "__main__":
    main()
