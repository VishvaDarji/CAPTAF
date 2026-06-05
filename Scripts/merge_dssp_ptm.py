#!/usr/bin/env python3
"""
Merge PTM predictions with DSSP structural data.
No filtering - just merge on residue position.

Usage:
    python merge_dssp_ptm.py ptm_file dssp_file output_file [--verbose]
"""

import pandas as pd
import sys
from pathlib import Path

def merge_ptm_dssp(ptm_file, dssp_file, output_file, verbose=False):
    """
    Merge PTM predictions with DSSP structural features.
    
    Parameters:
    -----------
    ptm_file : str
        Path to PTM predictions TSV
    dssp_file : str
        Path to DSSP predictions TSV
    output_file : str
        Path to save merged output
    verbose : bool
        Print progress information
    """
    try:
        # Load PTM file
        if verbose:
            print(f"✅ Loading PTM file: {ptm_file}")
        
        ptm_df = pd.read_csv(ptm_file, sep='\t')
        
        if verbose:
            print(f"   Columns: {list(ptm_df.columns)}")
            print(f"   Total PTM predictions: {len(ptm_df)}")
        
        # Load DSSP file
        if verbose:
            print(f"✅ Loaded DSSP file: {dssp_file}")
        
        dssp_df = pd.read_csv(dssp_file, sep='\t')
        
        if verbose:
            print(f"   Columns: {list(dssp_df.columns)}")
            print(f"   Total DSSP entries: {len(dssp_df)}")
        
        # Identify position columns
        # PTM uses "Position", DSSP uses "ResID"
        ptm_pos_col = None
        dssp_pos_col = None
        
        # Find position column in PTM file
        for col in ['Position', 'ResID', 'Pos', 'Residue']:
            if col in ptm_df.columns:
                ptm_pos_col = col
                break
        
        # Find position column in DSSP file
        for col in ['ResID', 'Position', 'Pos', 'Residue']:
            if col in dssp_df.columns:
                dssp_pos_col = col
                break
        
        if ptm_pos_col is None or dssp_pos_col is None:
            print(f"❌ Could not find position columns to merge on", file=sys.stderr)
            print(f"   PTM columns: {list(ptm_df.columns)}", file=sys.stderr)
            print(f"   DSSP columns: {list(dssp_df.columns)}", file=sys.stderr)
            return False
        
        if verbose:
            print(f"✅ Found position columns:")
            print(f"   PTM: {ptm_pos_col}")
            print(f"   DSSP: {dssp_pos_col}")
        
        # Merge on position
        if verbose:
            print(f"🔄 Merging on position...")
        
        # Get DSSP columns excluding AminoAcid only
        dssp_cols_to_use = [col for col in dssp_df.columns if col != 'AminoAcid']
        
        # Rename dssp_pos_col to Position if it's different
        dssp_temp = dssp_df[dssp_cols_to_use].copy()
        if dssp_pos_col != ptm_pos_col and dssp_pos_col in dssp_temp.columns:
            dssp_temp = dssp_temp.rename(columns={dssp_pos_col: ptm_pos_col})
        
        merged_df = pd.merge(
            ptm_df,
            dssp_temp,
            left_on=ptm_pos_col,
            right_on=ptm_pos_col,
            how='left'
        )
        
        if verbose:
            print(f"✅ Merged successfully!")
            print(f"   Total merged rows: {len(merged_df)}")
            print(f"   Columns: {list(merged_df.columns)}")
        
        # Save output
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        merged_df.to_csv(output_file, sep='\t', index=False)
        
        if verbose:
            print(f"✅ Saved merged file: {output_file}")
            print(f"   First few rows:")
            print(merged_df.head(3).to_string())
        
        return True
        
    except FileNotFoundError as e:
        print(f"❌ File not found: {e}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        return False

def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Merge PTM predictions with DSSP structural data"
    )
    parser.add_argument("ptm_file", help="PTM predictions TSV file")
    parser.add_argument("dssp_file", help="DSSP predictions TSV file")
    parser.add_argument("output_file", help="Output merged TSV file")
    parser.add_argument("--verbose", action="store_true", help="Print progress information")
    
    args = parser.parse_args()
    
    success = merge_ptm_dssp(args.ptm_file, args.dssp_file, args.output_file, args.verbose)
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()

