#!/usr/bin/env python3
"""
Create master table by merging PTM predictions with structural and evolutionary data.

This script merges:
1. PTM predictions (from MusiteDeep)
2. DSSP structural data (RSA and Secondary Structure)
3. Conservation scores (from MMseqs2 + JSD)
4. Disorder predictions (from IUPred2A)

Output: Master table with all metrics for each PTM site
This table will be used for parallel filtering to create variants.

Usage: python create_master_table.py --protein_id <ID> --output_dir <DIR>
"""

import sys
import pandas as pd
import argparse
from pathlib import Path

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Create master table merging PTM, DSSP, conservation, and disorder data'
    )
    parser.add_argument('--protein_id', required=True,
                       help='Protein ID (e.g., P27348)')
    parser.add_argument('--ptm_dssp_file', required=True,
                       help='Input PTM+DSSP merged file')
    parser.add_argument('--conservation_file', required=True,
                       help='Input conservation scores file')
    parser.add_argument('--disorder_file', required=True,
                       help='Input disorder predictions file')
    parser.add_argument('--output_dir', default='master_tables_poi',
                       help='Output directory for master tables')
    parser.add_argument('--verbose', action='store_true',
                       help='Print detailed progress')
    return parser.parse_args()


def load_ptm_dssp_data(ptm_dssp_file):
    """Load PTM predictions merged with DSSP data."""
    try:
        df = pd.read_csv(ptm_dssp_file, sep='\t')
        print(f"✅ Loaded PTM+DSSP file: {ptm_dssp_file}")
        print(f"   Columns: {list(df.columns)}")
        print(f"   Entries: {len(df)}")
        return df
    except Exception as e:
        print(f"❌ Error loading PTM+DSSP file: {e}")
        sys.exit(1)


def load_conservation_data(conservation_file):
    """Load conservation scores."""
    try:
        df = pd.read_csv(conservation_file, sep='\t')
        print(f"✅ Loaded conservation file: {conservation_file}")
        print(f"   Columns: {list(df.columns)}")
        print(f"   Entries: {len(df)}")
        return df
    except Exception as e:
        print(f"❌ Error loading conservation file: {e}")
        sys.exit(1)


def load_disorder_data(disorder_file):
    """Load disorder predictions."""
    try:
        df = pd.read_csv(disorder_file, sep='\t')
        print(f"✅ Loaded disorder file: {disorder_file}")
        print(f"   Columns: {list(df.columns)}")
        print(f"   Entries: {len(df)}")
        return df
    except Exception as e:
        print(f"❌ Error loading disorder file: {e}")
        sys.exit(1)


def merge_all_data(ptm_dssp_df, conservation_df, disorder_df, verbose=False):
    """
    Merge all datasets on residue position.

    Returns master table with columns:
    - Position, AA, PTM_Type, PTM_Score (from PTM predictions)
    - SS, RSA (from DSSP)
    - conservation_score, conservation_grade (from conservation)
    - iupred2_score, disorder_status (from disorder)
    """
    print(f"\n🔗 Merging all datasets...")

    # Step 1: Merge PTM+DSSP with conservation
    # Match on Position (PTM) with residue_number (conservation)
    position_col = 'Position' if 'Position' in ptm_dssp_df.columns else 'position'

    merged_df = pd.merge(
        ptm_dssp_df,
        conservation_df[['residue_number', 'conservation_score', 'conservation_grade']],
        left_on=position_col,
        right_on='residue_number',
        how='left'  # Keep all PTMs even if no conservation data
    )

    print(f"   After conservation merge: {len(merged_df)} entries")

    # Step 2: Merge with disorder predictions
    merged_df = pd.merge(
        merged_df,
        disorder_df[['residue_number', 'iupred2_score', 'disorder_status']],
        left_on=position_col,
        right_on='residue_number',
        how='left',  # Keep all PTMs even if no disorder data
        suffixes=('', '_disorder')
    )

    print(f"   After disorder merge: {len(merged_df)} entries")

    # Step 3: Clean up redundant columns
    # Remove duplicate residue_number columns from merges
    cols_to_remove = ['residue_number', 'residue_number_disorder']
    for col in cols_to_remove:
        if col in merged_df.columns:
            merged_df = merged_df.drop(columns=[col])

    # Step 4: Handle missing values
    # Fill NaN conservation scores with 0 (assume variable if no data)
    if 'conservation_score' in merged_df.columns:
        merged_df = merged_df.assign(conservation_score=merged_df['conservation_score'].fillna(0.0))
        merged_df = merged_df.assign(conservation_grade=merged_df['conservation_grade'].fillna(1))

    # Fill NaN disorder scores with 0.5 (neutral if no data)
    if 'iupred2_score' in merged_df.columns:
        merged_df = merged_df.assign(iupred2_score=merged_df['iupred2_score'].fillna(0.5))
        merged_df = merged_df.assign(disorder_status=merged_df['disorder_status'].fillna('Unknown'))

    print(f"\n✅ Master table created with {len(merged_df)} PTM sites")

    if verbose:
        print(f"\n📋 Master table columns:")
        for col in merged_df.columns:
            print(f"   - {col}")

        print(f"\n📊 Data completeness:")
        print(f"   PTMs with RSA data: {merged_df['RSA'].notna().sum()}")
        print(f"   PTMs with SS data: {merged_df['SS'].notna().sum()}")
        print(f"   PTMs with conservation: {(merged_df['conservation_score'] > 0).sum()}")
        print(f"   PTMs with disorder data: {merged_df['iupred2_score'].notna().sum()}")

    return merged_df


def validate_master_table(master_df):
    """Validate that master table has all required columns."""
    required_cols = ['Position', 'AA', 'PTM_Type', 'PTM_Score', 'RSA', 'SS',
                     'conservation_score', 'conservation_grade', 'iupred2_score', 'disorder_status']

    missing_cols = [col for col in required_cols if col not in master_df.columns]

    if missing_cols:
        print(f"⚠️  Warning: Missing columns in master table: {missing_cols}")
        print(f"   Available columns: {list(master_df.columns)}")
        return False

    print(f"✅ Master table validation passed")
    return True


def remove_unwanted_columns(master_df):
    """Remove ResID and AminoAcid columns if they exist."""
    cols_to_remove = ['ResID', 'AminoAcid']
    for col in cols_to_remove:
        if col in master_df.columns:
            master_df = master_df.drop(columns=[col])
            print(f"   Removed column: {col}")
    return master_df


def main():
    args = parse_arguments()

    print("="*60)
    print(f"🧬 Creating Master Table for Protein: {args.protein_id}")
    print("="*60)
    print(f"Input files:")
    print(f"  PTM+DSSP: {args.ptm_dssp_file}")
    print(f"  Conservation: {args.conservation_file}")
    print(f"  Disorder: {args.disorder_file}")
    print(f"Output directory: {args.output_dir}")

    # Load all data
    ptm_dssp_df = load_ptm_dssp_data(args.ptm_dssp_file)
    conservation_df = load_conservation_data(args.conservation_file)
    disorder_df = load_disorder_data(args.disorder_file)

    # Merge all data
    master_df = merge_all_data(ptm_dssp_df, conservation_df, disorder_df, args.verbose)

    # Validate
    validate_master_table(master_df)

    # Save master table
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"master_table_{args.protein_id}.tsv"
    master_df.to_csv(output_file, sep='\t', index=False)

    print(f"\n✅ Master table saved: {output_file}")
    print(f"   Total PTM sites: {len(master_df)}")
    print(f"   Total columns: {len(master_df.columns)}")

    # Summary statistics
    if args.verbose:
        print(f"\n📊 Summary Statistics:")
        print(f"   Mean RSA: {master_df['RSA'].mean():.3f}")
        print(f"   Mean conservation: {master_df['conservation_score'].mean():.3f}")
        print(f"   Mean disorder: {master_df['iupred2_score'].mean():.3f}")

        print(f"\n   Secondary structure distribution:")
        print(master_df['SS'].value_counts())

        print(f"\n   Disorder status distribution:")
        print(master_df['disorder_status'].value_counts())

        print(f"\n📋 Sample of master table:")
        display_cols = ['Position', 'PTM_Type', 'RSA', 'SS',
                       'conservation_score', 'iupred2_score', 'disorder_status']
        available_cols = [col for col in display_cols if col in master_df.columns]
        print(master_df[available_cols].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
