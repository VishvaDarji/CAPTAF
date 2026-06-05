#!/usr/bin/env python3
"""
Generate 8 PTM variants from master table using parallel filtering.

Variant Definitions:
v1: buried_conserved         - RSA<20%, SS=H/E, Disorder<0.5, Conservation>0.7
v2: buried_variable          - RSA<20%, SS=H/E, Disorder<0.5, Conservation<0.4
v3: interface_conserved      - RSA 20-50%, SS=H/E, Disorder<0.5, Conservation>0.7
v4: interface_variable       - RSA 20-50%, SS=H/E, Disorder<0.5, Conservation<0.4
v5: exposed_ordered          - RSA>50%, Disorder<0.5
v6: disordered_conserved     - Disorder>0.5, Conservation>0.7
v7: disordered_variable      - Disorder>0.5, Conservation<0.4
v8: high_confidence          - Score ≥0.1 OR phosphorylation (any score)

Each variant represents a different structural/evolutionary context for PTMs.
Uses parallel filtering: all filters applied to same master table independently.

Confidence Filtering (Applied to ALL variants):
- Remove PTMs with PTM_Score < 0.1 (low confidence predictions)
- EXCEPTION: Keep ALL phosphorylation PTMs regardless of score
  Rationale: Phosphorylation PTMs with low scores match validated dbPTM data,
  indicating they are biologically relevant despite lower prediction confidence.

Usage: python generate_variants.py --protein_id <ID> --master_table <FILE>
"""

import sys
import pandas as pd
import argparse
from pathlib import Path

# Variant definitions with filtering criteria
VARIANT_DEFINITIONS = {
    'v1_buried_conserved': {
        'description': 'PTMs in protein core, highly conserved',
        'filters': {
            'rsa': ('lt', 0.20),
            'ss': ('in', ['H', 'E']),
            'disorder': ('lt', 0.5),
            'conservation': ('gt', 0.7)
        }
    },
    'v2_buried_variable': {
        'description': 'PTMs in protein core, poorly conserved',
        'filters': {
            'rsa': ('lt', 0.20),
            'ss': ('in', ['H', 'E']),
            'disorder': ('lt', 0.5),
            'conservation': ('lt', 0.4)
        }
    },
    'v3_interface_conserved': {
        'description': 'PTMs at protein interfaces, highly conserved',
        'filters': {
            'rsa': ('between', 0.20, 0.50),
            'ss': ('in', ['H', 'E']),
            'disorder': ('lt', 0.5),
            'conservation': ('gt', 0.7)
        }
    },
    'v4_interface_variable': {
        'description': 'PTMs at protein interfaces, poorly conserved',
        'filters': {
            'rsa': ('between', 0.20, 0.50),
            'ss': ('in', ['H', 'E']),
            'disorder': ('lt', 0.5),
            'conservation': ('lt', 0.4)
        }
    },
    'v5_exposed_ordered': {
        'description': 'PTMs on protein surface, ordered structure',
        'filters': {
            'rsa': ('gt', 0.50),
            'disorder': ('lt', 0.5)
        }
    },
    'v6_disordered_conserved': {
        'description': 'PTMs in disordered regions, highly conserved',
        'filters': {
            'disorder': ('gt', 0.5),
            'conservation': ('gt', 0.7)
        }
    },
    'v7_disordered_variable': {
        'description': 'PTMs in disordered regions, poorly conserved',
        'filters': {
            'disorder': ('gt', 0.5),
            'conservation': ('lt', 0.4)
        }
    },
    'v8_high_confidence': {
        'description': 'High-confidence PTMs: Score ≥0.1 OR phosphorylation (any score)',
        'filters': {
            # v8 has NO structural filters - only confidence filter applied later
        }
    }
}


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Generate PTM variants using parallel filtering'
    )
    parser.add_argument('--protein_id', required=True,
                       help='Protein ID (e.g., P27348)')
    parser.add_argument('--master_table', required=True,
                       help='Input master table file')
    parser.add_argument('--output_dir', default='poi_variants',
                       help='Output directory for variant files')
    parser.add_argument('--verbose', action='store_true',
                       help='Print detailed progress')
    return parser.parse_args()


def load_master_table(master_table_file):
    """Load master table."""
    try:
        df = pd.read_csv(master_table_file, sep='\t')
        print(f"✅ Loaded master table: {master_table_file}")
        print(f"   Total PTM sites: {len(df)}")
        print(f"   Columns: {list(df.columns)}")
        return df
    except Exception as e:
        print(f"❌ Error loading master table: {e}")
        sys.exit(1)


def apply_filter(df, column, operation, value1, value2=None):
    """
    Apply a single filter to dataframe.

    Args:
        df: DataFrame to filter
        column: Column name to filter on
        operation: Filter operation ('lt', 'gt', 'between', 'in', 'gte')
        value1: First value for comparison
        value2: Second value for 'between' operation
    """
    if column not in df.columns:
        print(f"⚠️  Warning: Column '{column}' not found in master table")
        return pd.Series([True] * len(df))

    if operation == 'lt':
        return df[column] < value1
    elif operation == 'gt':
        return df[column] > value1
    elif operation == 'between':
        return (df[column] >= value1) & (df[column] <= value2)
    elif operation == 'in':
        return df[column].isin(value1)
    elif operation == 'gte':
        return df[column] >= value1
    else:
        print(f"⚠️  Warning: Unknown operation '{operation}'")
        return pd.Series([True] * len(df))


def generate_variant(master_df, variant_name, variant_def, verbose=False):
    """
    Generate a single variant by applying filters to master table.

    Returns filtered DataFrame for this variant.
    """
    print(f"\n🔬 Generating {variant_name}")
    print(f"   Description: {variant_def['description']}")

    # Start with all PTMs
    mask = pd.Series([True] * len(master_df))

    # Apply each filter
    filters = variant_def['filters']

    for filter_name, filter_spec in filters.items():
        operation = filter_spec[0]

        if filter_name == 'rsa':
            column = 'rsa'
        elif filter_name == 'ss':
            column = 'ss'
        elif filter_name == 'disorder':
            column = 'iupred2_score'
        elif filter_name == 'conservation':
            column = 'conservation_score'
        elif filter_name == 'ptm_score':
            column = 'PTM_Score'
        else:
            print(f"⚠️  Unknown filter: {filter_name}")
            continue

        if operation == 'between':
            filter_mask = apply_filter(master_df, column, operation,
                                      filter_spec[1], filter_spec[2])
            if verbose:
                print(f"   Filter {filter_name}: {filter_spec[1]} ≤ {column} ≤ {filter_spec[2]} → {filter_mask.sum()} PTMs")
        elif operation == 'in':
            filter_mask = apply_filter(master_df, column, operation, filter_spec[1])
            if verbose:
                print(f"   Filter {filter_name}: {column} in {filter_spec[1]} → {filter_mask.sum()} PTMs")
        else:
            filter_mask = apply_filter(master_df, column, operation, filter_spec[1])
            if verbose:
                op_symbol = '<' if operation == 'lt' else '>' if operation == 'gt' else '>='
                print(f"   Filter {filter_name}: {column} {op_symbol} {filter_spec[1]} → {filter_mask.sum()} PTMs")

        mask = mask & filter_mask

    # Apply combined mask
    variant_df = master_df[mask].copy()

    print(f"   ✅ Final {variant_name}: {len(variant_df)} PTMs")

    return variant_df


def apply_confidence_filter(variant_df, variant_name):
    """
    Apply PTM confidence filter: PTM_Score >= 0.1
    Exception: Keep phosphorylation PTMs regardless of score (they match validated data)
    """
    if 'PTM_Score' not in variant_df.columns or 'PTM_Type' not in variant_df.columns:
        print(f"   ⚠️  Warning: PTM_Score or PTM_Type column not found")
        return variant_df
    
    before = len(variant_df)
    
    # Create mask for PTMs to keep:
    # 1. PTM_Score >= 0.1 (high confidence)
    # 2. OR PTM_Type starts with "Phospho" (phosphorylation exception)
    keep_mask = (
        (variant_df['PTM_Score'] >= 0.1) |
        (variant_df['PTM_Type'].str.lower().str.startswith('phospho', na=False))
    )
    
    variant_df = variant_df[keep_mask].copy()
    removed = before - len(variant_df)
    
    # Count phospho PTMs with low scores that were kept
    if removed > 0:
        phospho_exception = (
            (variant_df['PTM_Score'] < 0.1) &
            (variant_df['PTM_Type'].str.lower().str.startswith('phospho', na=False))
        )
        n_phospho_exceptions = phospho_exception.sum()
        
        if n_phospho_exceptions > 0:
            print(f"   🔥 Removed {removed} low-confidence PTMs (PTM_Score < 0.1)")
            print(f"   ✅ Kept {n_phospho_exceptions} phosphorylation PTMs despite low scores (validated matches)")
        else:
            print(f"   🔥 Removed {removed} low-confidence PTMs (PTM_Score < 0.1)")
    
    return variant_df


def save_variant(variant_df, protein_id, variant_name, output_dir):
    """Save variant to TSV file."""
    output_file = Path(output_dir) / f"{protein_id}_{variant_name}.tsv"

    if len(variant_df) == 0:
        print(f"   ⚠️  Warning: {variant_name} has 0 PTMs - saving empty file")

    variant_df.to_csv(output_file, sep='\t', index=False)
    print(f"   💾 Saved: {output_file}")

    return output_file


def print_variant_summary(all_variants):
    """Print summary of all generated variants."""
    print("\n" + "="*60)
    print("📊 VARIANT GENERATION SUMMARY")
    print("="*60)

    total_ptms = 0
    for variant_name, variant_df in all_variants.items():
        n_ptms = len(variant_df)
        total_ptms += n_ptms
        print(f"{variant_name:30s}: {n_ptms:4d} PTMs")

    print(f"\n{'TOTAL (with overlaps)':30s}: {total_ptms:4d} PTM assignments")

    # Note about overlaps
    print("\n⚠️  Note: PTMs may appear in multiple variants (by design)")
    print("   Example: A highly conserved PTM at RSA=25% appears in both")
    print("            v3 (interface_conserved) AND possibly others")


def main():
    args = parse_arguments()

    print("="*60)
    print(f"🔬 Generating PTM Variants for: {args.protein_id}")
    print("="*60)
    print(f"Master table: {args.master_table}")
    print(f"Output directory: {args.output_dir}")
    print(f"\nGenerating {len(VARIANT_DEFINITIONS)} variants...")

    # Load master table
    master_df = load_master_table(args.master_table)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate all variants
    all_variants = {}

    for variant_name, variant_def in VARIANT_DEFINITIONS.items():
        variant_df = generate_variant(master_df, variant_name, variant_def, args.verbose)
        
        # Apply confidence filter with phosphorylation exception
        variant_df = apply_confidence_filter(variant_df, variant_name)
        
        all_variants[variant_name] = variant_df
        save_variant(variant_df, args.protein_id, variant_name, output_dir)

    # Print summary
    print_variant_summary(all_variants)

    print(f"\n✅ All variants generated successfully!")
    print(f"   Output directory: {output_dir}")
    print(f"   Files: {args.protein_id}_v1_buried_conserved.tsv ... v8_high_confidence.tsv")


if __name__ == "__main__":
    main()
