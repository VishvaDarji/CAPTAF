#!/usr/bin/env python3
"""
Predict disorder for a SINGLE protein using IUPred2A.

Usage: python predict_disorder.py --fasta <file> --protein_id <id> --output_dir <dir>
"""

import sys
import subprocess
import argparse
from pathlib import Path
import pandas as pd
from Bio import SeqIO

DEFAULT_IUPRED_PATH = "iupred2a/iupred2a.py"
DEFAULT_THRESHOLD = 0.5

# Amino acid code conversion (1-letter to 3-letter)
AA_1_TO_3 = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
    'U': 'SEC', 'O': 'PYL', 'X': 'UNK'
}

def convert_aa_code(aa_1letter):
    """Convert 1-letter amino acid code to 3-letter code."""
    return AA_1_TO_3.get(aa_1letter.upper(), 'UNK')


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Predict disorder for single protein'
    )
    parser.add_argument('--fasta', required=True, help='Input FASTA file')
    parser.add_argument('--protein_id', required=True, help='Protein ID')
    parser.add_argument('--output_dir', default='disordered_regions_outputs', help='Output directory')
    parser.add_argument('--iupred_path', default=DEFAULT_IUPRED_PATH, help='Path to iupred2a.py')
    parser.add_argument('--mode', default='long', choices=['long', 'short', 'glob'], help='IUPred mode')
    parser.add_argument('--threshold', type=float, default=DEFAULT_THRESHOLD, help='Disorder threshold')
    parser.add_argument('--verbose', action='store_true', help='Print detailed progress')
    return parser.parse_args()


def run_iupred2a(fasta_file, iupred_path, mode='long'):
    """Run IUPred2A on a FASTA file."""
    print(f"🔮 Running IUPred2A in '{mode}' mode...")
    
    try:
        cmd = ['python', str(iupred_path), str(fasta_file), mode]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        disorder_data = []
        for line in result.stdout.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) >= 3:
                try:
                    pos = int(parts[0])
                    residue = parts[1]
                    score = float(parts[2])
                    disorder_data.append((pos, residue, score))
                except (ValueError, IndexError):
                    continue
        
        print(f"✅ IUPred2A completed: {len(disorder_data)} residues analyzed")
        return disorder_data
        
    except subprocess.CalledProcessError as e:
        print(f"❌ IUPred2A error: {e}")
        print(f"   stdout: {e.stdout}")
        print(f"   stderr: {e.stderr}")
        return None
    except Exception as e:
        print(f"❌ Error running IUPred2A: {e}")
        return None


def classify_disorder(score, threshold=0.5):
    """Classify residue as Ordered or Disordered."""
    return "Disordered" if score >= threshold else "Ordered"


def main():
    args = parse_arguments()
    
    print("="*60)
    print(f"🔮 Disorder Prediction: {args.protein_id}")
    print("="*60)
    print(f"Input FASTA: {args.fasta}")
    print(f"Output directory: {args.output_dir}")
    print(f"Mode: {args.mode}")
    print(f"Threshold: {args.threshold}")
    
    if not Path(args.iupred_path).exists():
        print(f"❌ IUPred2A not found: {args.iupred_path}")
        sys.exit(1)
    
    try:
        record = SeqIO.read(args.fasta, 'fasta')
        query_seq = str(record.seq)
        print(f"   Sequence length: {len(query_seq)} residues")
    except Exception as e:
        print(f"❌ Error reading FASTA: {e}")
        sys.exit(1)
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    disorder_data = run_iupred2a(args.fasta, args.iupred_path, args.mode)
    
    if not disorder_data:
        print(f"❌ Failed to get disorder predictions")
        sys.exit(1)
    
    if len(disorder_data) != len(query_seq):
        print(f"⚠️  Warning: Disorder data length ({len(disorder_data)}) != sequence length ({len(query_seq)})")
    
    output_data = []
    for pos, residue, score in disorder_data:
        status = classify_disorder(score, args.threshold)
        output_data.append({
            'residue_number': pos,
            'aa': convert_aa_code(residue),  # ← CONVERSION APPLIED HERE
            'iupred2_score': round(score, 4),
            'disorder_status': status
        })
    
    df = pd.DataFrame(output_data)
    
    output_file = output_dir / f"disorder_predictions_{args.protein_id}.tsv"
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"✅ Disorder predictions saved: {output_file}")
    
    n_disordered = (df['disorder_status'] == 'Disordered').sum()
    n_ordered = (df['disorder_status'] == 'Ordered').sum()
    pct_disordered = 100 * n_disordered / len(df)
    mean_score = df['iupred2_score'].mean()
    
    print(f"   Mean disorder score: {mean_score:.3f}")
    print(f"   Disordered residues (≥{args.threshold}): {n_disordered} ({pct_disordered:.1f}%)")
    print(f"   Ordered residues (<{args.threshold}): {n_ordered} ({100-pct_disordered:.1f}%)")
    
    if args.verbose:
        print(f"\n📋 Sample disorder predictions:")
        print(df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
