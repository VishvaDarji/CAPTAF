#!/usr/bin/env python3
"""
Biotite-based DSSP prediction: Secondary structure and relative solvent accessibility.
Snakemake-compatible version for single PDB file processing.

Usage:
    python predict_dssp.py --pdb_file input.pdb --output_file output.tsv [--verbose]
"""

import biotite.structure as struc
import biotite.structure.io as strucio
import numpy as np
from pathlib import Path
import argparse
import sys
import urllib.request
import json

def download_alphafold_pdb(uniprot_id, output_path, verbose=False):
    """
    Download AlphaFold structure from AlphaFold DB.
    
    Parameters:
    -----------
    uniprot_id : str
        UniProt ID (e.g., "Q04917")
    output_path : str or Path
        Where to save the PDB file
    verbose : bool
        Print progress information
        
    Returns:
    --------
    bool : True if download successful, False otherwise
    """
    try:
        if verbose:
            print(f"  🔽 Downloading AlphaFold structure for {uniprot_id}...")
        
        # Query AlphaFold API to get prediction metadata
        api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
        
        if verbose:
            print(f"    Querying: {api_url}")
        
        with urllib.request.urlopen(api_url) as response:
            data = json.loads(response.read().decode())
        
        if not data:
            print(f"  ❌ No AlphaFold prediction found for {uniprot_id}", file=sys.stderr)
            return False
        
        # Get PDB download URL from API response
        pdb_url = data[0]["pdbUrl"]
        af_id = data[0]["entryId"]
        
        if verbose:
            print(f"    Found: {af_id}")
            print(f"    Downloading from: {pdb_url}")
        
        # Create parent directory if needed
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Download PDB file
        urllib.request.urlretrieve(pdb_url, output_path)
        
        if verbose:
            print(f"  ✅ Downloaded successfully to {output_path}")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Download failed for {uniprot_id}: {e}", file=sys.stderr)
        return False

def calculate_ss_rsa(pdb_file, output_file, verbose=False):
    """
    Calculate secondary structure and relative solvent accessibility for a PDB file.
    Downloads from AlphaFold DB if PDB file doesn't exist.
    
    Parameters:
    -----------
    pdb_file : str
        Path to the PDB file
    output_file : str
        Path to save results in TSV format
    verbose : bool
        Print detailed progress information
    """
    try:
        # Convert path to use forward slashes (cross-platform)
        pdb_file = pdb_file.replace('\\', '/')
        pdb_path = Path(pdb_file)
        
        # Check if PDB file exists, if not download it
        if not pdb_path.exists():
            if verbose:
                print(f"  ⚠️  PDB file not found: {pdb_file}")
            
            # Extract UniProt ID from filename
            uniprot_id = pdb_path.stem
            
            if verbose:
                print(f"  Attempting to download from AlphaFold DB...")
            
            # Try to download
            success = download_alphafold_pdb(uniprot_id, pdb_file, verbose)
            if not success:
                raise FileNotFoundError(f"Could not obtain PDB file for {uniprot_id}")
        
        # Load structure
        if verbose:
            print(f"  📂 Loading structure from {pdb_file}...")
        
        if pdb_path.stat().st_size == 0:
            raise ValueError(f"PDB file is empty: {pdb_file}")
        
        structure = strucio.load_structure(pdb_file)
        
        # Handle stacked structures (multiple models)
        if isinstance(structure, struc.AtomArrayStack):
            if verbose:
                print(f"  📊 Found {len(structure)} models, using first model...")
            structure = structure[0]
        
        if verbose:
            print(f"  🧬 Found {len(structure)} atoms")
        
        # 8-class secondary structure (B,E,G,H,I,S,T,C)
        if verbose:
            print(f"  🔄 Calculating secondary structure...")
        ss8 = struc.annotate_sse(structure)
        
        # Get CA atoms first
        ca_mask = structure.atom_name == 'CA'
        ca_atoms = structure[ca_mask]
        
        # Absolute ASA at CA positions
        if verbose:
            print(f"  💧 Calculating solvent accessible surface area...")
        asa = struc.sasa(structure)
        asa_ca = asa[ca_mask]
        res_ids = ca_atoms.res_id
        res_names = ca_atoms.res_name
        
        # Get secondary structure per residue (ss8 is per residue, not per atom)
        ss8_ca = struc.annotate_sse(ca_atoms)
        
        # RSA = ASA / max ASA (per residue type)
        # Standard maximum ASA values (Tien et al. 2013)
        max_asa = {
            'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0, 'CYS': 167.0,
            'GLN': 225.0, 'GLU': 214.0, 'GLY': 104.0, 'HIS': 216.0, 'ILE': 197.0,
            'LEU': 201.0, 'LYS': 230.0, 'MET': 224.0, 'PHE': 228.0, 'PRO': 154.0,
            'SER': 155.0, 'THR': 163.0, 'TRP': 281.0, 'TYR': 255.0, 'VAL': 165.0
        }
        max_vals = np.array([max_asa.get(rn, 165.0) for rn in res_names])
        rsa_vals = np.clip(asa_ca / max_vals, 0, 1)
        
        # Build TSV output
        if verbose:
            print(f"  📝 Building output...")
        
        output_lines = ["ResID\tAminoAcid\tSS\tRSA"]
        for rid, rn, ss, rsa_val in zip(res_ids, res_names, ss8_ca, rsa_vals):
            line = f"{rid}\t{rn}\t{ss}\t{rsa_val:.4f}"
            output_lines.append(line)
        
        output_text = "\n".join(output_lines)
        
        # Ensure output directory exists
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Save results
        with open(output_file, 'w') as f:
            f.write(output_text)
        
        if verbose:
            print(f"  ✅ Done! Total residues: {len(res_ids)}")
        
        return True
        
    except FileNotFoundError as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.stderr.flush()
        return False
    except ValueError as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.stderr.flush()
        return False
    except Exception as e:
        print(f"❌ Unexpected error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.stderr.flush()
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Predict secondary structure and RSA using Biotite"
    )
    parser.add_argument(
        "--pdb_file",
        required=True,
        help="Input PDB file path"
    )
    parser.add_argument(
        "--output_file",
        required=True,
        help="Output TSV file path"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed progress information"
    )
    
    args = parser.parse_args()
    
    success = calculate_ss_rsa(args.pdb_file, args.output_file, args.verbose)
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
