import biotite.structure as struc
import biotite.structure.io as strucio
import numpy as np
from pathlib import Path
import urllib.request
import json
import sys

def calculate_ss_rsa(pdb_file, output_file):
    """
    Calculate secondary structure and relative solvent accessibility for a PDB file.
    Uses Biotite's built-in SSE annotation and SASA calculation.
    """
    try:
        # Convert path to use forward slashes (cross-platform)
        pdb_file = pdb_file.replace('\\', '/')

        # Load structure
        print(f"  Loading structure from {pdb_file}...")

        pdb_path = Path(pdb_file)
        if not pdb_path.exists():
            raise FileNotFoundError(f"File {pdb_file} does not exist")
        if pdb_path.stat().st_size == 0:
            raise ValueError(f"PDB file is empty: {pdb_file}")

        structure = strucio.load_structure(pdb_file)
        if isinstance(structure, struc.AtomArrayStack):
            print(f"  Found {len(structure)} models, using first model...")
            structure = structure[0]

        print(f"  Found {len(structure)} atoms")

        # 8-class secondary structure (B,E,G,H,I,S,T,C)
        print(f"  Calculating secondary structure...")
        ss8 = struc.annotate_sse(structure)               # 1 value / residue

        # Absolute ASA at CA positions
        print(f"  Calculating solvent accessible surface area...")
        ca_mask   = structure.atom_name == 'CA'
        asa_ca    = struc.sasa(structure)[ca_mask]        # ASA at CA atoms
        res_ids   = structure.res_id[ca_mask]
        res_names = structure.res_name[ca_mask]

        # RSA = ASA / max ASA (per residue)
        max_asa = {
            'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0, 'CYS': 167.0,
            'GLN': 225.0, 'GLU': 214.0, 'GLY': 104.0, 'HIS': 216.0, 'ILE': 197.0,
            'LEU': 201.0, 'LYS': 230.0, 'MET': 224.0, 'PHE': 228.0, 'PRO': 154.0,
            'SER': 155.0, 'THR': 163.0, 'TRP': 281.0, 'TYR': 255.0, 'VAL': 165.0
        }
        max_vals = np.array([max_asa.get(rn, 165.0) for rn in res_names])
        rsa_vals = np.clip(asa_ca / max_vals, 0, 1)

        # Build TSV output
        ss8_names = {
            'B': 'Beta-bridge', 'E': 'Beta-sheet', 'G': '3-helix',
            'H': 'Alpha-helix', 'I': 'Pi-helix',   'S': 'Bend',
            'T': 'Turn',        'C': 'Coil'
        }
        output_lines = ["ResID\tAminoAcid\tSS\tRSA"]
        for rid, rn, ss, rsa_val in zip(res_ids, res_names, ss8, rsa_vals):
            line = f"{rid}\t{rn}\t{ss}\t{rsa_val:.4f}"
            output_lines.append(line)

        output_text = "\n".join(output_lines)
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, 'w') as f:
            f.write(output_text)

        print(f"  Done! Total residues: {len(res_ids)}")

    except FileNotFoundError:
        print(f"  Error: File {pdb_file} not found")
        sys.exit(1)
    except ValueError as e:
        print(f"  Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"  Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

def fetch_alphafold_structures(uniprot_ids, output_dir="pdb_files_poi"):
    """
    Download AlphaFold predicted structures from AlphaFold DB using REST API.
    Only downloads if file doesn't already exist.
    
    Parameters:
    -----------
    uniprot_ids : str or list
        Single UniProt ID or list of UniProt IDs (e.g., "Q04917" or ["Q04917", "P27348"])
    output_dir : str
        Directory to save PDB files
    """
    
    # Ensure it's a list
    if isinstance(uniprot_ids, str):
        uniprot_ids = [uniprot_ids]
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    print(f"Checking {len(uniprot_ids)} AlphaFold structure(s) in {output_dir}/...")
    
    for uniprot_id in uniprot_ids:
        output_file = Path(output_dir) / f"{uniprot_id}.pdb"
        
        # Skip if file already exists
        if output_file.exists():
            print(f"  {uniprot_id}: Already exists, skipping download")
            continue
        
        try:
            print(f"  Querying AlphaFold DB for {uniprot_id}...", end=" ")
            
            # Query AlphaFold API to get prediction metadata
            api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
            with urllib.request.urlopen(api_url) as response:
                data = json.loads(response.read().decode())
            
            if not data:
                print(f"No prediction found!")
                continue
            
            # Get PDB download URL from API response
            pdb_url = data[0]["pdbUrl"]
            af_id = data[0]["entryId"]
            
            print(f"Downloading {af_id}...", end=" ")
            
            # Download PDB file
            urllib.request.urlretrieve(pdb_url, output_file)
            print("Done!")
            
        except Exception as e:
            print(f"Failed! Error: {e}")
    
    print("Download check complete!")

def main():
    # Download AlphaFold structures using UniProt IDs (only if they don't exist)
    fetch_alphafold_structures(["Q04917", "P27348", "P61981"], "pdb_files_poi")
    
    # Create output directory
    output_dir = Path("dssp_outputs")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Analyze all PDB files in the directory
    pdb_dir = Path("pdb_files_poi")
    pdb_files = sorted(pdb_dir.glob("*.pdb"))
    
    print(f"\n{'='*60}")
    print(f"Analyzing {len(pdb_files)} PDB file(s)...")
    print(f"{'='*60}\n")
    
    for pdb_file in pdb_files:
        uniprot_id = pdb_file.stem
        output_file = output_dir / f"{uniprot_id}_ss_rsa.tsv"
        print(f"Processing: {pdb_file.name}")
        calculate_ss_rsa(str(pdb_file), str(output_file))
    
    print(f"\n{'='*60}")
    print(f"All files processed! Results saved to {output_dir}/")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
