#!/usr/bin/env python3
"""
Add CCD-code column to every filtered_ptm_*.tsv inside
~/Thesis/thesis_files/dssp_outputs/
PTM_Type is mapped residue-by-residue using the official
PDB Chemical Component Dictionary (components.cif).
"""

import os, sys, csv, glob, gemmi

# ---------- paths ----------
ROOT      = os.path.expanduser("~/Thesis/thesis_files")
CCD_FILE  = os.path.join(ROOT, "components.cif")
IN_DIR    = os.path.join(ROOT, "dssp_outputs")
SCRIPT_DIR= os.path.dirname(os.path.abspath(__file__))

# ---------- build lookup table ----------
def build_exact_map(cif_path):
    """Return dict: (residue, normalised_ptm_type) -> ccd_code."""
    exact = {}
    doc = gemmi.cif.read(cif_path)
    for block in doc:
        if block.find_value("_chem_comp.pdbx_pcm") != "Y":
            continue
        ccd  = block.name
        res  = (block.find_value("_pdbx_chem_comp_pcm.target_residue") or "").strip().upper()
        ptm  = (block.find_value("_pdbx_chem_comp_pcm.type") or "").strip().lower()
        if res and ptm:
            exact[(res, ptm)] = ccd
    return exact

EXACT = build_exact_map(CCD_FILE)

# ---------- normaliser ----------
def normalise(name: str) -> str:
    """Canonical key for fuzzy matching."""
    return (name.lower()
            .replace("n6-", "n-")
            .replace("sumoylation", "ubiquitination")
            .replace(" ", "_")
            .replace("-", "_"))

# ---------- mapper ----------
def map_ptm(residue: str, musite_name: str) -> str:
    """Return CCD code or empty string."""
    key = (residue.upper(), normalise(musite_name))
    return EXACT.get(key, "")

# ---------- single-file processor ----------
def add_ccd_to_tsv(tsv_path):
    """Read TSV, add CCD_code column, overwrite file."""
    tmp_path = tsv_path + ".tmp"
    with open(tsv_path, newline='') as fin, open(tmp_path, 'w', newline='') as fout:
        rdr = csv.DictReader(fin, delimiter='\t')
        flds = rdr.fieldnames + ["CCD_code"]
        wrt  = csv.DictWriter(fout, fieldnames=flds, delimiter='\t')
        wrt.writeheader()

        for row in rdr:
            residue = row.get("AA", "")          # one-letter code
            ptm     = row.get("PTM_Type", "")
            row["CCD_code"] = map_ptm(residue, ptm)
            wrt.writerow(row)

    os.replace(tmp_path, tsv_path)
    print(f"Updated  {os.path.basename(tsv_path)}")

# ---------- main loop ----------
def main():
    pattern = os.path.join(IN_DIR, "filtered_ptm_*.tsv")
    files = glob.glob(pattern)
    if not files:
        print(f"No files match {pattern}")
        sys.exit(1)
    for f in files:
        add_ccd_to_tsv(f)
    print("All done.")

if __name__ == "__main__":
    main()
