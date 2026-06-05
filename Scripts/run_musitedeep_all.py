import os
import subprocess
import sys

def run_musitedeep_all(input_fasta, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Path to MusiteDeep repo (relative to thesis_files/)
    repo_dir = os.path.join(os.path.dirname(__file__), "..", "MusiteDeep_web", "MusiteDeep")

    # Run all 13 PTM models in one go using the exact model-prefix string
    model_prefix = (
        "models/Phosphoserine_Phosphothreonine;"
        "Phosphotyrosine;"
        "N6-acetyllysine;"
        "Methyllysine;"
        "Methylarginine;"
        "Ubiquitination;"
        "SUMOylation;"
        "N-linked_glycosylation;"
        "O-linked_glycosylation;"
        "Hydroxyproline;"
        "Hydroxylysine;"
        "S-palmitoyl_cysteine;"
        "Pyrrolidone_carboxylic_acid"
    )

    output_path = os.path.join(output_dir, os.path.basename(input_fasta).replace(".fasta", "_all"))

    cmd = [
        sys.executable,
        "predict_multi_batch.py",
        "-input", os.path.abspath(input_fasta),
        "-output", os.path.abspath(output_path),
        "-model-prefix", model_prefix
    ]

    print("Running in:", repo_dir)
    print("Command:", " ".join(cmd))

    subprocess.run(cmd, check=True, cwd=repo_dir)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python run_musitedeep_all.py <input_fasta> <output_dir>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_dir = sys.argv[2]

    run_musitedeep_all(input_fasta, output_dir)

