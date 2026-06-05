# CAPTAF - Context-Aware Pipeline for PTM-Driven Binding Optimization with AlphaFold3

**CAPTAF** is an automated pipeline that systematically explores how post-translational modifications (PTMs) affect protein-protein binding by generating context-aware PTM variants and predicting their structures using AlphaFold3.

---

## Overview

Post-translational modifications (PTMs) play critical roles in regulating protein-protein interactions, but predicting their structural impact remains challenging. CAPTAF addresses this by:

1. **Fetching validated PTM sites** from dbPTM for the Protein of Interest [CONFIRM: also MusiteDeep?]
2. **Analyzing structural context** (buried/interface/exposed, ordered/disordered, conserved/variable)
3. **Generating 8 strategic variants** that sample different PTM contexts
4. **Predicting structures** with AlphaFold3 for each variant plus baseline and positive controls
5. **Ranking variants** by predicted binding quality (ipTM scores)

This reduces the combinatorial explosion from 2^N possible PTM combinations to just 8 strategically selected variants, making PTM exploration computationally tractable while maintaining biological relevance.

---

## Key Features

✅ **Automated PTM retrieval** - Fetches validated PTMs from dbPTM  
✅ **Context-aware variant generation** - Selects PTMs based on structural location, conservation, and disorder  
✅ **AlphaFold3 integration** - Predicts structures for all variants automatically  
✅ **Intelligent ranking** - Ranks variants by ipTM, interface quality, and contact metrics  
✅ **Validation controls** - Baseline (no PTMs) and positive control (all PTMs) comparisons  
✅ **Rich outputs** - TSV rankings, HTML reports, and text recommendations  
✅ **One-command execution** - Simple wrapper script handles entire pipeline  

> ⚠️ **AlphaFold3 PTM support limitation:** Not all PTM types are supported. Supported types include Phosphorylation (pSer, pThr, pTyr), Acetylation, and Methylation. **Ubiquitination and SUMOylation are NOT supported** and cannot be modelled. Check `PTM_to_CCD_mapping.csv` for the full list of supported modifications before running.

---

## System Requirements

**Important:** CAPTAF requires significant computational resources. **Use a server with both CPU and GPU**, not a laptop.

### Minimum Requirements
- **CPU**: 8+ cores
- **RAM**: 32 GB minimum, 64 GB recommended
- **GPU**: NVIDIA GPU with 16+ GB VRAM (for AlphaFold3)
- **Storage**: 2+ TB for databases and outputs
- **OS**: Linux (Ubuntu 20.04+ or CentOS 7+)
- **CUDA**: 11.8+

### Software Prerequisites
- Conda (Miniconda or Anaconda)
- Git
- Wget

---

## Installation

### Step 1: Clone CAPTAF Repository

```bash
git clone https://github.com/yourusername/captaf.git
cd captaf
chmod +x captaf.sh setup_captaf.sh verify_captaf.sh envs/install_envs.sh
```

---

### Step 2: Install Conda Environments

CAPTAF requires three conda environments:

```bash
# Install all three environments (takes 10-20 minutes)
./envs/install_envs.sh
```

This creates:
- **new_thesis_env**: Main pipeline environment (Snakemake, BioPython, Pandas, etc.)
- **musite_env**: MusiteDeep PTM prediction (TensorFlow 1.x, Keras)
- **af3**: AlphaFold3 environment (JAX, Haiku, etc.)

---

### Step 3: Install AlphaFold3

```bash
# Clone AlphaFold3 repository
git clone https://github.com/google-deepmind/alphafold3.git
cd alphafold3/

# Create and activate conda environment
conda create -n af3 python=3.11.0
conda activate af3

# Install dependencies
conda install hmmer=3.4
conda install cmake=3.30.2
conda install -c conda-forge gcc=12.4.0 libgcc=12.4.0
conda install -c conda-forge gxx_linux-64=12.4.0
conda install -c conda-forge zlib

# Install Python packages
pip install jmp==0.0.4 ml-dtypes==0.5.0 opt-einsum==3.4.0
pip3 install -r dev-requirements.txt
pip3 install --no-deps .

# Build data processing tools
build_data

# Download AlphaFold3 databases (this takes several hours and ~2TB)
mkdir database
cd database
nohup ../fetch_databases.sh . &

# Monitor download progress
tail -f nohup.out

cd ../..
```

**Verify AlphaFold3 installation:**

```bash
conda activate af3
python alphafold3/run_alphafold.py --help
conda deactivate
```

If the help message displays, AlphaFold3 is correctly installed.

---

### Step 4: Download MMseqs2 Database

```bash
# Create database directory
mkdir -p alphafold3/database
cd alphafold3/database

# Download UniRef50 (~30 GB)
wget https://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/UniRef50_2020_06_hhsuite.tar.gz
tar -xzf UniRef50_2020_06_hhsuite.tar.gz

# Create MMseqs2 database
mmseqs createdb UniRef50_2020_06_consensus.fasta uniref50DB

# Create index (required for searches)
mmseqs createindex uniref50DB tmp --threads 8

cd ../..
```

---

### Step 5: Verify Installation

```bash
# Verify all components
./verify_captaf.sh
# Expected: ✓ All required components are present!

# Check setup
./setup_captaf.sh
# Expected: ✓ CAPTAF is ready to use!
```

---

## Before Running

### Download PDB Structure

CAPTAF requires a PDB structure file for the Protein of Interest:

```bash
# Download from AlphaFold Database
wget https://alphafold.ebi.ac.uk/files/AF-P27348-F1-model_v4.pdb \
     -O pdb_files_poi/P27348.pdb

# Or from PDB (for experimental structures)
wget https://files.rcsb.org/download/1YCR.pdb \
     -O pdb_files_poi/P04637.pdb
```

**Important:** Place PDB files in `pdb_files_poi/` directory within CAPTAF installation, named as `PROTEIN_ID.pdb`

---

## Working with Partial Sequences and Peptide Targets

> ⚠️ **This section is critical if your target is a peptide fragment or partial sequence.**

CAPTAF handles two types of targets differently:

| Target Type | PTM Handling | Action Required |
|-------------|-------------|-----------------|
| Full protein in dbPTM | Automatic | None |
| Peptide fragment / partial sequence | Manual | Create PTM files (see below) |
| Protein not in dbPTM | Manual | Create PTM files (see below) |

When a peptide or partial sequence is used as the target, dbPTM cannot find a match (it stores only full protein sequences). The pipeline will **stop with an error** at the PTM fetch step. This is by design — the pipeline exits cleanly rather than continuing with empty PTM data.

**Key rule for partial sequences:** PTM positions in the manual files must be **relative to your peptide/fragment** (1-indexed from position 1 of your FASTA sequence), NOT the full protein positions from UniProt.

Example: If your peptide is residues 255–264 of a full protein, and the relevant PTM is at full-protein position 259, then the peptide-relative position is 259 − 255 + 1 = **5**.

---

## Usage

### Basic Command

```bash
./captaf.sh --poi POI.fasta --target TARGET.fasta [OPTIONS]
```

### Required Arguments

| Argument | Description |
|----------|-------------|
| `--poi FILE` | Protein of Interest FASTA file (can be anywhere) |
| `--target FILE` | Target/Partner protein FASTA file (can be anywhere) |

### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--output DIR` | `captaf_output` | Output directory |
| `--threads N` | `4` | Number of CPU threads |
| `--use-controls` | `false` | Include baseline and positive controls |
| `--keep-intermediates` | `true` | Keep all intermediate files |
| `--dry-run` | `false` | Show pipeline steps without execution |
| `--verbose` | `false` | Show detailed Snakemake output |
| `-h, --help` | - | Show help message |
| `-v, --version` | - | Show version |

> **Recommendation:** Always use `--use-controls` when testing known interactions. Controls (baseline = no PTMs, positive = all PTMs) are essential for interpreting whether variant scores represent genuine PTM effects.

### Example Usage

**Basic analysis (8 variants only):**
```bash
./captaf.sh --poi /data/protein.fasta --target /data/partner.fasta
```

**With validation controls (recommended):**
```bash
./captaf.sh --poi 14-3-3z.fasta --target RAF1_peptide.fasta \
    --output 14-3-3-RAF1-validation \
    --use-controls
```

**Custom output directory with more threads:**
```bash
./captaf.sh --poi /path/to/POI.fasta --target /path/to/TARGET.fasta \
    --output /results/my_analysis \
    --threads 16
```

**Test pipeline without execution:**
```bash
./captaf.sh --poi protein.fasta --target partner.fasta --dry-run
```

### Resuming a Stopped Run

If the pipeline stops (e.g., for manual PTM file creation), always resume with `--rerun-triggers mtime` to avoid Snakemake re-running steps due to code change detection:

```bash
snakemake --configfile {output_dir}/config.yaml \
          --rerun-triggers mtime \
          --cores 8 \
          --use-conda
```

If Snakemake still complains about changed code, clear the metadata for the affected file:
```bash
snakemake --configfile {output_dir}/config.yaml \
          --cleanup-metadata {output_dir}/target_ptm_dir/PROTEIN_Target_dbptm.csv
```

---

# Manual PTM File Creation for Target Proteins

> ⚠️ **IMPORTANT:** This process is **ONLY for TARGET proteins**, not POI. The pipeline automatically handles POI PTMs from dbPTM.

---

## When Manual Target PTM Files Are Needed

The pipeline automatically fetches POI PTMs from dbPTM and generates variants. However, you must create PTM files manually for the **TARGET** in these situations:

1. **Peptide fragments as target** — dbPTM only contains full protein sequences
2. **Partial sequences as target** — only a domain or region of the full protein used
3. **Proteins not in dbPTM** — novel or poorly characterized targets
4. **Custom PTM testing** — testing specific PTM sites from literature

---

## File Format Requirements

### Required Files:

```
{output_dir}/target_ptm_dir/{target_uniprot_id}_Target_dbptm.csv
{output_dir}/target_ptm_dir/{target_uniprot_id}_Target_processed.csv
```

### File Formats:

**1. Raw dbPTM format** (`{target_id}_Target_dbptm.csv`):
```csv
UniProt_ID,Position,PTM_Type,Source
P04049,3,Phosphorylation,dbPTM
P04049,5,Phosphorylation,dbPTM
```

**2. Processed format** (`{target_id}_Target_processed.csv`):
```csv
UniProt_ID,Position,AA,PTM_Type
P04049,3,S,Phosphorylation
P04049,5,S,Phosphorylation
```

---

## Complete Example: RAF1 Peptide Validation

**Scenario:** Testing 14-3-3ζ binding to RAF1 phosphopeptide (residues 255-264: QRSTSTPNVH)

### Step 1: Create Target FASTA
```bash
cat > P04049_Target.fasta << 'EOF'
>sp|P04049|RAF1_HUMAN_255-264
QRSTSTPNVH
EOF
```

### Step 2: Run Pipeline (Will Stop at PTM Fetch)
```bash
./captaf.sh --poi P63104.fasta \
            --target P04049_Target.fasta \
            --output 14-3-3-RAF1-validation \
            --use-controls
```

The pipeline will stop with an error message indicating no PTMs found in dbPTM. This is expected behaviour for peptide targets.

### Step 3: Create PTM Directory
```bash
mkdir -p 14-3-3-RAF1-validation/target_ptm_dir
```

### Step 4: Create Raw dbPTM File

**Position mapping:** Map full protein positions to peptide-relative positions.
The peptide QRSTSTPNVH corresponds to full protein residues 255–264.

- Full protein Ser257 → Peptide position 3
- Full protein Thr258 → Peptide position 4
- Full protein Ser259 → Peptide position 5 ← key 14-3-3 binding site
- Full protein Thr260 → Peptide position 6

```bash
cat > 14-3-3-RAF1-validation/target_ptm_dir/P04049_Target_dbptm.csv << 'EOF'
UniProt_ID,Position,PTM_Type,Source
P04049,3,Phosphorylation,dbPTM
P04049,4,Phosphorylation,dbPTM
P04049,5,Phosphorylation,dbPTM
P04049,6,Phosphorylation,dbPTM
EOF
```

### Step 5: Create Processed File

```bash
cat > 14-3-3-RAF1-validation/target_ptm_dir/P04049_Target_processed.csv << 'EOF'
UniProt_ID,Position,AA,PTM_Type
P04049,3,S,Phosphorylation
P04049,4,T,Phosphorylation
P04049,5,S,Phosphorylation
P04049,6,T,Phosphorylation
EOF
```

### Step 6: Resume Pipeline
```bash
snakemake --configfile 14-3-3-RAF1-validation/config.yaml \
          --rerun-triggers mtime \
          --cores 8 \
          --use-conda
```

---

## Important Guidelines

### Position Numbering:
- Use positions **relative to your peptide/fragment** (1-indexed from position 1 of your FASTA)
- **NOT** the full protein positions from UniProt
- Formula: peptide_position = full_protein_position − fragment_start + 1

### Amino Acid Matching:
- Must exactly match residues in your FASTA sequence
- Pipeline validates AA matches sequence

### PTM Types:
- Must match supported types in `PTM_to_CCD_mapping.csv`
- Common supported types: Phosphorylation, Acetylation, Methylation
- **Not supported by AlphaFold3:** Ubiquitination, SUMOylation — do not include these

### File Naming:
- Must use exact UniProt ID from FASTA header
- Format: `{ID}_Target_dbptm.csv` and `{ID}_Target_processed.csv`
- Case-sensitive

---

## POI vs Target PTM Handling

| Aspect | POI (Protein of Interest) | Target Protein |
|--------|--------------------------|----------------|
| **PTM Source** | Automatic from dbPTM | Manual (if peptide/partial) or automatic |
| **Filtering** | Yes (8 variant strategies) | No filtering |
| **Usage** | Generate 8 test variants | Used in positive control |
| **Manual Creation** | ❌ Never needed | ✅ Required for peptides and partial sequences |

---

## Output Files

### Directory Structure

```
my_analysis/                         # Your output directory
├── config.yaml                      # Generated configuration
├── ptm_predictions/                 # PTM predictions
├── conservation_outputs/            # Conservation scores
├── disordered_regions_outputs/      # Disorder predictions
├── dssp_outputs/                    # Secondary structure predictions
├── master_tables_poi/               # Merged feature tables
├── target_ptm_dir/                  # Target PTM files (manual or auto)
├── poi_variants/                    # Generated variants (v1-v8)
├── alphafold_inputs/                # JSON files for AlphaFold3
├── alphafold_outputs/               # AlphaFold3 structure predictions
│   ├── POI_TARGET_v1/
│   ├── POI_TARGET_v2/
│   ├── ...
│   ├── POI_TARGET_v8/
│   ├── POI_TARGET_baseline/         # Only if --use-controls
│   └── POI_TARGET_positive_control/ # Only if --use-controls
├── POI_TARGET_ranking.tsv           # ← Final ranking table
├── POI_TARGET_summary.html          # ← Interactive HTML report
└── POI_TARGET_recommendations.txt   # ← Top recommendations
```

### Main Results (Top Level)

The three main result files are placed directly in the output directory for easy access:

**1. `POI_TARGET_ranking.tsv`** - Tab-separated ranking table

| Column | Good Value | Description |
|--------|------------|-------------|
| `ipTM` | > 0.6 | Interface predicted TM-score (binding confidence) |
| `pTM` | > 0.5 | Overall predicted TM-score (structure confidence) |
| `pLDDT_mean` | > 70 | Mean predicted LDDT score (structure quality) |
| `interface_PAE` | < 8 Å | Interface position error |
| `contacts` | > 20 | Number of interface residue contacts |
| `variant` | - | Variant name (v1-v8, baseline, positive) |
| `type` | - | variant or control |
| `composite_score` | > 0.7 | Combined ranking score |

**2. `POI_TARGET_summary.html`** - Interactive HTML report with sortable tables and visualizations

**3. `POI_TARGET_recommendations.txt`** - Human-readable recommendations

---

## Understanding Results

### The 8 Variant Strategies

CAPTAF generates 8 variants based on PTM structural context:

| Variant | RSA | Conservation | Secondary Structure | Disorder | Context | Biological Rationale |
|---------|-----|-------------|---------------------|----------|---------|---------------------|
| **v1** | < 0.20 | > 0.7 | H or E | < 0.5 | Buried Conserved | Allosteric regulation in protein core |
| **v2** | < 0.20 | < 0.4 | H or E | < 0.5 | Buried Variable | Species-specific core modifications |
| **v3** | 0.20–0.50 | > 0.7 | H or E | < 0.5 | Interface Conserved | Universal interface binding modulators |
| **v4** | 0.20–0.50 | < 0.5 | H or E | < 0.5 | Interface Variable | Adaptive interface plasticity |
| **v5** | > 0.50 | — | — | < 0.5 | Exposed Ordered | Surface recognition signals |
| **v6** | — | > 0.7 | — | > 0.5 | Disordered Conserved | Conserved disorder-mediated regulation |
| **v7** | — | < 0.4 | — | > 0.5 | Disordered Variable | Lineage-specific flexible regulation |
| **v8** | — | — | — | — | High Confidence | All high-scoring PTMs (score ≥ 0.1) + all phosphorylation |

> **Key insight from validation:** No single strategy is universally optimal. Disordered/conserved-site selection (v6) performs best for compact reader-peptide interactions (e.g. 14-3-3 complexes). Exposed/high-confidence strategies (v5, v8) may perform better for larger folded-domain complexes. Always compare across variants rather than relying on one.

### Interpreting ipTM Scores

The interface predicted TM-score (ipTM) is the primary metric for ranking variants:

| ipTM Range | Interpretation | Recommendation |
|------------|----------------|----------------|
| **> 0.8** | High confidence binding | Strong candidate for experimental validation |
| **0.6 – 0.8** | Moderate confidence | Consider for further computational analysis |
| **0.4 – 0.6** | Low confidence | May indicate weak or transient binding |
| **< 0.4** | Very low confidence | Likely not a strong direct interaction |

### Interpreting Δ ipTM (with controls enabled)

When controls are included, compare the best variant score against the baseline:

- **Δ ipTM > 0** — PTMs enhance predicted binding (positive regulation)
- **Δ ipTM ≈ 0** — PTMs have no effect on binding (constitutive interaction)
- **Δ ipTM < 0** — PTMs disrupt predicted binding (negative regulation)

### Quality Thresholds

A good prediction typically meets these criteria:
- ipTM > 0.6
- pTM > 0.5
- pLDDT > 70
- interface_PAE < 8 Å
- contacts > 20

---

## Validation Controls

For well-studied proteins with known PTM effects, use validation controls:

```bash
./captaf.sh --poi 14-3-3z.fasta --target RAF1_peptide.fasta \
    --output my_analysis --use-controls
```

This adds two control conditions:
- **Baseline control**: Unmodified proteins (no PTMs)
- **Positive control**: All validated PTMs from dbPTM database

Compare variant performance against these controls to understand PTM-dependent binding changes.

### Recommended Validation Cases

Test CAPTAF with these well-characterised PTM-regulated interactions:

1. **14-3-3ζ + RAF1 peptide (pS259)** — PTMs required for binding
   - POI: P63104 (14-3-3ζ, 245 aa, full sequence)
   - Target: P04049 (RAF1, peptide residues 255–264: QRSTSTPNVH)
   - Expected: Baseline lower ipTM, phosphorylated variant (v6) high ipTM
   - Key PTM: S259 phosphorylation (peptide position 5)

2. **INSR + IRS1** — PTMs enhance binding
   - POI: P06213 (INSR kinase domain, residues 1023–1298)
   - Target: P35568 (IRS1 PTB + pTyr region, residues 160–500)
   - Expected: Baseline lower ipTM, phosphorylated variants higher ipTM

3. **CDK2 + p27Kip1** — PTMs disrupt binding
   - POI: P24941 (CDK2, 298 aa, full sequence)
   - Target: P46527 (p27Kip1, 198 aa, full sequence)
   - Expected: Baseline high ipTM (strong inhibitory complex), PTM variants lower ipTM
   - Key PTM: Thr187 phosphorylation on p27

---

## Troubleshooting

### Common Issues

**1. "PDB file not found: pdb_files_poi/PROTEIN.pdb"**

Download PDB structure for your protein:
```bash
# From AlphaFold Database
wget https://alphafold.ebi.ac.uk/files/AF-PROTEIN_ID-F1-model_v4.pdb \
     -O pdb_files_poi/PROTEIN_ID.pdb
```

**2. "No PTMs found in dbPTM" — pipeline stops**

Expected behaviour for peptide/partial-sequence targets. Create manual PTM files and resume:
```bash
mkdir -p {output_dir}/target_ptm_dir
# Create files as described in the Manual PTM section above
snakemake --configfile {output_dir}/config.yaml \
          --rerun-triggers mtime --cores 8 --use-conda
```

**3. "Snakemake not found in new_thesis_env"**

```bash
conda activate new_thesis_env
conda install -c conda-forge -c bioconda snakemake
conda deactivate
```

**4. "MMseqs2 database index not found"**

```bash
cd alphafold3/database
mmseqs createindex uniref50DB tmp --threads 8
cd ../..
```

**5. "Could not extract protein ID from FASTA header"**

Ensure FASTA header follows UniProt format:
```
>sp|P27348|1433T_HUMAN 14-3-3 protein theta OS=Homo sapiens
```

**6. AlphaFold3 out of GPU memory**

Reduce the number of parallel jobs or use a GPU with more memory.

**7. "No space left on device"**

```bash
df -h
sudo find /tmp -type f -atime +7 -delete
```

**8. Low ipTM scores for all variants**

This can occur when:
- The proteins do not directly interact
- PTMs do not significantly affect this interaction
- AlphaFold3 struggles with this particular complex size or type

Validate pipeline functionality using the recommended test cases above.

**9. Snakemake re-runs completed steps after code change**

Use `--rerun-triggers mtime` to prevent this, or clear metadata for specific files:
```bash
snakemake --configfile {output_dir}/config.yaml \
          --cleanup-metadata {output_dir}/target_ptm_dir/FILE.csv
snakemake --configfile {output_dir}/config.yaml \
          --rerun-triggers mtime --cores 8 --use-conda
```

---

## Pipeline Workflow

```
Input FASTA files (POI + Target, from any location)
    ↓
PTM Retrieval (dbPTM) — validated PTMs for POI
    ↓
Structure Analysis (DSSP) — RSA, secondary structure
    ↓
Conservation (MMseqs2) + Disorder (IUPred2A/3)
    ↓
Master Table — merge all features per residue
    ↓
Variant Generation — 8 strategic PTM combinations (v1–v8)
    ↓
AlphaFold3 JSON Generation
    ↓
Structure Prediction (AlphaFold3) — variants + baseline + positive control
    ↓
Analysis and Ranking — ipTM, pTM, pLDDT, interface PAE, contacts
    ↓
Results — TSV ranking table, HTML report, recommendations
```

---

## Getting Help

```bash
# Verify installation
./verify_captaf.sh

# Check setup
./setup_captaf.sh

# Test pipeline (no execution)
./captaf.sh --poi POI.fasta --target TARGET.fasta --dry-run

# View help
./captaf.sh --help
```

For bug reports and questions, open an issue on GitHub:
https://github.com/yourusername/captaf/issues

---

**CAPTAF v1.0.0** — Universität des Saarlandes
