# CAPTAF Complete Component Checklist
# Context-Aware Snakemake Pipeline for PTM-Driven Binding Optimization with AlphaFold3
# Use this to verify you have all required files

## ═══════════════════════════════════════════════════════════
## CORE PIPELINE FILES
## ═══════════════════════════════════════════════════════════

[ ] Snakefile                           # Main pipeline definition
[ ] config.yaml.template                # Configuration template
[ ] PTM_to_CCD_mapping.csv             # PTM to CCD code mapping

## ═══════════════════════════════════════════════════════════
## WRAPPER & AUTOMATION SCRIPTS
## ═══════════════════════════════════════════════════════════

[ ] captaf.sh                           # Main wrapper script (executable)
[ ] setup_captaf.sh                     # Setup script (executable)
[ ] verify_captaf.sh                    # Verification script (executable)

## ═══════════════════════════════════════════════════════════
## PYTHON SCRIPTS (scripts/ directory)
## ═══════════════════════════════════════════════════════════

### PTM Prediction & Processing
[ ] scripts/run_musitedeep_all.py      # Run MusiteDeep for PTM prediction
[ ] scripts/convert_musitedeep_to_tsv.py  # Convert MusiteDeep output to TSV

### Structural Features
[ ] scripts/predict_dssp.py            # Predict DSSP features (RSA, SS)
[ ] scripts/merge_dssp_ptm.py          # Merge PTM with DSSP data

### Conservation & Disorder
[ ] scripts/predict_conservation.py    # Predict conservation with MMseqs2
[ ] scripts/predict_disorder.py        # Predict disorder with IUPred2A

### Master Table & Variants
[ ] scripts/create_master_table.py     # Merge all data sources
[ ] scripts/generate_variants.py       # Generate 8 PTM variants

### Target PTM Data
[ ] scripts/fetch_dbptm_ptms.py        # Fetch PTMs from dbPTM
[ ] scripts/process_dbptm_ptms.py      # Process dbPTM data

### AlphaFold3 JSON Generation
[ ] scripts/generate_af3_json.py       # Generate AF3 JSON for variants
[ ] scripts/generate_baseline_control.py    # Generate baseline JSON
[ ] scripts/generate_positive_control.py    # Generate positive control JSON

### Results Analysis
[ ] scripts/analyze_af3_results.py     # Analyze & rank AF3 outputs

## ═══════════════════════════════════════════════════════════
## CONDA ENVIRONMENTS
## ═══════════════════════════════════════════════════════════

### Environment Definition Files
[ ] envs/new_thesis_env.yml            # Main pipeline environment
[ ] envs/musite_env.yml                # MusiteDeep environment
[ ] envs/af3.yml                       # AlphaFold3 environment

### Environment Installation
[ ] envs/install_envs.sh               # Install all environments (executable)

## ═══════════════════════════════════════════════════════════
## EXTERNAL TOOLS (must be present)
## ═══════════════════════════════════════════════════════════

[ ] iupred2a/                          # IUPred2A directory
    [ ] iupred2a/iupred2a.py          # IUPred2A main script
    [ ] iupred2a/iupred2a_lib/        # IUPred2A libraries

[ ] MusiteDeep_web/                    # MusiteDeep directory
    [ ] MusiteDeep_web/MusiteDeep/models/  # PTM models

## ═══════════════════════════════════════════════════════════
## DATABASES (must be configured)
## ═══════════════════════════════════════════════════════════

[ ] alphafold3/database/uniref50DB     # MMseqs2 database
[ ] alphafold3/database/uniref50DB.index  # MMseqs2 index
[ ] af3_models/                        # AlphaFold3 model weights
[ ] alphafold3/database/               # AF3 genetic databases

## ═══════════════════════════════════════════════════════════
## DOCUMENTATION (for distribution)
## ═══════════════════════════════════════════════════════════

[ ] README.md                          # Main documentation
[ ] LICENSE                            # Software license
[ ] CITATION.cff                       # Citation information
[ ] INSTALL.md                         # Installation guide
[ ] USAGE.md                           # Usage guide

## ═══════════════════════════════════════════════════════════
## VERIFICATION COMMANDS
## ═══════════════════════════════════════════════════════════

Run these commands to verify your installation:

# Check all required Python scripts exist
ls scripts/*.py | wc -l
# Should show: 14 files

# Check environment files
ls envs/*.yml | wc -l
# Should show: 3 files

# Check executables are executable
ls -l captaf.sh setup_captaf.sh verify_captaf.sh envs/install_envs.sh
# Should show: -rwxr-xr-x (executable)

# Check Snakefile exists
ls -l Snakefile
# Should exist

# Check IUPred2A exists
ls iupred2a/iupred2a.py
# Should exist

# Check MusiteDeep exists
ls MusiteDeep_web/MusiteDeep/models/ | wc -l
# Should show multiple PTM models

# Check conda environments
conda env list | grep -E "new_thesis_env|musite_env|af3"
# Should show all 3 environments

## ═══════════════════════════════════════════════════════════
## SUMMARY COUNTS
## ═══════════════════════════════════════════════════════════

Total Python scripts: 14
Total shell scripts: 4 (captaf.sh, setup_captaf.sh, verify_captaf.sh, install_envs.sh)
Total environment files: 3
Total conda environments: 3
External tools: 2 (IUPred2A, MusiteDeep)
Databases: 3 locations (MMseqs2, AF3 models, AF3 database)

## ═══════════════════════════════════════════════════════════
## NOTES
## ═══════════════════════════════════════════════════════════

1. All paths should be RELATIVE (no absolute paths to user directories)
2. All shell scripts should be executable (chmod +x)
3. Environment YAML files should be exported with --no-builds
4. Run ./verify_captaf.sh to automatically check all components
