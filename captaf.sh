#!/bin/bash
###############################################################################
# CAPTAF: Context-Aware Pipeline for PTM-Driven Binding Optimization
#         with AlphaFold3
#
# Main wrapper script for automated analysis
#
# Author: Vishva Darji
# Institution: Universität des Saarlandes
###############################################################################

set -e  # Exit on error
set -o pipefail

# Version
VERSION="1.0.0"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Default values
OUTPUT_DIR="captaf_output"
THREADS=4
USE_CONTROLS=false
DRY_RUN=false
VERBOSE=false
KEEP_INTERMEDIATES=true

###############################################################################
# Functions
###############################################################################

print_banner() {
    echo -e "${CYAN}"
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║                       CAPTAF v${VERSION}                        ║"
    echo "║  Context-Aware PTM Analysis Pipeline with AlphaFold3      ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

print_usage() {
    cat << EOF
${CYAN}CAPTAF${NC} - Context-Aware Pipeline for PTM-Driven Binding Optimization

${YELLOW}Usage:${NC}
    ./captaf.sh --poi POI.fasta --target TARGET.fasta [OPTIONS]

${YELLOW}Required Arguments:${NC}
    --poi FILE          Protein of Interest FASTA file (any path)
    --target FILE       Target/Partner protein FASTA file (any path)

${YELLOW}Optional Arguments:${NC}
    --output DIR        Output directory (default: captaf_output)
    --threads N         Number of CPU threads (default: 4)
    --use-controls      Include baseline and positive controls
    --keep-intermediates Keep all intermediate files (default: true)
    --dry-run           Show what would be executed without running
    --verbose           Show detailed Snakemake output
    -h, --help          Show this help message
    -v, --version       Show version

${YELLOW}Examples:${NC}
    ${GREEN}# Basic run${NC}
    ./captaf.sh --poi /data/protein.fasta --target /data/partner.fasta

    ${GREEN}# With controls${NC}
    ./captaf.sh --poi p53.fasta --target MDM2.fasta --use-controls

${YELLOW}Requirements:${NC}
    - PDB structure for POI in: pdb_files_poi/PROTEIN_ID.pdb
    - Run from CAPTAF installation directory

EOF
}

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

check_dependencies() {
    log_info "Checking dependencies..."
    
    if ! command -v conda &> /dev/null; then
        log_error "Conda not found"
        exit 1
    fi
    
    if ! conda run -n new_thesis_env snakemake --version &> /dev/null; then
        log_error "Snakemake not found in new_thesis_env"
        exit 1
    fi
    
    local missing_envs=false
    for env in new_thesis_env musite_env af3; do
        if conda env list | grep -q "^${env} "; then
            log_success "Found environment: ${env}"
        else
            log_error "Missing environment: ${env}"
            missing_envs=true
        fi
    done
    
    if [ "$missing_envs" = true ]; then
        log_error "Run: ./envs/install_envs.sh"
        exit 1
    fi
    
    # Check required files
    for item in Snakefile scripts MusiteDeep_web alphafold3 iupred2a PTM_to_CCD_mapping.csv; do
        if [ ! -e "$item" ]; then
            log_error "$item not found - run from CAPTAF directory"
            exit 1
        fi
    done
    
    log_success "All dependencies found"
}

validate_fasta() {
    local fasta_file=$1
    local label=$2
    
    if [ ! -f "$fasta_file" ]; then
        log_error "$label FASTA not found: $fasta_file"
        exit 1
    fi
    
    if ! grep -q "^>" "$fasta_file"; then
        log_error "$label is not valid FASTA format"
        exit 1
    fi
}

extract_protein_id() {
    local fasta_file=$1
    local header=$(grep "^>" "$fasta_file" | head -1)
    local protein_id=""
    
    if echo "$header" | grep -q '^>sp|'; then
        protein_id=$(echo "$header" | sed 's/^>sp|\([^|]*\)|.*/\1/')
    elif echo "$header" | grep -q '^>tr|'; then
        protein_id=$(echo "$header" | sed 's/^>tr|\([^|]*\)|.*/\1/')
    else
        protein_id=$(echo "$header" | sed 's/^>//;s/ .*//')
    fi
    
    echo "$protein_id"
}

setup_directories() {
    log_info "Creating output directory: ${OUTPUT_DIR}"
    
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR/tmp"
    
    log_success "Output directory ready"
}

check_pdb_file() {
    log_info "Checking for PDB structure..."
    
    local pdb_file="pdb_files_poi/${POI_ID}.pdb"
    
    if [ ! -f "$pdb_file" ]; then
        log_error "PDB file not found: $pdb_file"
        echo ""
        log_info "Download structure for ${POI_ID}:"
        log_info "  wget https://alphafold.ebi.ac.uk/entry/AF-${POI_ID}-F1-model_v4.pdb -O pdb_files_poi/${POI_ID}.pdb"
        exit 1
    fi
    
    log_success "Found PDB: ${pdb_file}"
}

generate_config() {
    log_info "Generating configuration..."
    
    local captaf_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    local config_file="$OUTPUT_DIR/config.yaml"
    
    # Get absolute paths to input files
    local poi_fasta_abs="$(cd "$(dirname "$POI_FASTA")" && pwd)/$(basename "$POI_FASTA")"
    local target_fasta_abs="$(cd "$(dirname "$TARGET_FASTA")" && pwd)/$(basename "$TARGET_FASTA")"
    local pdb_file_abs="${captaf_dir}/pdb_files_poi/${POI_ID}.pdb"
    
    cat > "$config_file" << EOF
# CAPTAF Pipeline Configuration
# Generated: $(date)
# POI: ${POI_ID}
# Target: ${TARGET_ID}

# Output directory
output_dir: "${OUTPUT_DIR}"

samples:
  - ${POI_ID}

targets:
  - ${TARGET_ID}_Target

# Input file paths (absolute paths - no copying)
poi_fasta: "${poi_fasta_abs}"
target_fasta: "${target_fasta_abs}"
pdb_file: "${pdb_file_abs}"

# Output directories
output_dirs:
  ptm_predictions: "ptm_predictions"
  rsa: "rsa_outputs"
  dssp: "dssp_outputs"
  filtered: "filtered_ptms"
  jsons: "alphafold_inputs"
  target_ptms: "target_ptm_dir"

# CAPTAF installation paths
captaf_dir: "${captaf_dir}"
scripts_dir: "${captaf_dir}/scripts"
musitedeep_dir: "${captaf_dir}/MusiteDeep_web"
alphafold3_dir: "${captaf_dir}/alphafold3"
iupred2a_dir: "${captaf_dir}/iupred2a"
ptm_mapping: "${captaf_dir}/PTM_to_CCD_mapping.csv"

# Database paths
mmseqs_db: "${captaf_dir}/alphafold3/database/uniref50DB"
af3_models: "${captaf_dir}/af3_models"
af3_db: "${captaf_dir}/alphafold3/database"

# Pipeline options
use_controls: ${USE_CONTROLS}
EOF
    
    log_success "Config created: $config_file"
}

run_pipeline() {
    log_info "Starting CAPTAF pipeline..."
    echo ""
    log_info "POI: ${POI_ID}"
    log_info "Target: ${TARGET_ID}"
    log_info "POI FASTA: ${POI_FASTA}"
    log_info "Target FASTA: ${TARGET_FASTA}"
    log_info "Output: $(cd "$OUTPUT_DIR" && pwd)"
    log_info "Threads: ${THREADS}"
    echo ""
    
    local captaf_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    local output_abs="$(cd "$OUTPUT_DIR" && pwd)"
    local snakefile="${captaf_dir}/Snakefile"
    local config_file="${output_abs}/config.yaml"
    
    # Set temp directory
    export TMPDIR="${output_abs}/tmp"
    log_info "Temp directory: $TMPDIR"
    
    # Build snakemake command
    local snake_cmd="snakemake"
    snake_cmd="$snake_cmd --snakefile ${snakefile}"
# REMOVED:     snake_cmd="$snake_cmd --directory ${output_abs}"
    snake_cmd="$snake_cmd --configfile ${config_file}"
    snake_cmd="$snake_cmd --cores $THREADS"
# REMOVED:     snake_cmd="$snake_cmd --rerun-incomplete"
    
    if [ "$DRY_RUN" = true ]; then
        snake_cmd="$snake_cmd --dry-run"
        log_info "DRY RUN MODE"
    fi
    
    if [ "$VERBOSE" = true ]; then
        snake_cmd="$snake_cmd --verbose --printshellcmds"
    fi
    
    log_info "Command: $snake_cmd"
    echo ""
    echo "═══════════════════════════════════════════════════════════"
    
    eval $snake_cmd
    
    local exit_code=$?
    
    echo "═══════════════════════════════════════════════════════════"
    echo ""
    
    return $exit_code
}

show_results() {
    log_success "Pipeline completed!"
    echo ""
    echo "═══════════════════════════════════════════════════════════"
    echo -e "${CYAN}Results Summary${NC}"
    echo "═══════════════════════════════════════════════════════════"
    
    local ranking="$OUTPUT_DIR/${POI_ID}_${TARGET_ID}_ranking.tsv"
    local summary="$OUTPUT_DIR/${POI_ID}_${TARGET_ID}_summary.html"
    local recommendations="$OUTPUT_DIR/${POI_ID}_${TARGET_ID}_recommendations.txt"
    
    if [ -f "$ranking" ]; then
        echo ""
        echo -e "${YELLOW}📊 Top 5 Variants:${NC}"
        echo "───────────────────────────────────────────────────────────"
        head -6 "$ranking" | column -t -s $'\t'
        echo ""
    fi
    
    if [ -f "$recommendations" ]; then
        echo ""
        echo -e "${YELLOW}🎯 Recommendations:${NC}"
        echo "───────────────────────────────────────────────────────────"
        head -30 "$recommendations"
        echo ""
    fi
    
    echo "═══════════════════════════════════════════════════════════"
    echo -e "${GREEN}Output Files:${NC}"
    echo "  • Ranking:        ${ranking}"
    echo "  • HTML summary:   ${summary}"
    echo "  • Recommendations: ${recommendations}"
    echo ""
    echo -e "${CYAN}View report:${NC}"
    echo "  firefox ${summary}"
    echo "═══════════════════════════════════════════════════════════"
}

cleanup_intermediates() {
    if [ "$KEEP_INTERMEDIATES" = false ]; then
        log_info "Cleaning intermediates..."
        
        cd "$OUTPUT_DIR"
        rm -rf ptm_predictions/ conservation_outputs/ disordered_regions_outputs/
        rm -rf dssp_outputs/ master_tables_poi/ poi_variants/
        rm -rf target_ptm_dir/ poi_validated_ptms/
        rm -rf baseline_control_case/ positive_control_case/
        rm -rf alphafold_inputs/ tmp/
        cd - > /dev/null
        
        log_success "Cleaned up"
    fi
}

###############################################################################
# Main
###############################################################################

while [[ $# -gt 0 ]]; do
    case $1 in
        --poi) POI_FASTA="$2"; shift 2 ;;
        --target) TARGET_FASTA="$2"; shift 2 ;;
        --output) OUTPUT_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --use-controls) USE_CONTROLS=true; shift ;;
        --keep-intermediates) KEEP_INTERMEDIATES=true; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        --verbose) VERBOSE=true; shift ;;
        -v|--version) echo "CAPTAF v${VERSION}"; exit 0 ;;
        -h|--help) print_banner; print_usage; exit 0 ;;
        *) log_error "Unknown: $1"; print_usage; exit 1 ;;
    esac
done

print_banner

if [ -z "$POI_FASTA" ] || [ -z "$TARGET_FASTA" ]; then
    log_error "Missing --poi and --target"
    print_usage
    exit 1
fi

log_info "Validating inputs..."
validate_fasta "$POI_FASTA" "POI"
validate_fasta "$TARGET_FASTA" "Target"
log_success "Inputs validated"
echo ""

POI_ID=$(extract_protein_id "$POI_FASTA")
TARGET_ID=$(extract_protein_id "$TARGET_FASTA")

if [ -z "$POI_ID" ] || [ -z "$TARGET_ID" ]; then
    log_error "Could not extract protein IDs"
    exit 1
fi

log_info "POI ID: ${POI_ID}"
log_info "Target ID: ${TARGET_ID}"
echo ""

check_dependencies
echo ""

check_pdb_file
echo ""

setup_directories
echo ""

generate_config
echo ""

if run_pipeline; then
    if [ "$DRY_RUN" = false ]; then
        echo ""
        show_results
        echo ""
        cleanup_intermediates
        echo ""
        log_success "CAPTAF complete!"
    else
        log_success "Dry run successful"
    fi
    exit 0
else
    echo ""
    log_error "Pipeline failed"
    exit 1
fi
