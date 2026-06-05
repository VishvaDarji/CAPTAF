#!/bin/bash
###############################################################################
# CAPTAF: Context-Aware Snakemake Pipeline for PTM-Driven Binding
#         Optimization with AlphaFold3
# 
# Setup Script - Run once to verify installation and configure databases
###############################################################################

set -e

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║                      CAPTAF Installation Setup                   ║"
echo "║   Context-Aware Snakemake Pipeline for PTM-Driven Binding       ║"
echo "║              Optimization with AlphaFold3                        ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""

# Check conda
echo -e "${BLUE}[1/5]${NC} Checking Conda installation..."
if ! command -v conda &> /dev/null; then
    echo -e "${RED}✗${NC} Conda not found"
    echo "    Please install Miniconda/Anaconda first"
    echo "    Download: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
else
    echo -e "${GREEN}✓${NC} Conda found: $(conda --version)"
fi
echo ""

# Check Snakemake
echo -e "${BLUE}[2/5]${NC} Checking Snakemake..."
if conda run -n new_thesis_env snakemake --version &> /dev/null; then
    VERSION=$(conda run -n new_thesis_env snakemake --version 2>/dev/null)
    echo -e "${GREEN}✓${NC} Snakemake found: $VERSION"
else
    echo -e "${RED}✗${NC} Snakemake not found in new_thesis_env"
    echo "    Install with: conda install -n new_thesis_env snakemake"
    exit 1
fi
echo ""

# Check conda environments
echo -e "${BLUE}[3/5]${NC} Checking conda environments..."
ALL_PRESENT=true

for env in new_thesis_env musite_env af3; do
    if conda env list | grep -q "^${env} "; then
        echo -e "${GREEN}✓${NC} ${env}"
    else
        echo -e "${RED}✗${NC} ${env} ${YELLOW}(NOT FOUND)${NC}"
        echo "    Run: ./envs/install_envs.sh"
        ALL_PRESENT=false
    fi
done

if [ "$ALL_PRESENT" = false ]; then
    exit 1
fi
echo ""

# Check database paths
echo -e "${BLUE}[4/5]${NC} Checking database locations..."

echo "Database paths (can be configured in config.yaml):"
echo ""
echo "  ${YELLOW}MMseqs2 database:${NC}"
if [ -f "alphafold3/database/uniref50DB" ]; then
    echo -e "    ${GREEN}✓${NC} alphafold3/database/uniref50DB"
else
    echo -e "    ${YELLOW}⚠${NC} alphafold3/database/uniref50DB ${YELLOW}(not found)${NC}"
    echo "      You can specify a different path in config.yaml"
fi

if [ -f "alphafold3/database/uniref50DB.index" ]; then
    echo -e "    ${GREEN}✓${NC} alphafold3/database/uniref50DB.index"
else
    echo -e "    ${YELLOW}⚠${NC} alphafold3/database/uniref50DB.index ${YELLOW}(not found)${NC}"
    echo "      Create index with:"
    echo "      cd alphafold3/database && mmseqs createindex uniref50DB tmp"
fi

echo ""
echo "  ${YELLOW}AlphaFold3 models:${NC}"
if [ -d "af3_models" ]; then
    echo -e "    ${GREEN}✓${NC} af3_models/"
else
    echo -e "    ${YELLOW}⚠${NC} af3_models/ ${YELLOW}(not found)${NC}"
    echo "      Download AlphaFold3 models as per AF3 installation guide"
fi

echo ""
echo "  ${YELLOW}AlphaFold3 databases:${NC}"
if [ -d "alphafold3/database" ]; then
    echo -e "    ${GREEN}✓${NC} alphafold3/database/"
else
    echo -e "    ${YELLOW}⚠${NC} alphafold3/database/ ${YELLOW}(not found)${NC}"
fi
echo ""

# Check external tools
echo -e "${BLUE}[5/5]${NC} Checking external tools..."

if [ -f "iupred2a/iupred2a.py" ]; then
    echo -e "${GREEN}✓${NC} IUPred2A"
else
    echo -e "${RED}✗${NC} IUPred2A not found"
    echo "    Expected: iupred2a/iupred2a.py"
    exit 1
fi

if [ -d "MusiteDeep_web/MusiteDeep/models" ]; then
    MODEL_COUNT=$(ls MusiteDeep_web/MusiteDeep/models/*.pt 2>/dev/null | wc -l)
    echo -e "${GREEN}✓${NC} MusiteDeep (${MODEL_COUNT} models found)"
else
    echo -e "${RED}✗${NC} MusiteDeep not found"
    echo "    Expected: MusiteDeep_web/MusiteDeep/models/"
    exit 1
fi
echo ""

# Summary
echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║                         Setup Summary                            ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""
echo -e "${GREEN}✓${NC} CAPTAF is ready to use!"
echo ""
echo "Quick start:"
echo "  1. Place your PDB structure in: pdb_files_poi/"
echo "  2. Run analysis:"
echo "     ./captaf.sh --poi protein.fasta --target partner.fasta"
echo ""
echo "For help:"
echo "  ./captaf.sh --help"
echo ""
echo "To verify all components:"
echo "  ./verify_captaf.sh"
echo ""

