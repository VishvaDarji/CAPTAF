#!/bin/bash
###############################################################################
# CAPTAF: Context-Aware Snakemake Pipeline for PTM-Driven Binding
#         Optimization with AlphaFold3
# 
# Component Verification Script - Check if all required files are present
###############################################################################

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

MISSING=0
TOTAL=0

check_file() {
    TOTAL=$((TOTAL + 1))
    if [ -f "$1" ]; then
        echo -e "${GREEN}✓${NC} $1"
    else
        echo -e "${RED}✗${NC} $1 ${YELLOW}(MISSING)${NC}"
        MISSING=$((MISSING + 1))
    fi
}

check_dir() {
    TOTAL=$((TOTAL + 1))
    if [ -d "$1" ]; then
        echo -e "${GREEN}✓${NC} $1/"
    else
        echo -e "${RED}✗${NC} $1/ ${YELLOW}(MISSING)${NC}"
        MISSING=$((MISSING + 1))
    fi
}

check_executable() {
    TOTAL=$((TOTAL + 1))
    if [ -x "$1" ]; then
        echo -e "${GREEN}✓${NC} $1 ${GREEN}(executable)${NC}"
    elif [ -f "$1" ]; then
        echo -e "${YELLOW}⚠${NC} $1 ${YELLOW}(exists but not executable)${NC}"
        echo "     Run: chmod +x $1"
        MISSING=$((MISSING + 1))
    else
        echo -e "${RED}✗${NC} $1 ${YELLOW}(MISSING)${NC}"
        MISSING=$((MISSING + 1))
    fi
}

echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║                    CAPTAF Component Verification                 ║"
echo "║   Context-Aware Snakemake Pipeline for PTM-Driven Binding       ║"
echo "║              Optimization with AlphaFold3                        ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""

echo "Core Pipeline Files:"
echo "───────────────────────────────────────────────────────────"
check_file "Snakefile"
check_file "config.yaml.template"
check_file "PTM_to_CCD_mapping.csv"
echo ""

echo "Wrapper Scripts:"
echo "───────────────────────────────────────────────────────────"
check_executable "captaf.sh"
check_executable "setup_captaf.sh"
check_executable "verify_captaf.sh"
echo ""

echo "Python Scripts (scripts/):"
echo "───────────────────────────────────────────────────────────"
check_file "scripts/run_musitedeep_all.py"
check_file "scripts/convert_musitedeep_to_tsv.py"
check_file "scripts/predict_dssp.py"
check_file "scripts/merge_dssp_ptm.py"
check_file "scripts/predict_conservation.py"
check_file "scripts/predict_disorder.py"
check_file "scripts/create_master_table.py"
check_file "scripts/generate_variants.py"
check_file "scripts/fetch_dbptm_ptms.py"
check_file "scripts/process_dbptm_ptms.py"
check_file "scripts/generate_af3_json.py"
check_file "scripts/generate_baseline_control.py"
check_file "scripts/generate_positive_control.py"
check_file "scripts/analyze_af3_results.py"
echo ""

echo "Conda Environment Files:"
echo "───────────────────────────────────────────────────────────"
check_file "envs/new_thesis_env.yml"
check_file "envs/musite_env.yml"
check_file "envs/af3.yml"
check_executable "envs/install_envs.sh"
echo ""

echo "External Tools:"
echo "───────────────────────────────────────────────────────────"
check_dir "iupred2a"
check_file "iupred2a/iupred2a.py"
check_dir "MusiteDeep_web"
check_dir "MusiteDeep_web/MusiteDeep/models"
echo ""

echo "Databases (checking existence):"
echo "───────────────────────────────────────────────────────────"
check_file "alphafold3/database/uniref50DB"
check_file "alphafold3/database/uniref50DB.index"
check_dir "af3_models"
check_dir "alphafold3/database"
echo ""

echo "Conda Environments (checking if installed):"
echo "───────────────────────────────────────────────────────────"
if conda env list | grep -q "^new_thesis_env "; then
    echo -e "${GREEN}✓${NC} new_thesis_env"
else
    echo -e "${RED}✗${NC} new_thesis_env ${YELLOW}(NOT INSTALLED)${NC}"
    echo "     Run: ./envs/install_envs.sh"
    MISSING=$((MISSING + 1))
fi
TOTAL=$((TOTAL + 1))

if conda env list | grep -q "^musite_env "; then
    echo -e "${GREEN}✓${NC} musite_env"
else
    echo -e "${RED}✗${NC} musite_env ${YELLOW}(NOT INSTALLED)${NC}"
    echo "     Run: ./envs/install_envs.sh"
    MISSING=$((MISSING + 1))
fi
TOTAL=$((TOTAL + 1))

if conda env list | grep -q "^af3 "; then
    echo -e "${GREEN}✓${NC} af3"
else
    echo -e "${RED}✗${NC} af3 ${YELLOW}(NOT INSTALLED)${NC}"
    echo "     Run: ./envs/install_envs.sh"
    MISSING=$((MISSING + 1))
fi
TOTAL=$((TOTAL + 1))

echo ""
echo "═══════════════════════════════════════════════════════════"
echo "Summary:"
echo "═══════════════════════════════════════════════════════════"
echo "Total items checked: $TOTAL"
echo "Items present: $((TOTAL - MISSING))"
echo "Items missing: $MISSING"
echo ""

if [ $MISSING -eq 0 ]; then
    echo -e "${GREEN}✅ All required components are present!${NC}"
    echo ""
    echo "CAPTAF is ready to use. Try:"
    echo "  ./captaf.sh --help"
    exit 0
else
    echo -e "${RED}⚠️  Missing $MISSING required components${NC}"
    echo ""
    echo "Please install missing components before using CAPTAF."
    echo "See CAPTAF_CHECKLIST.md for details."
    exit 1
fi
