#!/bin/bash
###############################################################################
# Install CAPTAF Conda Environments
###############################################################################

set -e

echo "=================================================="
echo "Installing CAPTAF Conda Environments"
echo "=================================================="
echo ""

if ! command -v conda &> /dev/null; then
    echo "❌ Conda not found. Please install Miniconda/Anaconda first."
    exit 1
fi

echo "Creating environments (this may take 10-30 minutes)..."
echo ""

# 1. new_thesis_env (main pipeline)
if conda env list | grep -q "^new_thesis_env "; then
    echo "⚠️  new_thesis_env already exists. Skipping..."
else
    echo "Creating new_thesis_env..."
    conda env create -f envs/new_thesis_env.yml
    echo "✅ new_thesis_env created"
fi

# 2. musite_env (PTM prediction)
if conda env list | grep -q "^musite_env "; then
    echo "⚠️  musite_env already exists. Skipping..."
else
    echo "Creating musite_env..."
    conda env create -f envs/musite_env.yml
    echo "✅ musite_env created"
fi

# 3. scraping_env (dbPTM fetching)  ← ADDED
if conda env list | grep -q "^scraping_env "; then
    echo "⚠️  scraping_env already exists. Skipping..."
else
    echo "Creating scraping_env..."
    conda env create -f envs/scraping_env.yml
    echo "✅ scraping_env created"
fi

# 4. af3 (AlphaFold3)
if conda env list | grep -q "^af3 "; then
    echo "⚠️  af3 already exists. Skipping..."
else
    echo "Creating af3..."
    conda env create -f envs/af3.yml
    echo "✅ af3 created"
fi

echo ""
echo "=================================================="
echo "All 4 environments installed successfully!"
echo "=================================================="
echo ""
echo "Verify installation:"
echo "  conda env list"
echo ""
