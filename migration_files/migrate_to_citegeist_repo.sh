#!/bin/bash

# CITEgeist Documentation Migration Script
# This script helps migrate documentation from dal413 to the dedicated CITEgeist repository

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_status "CITEgeist Documentation Migration Guide"
print_status "========================================"
echo ""

print_status "This script will help you migrate your documentation to the CITEgeist repository."
print_status "You'll need to run these commands manually on the target repository."
echo ""

print_status "STEP 1: Clone the CITEgeist repository and switch to dev branch"
echo "git clone https://github.com/leeoesterreich/CITEgeist.git"
echo "cd CITEgeist"
echo "git checkout dev"
echo ""

print_status "STEP 2: Create the new directory structure"
echo "mkdir -p core_scripts/model"
echo "mkdir -p core_scripts/Jupyter"
echo "mkdir -p docs/source/notebooks"
echo "mkdir -p docs/source/_static"
echo "mkdir -p docs/source/_templates"
echo ""

print_status "STEP 3: Copy Python source files (rename CITEgeist -> core_scripts)"
echo "# Copy from your current dal413 repository:"
echo "cp -r /path/to/dal413/CITEgeist/CITEgeist/model/* core_scripts/model/"
echo "cp -r /path/to/dal413/CITEgeist/CITEgeist/Jupyter/* core_scripts/Jupyter/"
echo ""

print_status "STEP 4: Copy documentation files"
echo "# Copy updated configuration files:"
echo "cp /path/to/migration_files/conf.py docs/source/"
echo "cp /path/to/migration_files/requirements.txt ."
echo "cp /path/to/migration_files/postBuild ."
echo ""
echo "# Copy all RST files:"
echo "cp -r /path/to/dal413/CITEgeist/docs/source/*.rst docs/source/"
echo "cp /path/to/migration_files/notebooks.rst docs/source/"
echo "cp /path/to/migration_files/vignette_1_biopsy_heterogeneity.rst docs/source/notebooks/"
echo ""
echo "# Copy remaining vignette RST files and update their Binder links"
echo ""

print_status "STEP 5: Update remaining vignette files"
echo "# You'll need to update the Binder/Colab links in:"
echo "# - vignette_2_surgical_d538g.rst"
echo "# - vignette_3_responder_macrophages.rst" 
echo "# - vignette_4_external_tools.rst"
echo ""
echo "# Replace URLs from:"
echo "# https://mybinder.org/v2/gh/leeoesterreich/dal413/HEAD?labpath=CITEgeist%2FCITEgeist%2FJupyter%2F"
echo "# TO:"
echo "# https://mybinder.org/v2/gh/leeoesterreich/CITEgeist/dev?labpath=core_scripts%2FJupyter%2F"
echo ""

print_status "STEP 6: Test the documentation build"
echo "cd docs"
echo "make html"
echo "# Check that build/html/ contains your documentation"
echo ""

print_status "STEP 7: Set up GitHub Pages (if needed)"
echo "# In GitHub repository settings, enable GitHub Pages from docs/ folder"
echo "# Or set up GitHub Actions for automated deployment"
echo ""

print_status "STEP 8: Commit and push"
echo "git add ."
echo "git commit -m \"Migrate documentation from dal413 repository"
echo ""
echo "- Add core_scripts/ with CITEgeist Python modules"
echo "- Add comprehensive Sphinx documentation"
echo "- Include interactive Jupyter notebooks with Binder/Colab integration"
echo "- Set up proper repository structure for documentation\""
echo ""
echo "git push origin dev"
echo ""

print_status "STEP 9: Create Pull Request"
echo "# Go to GitHub and create a pull request from dev to main branch"
echo "# Include description of the migration and new features"
echo ""

print_success "Migration preparation complete!"
print_status "All necessary files have been created in the migration_files/ directory"
print_status "Follow the steps above to complete the migration to the CITEgeist repository"

