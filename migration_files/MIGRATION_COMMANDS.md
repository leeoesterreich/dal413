# CITEgeist Documentation Migration Commands

## Phase 2: Exact Git Commands

Follow these commands step by step to migrate your documentation to the CITEgeist repository.

### Step 2.1: Prepare Local Environment

```bash
# Navigate to a clean directory for cloning
cd ~/Desktop  # or wherever you want to work

# Clone the CITEgeist repository
git clone https://github.com/leeoesterreich/CITEgeist.git
cd CITEgeist

# Switch to dev branch
git checkout dev

# Verify you're on the right branch
git branch
```

### Step 2.2: Create Directory Structure

```bash
# Create the new directory structure
mkdir -p docs/core_scripts/model
mkdir -p docs/core_scripts/Jupyter  
mkdir -p docs/source/notebooks
mkdir -p docs/source/_static
mkdir -p docs/source/_templates
mkdir -p docs/build

# Verify structure
tree -L 3  # or ls -la to check
```

### Step 2.3: Copy Python Source Files

```bash
# Copy Python modules (rename CITEgeist -> core_scripts)
cp /ix1/alee/LO_LAB/Personal/Daisong/Projects/Github/dal413/CITEgeist/CITEgeist/model/*.py docs/core_scripts/model/

# Copy Jupyter notebooks
cp /ix1/alee/LO_LAB/Personal/Daisong/Projects/Github/dal413/CITEgeist/CITEgeist/Jupyter/*.ipynb docs/core_scripts/Jupyter/

# Verify files copied
ls -la docs/core_scripts/model/
ls -la docs/core_scripts/Jupyter/
```

### Step 2.4: Copy Documentation Files

```bash
# Copy updated configuration files from migration_files
cp /ix1/alee/LO_LAB/Personal/Daisong/Projects/Github/dal413/migration_files/conf.py docs/source/
cp /ix1/alee/LO_LAB/Personal/Daisong/Projects/Github/dal413/migration_files/requirements.txt .
cp /ix1/alee/LO_LAB/Personal/Daisong/Projects/Github/dal413/migration_files/postBuild .

# Copy documentation RST files
cp /ix1/alee/LO_LAB/Personal/Daisong/Projects/Github/dal413/CITEgeist/docs/source/*.rst docs/source/
cp /ix1/alee/LO_LAB/Personal/Daisong/Projects/Github/dal413/migration_files/notebooks.rst docs/source/

# Copy notebook documentation  
cp /ix1/alee/LO_LAB/Personal/Daisong/Projects/Github/dal413/CITEgeist/docs/source/notebooks/*.ipynb docs/source/notebooks/

# Copy ALL updated vignette RST files with corrected Binder links
cp /ix1/alee/LO_LAB/Personal/Daisong/Projects/Github/dal413/migration_files/vignette_1_biopsy_heterogeneity.rst docs/source/notebooks/
cp /ix1/alee/LO_LAB/Personal/Daisong/Projects/Github/dal413/migration_files/vignette_2_surgical_d538g.rst docs/source/notebooks/
cp /ix1/alee/LO_LAB/Personal/Daisong/Projects/Github/dal413/migration_files/vignette_3_responder_macrophages.rst docs/source/notebooks/
cp /ix1/alee/LO_LAB/Personal/Daisong/Projects/Github/dal413/migration_files/vignette_4_external_tools.rst docs/source/notebooks/

# Copy the enhanced API documentation
cp /ix1/alee/LO_LAB/Personal/Daisong/Projects/Github/dal413/CITEgeist/docs/source/api.rst docs/source/
```

### Step 2.5: Verify All Binder Links Are Updated

```bash
# All vignette files now have correct Binder links! 
# The migration_files contain updated RST files with proper links:
# - vignette_1_biopsy_heterogeneity.rst ✅
# - vignette_2_surgical_d538g.rst ✅  
# - vignette_3_responder_macrophages.rst ✅
# - vignette_4_external_tools.rst ✅

# Verify the links are correct
echo "Checking Binder links in all vignettes..."
grep -r "mybinder.org" docs/source/notebooks/
grep -r "colab.research.google.com" docs/source/notebooks/

# Should show links pointing to:
# - Repository: leeoesterreich/CITEgeist 
# - Branch: dev
# - Path: docs/core_scripts/Jupyter/
```

### Step 2.6: Test Documentation Build

```bash
# Try building the documentation (may need to install Sphinx)
cd docs
make html

# Check if build succeeded
ls -la build/html/

# If build fails, check error messages and fix paths
```

### Step 2.7: Commit and Push Changes

```bash
# Return to repository root
cd ..

# Check what files were added
git status

# Add all new files
git add .

# Commit with descriptive message
git commit -m "Migrate comprehensive documentation from dal413 repository

Features:
- Add docs/core_scripts/ containing CITEgeist Python modules
- Include enhanced API documentation with workflow organization  
- Add interactive Jupyter notebooks with Binder/Colab integration
- Set up complete Sphinx documentation system
- Configure proper repository structure for documentation
- Update all links to point to CITEgeist repository

This migration consolidates all CITEgeist code and documentation 
into a single, well-organized repository following best practices."

# Push to dev branch
git push origin dev
```

### Step 2.8: Verify Migration

```bash
# Check that everything was pushed
git log --oneline -3

# Verify repository structure
tree -L 3

# Test that Python imports work (if Python available)
cd docs/core_scripts
python -c "from model.citegeist_model import CitegeistModel; print('Import successful!')"
```

## Next: Create Pull Request

After completing these commands, go to GitHub and create a pull request from the `dev` branch to `main` branch.
