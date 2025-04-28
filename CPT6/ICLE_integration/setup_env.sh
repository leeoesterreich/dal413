#!/bin/bash

# Script to set up the environment using system modules

# Change to the correct directory
cd /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/

# Load Python module
echo "Loading Python module..."
module load python/3.7.0

# Install required Python packages
echo "Installing Python packages..."
pip install --user pandas numpy matplotlib seaborn

# Load modules for R
echo "Loading modules for R..."
module load gcc/12.2.0
module load r/4.4.0

# Create results directory if it doesn't exist
mkdir -p results

# Install R packages required for maftools
echo "Installing R packages..."
R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='http://cran.us.r-project.org'); BiocManager::install('maftools')"
R -e "if (!requireNamespace('BSgenome.Hsapiens.UCSC.hg38', quietly = TRUE)) BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')"
R -e "if (!requireNamespace('ggplot2', quietly = TRUE)) install.packages('ggplot2', repos='http://cran.us.r-project.org')"

echo "Environment setup complete."
echo "To run the analysis, execute: ./run_analysis.sh" 