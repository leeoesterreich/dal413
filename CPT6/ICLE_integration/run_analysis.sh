#!/bin/bash

# Script to run the analysis using system Python module

# Change to the correct directory
cd /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/

# Load Python module
echo "Loading Python module..."
module load python/3.7.0

# Load modules for R
echo "Loading modules for R..."
module load gcc/12.2.0
module load r/4.4.0

# Create results directory if it doesn't exist
mkdir -p results

# Set matplotlib backend to non-interactive
export MPLBACKEND=Agg

# Install required Python packages if needed
echo "Checking for required Python packages..."
pip install --user pandas numpy matplotlib seaborn

# Run the Python analysis script
echo "Running Python analysis script..."
python analyze_mutations.py

# Run the R analysis script
echo "Running R analysis script..."
Rscript analyze_signatures.R

echo "Analysis complete. Check the 'results' directory for output files." 