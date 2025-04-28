#!/bin/bash

# Script to set up the conda environment and run the analysis

# Change to the correct directory
cd /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/

# Define conda path and environment
CONDA_PATH="/bgfs/alee/LO_LAB/Personal/Daisong/miniconda3"
ENV_PREFIX="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/env_ICLE"

# Create or activate the environment
if [ ! -d "$ENV_PREFIX" ]; then
    echo "Creating environment at $ENV_PREFIX..."
    source "${CONDA_PATH}/etc/profile.d/conda.sh"
    conda create -y --prefix $ENV_PREFIX python=3.9 pandas numpy matplotlib seaborn pip
    $ENV_PREFIX/bin/pip install gseapy
else
    echo "Environment exists at $ENV_PREFIX"
fi

# Activate the environment
source "${CONDA_PATH}/etc/profile.d/conda.sh"
conda activate $ENV_PREFIX

# Load modules for R
module load gcc/12.2.0
module load r/4.4.0

# Run the Python analysis script
echo "Running Python analysis script..."
$ENV_PREFIX/bin/python analyze_mutations.py

# Run the R analysis script
echo "Running R analysis script..."
Rscript analyze_signatures.R

echo "Analysis complete. Check the 'results' directory for output files."