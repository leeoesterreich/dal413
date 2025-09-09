#!/bin/bash

# Survival Analysis for Top 10 CNA Signatures
# This script runs the survival analysis for the top 10 signatures

# Set up environment
echo "Setting up environment..."

# Load required modules
echo "Loading required modules..."
module load gcc/12.2.0
module load r/4.4.0

# Run the survival analysis
echo "Running survival analysis..."
Rscript survival_analysis.R

echo "Analysis complete! Results are available in the current directory." 