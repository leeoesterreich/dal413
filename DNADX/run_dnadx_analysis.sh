#!/bin/bash

# DNADX Analysis Pipeline for Liquid Biopsy Shallow WGS Data
# This script runs the entire analysis pipeline

# Set up environment
echo "Setting up environment..."
mkdir -p output

# Load required modules
echo "Loading required modules..."
module load gcc/12.2.0
module load r/4.4.0

# Step 1: Prepare CNR data
echo "Step 1: Preparing CNR data..."
Rscript prepare_cnr_data.R

# Step 2: Analyze signatures
echo "Step 2: Analyzing CNA signatures..."
Rscript analyze_signatures.R

# Step 3: Generate report
echo "Step 3: Generating comprehensive report..."
Rscript generate_report.R

echo "Analysis complete! Results are available in the output directory."
echo "To view the HTML report, open output/cna_signature_report.html"
