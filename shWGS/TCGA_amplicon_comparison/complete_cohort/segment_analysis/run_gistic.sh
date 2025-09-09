#!/bin/bash

#SBATCH --job-name=run_GISTIC
#SBATCH --output=GISTIC/run_GISTIC_%j.out
#SBATCH --error=GISTIC/run_GISTIC_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu

# This script runs GISTIC2.0 to identify significant copy number alterations.

# Load required modules
module purge
module load singularity/3.9.6
module load gistic/2.0.23

# --- Configuration ---
# Define paths to the files we generated
base_dir=$(pwd)
idc_seg_file="${base_dir}/idc_gistic.seg"
ilc_seg_file="${base_dir}/ilc_gistic.seg"
markers_file="${base_dir}/markers.txt"
ref_gene_file="/ix1/alee/LO_LAB/Personal/Daisong/Test_sample/GISTIC2.0/gistic2-2.0.23/refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat"

# Define output directories
idc_output_dir="${base_dir}/GISTIC/idc"
ilc_output_dir="${base_dir}/GISTIC/ilc"

# Create output directories
mkdir -p "$idc_output_dir"
mkdir -p "$ilc_output_dir"

# --- Run GISTIC2 for IDC ---
echo "Starting GISTIC2 analysis for IDC cohort..."

gistic2 -b "$idc_output_dir" \
            -seg "$idc_seg_file" \
            -refgene "$ref_gene_file" \
            -ta 0.3 \
            -td 0.3 \
            -conf 0.95 \
            -qvt 0.25 \
            -brlen 0.98 \
            -genegistic 1 \
            -armpeel 1 \
            -broad 1 \
            -savegene 1 \
            -gcm extreme \
            -v 30

echo "GISTIC2 IDC analysis completed. Results are in ${idc_output_dir}"
echo "--------------------------------------------------"

# --- Run GISTIC2 for ILC ---
echo "Starting GISTIC2 analysis for ILC cohort..."

gistic2 -b "$ilc_output_dir" \
            -seg "$ilc_seg_file" \
            -refgene "$ref_gene_file" \
            -ta 0.3 \
            -td 0.3 \
            -conf 0.95 \
            -qvt 0.25 \
            -brlen 0.98 \
            -genegistic 1 \
            -armpeel 1 \
            -broad 1 \
            -savegene 1 \
            -gcm extreme \
            -v 30

echo "GISTIC2 ILC analysis completed. Results are in ${ilc_output_dir}" 