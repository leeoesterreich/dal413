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
# Define file paths
seg_file="/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb_Broad/results/all_segments_seg.txt"
ref_gene_file="/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb_Broad/GISTIC/hg38.UCSC.add_miR.160920.refgene.mat"

# Create output directory
output_dir_base="GISTIC/gistic_run"
mkdir -p $output_dir_base

# --- Run GISTIC2 ---
echo "Starting GISTIC2 analysis..."

gistic2 -b $output_dir_base \
            -seg $seg_file \
            -refgene $ref_gene_file \
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

# Print when gistic is finished
echo "GISTIC2 analysis completed. Results are in ${output_dir_base}" 