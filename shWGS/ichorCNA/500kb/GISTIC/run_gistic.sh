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
seg_file="/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/500kb/results/all_segments.txt"
ref_gene_file="/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/500kb/GISTIC/hg38.UCSC.add_miR.160920.refgene.mat"
output_dir="GISTIC/gistic_output"

# Create output directory
mkdir -p $output_dir

# Prepare the GISTIC segmentation file
# GISTIC requires a header and specific column names. It also requires numeric chromosome values (1-22, X=23, Y=24).
# This awk command will:
# 1. Skip the original header (NR>1).
# 2. Reformat the columns to match GISTIC's requirements.
# 3. Convert 'chrX' to '23' and 'chrY' to '24' and remove 'chr' prefix.
# 4. Save to a new file in the output directory.
gistic_seg_file="${output_dir}/input.seg"

echo -e "Sample\tChromosome\tStart Position\tEnd Position\tNum Markers\tSeg.CN" > $gistic_seg_file
tail -n +2 "$seg_file" | awk 'BEGIN{FS=OFS="\t"} {gsub("chr", "", $2); if ($2 == "X") $2 = 23; if ($2 == "Y") $2 = 24; print $1, $2, $3, $4, $5, $6}' >> $gistic_seg_file


# GISTIC Parameters
# The -b parameter specifies the output directory.
output_dir_base="GISTIC/gistic_run"

# Create the output directory
mkdir -p $output_dir_base

# --- Run GISTIC2 ---
echo "Starting GISTIC2 analysis..."

gistic2 -b $output_dir_base \
            -seg $gistic_seg_file \
            -refgene $ref_gene_file \
            -ta 0.1 \
            -td 0.1 \
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