#!/bin/bash

# This is the quality control script for ctDNA NGS raw data

# Author=Daisong_Liu/LeeOesterreich_lab

#SBATCH --job-name=FastQC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --time=48:00:00

# Loading modules
module purge
module load fastqc/0.11.7
module load multiqc/1.19

# Input directory
input_dir="/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/Raw_file"

# Output directory
output_dir="/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/1_fastqc_result"

# Traversing fastq files in input directory and running FastQC in parallel
for fastq_file in $input_dir/*.fastq.gz; do
    # Extract the filename without extension
    sample=$(basename "$fastq_file" .fastq.gz)

    # Check if the FastQC report already exists
    if [ -f "$output_dir/${sample}_fastqc.html" ]; then
        echo "FastQC report already exists for $sample, skipping..."
        continue
    fi

    # Run FastQC in the background
    fastqc -o "$output_dir" "$fastq_file" > "$output_dir/${sample}_fastqc.log" &

    echo "FastQC analysis launched for $sample"
done

# Wait for all background jobs to finish
wait

echo "All FastQC analyses completed"
   
# MultiQC analysis
# Define input directory containing analysis results
input_dir1="/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/1_fastqc_result"

# Define output directory where MultiQC report will be saved
output_dir1="/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/1_fastqc_result"

# Run MultiQC on the input directory
multiqc "$input_dir1" -o "$output_dir1"

echo "MultiQC analysis completed"
