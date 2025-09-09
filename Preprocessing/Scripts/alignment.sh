#!/bin/bash
# ==============================================================================
#
#            APOLLO Project - Preprocessing Pipeline: Alignment
#
# This script performs BWA-MEM alignment of paired-end FASTQ files to a 
# reference genome, followed by SAM to BAM conversion, sorting, indexing,
# and quality control using Qualimap and MultiQC.
#
# Author: Daisong Liu / Lee Oestereich Lab
# Created: [Date]
# Last Modified: [Date]
# Version: 1.0
#
# Dependencies:
#   - BWA (v0.7.17)
#   - SAMtools (v1.14)
#   - Qualimap (v2.2.2)
#   - MultiQC (v1.19)
#
# Usage:
#   sbatch alignment.sh
#
# Git Workflow:
#   This script is part of the APOLLO project preprocessing pipeline.
#   Before making changes:
#   1. Create a feature branch: git checkout -b feature/update-alignment
#   2. Make your changes and test thoroughly
#   3. Commit: git commit -m "feat: improve alignment pipeline error handling"
#   4. Push: git push -u origin feature/update-alignment
#   5. Create a pull request for review
#
# ==============================================================================

# BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome
# BWA-MEM, for high-quality queries and better performance, as it is faster and more accurate
#SBATCH --job-name=Alignment
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 96:00:00
#SBATCH --mem=256G

# ==============================================================================
# CONFIGURATION SECTION
# ==============================================================================

# Error handling: Exit on any error
set -e
set -u
set -o pipefail

# Enable debugging (uncomment for verbose output)
# set -x

echo "=== APOLLO Alignment Pipeline Started at $(date) ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODEID"

# Module required
echo "Loading required modules..."
module purge
module load gcc/8.2.0
module load bwa/0.7.17
module load samtools/1.14
module load qualimap/2.2.2
module load multiqc/1.19

echo "Modules loaded successfully."

# ==============================================================================
# PATH CONFIGURATION
# ==============================================================================
# NOTE: Update these paths according to your project structure
# For production use, consider using environment variables or config files

# Base project directory
BASE_DIR="/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample"

# Reference genome file
REFERENCE_GENOME="/bgfs/alee/LO_LAB/Personal/Daisong/Reference_genome/Human_ref_hg38_p7/Human_genome_hg38_p7_chr_only.fna"

# Input directory for trimmed FASTQ files
INPUT_DIR="${BASE_DIR}/2_trimmed_file/paired"

# Output directories
ALIGNMENT_DIR="${BASE_DIR}/4_alignment_file/A5"
QUALIMAP_DIR="${ALIGNMENT_DIR}/Qualimap"

# Create output directories
echo "Creating output directories..."
mkdir -p "$ALIGNMENT_DIR"
mkdir -p "$QUALIMAP_DIR"

# Verify reference genome exists
if [ ! -f "$REFERENCE_GENOME" ]; then
    echo "ERROR: Reference genome file not found: $REFERENCE_GENOME"
    exit 1
fi

# Verify input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

echo "Directory structure verified successfully."

# ==============================================================================
# SAMPLE CONFIGURATION
# ==============================================================================
# TODO: Update sample list or read from external file for production use
samples=( "TP19-M480_FOL6151A5_S12" )

echo "Processing ${#samples[@]} samples: ${samples[*]}"

# ==============================================================================
# ALIGNMENT PROCESSING
# ==============================================================================

for sample in "${samples[@]}"; do
    echo ""
    echo "=== Processing sample: $sample ==="

    # Input FASTQ files for the sample
    forward_read="${INPUT_DIR}/${sample}_R1_paired.fastq.gz"
    reverse_read="${INPUT_DIR}/${sample}_R2_paired.fastq.gz"
    
    # Verify input files exist
    if [ ! -f "$forward_read" ]; then
        echo "ERROR: Forward read file not found: $forward_read"
        exit 1
    fi
    
    if [ ! -f "$reverse_read" ]; then
        echo "ERROR: Reverse read file not found: $reverse_read"
        exit 1
    fi

    # Output SAM file for the sample
    output_sam="${ALIGNMENT_DIR}/${sample}.sam"  

    # Check if SAM file already exists
    if [ -f "$output_sam" ]; then
        echo "SAM file for $sample already exists, skipping alignment..."
    else     
        # Perform mapping using bwa-mem
        echo "Aligning sample: $sample"
        echo "Forward read: $forward_read"
        echo "Reverse read: $reverse_read"
        echo "Output SAM: $output_sam"
        
        bwa mem -K 100000000 -t 64 "$REFERENCE_GENOME" "$forward_read" "$reverse_read" > "$output_sam"
        
        if [ $? -ne 0 ]; then
            echo "ERROR: BWA alignment failed for $sample"
            exit 1
        fi      
        echo "Alignment completed successfully for $sample"
    fi  
done

echo "All samples were aligned successfully."

# ==============================================================================
# SAM TO BAM CONVERSION
# ==============================================================================
# samtools is a suite of programs for interacting with high-throughput sequencing data. 
# samtools view used to convert a sam file to a bam file
# samtools sort used to sort the sam/bam file by coordinates.
# samtools index used to create an index file for the sorted bam file.  

echo ""
echo "=== Converting SAM files to BAM format ==="

# Traverse each SAM file in the alignment directory
for sam_file in "$ALIGNMENT_DIR"/*.sam; do
    # Check if any SAM files exist
    if [ ! -f "$sam_file" ]; then
        echo "No SAM files found in $ALIGNMENT_DIR"
        break
    fi

    # Extract sample name from the SAM file
    sample=$(basename "$sam_file" .sam)
    
    echo "Processing sample: $sample"

    # Output BAM file
    bam_file="${ALIGNMENT_DIR}/${sample}.bam"

    # Convert SAM to BAM
    echo "Converting $sam_file to BAM format..."
    samtools view -bS "$sam_file" > "$bam_file"
    
    if [ $? -ne 0 ]; then
        echo "ERROR: SAM to BAM conversion failed for $sample"
        exit 1
    fi

    echo "BAM file created: $bam_file"
done

echo "BAM files generated for all SAM files."

# ==============================================================================
# BAM SORTING
# ==============================================================================

echo ""
echo "=== Sorting BAM files ==="

# Traverse each BAM file in the alignment directory
for bam_file in "$ALIGNMENT_DIR"/*.bam; do
    # Check if any BAM files exist
    if [ ! -f "$bam_file" ]; then
        echo "No BAM files found in $ALIGNMENT_DIR"
        break
    fi

    # Extract sample name from the BAM file
    sample=$(basename "$bam_file" .bam)
    
    echo "Processing sample: $sample"

    # Output Sorted BAM file for the sample
    sorted_bam="${ALIGNMENT_DIR}/${sample}_sorted.bam"
    
    # Sort BAM file
    echo "Sorting $bam_file..."
    samtools sort "$bam_file" -o "$sorted_bam"
    
    if [ $? -ne 0 ]; then
        echo "ERROR: BAM sorting failed for $sample"
        exit 1
    fi

    echo "Sorted BAM file created: $sorted_bam"
done

echo "All samples were sorted successfully."

# ==============================================================================
# BAM INDEXING
# ==============================================================================

echo ""
echo "=== Indexing sorted BAM files ==="

# Traverse each sorted BAM file in the alignment directory
for sorted_bam_file in "$ALIGNMENT_DIR"/*_sorted.bam; do
    # Check if any sorted BAM files exist
    if [ ! -f "$sorted_bam_file" ]; then
        echo "No sorted BAM files found in $ALIGNMENT_DIR"
        break
    fi

    # Extract sample name from the sorted BAM file
    sample=$(basename "$sorted_bam_file" _sorted.bam)
    
    echo "Processing sample: $sample"

    # Index BAM file
    echo "Indexing $sorted_bam_file..."
    samtools index "$sorted_bam_file"
    
    if [ $? -ne 0 ]; then
        echo "ERROR: BAM indexing failed for $sample"
        exit 1
    fi

    echo "Index file created for $sample"
done

echo "All samples were converted, sorted, and indexed successfully."

# ==============================================================================
# QUALITY CONTROL WITH QUALIMAP
# ==============================================================================

echo ""
echo "=== Running Qualimap quality control ==="

# Traverse each sorted BAM file for quality control
for bam_file in "$ALIGNMENT_DIR"/*_sorted.bam; do
    # Check if any sorted BAM files exist
    if [ ! -f "$bam_file" ]; then
        echo "No sorted BAM files found for quality control"
        break
    fi

    # Extract sample name from the BAM file
    sample=$(basename "$bam_file" _sorted.bam)
    
    echo "Running quality control for sample: $sample"

    # Create an output directory for the sample
    sample_output_dir="${QUALIMAP_DIR}/${sample}"
    mkdir -p "$sample_output_dir"

    # Run Qualimap
    qualimap bamqc \
        -bam "$bam_file" \
        -outdir "$sample_output_dir" \
        -outfile "${sample}_qualimap_report.pdf" \
        -outformat pdf \
        --java-mem-size=4G

    if [ $? -ne 0 ]; then
        echo "ERROR: Qualimap analysis failed for $sample"
        exit 1
    fi

    echo "Qualimap report generated for $sample"
done

echo "Qualimap analysis completed for all samples."

# ==============================================================================
# MULTIQC SUMMARY REPORT
# ==============================================================================

echo ""
echo "=== Generating MultiQC summary report ==="

# Run MultiQC on the Qualimap output directory
multiqc "$QUALIMAP_DIR" -o "$QUALIMAP_DIR"

if [ $? -ne 0 ]; then
    echo "ERROR: MultiQC analysis failed"
    exit 1
fi

echo "MultiQC analysis completed successfully."

# ==============================================================================
# PIPELINE COMPLETION
# ==============================================================================

echo ""
echo "=== APOLLO Alignment Pipeline Completed Successfully at $(date) ==="
echo "Results are available in: $ALIGNMENT_DIR"
echo "Quality control reports are available in: $QUALIMAP_DIR"
echo ""


