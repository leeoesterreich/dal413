#!/bin/bash

# This script is used to mark duplicate reads and recalibrate the base score.

#SBATCH --job-name=MarkDuplicates_BaseRecalibration
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 1-00:00:00


# Delete all existing modules
module purge

# Load required modules
module load gcc/8.2.0
module load gatk/4.5.0.0
Picard=/ihome/crc/install/picard/2.18.12
java=/ihome/crc/install/java/jdk-21.0.2/bin/java

# Create a new directory
designated_path="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6"
mkdir -p "$designated_path/6_recalibration"

# Set input  directory
input_dir1="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/5_mkdup"
input_dir2="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/6_recalibration"

# Set output directory
output_dir="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/6_recalibration"

# Set reference genome path
reference="/bgfs/alee/LO_LAB/Personal/Daisong/Reference_genome/Mouse/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna"

# Set known sites for BQSR
known_sites="/bgfs/alee/LO_LAB/Personal/Daisong/Reference_genome/Mouse"


# Iterate over each rg_mkdp_bam file in the input directory
for rg_mkdp_bam in "$input_dir1"/*rg.mkdp.bam; do

    # Extract sample name from the rg_mkdp_bam file name
    sample=$(basename "$rg_mkdp_bam" rg.mkdp.bam)
    
    # Output recalibration table for the sample
    recal_table="${output_dir}/${sample}.recal.table"
    
    # Output recalibrated BAM file
    recal_bam="${output_dir}/${sample}.recal.bam"
    
    # Check if recalibrated BAM file already exists
    if [[ -f "$recal_bam" && -f "$recal_table" ]]; then
        echo "Recalibrated BAM and table already exist for $sample. Skipping BQSR."
        continue
    fi
    
    # Run BaseRecalibrator utility
    java -jar /ihome/crc/install/gatk/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar BaseRecalibrator \
        -I "$rg_mkdp_bam" \
        -R "$reference" \
        --known-sites "$known_sites"/mgp_REL2021_indels.rsID.vcf.gz \
        -O "$recal_table" \
    
    # Check if BaseRecalibrator succeeded
    if [ $? -ne 0 ]; then
        echo "Error: BaseRecalibrator failed for $sample"
        exit 1  # Exit the script on error
    else
        echo "Recalibration table created for $sample"
    fi
    
    # Apply BQSR utility
    java -jar /ihome/crc/install/gatk/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar ApplyBQSR \
        -I "$rg_mkdp_bam" \
        -R "$reference" \
        --bqsr-recal-file "$recal_table" \
        -O "$recal_bam" \
     
     # Check if ApplyBQSR succeeded
    if [ $? -ne 0 ]; then
        echo "Error: ApplyBQSR failed for $sample"
        exit 1  # Exit the script on error
    else
        echo "Recalibrated BAM file created for $sample"
    fi
done

echo "Recalibrated BAM files generated for all rg_mkdup_bam files."

# BuildBamIndex, used to generate a BAM index ".bai" file.

# Iterate over each recal_bam file in the input directory
for recal_bam in "$input_dir2"/*.recal.bam; do

    # Extract sample name from the recal_bam file name
    sample=$(basename "$recal_bam" .recal.bam)

    # Output rg mkdp bam bai file for the sample
     recal_bam_bai="${output_dir}/${sample}.recal.bam.bai"
   
    # Run BuildBamIndex utility
    java -jar $Picard/picard.jar BuildBamIndex \
        I="$recal_bam" \
        O="$recal_bam_bai" 
    
    echo "Index file created for recalibrated BAM: $recal_bam_bai"
		
done

echo "Indexing completed for all recal_bam files."