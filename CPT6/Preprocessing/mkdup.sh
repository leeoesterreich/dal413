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
mkdir -p "$designated_path/5_mkdup"

# Set input  directory
input_dir1="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/4_alignment_file"
input_dir2="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/5_mkdup"

# Set output directory
output_dir="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/5_mkdup"

# Set reference genome path
reference="/bgfs/alee/LO_LAB/Personal/Daisong/Reference_genome/Mouse/GRCm38.p4.genome.fasta"

# Iterate over each sorted BAM file in the input directory
for sorted_bam in "$input_dir1"/*_sorted.bam; do

    # Extract sample name from the sorted BAM file name
    sample=$(basename "$sorted_bam" _sorted.bam)

    # Output RG BAM file for the sample
    rg_bam="${output_dir}/${sample}_rg.bam"
       
       # Check if the output file already exists
       if [ -f "$rg_bam" ]; then
           echo "Skipping AddOrReplaceReadGroups for $sample, file already exists: $rg_bam"
           continue
       fi

    # Run AddOrReplaceReadGroups tool
    java -jar $Picard/picard.jar AddOrReplaceReadGroups \
        I="$sorted_bam" \
        O="$rg_bam" \
        RGID=1 \
        RGLB=library \
        RGPL=illumina \
        RGPU=unit \
        RGSM="$sample"
    
    echo "Processed BAM file created: $rg_bam"
	
done

echo "Processing completed for all BAM files."

# MarkDuplicates by Picard, locates and tags duplicate reads in a BAM file.

# Iterate over each RG BAM file in the input directory
for rg_bam in "$input_dir2"/*_rg.bam; do

    # Extract sample name from the RG BAM file name
    sample=$(basename "$rg_bam" _rg.bam)

    # Output marked duplicates BAM file for the sample
    rg_mkdp_bam="${output_dir}/${sample}rg.mkdp.bam"
        
        # Check if the mkdup.bam file exists 
  
        if [ -f "$rg_mkdp_bam" ]; then
           echo "Skipping MarkDuplicates for $sample, file already exists: $rg_mkdp_bam"
           continue
        fi
        
    # Run Picard's MarkDuplicates tool
    java -jar $Picard/picard.jar MarkDuplicates \
        I="$rg_bam" \
        O="$rg_mkdp_bam" \
        M="$output_dir/${sample}rg_dup_metrics.txt" \
        REMOVE_DUPLICATES=false \
        VALIDATION_STRINGENCY=SILENT \
        MAX_RECORDS_IN_RAM=500000 > "$output_dir/${sample}_MarkDuplicates.log" 2>&1

    # Check if MarkDuplicates succeeded
    if [ $? -ne 0 ]; then
        echo "Error: MarkDuplicates failed for $sample. Check the log file: $output_dir/${sample}_MarkDuplicates.log"
        exit 1
    else
        echo "Marked duplicates BAM file created: $rg_mkdp_bam"
    fi

done
   
echo "Marking duplicates completed for all samples."

# BuildBamIndex, used to generate a BAM index ".bai" file.

# Iterate over each rg_mkdp_bam file in the input directory
for rg_mkdp_bam in "$input_dir2"/*rg.mkdp.bam; do

    # Extract sample name from the rg_mkdp_bam file name
    sample=$(basename "$rg_mkdp_bam" rg.mkdp.bam)

    # Output rg mkdp bam bai file for the sample
     rg_mkdp_bam_bai="${output_dir}/${sample}rg.mkdp.bam.bai"
   
    # Run BuildBamIndex utility
    java -jar $Picard/picard.jar BuildBamIndex \
        I="$rg_mkdp_bam" \
        O="$rg_mkdp_bam_bai" 
    
    echo "Index file created: $rg_mkdp_bam_bai"
		
done

echo "Indexing completed for all rg_mkdup_bam files."