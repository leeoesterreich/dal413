#!/bin/bash

# This script is used to align the samples with the reference genome and check the quality of BAM files.

# Author=Daisong_Liu/LeeOestereich lab

# BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome
# BWA-MEM, for high-quality queries and better performance, as it is faster and more accurate
#SBATCH --job-name=Alignment
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 96:00:00
#SBATCH --mem=256G


echo "The first parameter: $1"
echo "The second parameter: $2"
# echo "The third parameter: $3"


# Module required
module purge
module load gcc/8.2.0
module load bwa/0.7.17
module load samtools/1.14
module load qualimap/2.2.2
module load multiqc/1.19


len1=${#1}

# 去掉最后16位后的子字符串
name=${1:0:len1-19}
echo "This is the trimmed name!!!! $name"

# Create new directory
designated_path1="/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample"
mkdir -p "$designated_path1/4_alignment_file/$name"

echo "agfasdgdfhzdfgvzsfa"

designated_path2="/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/4_alignment_file/$name"
mkdir -p "$designated_path2/Qualimap"

# Reference genome file
reference="/bgfs/alee/LO_LAB/Personal/Daisong/Reference_genome/Human_ref_hg38_p7/Human_genome_hg38_p7_chr_only.fna"

# Input directory for trimmed files and quality check
input_dir1="/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/2_trimmed_file/paired"
input_dir2="/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/4_alignment_file/$name"
input_dir3="/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/4_alignment_file/$name/Qualimap"

# Output directory for aligned SAM files
output_dir1="/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/4_alignment_file/$name"
output_dir2="/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/4_alignment_file/$name/Qualimap"

# List of sample names
samples=( "$name" )

# Traverse each sample
for sample in "${samples[@]}"; do

    # Input FASTQ files for the sample
    Forward_read="${input_dir1}/${sample}_R1_paired.fastq.gz"
    Reverse_read="${input_dir1}/${sample}_R2_paired.fastq.gz"

    # Output SAM file for the sample
    output_sam="${output_dir1}/${sample}.sam"  

    # Check if SAM file already exists
    if [ -f "$output_sam" ]; then
        echo "SAM file for $sample already exists, skipping alignment..."
    else     

        # Perform mapping using bwa-mem
        echo "Aligning sample: $sample"
        bwa mem -K 100000000 -t 64 $reference $Forward_read $Reverse_read > $output_sam
        if [ $? -ne 0 ]; then
            echo "Error: BWA alignment failed for $sample"
            exit 1
        fi      
        echo "Alignment completed for $sample"
    fi  
done

echo "All the samples were aligned successfully."

# samtools is a suite of programs for interacting with high-throughput sequencing data. 
# samtools view used to convert a sam file to a bam file
# samtools sort used to sort the sam/bam file by coordinates.
# samtools index used to create an index file for the sorted bam file.  

# Traverse each SAM file in the input directory
for sam_file in $input_dir2/*.sam; do

    # Extract sample name from the SAM file
    sample=$(basename "$sam_file" .sam)

    # Output BAM file
    bam_file="${output_dir1}/${sample}.bam"

    # Convert SAM to BAM
    echo "Converting $sam_file to BAM format..."
    samtools view -bS "$sam_file" > "$bam_file"

    echo "BAM file created: $bam_file"
done

echo "BAM files generated for all the SAM files."

# Traverse each BAM file in the input directory
for bam_file in $input_dir2/*.bam; do

    # Extract sample name from the BAM file
    sample=$(basename "$bam_file" .bam)             

    # Output Sorted BAM file for the sample
    sorted_bam="${output_dir1}/${sample}_sorted.bam"
	
	# Sort BAM file
    echo "Sorting $bam_file..."
    samtools sort "$bam_file" -o "$sorted_bam"

    echo "Sorted BAM file created: $sorted_bam"
done

echo "All samples were sorted successfully."

# Traverse each Sorted BAM file in the input directory
for sorted_bam_file in $input_dir2/*_sorted.bam; do

    # Extract sample name from the sorted BAM file
    sample=$(basename "$sorted_bam_file" _sorted.bam)

    # Output bai file for the sample
    output_bai="${output_dir1}/${sample}.bai"   
	
	# Index BAM file
    echo "Indexing $sorted_bam_file..."
    samtools index "$sorted_bam_file"

    echo "Indexing bai file created: $output_bai"
done

echo "All samples were converted sorted and indexed successfully."

# BamQC, used to perform quality control on BAM files

# Traverse each BAM file in the input directory
for bam_file in "$input_dir2"/*_sorted.bam; do

    # Extract sample name from the BAM file
    sample=$(basename "$bam_file" _sorted.bam)

     # Create an output directory for the sample
    sample_output_dir="${output_dir2}/${sample}"    


    # Run Qualimap
    qualimap bamqc \
        -bam "$bam_file" \
        -outdir "$sample_output_dir" \
        -outfile "${sample}_qualimap_report.pdf" \
        -outformat pdf \
        --java-mem-size=4G

    echo "Qualimap report generated for $sample"
done

echo "Qualimap analysis completed for all samples."

# Load MultiQC module
module load multiqc/1.19

# Run MultiQC on the output directory
multiqc "$input_dir3" -o "$output_dir2"

echo "MultiQC analysis completed."


