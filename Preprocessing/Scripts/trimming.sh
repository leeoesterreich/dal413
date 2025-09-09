#!/bin/bash
#SBATCH --job-name=trimming_ctDNA
#SBATCH --output=trimming_ctDNA_%j.out
#SBATCH --error=trimming_ctDNA_%j.err
#SBATCH --time=2:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu

#This script is used to trim low quality reads and perform QC analysis on trimmed reads.

#Author:Daisong_Liu/LeeOestereich lab

#Before running trimmomatic, have a look at the MultiQC Report of raw sequencing data.
#Trimmomatic, trimming tool to remove adapter and low-quality reads for Illumina NGS data.

# Load modules
module purge
module load trimmomatic/0.38
module load fastqc/0.11.7
module load multiqc/1.30

# Define designated paths
designated_path="/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/Preprocessing"
mkdir -p "$designated_path/2_trimmed_file"
mkdir -p "$designated_path/3_trimmed_fastqc_result"

# Array of multiple designated folders
designated_folders=("$designated_path/2_trimmed_file" "$designated_path/3_trimmed_fastqc_result")
directory_names=("paired" "unpaired")

# Create necessary folders
for folder in "${designated_folders[@]}"; do
    # Iterate over each directory name
    for dir_name in "${directory_names[@]}"; do
        # Create new directory in the current designated folder
        mkdir -p "$folder/$dir_name"
    done
done

# Input directories
input_dir1="$designated_path/Raw_files"
input_dir2="$designated_path/2_trimmed_file/paired"
input_dir3="$designated_path/3_trimmed_fastqc_result/paired"
input_dir4="$designated_path/2_trimmed_file/unpaired"
input_dir5="$designated_path/3_trimmed_fastqc_result/unpaired"

# Output directories
output_dir_paired="$designated_path/2_trimmed_file/paired"
output_dir_unpaired="$designated_path/2_trimmed_file/unpaired"
output_dir="$designated_path/3_trimmed_fastqc_result/paired"
output_dir2="$designated_path/3_trimmed_fastqc_result/unpaired"


# Define adapter file
adapter_file="$designated_path/adapter_file/adapter.fa"
# Echo if adapter file is not found
if [[ ! -f "$adapter_file" ]]; then
    echo "Adapter file not found at $adapter_file"
    exit 1
fi


#Loop through all paired-end FASTQ files in the input directory
echo "Trimming adapter and low quality bases for each sample..."

for forward_read in $input_dir1/*_R1.fq; do               
    # Extract the sample without extension
    sample=$(basename "$forward_read" _R1.fq)
	
    # Define reverse read
    reverse_read="${input_dir1}/${sample}_R2.fq" 
 
       # Skip trimming if the paired output already exists
    if [[ -f "${input_dir2}/${sample}_R1_paired.fastq.gz" && -f "${input_dir2}/${sample}_R2_paired.fastq.gz" ]]; then
        echo "Skipping trimming for $sample, already completed."
        continue
    fi  
	
    # Run Trimmomatic for paired-end reads 
    java -jar /ihome/crc/install/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 -phred33 \
        "$forward_read" "$reverse_read" \
        "${output_dir_paired}/${sample}_R1_paired.fastq.gz" "${output_dir_unpaired}/${sample}_R1_unpaired.fastq.gz" \
        "${output_dir_paired}/${sample}_R2_paired.fastq.gz" "${output_dir_unpaired}/${sample}_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:"$adapter_file":2:30:10 HEADCROP:0 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:15 MINLEN:75
    
    echo "Trimming completed for $sample"
done

echo "Trimming completed for all sample"


# FastQC for paired trimmed samples
echo "Analyzing fastq files for each paired trimmed sample..."
for fastq_file in $input_dir2/*paired.fastq.gz; do
    sample=$(basename "$fastq_file" paired.fastq.gz)

    # Skip FastQC if output already exists
    if [[ -f "$input_dir3/${sample}_fastqc.zip" ]]; then
        echo "Skipping FastQC for paired sample $sample, already completed."
        continue
    fi

    fastqc -t 4 -f fastq -o "$input_dir3" "$fastq_file"
    echo "FastQC analysis completed for paired sample $sample"
done

# FastQC for unpaired trimmed samples
echo "Analyzing fastq files for each unpaired trimmed sample..."
for fastq_file in $input_dir4/*unpaired.fastq.gz; do
    sample=$(basename "$fastq_file" unpaired.fastq.gz)

    # Skip FastQC if output already exists
    if [[ -f "$input_dir5/${sample}_fastqc.zip" ]]; then
        echo "Skipping FastQC for unpaired sample $sample, already completed."
        continue
    fi

    fastqc -t 4 -f fastq -o "$input_dir5" "$fastq_file"
    echo "FastQC analysis completed for unpaired sample $sample"
done

#Run MultiQC on the input directory
echo "Analyzing fastq files for all paired trimmed samples..."

multiqc "$input_dir3" -o "$output_dir"

echo "Analyzing fastq files for all unpaired trimmed samples..."

multiqc "$input_dir5" -o "$output_dir2"
echo "MultiQC analysis completed"
echo "Trimming and Trimmed_FastQC completed for all samples"