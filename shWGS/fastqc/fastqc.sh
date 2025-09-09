#!/bin/bash
#!/bin/bash
#SBATCH --job-name=FastQC_BAM
#SBATCH --output=fastqc_%j.out
#SBATCH --error=fastqc_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu
#SBATCH --mem-per-cpu=2000


# Loading modules
module purge
module load fastqc/0.11.7
module load multiqc/1.19

# Input directory
input_dir="/ix1/alee/LO_LAB/General/Lab_Data/20250610_ILC_ctDNA_shWGS_01_Daisong/bam"

# Output directory (in current location)
output_dir="fastqc_results"
mkdir -p "$output_dir"

# Check if there are any bam files
if ! ls "$input_dir"/*.bam 1> /dev/null 2>&1; then
    echo "No BAM files found in $input_dir"
    exit 1
fi

#Traversing BAM files in input directory
for bam_file in "$input_dir"/*.bam; do
    # Extract the filename without extension
    sample=$(basename "$bam_file" .bam)

    # Run FastQC and redirect both stdout and stderr to the log file
    fastqc -o "$output_dir" "$bam_file" 2>&1 | tee "$output_dir/${sample}_fastqc.log"

    echo "FastQC analysis completed for $sample"
done

echo "All FastQC analyses completed"

# MultiQC analysis
# Run MultiQC on the output directory
multiqc "$output_dir" -o "$output_dir"

echo "MultiQC analysis completed"
