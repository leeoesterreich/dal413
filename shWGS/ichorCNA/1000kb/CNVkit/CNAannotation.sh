#!/bin/bash
   
# This script is used to annotate cnv by cnvkit.

#SBATCH --job-name=cnvkit
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 1-00:00
 
# Modules required
module purge
module load cnvkit/0.9.5

# Use the Slurm submission directory as the base directory
BASE_DIR="/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/CNVkit"

# Input directory 
input_dir="$BASE_DIR/"
refFlatb38="$BASE_DIR/refFlat.txt"
cnvkit_dir="$BASE_DIR"
mkdir -p "$input_dir/temp_file"
temp_dir="$input_dir/temp_file"

# Output directory
output_dir="$BASE_DIR"

# Ensure the cnv_annotate.py script is executable
chmod +x $cnvkit_dir/cnv_annotate.py

# Process each .cna.seg file in the input directory
for cna_seg in "$input_dir"/*.cna.seg; do
    # Extract sample name from cna seg file name
    sample=$(basename "${cna_seg}" .cna.seg)

    file1="$cna_seg"
    file2="$temp_dir/${sample}_edit1.cnr"
    file3="$temp_dir/${sample}_edit2.cnr"
    file4="$temp_dir/${sample}_edit3.cnr"
    cnv_file="$temp_dir/${sample}_reheader.cnr"

    # Step 1: Append 0 to the third column
    awk '{$3=$3" "0; print }' OFS='\t' $file1 > $file2

    # Step 2: Replace "chr" with "chromosome", ".*.logR" with "log2" in the 7th column, "0" with "gene" in the 4th column for the header
    awk 'NR==1 {gsub("chr","chromosome");gsub(".*.logR","log2",$7);gsub("0","gene",$4);print};1' OFS='\t' $file2 > $file3

    # Step 3: Replace "0" with "-" in the 4th column
    awk '{gsub("0","-",$4);print}' OFS='\t' $file3 > $file4

    # Step 4: Remove the header from the file
    awk '{if (NR!=1) {print}}' OFS='\t' $file4 > $cnv_file

    # Annotate the CNV file
    $cnvkit_dir/cnv_annotate.py $refFlatb38 $cnv_file -o $output_dir/${sample}.annotated.cnr

    echo "Processing $sample..."
    
done

echo "Annotated CNR file generated for $sample"


