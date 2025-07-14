#!/bin/bash
   
# This script is used to annotate cnv by cnvkit.

#SBATCH --job-name=cnvkit
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 1-00:00
 
# Modules required
module purge
module load cnvkit/0.9.5

# Use the current directory as the base directory
BASE_DIR="/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb_Broad/CNVkit"

# Input directory 
input_dir="$BASE_DIR"
refFlatb38="$BASE_DIR/../annotation/genes.bed"
cnvkit_dir="$BASE_DIR"
mkdir -p "$input_dir/temp_file"
temp_dir="$input_dir/temp_file"

# Output directory
output_dir="$BASE_DIR/annotated"
mkdir -p "$output_dir"

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

    # Step 1: Get the logR column name
    logr_col=$(head -n 1 "$file1" | tr '\t' '\n' | grep -n "logR" | grep -v "Copy_Number" | cut -d: -f1)

    # Step 2: Create CNR format with required columns
    awk -v logr_col="$logr_col" 'BEGIN {OFS="\t"}
        NR==1 {
            print "chromosome", "start", "end", "gene", "log2", "depth", "weight", "probes"
            next
        }
        {
            print $1, $2, $3, "-", $logr_col, "0", "1", "1"
        }' "$file1" > "$file2"

    # Step 3: Replace "chr" with "chromosome" in the data rows
    sed 's/^chr/chromosome/' "$file2" > "$file3"

    # Step 4: Remove any NA values in log2 column
    awk 'BEGIN {OFS="\t"}
        NR==1 {print; next}
        {
            if ($5 == "NA") $5 = "0"
            print
        }' "$file3" > "$cnv_file"

    # Annotate the CNV file
    $cnvkit_dir/cnv_annotate.py $refFlatb38 $cnv_file -o $output_dir/${sample}.annotated.cnr

    echo "Processing $sample..."
    
done

echo "All samples have been processed. Annotated CNR files are in $output_dir"

# Clean up temporary files
rm -rf "$temp_dir"


