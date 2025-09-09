#!/bin/bash
   
# This script is used to annotate CNV segments with gene names.

#SBATCH --job-name=cnv_annotate
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH -t 1-00:00
#SBATCH --output=annotation_%j.out
#SBATCH --error=annotation_%j.err
 
# Base directory setup
BASE_DIR=$(pwd)
input_dir="$BASE_DIR/cna_seg_files"
bed_file="/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/annotation/genes.bed"

echo "Starting CNV annotation process"
echo "Input directory: $input_dir"
echo "BED file: $bed_file"

# Ensure the Python script is executable
chmod +x "$input_dir/cnv_annotate.py"

# Process each .cna.seg file in the input directory
for cna_seg in "$input_dir"/*.cna.seg; do
    if [ -f "$cna_seg" ]; then
        # Extract sample name from cna seg file name
        sample=$(basename "${cna_seg}" .cna.seg)
        echo "Processing $sample..."
        
        # Run the annotation
        python3 "$input_dir/cnv_annotate.py" \
            "$bed_file" \
            "$cna_seg" \
            -o "$input_dir/${sample}.annotated.cnr"
            
        echo "Completed processing $sample"
    fi
done

echo "All files have been processed."


