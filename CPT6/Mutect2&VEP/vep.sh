#!/bin/bash

#SBATCH --job-name=vep
#SBATCH --output=vep_mouse_%j.out
#SBATCH --error=vep_mouse_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu
#SBATCH --mem-per-cpu=2000


module load vep/95
folder=/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Mutect2/tumor_vs_liver_AF_5_DP_10/

# Loop through each subfolder in the specified folder
for subfolder in "$folder"/*/; do
  # Extract the subfolder name without the trailing slash
  subfolder_name=$(basename "$subfolder")
  
  # Construct the input and output file paths
  input_vcf="$subfolder/${subfolder_name}.mutect2.vcf"
  output_txt="$subfolder/${subfolder_name}.mutect2.txt"
  
# Run the VEP command
  vep -i "$input_vcf" \
      -o "$output_txt" \
      --species mouse --cache --dir_cache /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/VEP/ --force_overwrite
done

