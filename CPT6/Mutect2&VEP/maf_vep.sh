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

# Construct the input and output file paths
input_vcf="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Mutect2/tumor_vs_liver_maf/maf/tumor_vs_liver.vcf"
output_vcf="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Mutect2/tumor_vs_liver_maf/maf/tumor_vs_liver.vep.vcf"
  
# Run the VEP command
  vep -i "$input_vcf" \
      -o "$output_vcf" \
      --species mouse --cache --dir_cache /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/VEP/ --force_overwrite
done

