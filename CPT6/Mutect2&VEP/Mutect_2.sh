#!/bin/bash

#SBATCH --job-name=Mutect2
#SBATCH --output=mutect2_%j.out
#SBATCH --error=mutect2_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu
#SBATCH --mem-per-cpu=2000

module load gcc/8.2.0
module load gatk/4.5.0.0
Picard=/ihome/crc/install/picard/2.18.12
java=/ihome/crc/install/java/jdk-21.0.2/bin/java

tumor_bam=/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Preprocessing/5_bed_extract/B6ILC-cells.bam
normal_bam=/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Preprocessing/5_bed_extract/BL6-mouse-liver.bam
ref_dir=/bgfs/alee/LO_LAB/Personal/Daisong/Reference_genome/Mouse/GRCm38.p4.genome/GRCm38.p4.genome.fasta
out_dir=/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Mutect2/tumor_vs_liver_maf

gatk Mutect2 \
    -R "${ref_dir}" \
    -I "${tumor_bam}" \
    -I "${normal_bam}" \
    -normal BL6-mouse-liver \
    -O "${out_dir}/tumor_vs_liver.vcf" \
    # --normal-lod 1.5 \
    # --tumor-lod-to-emit 2.0 \
    # --max-population-af 0.1 \
    # --af-of-alleles-not-in-resource 0.001 \
    # --min-base-quality-score 8 \
    # --initial-tumor-lod 1.5 \
