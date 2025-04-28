#!/bin/bash

#SBATCH --job-name=vcf
#SBATCH --output=vcf_%j.out
#SBATCH --error=vcf_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu
#SBATCH --mem-per-cpu=2000

# module load gcc/8.2.0
# module load samtools/1.9
# module load varscan/2.4.2

#  samtools mpileup -d 1000 -q 15 -Q 15 -A -f /bgfs/alee/LO_LAB/Personal/Daisong/Reference_genome/Mouse/GRCm38.p4.genome/GRCm38.p4.genome.fasta  /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Preprocessing/4_mkdup/B6ILC-cellsrg.mkdp.bam > /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Preprocessing/6_vcf/tumor.mpileup 
#  java -jar VarScan.v2.4.2.jar mpileup2cns /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Preprocessing/6_vcf/tumor.mpileup --variants --strand-filter 0 --p-value 0.01 --min-var-freq 0.02 --output-vcf 1 > /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Preprocessing/6_vcf/tumor.vcf
#  samtools mpileup -d 1000 -q 15 -Q 15 -A -f /bgfs/alee/LO_LAB/Personal/Daisong/Reference_genome/Mouse/GRCm38.p4.genome/GRCm38.p4.genome.fasta  /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Preprocessing/4_mkdup/BL6-mouse-liverrg.mkdp.bam > /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Preprocessing/6_vcf/liver.mpileup 
#  java -jar VarScan.v2.4.2.jar mpileup2cns /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Preprocessing/6_vcf/liver.mpileup  --variants --strand-filter 0 --p-value 0.01 --min-var-freq 0.02 --output-vcf 1 > /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Preprocessing/6_vcf/liver.vcf


module load gcc/8.2.0
module load gatk/4.5.0.0
Picard=/ihome/crc/install/picard/2.18.12
java=/ihome/crc/install/java/jdk-21.0.2/bin/java

tumor_bam=/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Preprocessing/4_mkdup/B6ILC-cellsrg.mkdp.bam
normal_bam=/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Preprocessing/4_mkdup/BL6-mouse-liverrg.mkdp.bam
ref_dir=/bgfs/alee/LO_LAB/Personal/Daisong/Reference_genome/Mouse/GRCm38.p4.genome/GRCm38.p4.genome.fasta
out_dir=/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Preprocessing/6_vcf

gatk Mutect2 \
    -R "${ref_dir}" \
    -I "${normal_bam}" \
    -O "${out_dir}/liver.vcf" 