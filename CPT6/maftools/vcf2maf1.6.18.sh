#!/bin/bash

#SBATCH --job-name=vcf2maf
#SBATCH --output=vcf2maf_%j.out
#SBATCH --error=vcf2maf_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu
#SBATCH --mem-per-cpu=2000

module load gcc/8.2.0
module load vcf2maf/1.6.18
module load samtools/1.9
module load vep/95

VCF_FILE="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Mutect2/tumor_vs_liver_manual/tumor_vs_liver/filtered_tumor_vs_liver.vcf"
TUMOR_SAMPLE="B6ILC-cells"
NORMAL_SAMPLE="BL6-mouse-liver"
OUTPUT_MAF="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/maftools/filtered/tumor_vs_liver.maf"
REFERENCE_FASTA="/bgfs/alee/LO_LAB/Personal/Daisong/Reference_genome/Mouse/GRCm38.p4.genome/GRCm38.p4.genome.fasta"
VEP_PATH="/ihome/crc/install/vep/python3.7_vep95/bin"
VEP_CACHE="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/VEP"

export PERL5LIB=/ihome/crc/install/vep/python3.7_vep95/lib/site_perl/5.26.2:$PERL5LIB
/ihome/crc/install/vcf2maf/1.6.18/python3.7/bin/perl -MBio::Perl -e 'print "BioPerl is now accessible\n"'

/ihome/crc/install/vcf2maf/1.6.18/python3.7/bin/perl /ihome/crc/install/vcf2maf/1.6.18/python3.7/bin/vcf2maf.pl \
    --input-vcf $VCF_FILE \
    --output-maf $OUTPUT_MAF \
    --tumor-id $TUMOR_SAMPLE \
    --normal-id $NORMAL_SAMPLE \
    --vcf-tumor-id "B6ILC-cells" \
    --vcf-normal-id "BL6-mouse-liver" \
    --ref-fasta $REFERENCE_FASTA \
    --species mus_musculus \
    --vep-path $VEP_PATH \
    --vep-data $VEP_CACHE \
    --ncbi-build GRCm38 \
    --vep-forks 8
