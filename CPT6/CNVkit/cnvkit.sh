#!/bin/bash

#SBATCH --job-name=cnvkit
#SBATCH --output=cnvkit_%j.out
#SBATCH --error=cnvkit_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu
#SBATCH --mem-per-cpu=2000

module load cnvkit/0.9.5

cnvkit.py batch \
        /bgfs/alee/LO_LAB/Personal/Daisong/Sarek/FASTQ2SAREK/CPT6_bam/preprocessing/recalibrated/tumor/tumor.recal.bam \
        --normal /bgfs/alee/LO_LAB/Personal/Daisong/Sarek/FASTQ2SAREK/CPT6_bam/preprocessing/recalibrated/liver/liver.recal.bam \
        --targets /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/reference/Twist_Mouse_Exome_Target_Rev1_GRCm38_no_chr.bed \
        --annotate /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/reference/refFlat_no_chr.txt \
        --fasta /bgfs/alee/LO_LAB/Personal/Daisong/Reference_genome/Mouse/GRCm38.p4.genome/GRCm38.p4.genome.fasta \
        --output-reference /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/CNVkit/GRCm38_CKVkit.cnn \
        --drop-low-coverage \
        --output-dir /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/CNVkit/Sarek \
        --diagram \
        --scatter

