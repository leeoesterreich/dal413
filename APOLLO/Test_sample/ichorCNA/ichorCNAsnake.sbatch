#!/bin/bash
#
#SBATCH --job-name=ichorCNAsnakemake
#SBATCH -N 1
#SBATCH --cpus-per-task=8 # Request that ncpus be allocated per process.
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH --output=/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/ichorCNA/ichorScript/ichorCNAsnake.out



module purge
module load gcc/8.2.0 r/3.6.0
module load hmmcopy_utils/dev
module load snakemake/5.4.0
snakemake -s /bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/ichorCNA/ichorScript/ichorCNA.snakefile --cores 16

