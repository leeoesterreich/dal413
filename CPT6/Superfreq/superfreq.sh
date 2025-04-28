#!/bin/bash

#SBATCH --job-name=SuperFreq
#SBATCH --output=superfreq_%j.out
#SBATCH --error=superfreq_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu
#SBATCH --mem-per-cpu=2000

module load gcc/12.2.0
module load r/4.4.0
module load samtools/1.21
module load varscan/2.4.2
module load vep/95

chmod u+rx /ihome/alee/dal413/R/Superfreq/superfreq.R
Rscript /ihome/alee/dal413/R/Superfreq/superfreq.R


