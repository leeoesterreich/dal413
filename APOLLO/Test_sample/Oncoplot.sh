#!/bin/bash

#SBATCH --job-name=Oncoplot
#SBATCH --output=Oncoplot_j%.out
#SBATCH --error=Oncoplot_%j.err
#SBATCH --time=1-00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=64
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu


module purge
module load gcc/12.2.0
module load r/4.4.0

chmod +x /bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/Oncoplot/Oncoplot.r
R --no-save  /bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/Oncoplot/Oncoplot.r
