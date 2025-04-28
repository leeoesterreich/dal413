#!/bin/bash

#SBATCH --job-name=ControlFREEC
#SBATCH --output=ControlFREEC_%j.out
#SBATCH --error=ControlFREEC_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu
#SBATCH --mem-per-cpu=2000

module purge
module load gcc/8.2.0
module load samtools/1.9
module load control-freec/11.5

chmod +r /bgfs/alee/LO_LAB/Personal/Daisong/Reference_genome/Mouse/GRCm38
freec -conf /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ControlFREEC/config_CPT6.txt