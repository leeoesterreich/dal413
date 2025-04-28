#!/bin/bash

#SBATCH --job-name=downlaod
#SBATCH --output=download_%j.out
#SBATCH --error=download_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu
#SBATCH --mem-per-cpu=2000

cd /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/VEP

wget https://ftp.ensembl.org/pub/release-95/variation/indexed_vep_cache/mus_musculus_merged_vep_95_GRCm38.tar.gz
tar -xvzf /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/VEP/mus_musculus_merged_vep_95_GRCm38.tar.gz
