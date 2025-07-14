#!/bin/bash

#SBATCH --job-name=METABRICVAL
#SBATCH --output=metabricval_%j.out
#SBATCH --error=metabricval_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu
#SBATCH --mem-per-cpu=2000

# Load Python module
module load python/3.7.0

# Move to the script directory
cd /ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation

# Run the METABRIC validation script
python METABRIC_validation.py 