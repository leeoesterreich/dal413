#!/bin/bash

#SBATCH --job-name=ML
#SBATCH --output=ml_%j.out
#SBATCH --error=ml_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu
#SBATCH --mem-per-cpu=2000

# Load Python module
module load python/3.7.0

# Run the Python script
python association_analysis.py elasticnet_only 