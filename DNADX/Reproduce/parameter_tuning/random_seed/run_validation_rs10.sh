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

# Create results directory
mkdir -p results_rs10/models
mkdir -p results_rs10/validation_summary
mkdir -p results_rs10/plots

# Load Python module
module load python/3.7.0

# Run the validation script
python METABRIC_validation.py > results_rs10/validation.log 2>&1

# Print completion message
echo "Validation process completed. Check results_rs10/validation.log for results" 