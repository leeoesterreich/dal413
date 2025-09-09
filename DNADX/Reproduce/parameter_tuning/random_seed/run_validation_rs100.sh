#!/bin/bash
#SBATCH --job-name=METABRICVAL100
#SBATCH --output=metabricval100_%j.out
#SBATCH --error=metabricval100_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu
#SBATCH --mem-per-cpu=2000

# Create results directory
mkdir -p results_rs100/models
mkdir -p results_rs100/validation_summary
mkdir -p results_rs100/plots

# Load Python module
module load python/3.7.0

# Run the validation script
python METABRIC_validation_rs100.py > results_rs100/validation.log 2>&1

# Print completion message
echo "Validation process completed. Check results_rs100/validation.log for results" 