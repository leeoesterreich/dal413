#!/bin/bash
# Script to run the conda environment setup as a SLURM job

# General SLURM configurations
#SBATCH --job-name=ENV_SETUP
#SBATCH --output=env_setup_%j.out
#SBATCH --error=env_setup_%j.err
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4000

# Run the setup script
echo "Running conda environment setup script..."
bash /ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/METDX/00_setup_new_env.sh
echo "Setup script finished."
