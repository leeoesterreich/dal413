#!/bin/bash

#SBATCH --job-name=NN_ALL_TRAIN
#SBATCH --output=nn_all_train_%j.out
#SBATCH --error=nn_all_train_%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu
#SBATCH --mem-per-cpu=4000

# Load Python module
module load python/3.7.0

# Move to the script directory
cd /ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training

# Run the neural network training script for all signatures
echo "Starting neural network training for all signatures..."
echo "Start time: $(date)"
python3 model_training_nn_all.py
echo "End time: $(date)"
echo "Neural network training completed!" 