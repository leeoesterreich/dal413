#!/bin/bash
# Script to run XGBoost training

# General SLURM configurations
#SBATCH --job-name=XGB_GPU_TRAIN
#SBATCH --output=xgb_gpu_train_%j.out
#SBATCH --error=xgb_gpu_train_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu 
#SBATCH --mem-per-cpu=4000        

# GPU-specific configurations
#SBATCH --nodes=1
#SBATCH --time=24:00:00             # User-specified time
#SBATCH --ntasks-per-node=1         # Based on 1 GPUs per node
#SBATCH --gres=gpu:1                # User-specified 1 GPUs per node
#SBATCH --cluster=gpu
#SBATCH --partition=a100

# Script execution
echo "Loading Python 3.7.0 module..."
module load python/3.7.0

echo "Running XGBoost training script on GPU..."
# Using the absolute path to your script
python3.7 "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training/model_training_nn_all.py"

echo "Script execution finished." 